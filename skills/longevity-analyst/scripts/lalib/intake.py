"""Inventory a member's raw folder into state.json.

Only format sniffing lives here: file extension, header lines, id patterns,
value ranges, column counts. Anything that needs meaning (which assay made a
FASTQ, which platform measured a protein matrix, which column is the member,
what a PDF says, whether the files are the member's) becomes a pending
judgment the agent resolves with `la.py assign` / `la.py labs add`.

Every analysable table is rewritten into a normalised single-sample file under
<workspace>/work/derived/ (clean ids, one value column); methods read only
those derived files, never the raw upload.
"""
from __future__ import annotations

import csv
import gzip
import io
import os
import re
import statistics
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import (EXIT_INPUT, EXIT_USAGE, LAError, add_judgment, file_by_id, now_iso,
                     resolve_judgment, sha256_file)

INDEX_EXT = {".bai", ".crai", ".tbi", ".csi", ".fai", ".md5"}
DOC_EXT = {".pdf", ".png", ".jpg", ".jpeg", ".heic", ".webp", ".tif", ".tiff", ".docx", ".doc"}
READ_MODALITIES = ["wgs", "wes", "metagenome_shotgun", "amplicon_16s", "wgbs", "em_seq", "rnaseq", "scrnaseq", "unknown"]
PROTEIN_PLATFORMS = ["ms_dia", "ms_dda", "olink", "somascan", "unknown"]
VALUE_SCALES = ["raw_intensity", "log2_intensity", "npx", "rfu", "cohort_z", "unknown"]
TISSUES = ["whole_blood", "saliva", "buccal", "pbmc", "other", "unknown"]
ASSEMBLIES = ["GRCh38", "GRCh37", "unknown"]

UNIPROT = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")
CPG = re.compile(r"^cg\d{8}(?:_[A-Z]{2}\d+)?$")          # EPIC v2 adds _BC11/_TC21 replicate suffixes
PROBE = re.compile(r"^(?:cg\d{8}|ch\.\d+\.\d+[A-Z]?|ch\.[A-Za-z0-9.]+)(?:_[A-Z]{2}\d+)?$")
PVAL_COL = re.compile(r"detection|p[._ -]?val|pvalue|p\.value|detp", re.I)

NAME_COLUMNS = ("marker", "item", "name", "indicator", "test", "testname", "analyte", "项目", "指标", "检验项目", "名称",
                "检查项目", "检测项目", "项目名称", "检验名称", "中文名称", "检验项目名称", "检测指标")
VALUE_COLUMNS = ("value", "result", "results", "结果", "值", "检验结果", "数值", "测定值", "检测结果", "结果值", "测量值", "结果数值")
UNIT_COLUMNS = ("unit", "units", "单位", "参考单位", "单位名称", "计量单位")
REF_COLUMNS = ("referencerange", "reference", "refrange", "range", "参考范围", "参考区间", "参考值", "正常范围", "生物参考区间")


def _hnorm(h: str) -> str:
    """Header comparison form: lower-case, brackets and spaces removed ('结果(数值)' -> '结果数值', 'Test Name' -> 'testname')."""
    return re.sub(r"[\s()（）\[\]_\-]", "", h.strip().lower())
HEAD_BYTES = 2 * 1024 * 1024


# ------------------------------------------------------------ reading
def _raw_head(path: Path, nbytes: int = HEAD_BYTES) -> bytes:
    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return fh.read(nbytes)
    with open(path, "rb") as fh:
        return fh.read(nbytes)


def _decode(raw: bytes) -> str:
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16", errors="replace")
    for enc in ("utf-8-sig", "gb18030"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def _head_lines(path: Path, n: int = 5000) -> List[str]:
    text = _decode(_raw_head(path))
    lines = text.splitlines()
    if len(text) >= HEAD_BYTES - 1 and lines:
        lines = lines[:-1]            # last line may be cut
    return lines[:n]


def _full_text(path: Path) -> str:
    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return _decode(fh.read())
    return _decode(path.read_bytes())


def _sep(line: str) -> str:
    counts = {s: line.count(s) for s in ("\t", ",", ";")}
    return max(counts, key=counts.get) if max(counts.values()) else ","


def _split(line: str, sep: str) -> List[str]:
    return [c.strip().strip('"').strip("'") for c in next(csv.reader([line], delimiter=sep))]


def _stem_ext(path: Path) -> Tuple[str, str]:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".g.vcf"):
        return name[:-6], ".g.vcf"
    return (name.rsplit(".", 1)[0], "." + name.rsplit(".", 1)[1]) if "." in name else (name, "")


# ------------------------------------------------------------ sniffers
def _vcf_header(path: Path) -> Dict[str, Any]:
    info: Dict[str, Any] = {"samples": [], "assembly": None, "gvcf": False, "provenance_lines": []}
    for line in _head_lines(path, 5000):
        if re.match(r"##(\w*Command|source|SAMPLE|sample|bcftools_\w+|GATKCommandLine|DeepVariant)", line):
            info["provenance_lines"].append(line[:200])
        mc = re.match(r"##contig=<ID=(?:chr)?([12]),.*?length=(\d+)", line)
        if mc:
            builds = {("1", "248956422"): "GRCh38", ("1", "249250621"): "GRCh37",
                      ("2", "242193529"): "GRCh38", ("2", "243199373"): "GRCh37"}
            b = builds.get((mc.group(1), mc.group(2)))
            if b:
                if info.get("_build_seen") and info["_build_seen"] != b:
                    info["assembly"] = None          # contigs disagree: ask
                    info["_conflict"] = True
                elif not info.get("_conflict"):
                    info["assembly"] = info["_build_seen"] = b
        if "<NON_REF>" in line or line.startswith("##GVCFBlock"):
            info["gvcf"] = True
        if line.startswith("#CHROM"):
            info["samples"] = line.split("\t")[9:]
            break
    return info


def _sniff_table(lines: List[str], full_body=None) -> Dict[str, Any]:
    full_body = full_body or (lambda: [l for l in lines if l.strip() and not l.startswith("#")])
    lines = [l for l in lines if l.strip()]
    comments = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    text = "\n".join(lines[:80])
    res: Dict[str, Any] = {"comments": comments[:20]}
    mpa = re.search(r"#\s*(mpa_v\w+)", text)
    if "s__" in text:
        res["kind"] = "metaphlan_profile"
        res["mpa_version"] = mpa.group(1) if mpa else None
        hdr = next((l.lstrip("#").strip() for l in comments if "clade" in l.lower()), None)
        if hdr is None and body and not body[0].startswith(("k__", "s__", "UNCLASSIFIED", "UNKNOWN")):
            hdr = body[0]
        cols = _split(hdr, _sep(hdr)) if hdr else []
        res["columns"] = cols
        meta_cols = {"clade_name", "ncbi_tax_id", "clade_taxid", "additional_species", "#clade_name", "coverage",
                     "estimated_number_of_reads_from_the_clade", "number_of_reads", "relative_abundance_reads"}
        if {"ncbi_tax_id", "clade_taxid"} & {c.lower() for c in cols} or "|t__SGB" in text:
            res["mpa_version_from_content"] = "mpa_v3_or_later"
        sample_cols = [c for c in cols[1:] if c.lower() not in meta_cols]
        res["sample_columns"] = sample_cols
        res["n_samples"] = max(1, len(sample_cols)) if cols else 1
        return res
    if not body:
        res["kind"] = "unknown_table"
        return res
    # header row may sit under a title row
    hdr_i = 0
    for i, l in enumerate(body[:6]):
        low = [_hnorm(c) for c in _split(l, _sep(l))]
        if any(c in NAME_COLUMNS for c in low) and any(c in VALUE_COLUMNS for c in low):
            hdr_i = i
            break
    sep = _sep(body[hdr_i])
    header = _split(body[hdr_i], sep)
    low = [_hnorm(h) for h in header]
    rows = [_split(l, sep) for l in body[hdr_i + 1:hdr_i + 5000]]
    first = [r[0] for r in rows if r]
    if first and sum(bool(CPG.match(x)) for x in first) >= 0.8 * len(first):
        value_cols = [i for i in range(1, len(header)) if not PVAL_COL.search(header[i])]
        vals = []
        for r in rows:
            for i in value_cols:
                try:
                    vals.append(float(r[i]))
                except (ValueError, IndexError):
                    pass
        res.update(kind="methylation_beta", sep=sep, header_row=hdr_i,
                   sample_columns=[header[i] for i in value_cols],
                   pvalue_columns=[header[i] for i in range(1, len(header)) if PVAL_COL.search(header[i])],
                   value_range=[min(vals), max(vals)] if vals else None,
                   value_scale="beta" if vals and 0 <= min(vals) and max(vals) <= 1 else "not_beta",
                   suffixed_probes=any("_" in x for x in first))
        res["n_samples"] = len(value_cols)
        return res
    def _acc(x: str) -> str:
        x = x.split(";")[0].strip()
        parts = x.split("|")
        return parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else x
    sid = next((c for c in ("sampleid", "sample", "samplename", "sample_id") if c in low), None)
    if sid and "npx" in low and ({"olinkid", "uniprot", "assay"} & set(low)):
        si = low.index(sid)
        samples = sorted({r[si] for r in (_split(l, sep) for l in full_body()[hdr_i + 1:]) if len(r) > si and r[si]})
        res.update(kind="olink_long", sep=sep, header_row=hdr_i, sample_columns=samples, n_samples=len(samples),
                   id_column=header[low.index("uniprot")] if "uniprot" in low else header[low.index("olinkid" if "olinkid" in low else "assay")],
                   value_column=header[low.index("npx")], sample_id_column=header[si])
        return res
    uni_cols = []
    for ci in range(min(len(header), 4)):
        col = [_acc(r[ci]) for r in rows if len(r) > ci and r[ci].strip()]
        if col and sum(bool(UNIPROT.match(x)) for x in col) >= 0.8 * len(col):
            uni_cols.append(ci)
    if uni_cols:
        numeric_cols = []
        for i in range(len(header)):
            if i in uni_cols:
                continue
            col = [r[i] for r in rows if len(r) > i and r[i] and r[i].upper() not in ("NA", "NAN", "N/A", "NULL", "-")]
            ok = 0
            for x in col:
                try:
                    float(x)
                    ok += 1
                except ValueError:
                    pass
            if ok and ok >= 0.8 * len(col):          # missing values do not disqualify a sample column
                numeric_cols.append(header[i])
        res.update(kind="protein_matrix", sep=sep, header_row=hdr_i, id_column=header[uni_cols[0]],
                   sample_columns=numeric_cols, n_samples=len(numeric_cols),
                   npx_hint=any("npx" in h for h in low))
        return res
    name_c = next((i for i, h in enumerate(low) if h in NAME_COLUMNS), None)
    val_c = next((i for i, h in enumerate(low) if h in VALUE_COLUMNS), None)
    if name_c is not None and val_c is not None:
        res.update(kind="lab_table", sep=sep, header_row=hdr_i)
        return res
    res["kind"] = "unknown_table"
    return res


def _xlsx_to_lines(path: Path) -> List[str]:
    try:
        import openpyxl  # type: ignore
    except ImportError:
        raise LAError("reading .xlsx needs openpyxl (`pip install openpyxl`), or save the sheet as CSV", EXIT_INPUT)
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    out = []
    for row in ws.iter_rows(values_only=True):
        cells = ["" if v is None else str(v) for v in row]
        if any(cells):
            buf = io.StringIO()
            csv.writer(buf, delimiter="\t").writerow(cells)
            out.append(buf.getvalue().rstrip("\r\n"))
    return out


def classify(path: Path) -> Dict[str, Any]:
    stem, ext = _stem_ext(path)
    rec: Dict[str, Any] = {"ext": ext}
    if path.stat().st_size == 0:
        return {**rec, "kind": "empty"}
    if ext in INDEX_EXT:
        rec.update(kind="index")
    elif ext in (".fastq", ".fq"):
        mate = re.search(r"(?:_|\.)(r?[12])(?:_\d+)?$", stem)
        rec.update(kind="reads_fastq", mate=(mate.group(1)[-1] if mate else None),
                   pair_key=str(path.parent.resolve() / (re.sub(r"(?:_|\.)(r?[12])(?:_\d+)?$", "", stem) if mate else stem)))
        _raw_head(path, 4096)          # a truncated or corrupt gz fails here, at intake
    elif ext in (".bam", ".cram"):
        rec.update(kind="alignment", format=ext[1:])
    elif ext in (".vcf", ".g.vcf"):
        rec.update(kind="variants_vcf", **_vcf_header(path))
        if ext == ".g.vcf":
            rec["gvcf"] = True
    elif ext == ".bcf":
        rec.update(kind="unsupported_format", reason="BCF: convert with `bcftools view -Oz` to .vcf.gz")
    elif ext == ".idat":
        chan = "Grn" if "_grn" in stem else ("Red" if "_red" in stem else None)
        rec.update(kind="methylation_idat", channel=chan, pair_key=re.sub(r"_(grn|red)$", "", stem))
    elif ext == ".adat":
        rec.update(kind="somascan_adat")
    elif ext == ".biom":
        rec.update(kind="amplicon_table")
    elif ext in (".h5ad", ".h5", ".loom"):
        rec.update(kind="single_cell_matrix")
    elif ext in (".mzml", ".raw", ".wiff", ".mzxml"):
        rec.update(kind="ms_raw")
    elif ext in DOC_EXT:
        rec.update(kind="document")
    elif ext in (".dcm", ".nii"):
        rec.update(kind="imaging")
    elif ext in (".csv", ".tsv", ".txt", ".tab"):
        rec.update(_sniff_table(_head_lines(path), lambda: [l for l in _full_text(path).splitlines() if l.strip() and not l.startswith("#")]))
    elif ext in (".xlsx", ".xlsm"):
        rec.update(_sniff_table(_xlsx_to_lines(path)))
        rec["from_xlsx"] = True
    else:
        rec.update(kind="unknown")
    return rec


# ------------------------------------------------------------ derived files
def _table_lines(f: Dict[str, Any]) -> List[str]:
    p = Path(f["path"])
    if f.get("from_xlsx"):
        return _xlsx_to_lines(p)
    return _full_text(p).splitlines()


def _write_derived_beta(f: Dict[str, Any], column: str, dst: Path) -> Dict[str, Any]:
    lines = [l for l in _table_lines(f) if l.strip() and not l.startswith("#")]
    sep = f["sep"]
    header = _split(lines[f["header_row"]], sep)
    ci = header.index(column)
    acc: Dict[str, List[float]] = {}
    seen_raw = set()
    n_na = 0
    for l in lines[f["header_row"] + 1:]:
        r = _split(l, sep)
        if not r or not PROBE.match(r[0]):
            continue
        if r[0] in seen_raw:
            raise LAError(f"{f['name']}: probe {r[0]} appears twice; collapse replicates at the lab (SeSAMe) and re-export", EXIT_INPUT)
        seen_raw.add(r[0])
        probe = r[0].split("_")[0]          # EPIC v2 suffixed replicates of one probe are averaged
        try:
            v = float(r[ci])
        except (ValueError, IndexError):
            n_na += 1
            continue
        if not 0 <= v <= 1:
            raise LAError(f"{f['name']}: value {v} for {r[0]} is outside [0,1]; not a β column", EXIT_INPUT)
        acc.setdefault(probe, []).append(v)
    groups = [vs for vs in acc.values() if len(vs) > 1]
    clash = [k for k, vs in acc.items() if len(vs) > 1 and max(vs) - min(vs) > 0.3]
    # real EPIC v2 arrays have ~2% of replicate groups differing by >0.1; a merged or 1-beta copy disagrees broadly
    if groups and len(clash) > 0.1 * len(groups):
        raise LAError(f"{f['name']}: {len(clash)} of {len(groups)} replicate probe groups disagree by more than 0.3 "
                      f"(e.g. {clash[0]}); this is not one sample's array", EXIT_INPUT)
    vals = sorted(v for vs in acc.values() for v in vs)
    if vals:
        low = sum(1 for v in vals if v < 0.02) / len(vals)
        mid = vals[len(vals) // 2]
        if low > 0.8 or mid < 0.05:
            raise LAError(f"{f['name']} column {column!r}: {low:.0%} of values are below 0.02 (median {mid:.3g}); "
                          "that is the shape of detection p-values, not β values", EXIT_INPUT)
    dst.parent.mkdir(parents=True, exist_ok=True)
    collapsed = sum(1 for v in acc.values() if len(v) > 1)
    with gzip.open(dst, "wt") as fh:
        fh.write("cpg,beta\n")
        for k, vs in acc.items():
            fh.write(f"{k},{sum(vs) / len(vs):.6g}\n")
    return {"probes": len(acc), "missing_values": n_na, "replicate_probes_averaged": collapsed}


def _write_derived_metaphlan(f: Dict[str, Any], column: Optional[str], dst: Path) -> Dict[str, Any]:
    lines = _table_lines(f)
    cols = f.get("columns") or []
    ci = cols.index(column) if column and column in cols else None
    if ci is None:
        for name in ("relative_abundance",):
            if name in [c.lower() for c in cols]:
                ci = [c.lower() for c in cols].index(name)
    n = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as out:
        out.write(f"#{f.get('mpa_version') or 'mpa_unknown'}\n#clade_name\trelative_abundance\n")
        for l in lines:
            if l.startswith("#") or not l.strip():
                continue
            parts = _split(l, "\t" if "\t" in l else _sep(l))
            if not parts or ("s__" not in parts[0] and not parts[0].startswith("k__")):
                continue
            idx = ci if ci is not None else len(parts) - 1
            if idx >= len(parts):
                continue
            try:
                v = float(parts[idx])
            except ValueError:
                continue
            out.write(f"{parts[0]}\t{v}\n")
            n += 1
    return {"rows": n}


def _write_derived_protein(f: Dict[str, Any], column: str, dst: Path) -> Dict[str, Any]:
    lines = [l for l in _table_lines(f) if l.strip() and not l.startswith("#")]
    sep = f["sep"]
    header = _split(lines[f["header_row"]], sep)
    idc, vc = header.index(f["id_column"]), header.index(column)
    n = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as out:
        out.write("protein\tvalue\n")
        for l in lines[f["header_row"] + 1:]:
            r = _split(l, sep)
            if len(r) <= max(idc, vc):
                continue
            out.write(f"{r[idc]}\t{r[vc]}\n")
            n += 1
    return {"rows": n}


def _write_derived_vcf(f: Dict[str, Any], sample: str, dst: Path) -> Dict[str, Any]:
    src = Path(f["path"])
    samples = f.get("samples") or []
    si = samples.index(sample)
    opener = (lambda p: io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")) if src.name.endswith(".gz") else (lambda p: open(p, encoding="utf-8"))
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with opener(src) as fh, gzip.open(dst, "wt") as out:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                out.write("\t".join(parts[:9] + [sample]) + "\n")
                continue
            if len(parts) < 10 + si:
                continue
            out.write("\t".join(parts[:9] + [parts[9 + si]]) + "\n")
            n += 1
    return {"records": n}


SCALE_SHAPE = {  # (median low, median high) that values on each declared scale plausibly have
    "cohort_z": (-2.0, 2.0), "npx": (-6.0, 20.0), "log2_intensity": (8.0, 45.0), "raw_intensity": (100.0, float("inf")),
    "rfu": (20.0, float("inf"))}


def _check_scale(f: Dict[str, Any], derived: Path) -> None:
    """A declared scale that the numbers contradict is refused (MS log2 intensities declared as cohort z, etc.)."""
    vals = []
    for line in open(derived, encoding="utf-8").read().splitlines()[1:]:
        try:
            vals.append(float(line.split("\t")[1]))
        except (ValueError, IndexError):
            pass
    if len(vals) < 5:
        return
    med = statistics.median(vals)
    lo, hi = SCALE_SHAPE.get(f["value_scale"], (-float("inf"), float("inf")))
    if not lo <= med <= hi:
        derived.unlink(missing_ok=True)
        raise LAError(f"{f['name']}: declared value scale {f['value_scale']} but the median value is {med:.3g}; "
                      "check the lab's documentation before answering", EXIT_INPUT)


def _write_derived_olink_long(f: Dict[str, Any], sample: str, dst: Path) -> Dict[str, Any]:
    lines = [l for l in _table_lines(f) if l.strip() and not l.startswith("#")]
    sep = f["sep"]
    header = _split(lines[f["header_row"]], sep)
    si, ii, vi = header.index(f["sample_id_column"]), header.index(f["id_column"]), header.index(f["value_column"])
    n = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as out:
        out.write("protein\tvalue\n")
        for l in lines[f["header_row"] + 1:]:
            r = _split(l, sep)
            if len(r) > max(si, ii, vi) and r[si] == sample:
                out.write(f"{r[ii]}\t{r[vi]}\n")
                n += 1
    return {"rows": n}


def _derived_dir(st: Dict[str, Any]) -> Path:
    return Path(st["workspace"]) / "work" / "derived"


def _register(st: Dict[str, Any], modality: str, f: Dict[str, Any], path: str, **meta: Any) -> None:
    lst = st["processed"].setdefault(modality, [])
    lst[:] = [x for x in lst if x.get("file_id") != f["id"]]
    entry = {"file_id": f["id"], "source_file_ids": [f["id"]], "path": path, "provenance": {"origin": "delivered", **meta}}
    if "/work/derived/" in path:
        entry["sha256"] = sha256_file(Path(path), limit=None)      # methods refuse a derived file changed after intake
    if f.get("provenance_uncertain"):
        entry["provenance_uncertain"] = True
    lst.append(entry)
    if len(lst) > 1:
        add_judgment(st, "primary", modality,
                     f"{len(lst)} {modality} files: which one is this visit's measurement? The others stay on record.",
                     [x["file_id"] for x in lst])


def _finish_table(st: Dict[str, Any], f: Dict[str, Any]) -> None:
    """Write the derived single-sample file once the sample column (and platform/tissue) are known."""
    k = f["kind"]
    col = f.get("sample_column")
    dd = _derived_dir(st)
    if k == "methylation_beta" and col and f.get("tissue"):
        info = _write_derived_beta(f, col, dd / f"{f['id']}.beta.csv.gz")
        f["derived"] = info
        if f.get("tissue") == "unknown":
            return
        _register(st, "methylation", f, str(dd / f"{f['id']}.beta.csv.gz"), platform="beta_table", tissue=f["tissue"],
                  sample_column=col, preprocessing="delivered β table: normalisation and cell-composition handling not stated by the lab",
                  **info)
        f["status"] = "ready"
    elif k == "metaphlan_profile" and (col or f.get("n_samples", 1) <= 1):
        dst = dd / f"{f['id']}.metaphlan.tsv"
        info = _write_derived_metaphlan(f, col, dst)
        tot = 0.0
        for line in open(dst, encoding="utf-8"):
            if line.startswith("#"):
                continue
            name, v = line.rstrip("\n").split("\t")
            if name.split("|")[-1].startswith("s__") and "|t__" not in name:
                tot += float(v)
        if not 50 <= tot <= 101:
            dst.unlink(missing_ok=True)
            raise LAError(f"{f['name']}: species abundances in column {col!r} sum to {tot:.4g}, not a percentage profile "
                          "(read counts or coverage?); pick the relative-abundance column", EXIT_INPUT)
        _register(st, "gut_profile", f, str(dst), platform="metaphlan",
                  mpa_version=f.get("mpa_version_from_content") or f.get("mpa_version"), sample_column=col, **info)
        f["status"] = "ready"
    elif k in ("protein_matrix", "olink_long") and col and f.get("platform") and f.get("value_scale"):
        if "unknown" in (f["platform"], f["value_scale"]):
            f["status"] = "deferred_ask_lab"
            f["reason"] = "protein platform or value scale unknown"
            return
        if k == "olink_long":
            info = _write_derived_olink_long(f, col, dd / f"{f['id']}.protein.tsv")
        else:
            info = _write_derived_protein(f, col, dd / f"{f['id']}.protein.tsv")
        _check_scale(f, dd / f"{f['id']}.protein.tsv")
        _register(st, "proteomics", f, str(dd / f"{f['id']}.protein.tsv"), platform=f["platform"],
                  value_scale=f["value_scale"], sample_column=col, **info)
        f["status"] = "ready"
    elif k == "variants_vcf" and f.get("assembly") and (col or len(f.get("samples") or []) == 1):
        if f["assembly"] == "unknown":
            f["status"] = "deferred_ask_lab"
            f["reason"] = "genome assembly unknown"
            return
        path = f["path"]
        if col and len(f.get("samples") or []) > 1:
            dst = dd / f"{f['id']}.{re.sub(r'[^A-Za-z0-9_.-]', '_', col)}.vcf.gz"
            _write_derived_vcf(f, col, dst)
            path = str(dst)
        _register(st, "variants", f, path, platform="vcf", assembly=f["assembly"], gvcf=f.get("gvcf", False),
                  sample=col or (f.get("samples") or [None])[0])
        f["status"] = "ready"


# ------------------------------------------------------------ intake
def run_intake(st: Dict[str, Any], raw_dir: Path) -> Dict[str, Any]:
    known_real = {os.path.realpath(f["path"]) for f in st["files"]}
    added = 0
    for path in sorted(raw_dir.rglob("*")):
        rel = path.relative_to(raw_dir)
        if any(part.startswith(".") for part in rel.parts):
            continue                                   # hidden files and folders
        try:
            if not path.is_file():
                continue
            real = os.path.realpath(path)
        except OSError:
            continue
        if real in known_real:
            continue                                   # a symlink to a file already on record
        known_real.add(real)
        fid = f"F{len(st['files']) + 1:03d}"
        f: Dict[str, Any] = {"id": fid, "path": str(path.resolve()), "name": path.name, "added_at": now_iso()}
        try:
            f["size"] = path.stat().st_size
            f["sha256"] = sha256_file(path)
            f.update(classify(path))
        except (OSError, EOFError, zlib.error, gzip.BadGzipFile, UnicodeError, LAError, ValueError, csv.Error, IndexError, StopIteration) as e:
            f.setdefault("size", 0)
            f.setdefault("sha256", "unreadable")
            f.update(kind="unreadable", status="unreadable", reason=f"{type(e).__name__}: {str(e)[:200]}")
        st["files"].append(f)
        added += 1
        twin_of = next((x for x in st["files"] if x is not f and x.get("excluded") and x.get("sha256") == f.get("sha256")), None)
        if twin_of:
            f.update(excluded={"reason": f"same content as excluded {twin_of['id']}", "at": now_iso()},
                     status="excluded_not_member", reason=f"与已排除文件 {twin_of['id']} 内容相同")
            continue
        if f.get("status") != "unreadable":
            try:
                _route(st, f)
            except LAError as e:
                f.update(status="rejected", reason=str(e)[:300])
            except (OSError, EOFError, zlib.error, gzip.BadGzipFile, UnicodeError, ValueError, csv.Error, IndexError, StopIteration) as e:
                f.update(status="unreadable", reason=f"{type(e).__name__}: {str(e)[:200]}")
    _pair_reads(st)
    if added:
        hints = []
        for f in st["files"]:
            if f.get("kind") == "variants_vcf" and f.get("samples"):
                hints.append(f"{f['id']} {f['name']}: VCF sample id(s) {f['samples'][:10]}")
                hints += [f"{f['id']} {f['name']}: {l}" for l in f.get("provenance_lines", [])[:5]]
            for c in f.get("comments", []):
                hints.append(f"{f['id']} {f['name']}: header line {c[:160]!r}")
        add_judgment(st, "identity", "member",
                     "Do all files belong to this member? Compare VCF sample ids, names/sex/age printed on documents you "
                     "transcribe, and file headers with the member record. Files from another person or marked synthetic "
                     "must be excluded (`--files`) and disclosed in the report.",
                     ["consistent", "inconsistent", "cannot_tell"], hints=hints)
    return {"added": added, "files": len(st["files"])}


def _ask_column(st: Dict[str, Any], f: Dict[str, Any], options: List[str], what: str) -> None:
    add_judgment(st, "sample_column", f["id"],
                 f"{f['name']}: {len(options)} {what} columns. Which one is this member? Using another person's column "
                 "silently gives that person's results.", options)


def _route(st: Dict[str, Any], f: Dict[str, Any]) -> None:
    k = f["kind"]
    fid = f["id"]
    if k == "reads_fastq" or k == "alignment":
        f["status"] = "needs_judgment"
        add_judgment(st, "modality", fid,
                     f"{f['name']}: which assay produced these reads? Read the lab's delivery note or ask the user; do not guess from the file name alone.",
                     READ_MODALITIES)
    elif k == "protein_matrix":
        f["status"] = "needs_judgment"
        hints = list(f.get("comments", [])) + (["a column name contains 'NPX' (a hint, not proof)"] if f.get("npx_hint") else [])
        add_judgment(st, "platform", fid,
                     f"{f['name']}: which platform measured these proteins? Olink/SomaScan clocks refuse MS data, so this decides which methods may run.",
                     PROTEIN_PLATFORMS, hints=hints)
        add_judgment(st, "value_scale", fid,
                     f"{f['name']}: what scale are the values on? Published protein clocks need a specific scale (often cohort z-scores).",
                     VALUE_SCALES, hints=hints)
        if f.get("n_samples", 0) > 1:
            _ask_column(st, f, f["sample_columns"], "sample")
        elif f.get("n_samples") == 1:
            f["sample_column"] = f["sample_columns"][0]
        else:
            f.update(status="rejected", reason="no numeric value column found")
    elif k == "methylation_beta":
        if f.get("value_scale") != "beta":
            f["status"] = "rejected"
            f["reason"] = f"values {f.get('value_range')} are not β in [0,1] (M-values or percent are refused, not rescaled)"
            return
        f["status"] = "needs_judgment"
        add_judgment(st, "tissue", fid, f"{f['name']}: which tissue was sampled? The bundled clocks are whole-blood models.", TISSUES)
        if f.get("n_samples", 0) > 1:
            _ask_column(st, f, f["sample_columns"], "β")
        else:
            f["sample_column"] = (f.get("sample_columns") or [None])[0]
    elif k == "methylation_idat":
        f["status"] = "needs_judgment"
        add_judgment(st, "tissue", fid, f"{f['name']}: which tissue was sampled?", TISSUES)
    elif k == "variants_vcf":
        f["status"] = "needs_judgment"
        if not f.get("assembly"):
            add_judgment(st, "assembly", fid, f"{f['name']}: the header does not say GRCh38 or GRCh37. Which build? Positions differ between builds.", ASSEMBLIES)
        n = len(f.get("samples") or [])
        if n == 0:
            f.update(status="rejected", reason="VCF has no sample column (sites-only); no genotypes to read")
            return
        if n > 1:
            _ask_column(st, f, f["samples"], "VCF sample")
        _finish_table(st, f)
    elif k == "metaphlan_profile":
        f["status"] = "needs_judgment"
        if f.get("n_samples", 1) > 1:
            _ask_column(st, f, f["sample_columns"], "sample")
        _finish_table(st, f)
    elif k == "lab_table":
        rows = parse_lab_table(f)
        for r in rows:
            r.update(source_file=fid, source="table", transcribed_by="parser")
        st["labs"].extend(rows)
        f["status"] = "ready"
        f["lab_rows"] = len(rows)
    elif k == "document":
        f["status"] = "needs_judgment"
        add_judgment(st, "transcribe", fid,
                     f"{f['name']}: read every page with multimodal vision, write the lab rows to a CSV (marker,value,unit,ref_range,page) exactly as printed, then run `la.py labs add`.",
                     ["transcribed", "no_lab_values", "unreadable"])
    elif k == "olink_long":
        f.update(platform="olink", value_scale="npx", status="needs_judgment")
        if f.get("n_samples", 0) > 1:
            _ask_column(st, f, f["sample_columns"], "SampleID")
        elif f.get("n_samples") == 1:
            f["sample_column"] = f["sample_columns"][0]
            _finish_table(st, f)
    elif k == "unknown_table":
        f["status"] = "needs_judgment"
        add_judgment(st, "transcribe", fid,
                     f"{f['name']}: a table whose columns were not recognised. If it holds lab values, write them to a CSV "
                     "(marker,value,unit,ref_range,page) and `la.py labs add`; otherwise answer no_lab_values.",
                     ["transcribed", "no_lab_values", "unreadable"])
    elif k in ("somascan_adat",):
        f["status"] = "unsupported"
        f["reason"] = "ADAT parsing is not in this version; ask the lab for a CSV export"
    elif k in ("index",):
        f["status"] = "ignored"
    elif k == "empty":
        f["status"] = "rejected"
        f["reason"] = "empty file"
    else:
        f["status"] = "unsupported"
        f["reason"] = f.get("reason") or f"no handler for {k} in this version; kept in the record, not analysed"


def parse_lab_table(f: Dict[str, Any]) -> List[Dict[str, Any]]:
    lines = [l for l in _table_lines(f) if l.strip() and not l.startswith("#")]
    sep = f.get("sep", ",")
    hi = f.get("header_row", 0)
    header = [_hnorm(c) for c in _split(lines[hi], sep)]
    name_c = next(i for i, h in enumerate(header) if h in NAME_COLUMNS)
    val_c = next(i for i, h in enumerate(header) if h in VALUE_COLUMNS)
    unit_c = next((i for i, h in enumerate(header) if h in UNIT_COLUMNS), None)
    ref_c = next((i for i, h in enumerate(header) if h in REF_COLUMNS), None)
    out = []
    for line in lines[hi + 1:]:
        r = _split(line, sep)
        if len(r) <= max(name_c, val_c) or not r[name_c] or not r[val_c]:
            continue
        out.append({"marker": r[name_c], "value": r[val_c],
                    "unit": (r[unit_c] if unit_c is not None and len(r) > unit_c else ""),
                    "ref_range": (r[ref_c] if ref_c is not None and len(r) > ref_c else "")})
    return out


def _pair_reads(st: Dict[str, Any]) -> None:
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for f in st["files"]:
        if f.get("kind") == "reads_fastq":
            groups.setdefault(f["pair_key"], []).append(f)
    for fs in groups.values():
        for f in fs:
            f["pair_ids"] = sorted(x["id"] for x in fs)


def assign(st: Dict[str, Any], fid: str, what: str, value: str, reason: str) -> Dict[str, Any]:
    f = file_by_id(st, fid)
    if f.get("excluded"):
        raise LAError(f"{fid} was excluded as not this member's data", EXIT_USAGE)
    if f.get("status") == "rejected" and what in ("platform", "value_scale", "sample_column", "tissue"):
        # the data contradicted the earlier answer; a corrected answer may be given (the old one stays in history)
        prev = next((j for j in st["judgments"] if j["kind"] == what and j["target"] == fid), None)
        if prev and prev["status"] == "resolved":
            f.setdefault("judgment_history", []).append(dict(prev))
            prev.update(status="pending")
            f["status"] = "needs_judgment"
            f.pop("reason", None)
    resolve_judgment(st, what, fid, value, reason)
    f[what] = value
    f.setdefault("assigned", {})[what] = {"value": value, "reason": reason, "at": now_iso()}
    if what == "modality":
        mates = [o for o in st["files"] if o["id"] in f.get("pair_ids", []) and o["id"] != fid]
        for o in mates:
            if any(j["kind"] == "modality" and j["target"] == o["id"] and j["status"] == "pending" for j in st["judgments"]):
                resolve_judgment(st, "modality", o["id"], value, f"paired with {fid}: {reason}")
            o["modality"] = value
        for x in [f, *mates]:
            if value == "unknown":
                x.update(status="deferred_ask_lab", reason=f"modality unknown: {reason}")
            else:
                x["status"] = "needs_pipeline"
                if value in ("wgbs", "em_seq"):
                    add_judgment(st, "tissue", x["id"], f"{x['name']}: which tissue was sampled?", TISSUES)
        return f
    if what == "tissue" and f["kind"] in ("methylation_idat", "reads_fastq"):
        if value == "unknown":
            f.update(status="deferred_ask_lab", reason=f"tissue unknown: {reason}")
        elif f["kind"] == "methylation_idat":
            f["status"] = "needs_pipeline"
        return f
    if what == "tissue" and value == "unknown":
        f.update(status="deferred_ask_lab", reason=f"tissue unknown: {reason}")
        return f
    if what == "platform" and value in ("olink", "somascan"):
        ms_re = (r"(?<![A-Za-z0-9])(?:DIA|DDA|LFQ|iBAQ|MaxQuant|Spectronaut|DIA-NN|FragPipe|SWATH|Orbitrap|Q\s?Exactive|timsTOF|"
                 r"LC-?MS(?:/MS)?|mass\s*spec\w*|data[- ]independent|data[- ]dependent|proteinGroups)(?![A-Za-z0-9])|质谱")
        seen = list(f.get("comments", [])) + [f["name"]] + [str(c) for c in (f.get("columns") or [])[:200]]   # sample ids (DIA_001) are not platform evidence
        # a line that also names Olink/NPX/SomaScan is about the affinity assay ("validated vs mass spectrometry")
        ms_words = [c for c in seen if re.search(ms_re, c, re.I) and not re.search(r"(?i)olink|npx|soma\s*scan|soma\s*logic", c)]
        if ms_words:
            prev = next(j for j in st["judgments"] if j["kind"] == "platform" and j["target"] == fid)
            prev.update(status="pending")
            f.pop("platform", None)
            raise LAError(f"{f['name']}: the file header says mass spectrometry ({ms_words[0][:80]!r}); "
                          f"it cannot be {value}. Answer ms_dia/ms_dda, or unknown if the lab disagrees", EXIT_USAGE)
    if what == "platform" and value == "unknown":
        f.update(status="deferred_ask_lab", reason=f"platform unknown: {reason}")
        return f
    if what == "primary":
        lst = st["processed"].get(fid, [])
        for x in lst:
            x["primary"] = x["file_id"] == value
        return f
    try:
        _finish_table(st, f)
    except LAError as e:                      # the data contradicts what it claims to be: refuse, say why, move on
        f.update(status="rejected", reason=str(e)[:300])
    except (OSError, EOFError, zlib.error, gzip.BadGzipFile, UnicodeError, ValueError, csv.Error, IndexError) as e:
        f.update(status="unreadable", reason=f"{type(e).__name__}: {str(e)[:200]}")
    return f


def choose_primary(st: Dict[str, Any], modality: str, file_id: str, reason: str) -> None:
    lst = st["processed"].get(modality) or []
    keys = {x.get("file_id") or x["path"] for x in lst}
    if file_id not in keys:
        raise LAError(f"{file_id} is not one of the {modality} entries: {sorted(keys)}", EXIT_USAGE)
    if any(j["kind"] == "primary" and j["target"] == modality and j["status"] == "pending" for j in st["judgments"]):
        resolve_judgment(st, "primary", modality, file_id, reason)
    elif not reason.strip():
        raise LAError("--reason is required", EXIT_USAGE)
    for x in st["processed"].get(modality, []):
        x["primary"] = (x.get("file_id") or x["path"]) == file_id


def add_labs(st: Dict[str, Any], csv_path: Path, fid: str, reason: str) -> int:
    f = file_by_id(st, fid)
    if f.get("excluded"):
        raise LAError(f"{fid} was excluded as not this member's data; its values cannot be added", EXIT_USAGE)
    if f.get("kind") not in ("document", "lab_table", "unknown_table", "spreadsheet"):
        raise LAError(f"{fid} is a {f.get('kind')}, not a lab document; lab rows come from reports and lab tables only", EXIT_USAGE)
    if not csv_path.exists():
        raise LAError(f"{csv_path} not found", EXIT_USAGE)
    text = _decode(csv_path.read_bytes())
    rows = list(csv.DictReader(io.StringIO(text)))
    need = {"marker", "value", "unit"}
    if not rows or not need <= set(rows[0].keys()):
        raise LAError(f"{csv_path} needs columns marker,value,unit (+page); got {list(rows[0].keys()) if rows else []}", EXIT_INPUT)
    bad = [i + 2 for i, r in enumerate(rows) if not (r.get("marker") or "").strip() or not (r.get("value") or "").strip()]
    if bad:
        raise LAError(f"rows {bad[:10]} have an empty marker or value; transcribe exactly or drop the row", EXIT_INPUT)
    for r in rows:
        st["labs"].append({"marker": r["marker"].strip(), "value": r["value"].strip(),
                           "unit": (r.get("unit") or "").strip(), "page": r.get("page", ""),
                           "ref_range": (r.get("ref_range") or r.get("reference") or "").strip(),
                           "source_file": fid, "source": "transcription", "transcribed_by": "agent"})
    if any(j["kind"] == "transcribe" and j["target"] == fid and j["status"] == "pending" for j in st["judgments"]):
        resolve_judgment(st, "transcribe", fid, "transcribed", reason)
    f["status"] = "ready"
    f["lab_rows"] = f.get("lab_rows", 0) + len(rows)
    return len(rows)


def fix_lab(st: Dict[str, Any], marker: str, action: str, value: Optional[str], unit: Optional[str], reason: str,
            row: Optional[int] = None) -> int:
    """Correct a lab row that was read wrongly (e.g. '6.8↑' as printed), or drop it. Keeps the old value."""
    if not reason.strip():
        raise LAError("--reason is required", EXIT_USAGE)
    if action == "fix" and value is None and unit is None:
        raise LAError("labs fix needs --value and/or --unit", EXIT_USAGE)
    hit = [r for r in st["labs"] if r["marker"] == marker]
    if not hit:
        raise LAError(f"no lab row named {marker!r}", EXIT_USAGE)
    if len(hit) > 1 and row is None:
        raise LAError(f"{len(hit)} rows are named {marker!r}; pick one with --row 1..{len(hit)}", EXIT_USAGE)
    if row is not None:
        if not 1 <= row <= len(hit):
            raise LAError(f"--row must be 1..{len(hit)}", EXIT_USAGE)
        hit = [hit[row - 1]]
    for r in hit:
        r.setdefault("history", []).append({"value": r["value"], "unit": r.get("unit", ""), "at": now_iso(), "reason": reason})
        if action == "drop":
            r["dropped"] = True
        else:
            if value is not None:
                r["value"] = value
            if unit is not None:
                r["unit"] = unit
    if action == "drop":
        st["labs"] = [r for r in st["labs"] if not r.get("dropped")]
        st.setdefault("dropped_labs", []).extend(hit)
    return len(hit)


def exclude_files(st: Dict[str, Any], fids: List[str], reason: str) -> List[str]:
    """Files judged not to be this member's: kept on record, never analysed, and nothing derived from them survives."""
    done = []
    for fid in fids:
        f = file_by_id(st, fid)
        f["excluded"] = {"reason": reason, "at": now_iso()}
        f["status"] = "excluded_not_member"
        f["reason"] = f"不属于本人：{reason}"
        done.append(fid)
    ex = set(done)
    ddir = Path(st.get("workspace", ".")) / "work" / "derived"
    for p in ddir.glob("*"):
        if p.name.split(".")[0] in ex:
            p.unlink(missing_ok=True)
    from . import pipelines as _pl
    for run in st.get("pipelines", []):
        if set(run.get("file_ids", [])) & ex or set(run.get("source_file_ids", [])) & ex:
            if _pl.run_alive(run):
                _pl.stop(st, run["run_id"])          # a run on someone else's data is stopped, not left to finish
            run.update(status="void_excluded", blocked="input files were excluded as not this member's")
    for mod, lst in list(st["processed"].items()):
        st["processed"][mod] = [x for x in lst if not (set(x.get("source_file_ids") or [x.get("file_id")]) & ex)]
    st["labs"] = [l for l in st["labs"] if l.get("source_file") not in ex]
    for j in st["judgments"]:
        if j["target"] in ex and j["status"] == "pending":
            j.update(status="resolved", answer="excluded", reason=f"file excluded: {reason}", resolved_at=now_iso())
    return done


def mark_uncertain(st: Dict[str, Any], fids: List[str], reason: str) -> None:
    """Files whose ownership cannot be confirmed: analysed, but every readout from them is labelled as such."""
    for fid in fids:
        f = file_by_id(st, fid)
        f["provenance_uncertain"] = {"reason": reason, "at": now_iso()}
    for lst in st["processed"].values():
        for x in lst:
            if set(x.get("source_file_ids") or []) & set(fids):
                x["provenance_uncertain"] = True
