"""Genotype ↔ phenotype: the member's genotype at the loci that published GWAS tie to a lab analyte, a risk-allele
count placed in the East Asian distribution, and a ClinVar scan of the member's own variants in monogenic genes.

The harness gathers and computes; whether an abnormal value is "mainly inherited" is the agent's judgment
(workflows/04c-insights.md). Rules that keep a statement from being falsely reassuring:
- a site the VCF does not cover is 'not called', never reference, unless a gVCF reference block covers it or the
  agent recorded that absent sites are reference for this delivery (`--absent-as-ref "<reason>"`); a position that
  has any record (another allele, a failed call, a carried deletion over it) is never "absent";
- a monogenic gene counts as scanned only when the VCF covers its region (gVCF blocks/calls, or the recorded
  absent-as-reference judgment) and ClinVar was retrieved; otherwise it is 'not scanned', never 0;
- ClinVar matching is by position and alleles for SNVs and by NCBI canonical SPDI for indels (dup / ins / delins /
  unaligned calls of one change compare equal), never by an HGVS spelling;
- ClinVar P/LP counts only with a unanimous aggregate classification and review criteria; other pathogenic mentions
  are listed separately; a heterozygous hit in a recessive gene is a carrier finding.
"""
from __future__ import annotations

import bisect
import gzip
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from . import pubdata
from .common import EXIT_INPUT, LAError, data, now_iso, skillkit, write_json

CHR = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|M|MT)$", re.I)
END = re.compile(r"(?:^|;)END=(\d+)")
PLP_OK = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
GENE_PAD = 2000            # promoter / UTR margin around the Ensembl gene span
INDEL_WINDOW = 100         # an indel farther than this from every ClinVar indel record cannot be one of them
MIN_REGION_COVERAGE = 0.9


def _norm_chrom(c: str) -> str:
    m = CHR.match(str(c))
    return m.group(1).upper().replace("MT", "M") if m else str(c)


def _is_block(f: List[str]) -> bool:
    return bool(END.search(f[7])) and f[4] in ("<NON_REF>", "<*>", ".")


def _gt(fmt: List[str], val: List[str]) -> List[str]:
    return re.split(r"[/|]", dict(zip(fmt, val)).get("GT", "./."))


# ------------------------------------------------------------------ VCF access
class Vcf:
    """Genotype lookups on one sample: tabix when an index and the tabix binary both exist, otherwise one streaming
    pass that keeps only rows at the wanted positions/regions, deletions spanning them, and gVCF blocks overlapping them
    (a whole-genome gVCF is never held in memory)."""

    def __init__(self, path: Path, rules: Dict[str, Any], absent_is_ref: bool = False):
        self.path = Path(path)
        self.rules = rules
        self.absent_is_ref = absent_is_ref
        indexed = self.path.suffix == ".gz" and (Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists())
        self._tabix = indexed and shutil.which("tabix") is not None
        self._rows: Dict[Tuple[str, int], List[List[str]]] = {}
        self._spans: Dict[str, List[List[str]]] = {}
        self._blocks: Dict[str, List[Tuple[int, int, List[str], List[str], str]]] = {}
        self._bstarts: Dict[str, List[int]] = {}
        self._loaded = False

    def _open(self):
        return gzip.open(self.path, "rt") if self.path.suffix == ".gz" else open(self.path)

    def prefetch(self, positions: Set[Tuple[str, int]], regions: List[Tuple[str, int, int]]) -> None:
        if self._tabix or self._loaded:
            return
        pos_by: Dict[str, List[int]] = {}
        for c, p in positions:
            pos_by.setdefault(_norm_chrom(c), []).append(int(p))
        for v in pos_by.values():
            v.sort()
        reg: Dict[str, List[Tuple[int, int]]] = {}
        for c, s, e in regions:
            reg.setdefault(_norm_chrom(c), []).append((int(s), int(e)))

        def wanted_in(c: str, s: int, e: int) -> bool:
            ps = pos_by.get(c)
            if ps:
                i = bisect.bisect_left(ps, s)
                if i < len(ps) and ps[i] <= e:
                    return True
            return any(rs <= e and s <= re_ for rs, re_ in reg.get(c, ()))

        with self._open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                c, pos = _norm_chrom(f[0]), int(f[1])
                if _is_block(f):
                    e = int(END.search(f[7]).group(1))
                    if wanted_in(c, pos, e):
                        self._blocks.setdefault(c, []).append((pos, e, f[8].split(":"), f[9].split(":"), f[6]))
                    continue
                if wanted_in(c, pos, pos):
                    self._rows.setdefault((c, pos), []).append(f)
                if len(f[3]) > 1 and wanted_in(c, pos + 1, pos + len(f[3]) - 1):
                    self._spans.setdefault(c, []).append(f)
        for c, bl in self._blocks.items():
            bl.sort(key=lambda b: b[0])
            self._bstarts[c] = [b[0] for b in bl]
        self._loaded = True

    def _tabix_rows(self, chrom: str, start: int, end: int) -> List[List[str]]:
        for c in (f"chr{chrom}", chrom):
            try:
                r = subprocess.run(["tabix", str(self.path), f"{c}:{start}-{end}"], capture_output=True, text=True, timeout=120)
            except (OSError, subprocess.TimeoutExpired) as e:
                raise LAError(f"tabix failed on {self.path.name}: {e}", EXIT_INPUT)
            if r.returncode != 0:
                raise LAError(f"tabix failed on {self.path.name}: {r.stderr.strip()[:200]}", EXIT_INPUT)
            rows = [l.split("\t") for l in r.stdout.splitlines() if l]
            if rows:
                return rows
        return []

    def _need_load(self) -> None:
        if not self._tabix and not self._loaded:
            raise LAError("internal: Vcf.prefetch must run before lookups without tabix", EXIT_INPUT)

    def _blocks_over(self, chrom: str, s: int, e: int) -> List[Tuple[int, int, List[str], List[str], str]]:
        bl, starts = self._blocks.get(chrom, []), self._bstarts.get(chrom, [])
        i = bisect.bisect_right(starts, e)
        out = []
        for b in reversed(bl[:i]):
            if b[1] >= s:
                out.append(b)
            elif b[1] < s - 10_000_000:
                break
        return out

    def _site(self, chrom: str, pos: int):
        """(rows starting at pos, rows spanning pos from upstream, gVCF blocks covering pos)."""
        if self._tabix:
            got = self._tabix_rows(chrom, pos, pos)
            blocks = [(int(f[1]), int(END.search(f[7]).group(1)), f[8].split(":"), f[9].split(":"), f[6]) for f in got if _is_block(f)]
            rows = [f for f in got if not _is_block(f) and int(f[1]) == pos]
            spans = [f for f in got if not _is_block(f) and int(f[1]) < pos <= int(f[1]) + len(f[3]) - 1]
            return rows, spans, [b for b in blocks if b[0] <= pos <= b[1]]
        self._need_load()
        spans = [f for f in self._spans.get(chrom, []) if int(f[1]) < pos <= int(f[1]) + len(f[3]) - 1]
        return self._rows.get((chrom, pos), []), spans, self._blocks_over(chrom, pos, pos)

    def quality(self, fmt: List[str], val: List[str], filt: str) -> Tuple[bool, str]:
        if filt not in ("PASS", ".", ""):
            return False, f"filtered:{filt}"
        d = dict(zip(fmt, val))
        gq, dp = d.get("GQ"), d.get("DP")
        try:
            if gq not in (None, ".") and float(gq) < self.rules["min_gq"]:
                return False, "low_gq"
            if dp not in (None, ".") and float(dp) < self.rules["min_dp"]:
                return False, "low_dp"
        except ValueError:
            return False, "bad_quality_field"
        return True, "called" if (gq not in (None, ".") and dp not in (None, ".")) else "quality_unverified"

    def _upstream_deletion(self, spans: List[List[str]]) -> Optional[str]:
        """A deletion from upstream that the member carries (or whose call failed) makes this position unknown."""
        for f in spans:
            fmt, val = f[8].split(":"), f[9].split(":")
            gt = _gt(fmt, val)
            alts = f[4].split(",")
            ok, why = self.quality(fmt, val, f[6])
            carried = [g for g in gt if g not in (".", "0")]
            if "." in gt or not ok:
                return why if not ok else "no_call"
            if any(g.isdigit() and int(g) <= len(alts) and len(alts[int(g) - 1]) < len(f[3]) for g in carried):
                return "within_carried_deletion"
        return None

    def _absent(self, chrom: str, pos: int, rows, blocks, ref_copies: bool) -> Dict[str, Any]:
        full = 2 if ref_copies else 0
        if rows:                                              # a record exists here but not for this allele / REF
            return {"dosage": None, "status": "ref_mismatch"}
        for s, e, fmt, val, filt in blocks:
            ok, why = self.quality(fmt, val, filt)
            return {"dosage": full if ok else None, "status": "ref_block" if ok else why}
        if self.absent_is_ref:
            return {"dosage": full, "status": "assumed_ref"}
        return {"dosage": None, "status": "not_called"}

    def alt_dosage(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """Copies of this exact alt allele. Every record at the position is read (split multi-allelic rows too)."""
        chrom = _norm_chrom(chrom)
        rows, spans, blocks = self._site(chrom, pos)
        up = self._upstream_deletion(spans)
        if up:
            return {"dosage": None, "status": up}
        same = [f for f in rows if f[3].upper() == ref.upper()]
        other, fail = None, None
        for f in same:
            alts = [a.upper() for a in f[4].split(",")]
            fmt, val = f[8].split(":"), f[9].split(":")
            ok, why = self.quality(fmt, val, f[6])
            gt = _gt(fmt, val)
            if alt.upper() in alts:
                if not ok:
                    return {"dosage": None, "status": why}
                if "." in gt or len(gt) != 2:
                    return {"dosage": None, "status": "no_call"}
                idx = str(alts.index(alt.upper()) + 1)
                return {"dosage": sum(1 for g in gt if g == idx), "status": why}
            if not ok or "." in gt or len(gt) != 2:
                fail = fail or (why if not ok else "no_call")
            else:
                other = {"dosage": 0, "status": "called_other_alt"}
        if fail:
            return {"dosage": None, "status": fail}          # another row here failed: this allele is not known absent
        if other:
            return other
        return self._absent(chrom, pos, rows, blocks, ref_copies=False)

    def ref_dosage(self, chrom: str, pos: int, ref: str) -> Dict[str, Any]:
        """Copies of the reference base: 2 minus every alt allele called at the position, over all rows."""
        chrom = _norm_chrom(chrom)
        rows, spans, blocks = self._site(chrom, pos)
        up = self._upstream_deletion(spans)
        if up:
            return {"dosage": None, "status": up}
        same = [f for f in rows if f[3].upper() == ref.upper()]
        if same:
            alt_copies, status = 0, "called"
            for f in same:
                fmt, val = f[8].split(":"), f[9].split(":")
                ok, why = self.quality(fmt, val, f[6])
                gt = _gt(fmt, val)
                if not ok:
                    return {"dosage": None, "status": why}
                if "." in gt or len(gt) != 2:
                    return {"dosage": None, "status": "no_call"}
                alt_copies += sum(1 for g in gt if g != "0")
                status = why if why != "called" else status
            return {"dosage": max(0, 2 - alt_copies), "status": status}
        return self._absent(chrom, pos, rows, blocks, ref_copies=True)

    def _region_rows(self, chrom: str, start: int, end: int) -> List[List[str]]:
        if self._tabix:
            return self._tabix_rows(chrom, start, end)
        self._need_load()
        rows = [f for (c, p), fs in self._rows.items() if c == chrom and start <= p <= end for f in fs]
        return rows + [["", str(s), ".", "N", "<NON_REF>", ".", filt, f"END={e}", ":".join(fm), ":".join(vl)]
                       for s, e, fm, vl, filt in self._blocks_over(chrom, start, end)]

    def region_coverage(self, chrom: str, start: int, end: int) -> float:
        """Share of the region's bases the VCF speaks for: good gVCF reference blocks plus called records."""
        if self.absent_is_ref:
            return 1.0
        chrom = _norm_chrom(chrom)
        ivs = []
        for f in self._region_rows(chrom, start, end):
            fmt, val = f[8].split(":"), f[9].split(":")
            ok, _ = self.quality(fmt, val, f[6])
            if not ok:
                continue
            s = int(f[1])
            e = int(END.search(f[7]).group(1)) if _is_block(f) else s + max(len(f[3]), 1) - 1
            ivs.append((max(s, start), min(e, end)))
        ivs.sort()
        covered, cur_s, cur_e = 0, None, None
        for s, e in ivs:
            if cur_e is None or s > cur_e + 1:
                if cur_e is not None:
                    covered += cur_e - cur_s + 1
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        if cur_e is not None:
            covered += cur_e - cur_s + 1
        return covered / max(1, end - start + 1)

    def carried(self, chrom: str, start: int, end: int) -> List[Dict[str, Any]]:
        """Alt alleles the member carries in a region, each with zygosity (het / hom / hemi) and a quality verdict."""
        chrom = _norm_chrom(chrom)
        out = []
        for f in self._region_rows(chrom, start, end):
            if len(f) < 10 or _is_block(f):
                continue
            fmt, val = f[8].split(":"), f[9].split(":")
            d = dict(zip(fmt, val))
            gt = _gt(fmt, val)
            present = sorted({g for g in gt if g not in (".", "0")})
            if not present:
                continue
            ok, why = self.quality(fmt, val, f[6])
            if ok and "." in gt:
                ok, why = False, "half_call"
            alts = f[4].split(",")
            ad = d.get("AD", "")
            for g in present:
                if not g.isdigit():
                    continue
                i = int(g)
                if i > len(alts) or alts[i - 1] == "*" or alts[i - 1].startswith("<"):
                    continue
                reads = None
                if "," in ad:
                    try:
                        reads = int(ad.split(",")[i])
                    except (ValueError, IndexError):
                        reads = None
                good, verdict = ok, why
                if good and reads is not None and reads < self.rules["min_alt_reads"]:
                    good, verdict = False, "few_alt_reads"
                zyg = "hemi" if len(gt) == 1 else ("hom" if gt.count(g) == 2 else "het")
                out.append({"chrom": chrom, "pos": int(f[1]), "ref": f[3].upper(), "alt": alts[i - 1].upper(),
                            "zygosity": zyg, "quality_ok": good, "quality": verdict, "alt_reads": reads})
        return out


# ------------------------------------------------------------------ labs
def _ref_bounds(text: str) -> Tuple[Optional[float], Optional[float]]:
    t = str(text or "").replace("～", "-").replace("~", "-").replace("—", "-").replace("–", "-").strip()
    m = re.fullmatch(r"\s*([0-9.]+)\s*-\s*([0-9.]+)\s*", t)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.fullmatch(r"\s*[<≤]\s*=?\s*([0-9.]+)\s*", t)
    if m:
        return None, float(m.group(1))
    m = re.fullmatch(r"\s*[>≥]\s*=?\s*([0-9.]+)\s*", t)
    if m:
        return float(m.group(1)), None
    return None, None


def member_labs(st: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """analyte -> the member's confirmed value (as printed) with its printed range and a flag (high/low/in_range)."""
    from .labnames import clean_value, name_keys, usable_rows
    kit = skillkit()
    spec = data("trait_map.json")["analytes"]
    out: Dict[str, Dict[str, Any]] = {}
    rows = usable_rows(st)
    for key, a in spec.items():
        names = {kit.fold_name(n) for n in a["names"]}
        for r in rows:
            variants = name_keys(r["marker"]) | (name_keys(r["maps_to"]) if r.get("maps_to") else set())
            if not (variants & names):
                continue
            try:
                v = kit.parse_number(clean_value(r["value"]))
            except ValueError:
                continue
            lo, hi = _ref_bounds(r.get("ref_range", ""))
            flag = "high" if hi is not None and v > hi else "low" if lo is not None and v < lo else ("in_range" if (lo is not None or hi is not None) else "no_range")
            out[key] = {"analyte": key, "label_zh": a["label_zh"], "marker": r["marker"], "value": v, "unit": r.get("unit", ""),
                        "ref_range": r.get("ref_range", ""), "flag": flag}
            break
    return out


# ------------------------------------------------------------------ helpers
def _clump(loci: List[Dict[str, Any]], kb: int) -> List[Dict[str, Any]]:
    kept: List[Dict[str, Any]] = []
    for l in sorted(loci, key=lambda x: x["p_value"]):
        if all(not (k["chrom"] == l["chrom"] and abs(k["pos"] - l["pos"]) < kb * 1000) for k in kept):
            kept.append(l)
    return kept


def _assembly(st: Dict[str, Any]) -> str:
    for x in st.get("processed", {}).get("variants", []) or []:
        a = (x.get("provenance") or {}).get("assembly")
        if a:
            return "hg19" if "37" in str(a) or "19" in str(a) else "hg38"
    return "hg38"


PALINDROMIC = ({"A", "T"}, {"C", "G"})


def pick_hit(hits: List[Dict[str, Any]], effect: str) -> Tuple[Optional[Dict[str, Any]], str]:
    """The allele the GWAS effect refers to, with its East Asian frequency.
    - effect allele is one alt: that alt (its frequency);
    - effect allele is the reference base: the reference, whose frequency is 1 minus every alt's (a multi-allelic site
      is fine: the member's reference copies are 2 minus all alt copies);
    - an A/T or C/G SNP with a frequency near 50% is skipped: its strand cannot be told from frequency."""
    snv = [h for h in hits if h.get("pos") and len(h["ref_allele"]) == 1 and len(h["alt_allele"]) == 1]
    if not snv:
        return None, "no_snv_record"
    if len({(h["chrom"], h["pos"], h["ref_allele"]) for h in snv}) > 1:
        return None, "rsid_maps_to_several_positions"
    for h in snv:
        if {h["ref_allele"], h["alt_allele"]} in PALINDROMIC and h["af_eas"] is not None and 0.4 <= h["af_eas"] <= 0.6:
            return None, "palindromic_ambiguous"
    as_alt = [h for h in snv if h["alt_allele"] == effect]
    if len(as_alt) == 1:
        return (dict(as_alt[0], p_effect=as_alt[0]["af_eas"]) if as_alt[0]["af_eas"] is not None else None), \
            ("effect_is_alt" if as_alt[0]["af_eas"] is not None else "no_eas_frequency")
    if snv[0]["ref_allele"] == effect:
        if all(h["af_eas"] is None for h in snv):
            return None, "no_eas_frequency"
        p_ref = 1 - sum(h["af_eas"] or 0.0 for h in snv)      # an alt with no gnomAD EAS record is taken as absent there
        return dict(snv[0], alt_allele=None, p_effect=max(0.0, p_ref)), "effect_is_ref"
    return None, "effect_allele_not_in_record"


def clinvar_tier(rec: Dict[str, Any]) -> Optional[str]:
    """'plp': aggregate classification P/LP only, with review criteria (≥1 star).
    'plp_not_unanimous': mentions pathogenic otherwise (conflicting, with risk factor, no criteria). None: neither."""
    cls = str(rec.get("classification", "")).lower().strip()
    rev = str(rec.get("review_status", "")).lower().strip()
    if "pathogenic" not in cls:
        return None
    starred = rev.startswith("criteria provided") or rev.startswith("reviewed by expert panel") or rev == "practice guideline"
    if cls in PLP_OK and starred and "conflicting" not in rev:
        return "plp"
    return "plp_not_unanimous"


class ClinvarIndex:
    def __init__(self, records: List[Dict[str, Any]], assembly: str):
        name = "GRCh38" if assembly == "hg38" else "GRCh37"
        self.snv: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
        self.spdi: Dict[str, Dict[str, Any]] = {r["spdi"]: r for r in records}
        self.indel_pos: Dict[str, List[int]] = {}
        for r in records:
            loc = (r.get("loc") or {}).get(name)
            if not loc:
                continue
            if len(r["deleted"]) == 1 and len(r["inserted"]) == 1:
                self.snv[(_norm_chrom(loc["chr"]), loc["start"], r["deleted"].upper(), r["inserted"].upper())] = r
            else:
                self.indel_pos.setdefault(_norm_chrom(loc["chr"]), []).append(loc["start"])
        for v in self.indel_pos.values():
            v.sort()

    def near_indel(self, chrom: str, pos: int) -> bool:
        ps = self.indel_pos.get(chrom, [])
        i = bisect.bisect_left(ps, pos - INDEL_WINDOW)
        return i < len(ps) and ps[i] <= pos + INDEL_WINDOW


# ------------------------------------------------------------------ explain
def explain(st: Dict[str, Any], ws: Path, absent_as_ref: Optional[str] = None) -> Dict[str, Any]:
    if (st["member"].get("answers") or {}).get("genetic_disclosure") != "yes":
        raise LAError("genetic results need genetic_disclosure=yes (`la.py member <ws> genetic_disclosure=yes|no`)", EXIT_INPUT)
    if absent_as_ref is not None and len(re.sub(r"[\s\W_]", "", absent_as_ref)) < 6:
        raise LAError("--absent-as-ref needs the reason (what says this VCF lists every non-reference site)", EXIT_INPUT)
    from .methods import pick
    variants = pick(st, "variants")
    if not variants:
        raise LAError("no usable VCF for this member (intake, identity, primary)", EXIT_INPUT)
    tm = data("trait_map.json")
    rules, modes = tm["score_rules"], tm.get("gene_modes", {})
    asm = _assembly(st)
    male = st["member"].get("sex") == "male"
    vcf = Vcf(Path(variants["path"]), rules, absent_is_ref=bool(absent_as_ref))
    labs = member_labs(st)
    # 1. retrieve, isolating failures per analyte and per gene
    plan: Dict[str, Dict[str, Any]] = {}
    for key, lab in labs.items():
        a = tm["analytes"][key]
        it: Dict[str, Any] = {"a": a, "lab": lab, "loci": [], "skipped": {}, "error": None, "regions": {}, "unscanned": {}}
        try:
            g = pubdata.gwas_trait(ws, a["efo"], pmax=rules["gwas_pmax"])
            best: Dict[str, Dict[str, Any]] = {}
            for r in g["records"]:
                if r["rsid"] not in best or r["p_value"] < best[r["rsid"]]["p_value"]:
                    best[r["rsid"]] = r
            mv = pubdata.myvariant_rsids(ws, list(best)[:400], assembly=asm)
            cand = []
            for rs, r in best.items():
                h, how = pick_hit(mv.get(rs, []), r["effect_allele"])
                if not h:
                    it["skipped"][how] = it["skipped"].get(how, 0) + 1
                    continue
                cand.append({**r, "chrom": _norm_chrom(h["chrom"]), "pos": h["pos"], "ref": h["ref_allele"], "alt": h["alt_allele"],
                             "p_eff_eas": h["p_effect"], "variant_ref": h["ref"]})
            it["loci"] = _clump(cand, rules["clump_kb"])[:rules["max_loci"]]
        except LAError as e:
            it["error"] = str(e)[:200]
        for gene in a.get("monogenic_genes", []):
            reg = pubdata.gene_region(ws, gene, "GRCh38" if asm == "hg38" else "GRCh37")
            if reg and reg.get("start") and reg.get("end"):
                it["regions"][gene] = {**reg, "start": max(1, int(reg["start"]) - GENE_PAD), "end": int(reg["end"]) + GENE_PAD}
            else:
                it["unscanned"][gene] = "gene region not retrieved (Ensembl)"
        plan[key] = it
    vcf.prefetch({(l["chrom"], l["pos"]) for it in plan.values() for l in it["loci"]},
                 [(r["chrom"], r["start"], r["end"]) for it in plan.values() for r in it["regions"].values()])
    # 2. compute
    results: Dict[str, Any] = {}
    readouts: List[Dict[str, Any]] = []
    for key, it in plan.items():
        a, lab = it["a"], it["lab"]
        res: Dict[str, Any] = {"analyte": key, "label_zh": a["label_zh"], "member_lab": lab, "efo": a["efo"], "skipped_loci": it["skipped"]}
        if it["error"]:
            res.update(not_retrieved=it["error"], loci_tested=0, loci_informative=0, loci_called=0, coverage=0.0, loci=[])
        else:
            mean = var = score = 0.0
            informative = covered = 0
            detail = []
            for l in it["loci"]:
                raising = l["direction"] > 0
                # effect-allele dosage: an alt's own copies, or the reference base's copies (2 minus every alt)
                gt = vcf.alt_dosage(l["chrom"], l["pos"], l["ref"], l["alt"]) if l["alt"] else vcf.ref_dosage(l["chrom"], l["pos"], l["ref"])
                p_risk = l["p_eff_eas"] if raising else 1 - l["p_eff_eas"]
                risk_allele = l["effect_allele"] if raising else f"non-{l['effect_allele']}"
                useful = rules["min_informative_maf"] <= p_risk <= 1 - rules["min_informative_maf"]
                dos = None
                if gt["dosage"] is not None:
                    dos = gt["dosage"] if raising else 2 - gt["dosage"]
                if useful:
                    informative += 1
                    if dos is not None:
                        covered += 1
                        score += dos
                        mean += 2 * p_risk
                        var += 2 * p_risk * (1 - p_risk)
                detail.append({"rsid": l["rsid"], "genes": l["mapped_genes"], "risk_allele": risk_allele, "risk_allele_freq_eas": round(p_risk, 4),
                               "informative_in_eas": useful, "member_risk_alleles": dos, "call_status": gt["status"], "p_value": l["p_value"],
                               "effect_text": l["beta_text"], "gwas_ref": l["ref"], "pubmed_id": l["pubmed_id"]})
            coverage = covered / informative if informative else 0.0
            res.update(loci_tested=len(it["loci"]), loci_informative=informative, loci_called=covered, coverage=round(coverage, 3), loci=detail)
            n_absent = sum(1 for d in detail if d["informative_in_eas"] and d["call_status"] == "not_called")
            if informative < rules["min_informative_loci"]:
                res["score_not_computed"] = f"only {informative} lead loci vary in East Asians (need {rules['min_informative_loci']})"
            elif coverage < rules["min_coverage"]:
                res["score_not_computed"] = (
                    f"{n_absent} of {informative} loci are absent from the VCF; a variant-only VCF cannot tell reference from not "
                    "sequenced — use a gVCF, or record the judgment with --absent-as-ref \"<reason>\""
                    if n_absent > informative / 2 else f"coverage {coverage:.0%} is below {rules['min_coverage']:.0%}")
            else:
                z = (score - mean) / math.sqrt(var) if var > 0 else 0.0
                pct = 100 * 0.5 * (1 + math.erf(z / math.sqrt(2)))
                res.update(score=score, expected=round(mean, 2), z=round(z, 2), pct_eas=round(pct, 1))
                readouts.append({"id": f"gen.{key}.grs_pct", "label_zh": f"{a['label_zh']}遗传倾向（东亚人群百分位）", "value": round(pct, 1),
                                 "unit": "%", "kind": "genetic_score", "method": "native.genotype_phenotype", "analyte": key,
                                 "coverage_pct": round(100 * coverage, 1), "loci": covered, "absent_as_ref": bool(absent_as_ref),
                                 "note_zh": "在东亚人群中有变异的全基因组显著位点上，升高该指标的等位基因个数，按 gnomAD 东亚频率换算成百分位",
                                 "systems": []})
        res["monogenic_scan"] = _scan(ws, vcf, it, modes, asm, male)
        ms = res["monogenic_scan"]
        if a.get("monogenic_genes") and not ms["not_scanned"] and not ms["indels_unresolved"]:
            readouts.append({"id": f"gen.{key}.clinvar_plp", "label_zh": f"{a['label_zh']}相关单基因致病/可能致病变异（显性、纯合、半合或复合杂合）",
                             "value": sum(1 for x in ms["pathogenic_or_likely"] if not x["carrier_only"]), "unit": "个",
                             "kind": "genetic_finding", "method": "native.genotype_phenotype", "analyte": key, "systems": []})
        results[key] = res
    out = {"generated_at": now_iso(), "assembly": asm, "vcf": variants["path"], "absent_as_ref": absent_as_ref,
           "lookup": "tabix" if vcf._tabix else "stream", "analytes": results}
    p = ws / "work" / "insights" / "genotype_phenotype.json"
    write_json(p, out)
    return {"path": str(p), "absent_as_ref": absent_as_ref,
            "analytes": {k: {"flag": v["member_lab"]["flag"], "coverage": v.get("coverage"), "pct_eas": v.get("pct_eas"),
                             "not_computed": v.get("score_not_computed") or v.get("not_retrieved"),
                             "plp": len(v["monogenic_scan"]["pathogenic_or_likely"]),
                             "carrier_only": sum(1 for x in v["monogenic_scan"]["pathogenic_or_likely"] if x["carrier_only"]),
                             "plp_not_unanimous": len(v["monogenic_scan"]["plp_not_unanimous"]),
                             "not_scanned": sorted(v["monogenic_scan"]["not_scanned"])} for k, v in results.items()},
            "readouts": readouts}


def _scan(ws: Path, vcf: Vcf, it: Dict[str, Any], modes: Dict[str, str], asm: str, male: bool) -> Dict[str, Any]:
    """ClinVar scan of the member's own alleles in each monogenic gene of one analyte."""
    plp: List[Dict[str, Any]] = []
    other: List[Dict[str, Any]] = []
    lowq: List[Dict[str, Any]] = []
    unresolved = 0
    unscanned = dict(it["unscanned"])
    for gene, reg in it["regions"].items():
        try:
            cov = vcf.region_coverage(reg["chrom"], reg["start"], reg["end"])
            carried = vcf.carried(reg["chrom"], reg["start"], reg["end"])
            idx = ClinvarIndex(pubdata.clinvar_gene(ws, gene), asm)
        except (LAError, KeyError, ValueError, TypeError) as e:     # a source failure never becomes "nothing found"
            unscanned[gene] = f"{type(e).__name__}: {str(e)[:120]}"
            continue
        if cov < MIN_REGION_COVERAGE:
            unscanned[gene] = (f"the VCF speaks for {cov:.0%} of the gene region (variant-only VCF: absent sites are unknown; "
                               "a gVCF or a recorded --absent-as-ref judgment is needed); carried variants were still checked")
        hits = []
        for c in carried:
            rec = None
            if len(c["ref"]) == 1 and len(c["alt"]) == 1:
                rec = idx.snv.get((c["chrom"], c["pos"], c["ref"], c["alt"]))
            elif idx.near_indel(c["chrom"], c["pos"]):
                try:
                    sp = pubdata.canonical_spdi(ws, asm, c["chrom"], c["pos"], c["ref"], c["alt"])
                except (LAError, KeyError, ValueError, TypeError):
                    unresolved += 1
                    continue
                if not sp["spdi"] or (sp.get("warning") and "reference" in sp["warning"].lower()):
                    unresolved += 1
                    continue
                rec = idx.spdi.get(sp["spdi"])
            if not rec:
                continue
            tier = clinvar_tier(rec)
            if not tier:
                continue
            mode = modes.get(gene, "unknown")
            hits.append({"gene": gene, "variant": f"{c['chrom']}:{c['pos']}{c['ref']}>{c['alt']}", "spdi": rec["spdi"],
                         "title": rec["title"], "classification": rec["classification"], "review_status": rec["review_status"],
                         "conditions": rec["conditions"], "clinvar_ref": rec["ref"], "zygosity": c["zygosity"], "inheritance": mode,
                         "quality": c["quality"], "alt_reads": c["alt_reads"], "_ok": c["quality_ok"], "_tier": tier})
        good = [x for x in hits if x["_ok"] and x["_tier"] == "plp"]
        for x in hits:
            ok, tier = x.pop("_ok"), x.pop("_tier")
            if x["inheritance"] == "AR":
                x["carrier_only"] = x["zygosity"] == "het" and len(good) < 2
                x["possible_compound_het"] = x["zygosity"] == "het" and len(good) >= 2
            elif x["inheritance"] == "XLR":
                x["carrier_only"] = not male and x["zygosity"] == "het"
                x["possible_compound_het"] = False
            else:
                x["carrier_only"], x["possible_compound_het"] = False, False
            if not ok:
                lowq.append(x)
            elif tier == "plp":
                plp.append(x)
            else:
                other.append(x)
    return {"genes": it["a"].get("monogenic_genes", []), "not_scanned": unscanned, "pathogenic_or_likely": plp,
            "plp_not_unanimous": other, "low_quality_not_counted": lowq, "indels_unresolved": unresolved,
            "rule": "ClinVar aggregate classification P/LP with review criteria; SNVs by position+alleles, indels by canonical SPDI; "
                    "calls: FILTER, GQ/DP, alt reads; gene counted as scanned only when the VCF covers its region"}
