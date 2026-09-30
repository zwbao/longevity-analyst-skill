"""Genotype ↔ phenotype: the member's genotype at the loci that published GWAS tie to a lab analyte, a risk-allele
score placed in the East Asian distribution, and a ClinVar scan of the member's own variants in monogenic genes.

The harness gathers and computes; whether an abnormal value is "mainly inherited" is the agent's judgment
(workflows/04c-insights.md). A site the VCF does not cover is 'not called', never assumed reference, unless a
gVCF reference block covers it with adequate quality.
"""
from __future__ import annotations

import gzip
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import pubdata
from .common import EXIT_INPUT, LAError, data, now_iso, sha256_file, skillkit, write_json

CHR = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|M|MT)$", re.I)


def _norm_chrom(c: str) -> str:
    m = CHR.match(str(c))
    return m.group(1).upper().replace("MT", "M") if m else str(c)


# ------------------------------------------------------------------ VCF access
class Vcf:
    """Genotype lookups on one sample. Uses tabix when an index exists, otherwise one streaming pass."""

    def __init__(self, path: Path, rules: Dict[str, Any]):
        self.path = Path(path)
        self.rules = rules
        self._tabix = self.path.suffix == ".gz" and (Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists())
        self._records: Optional[Dict[Tuple[str, int], List[List[str]]]] = None
        self._blocks: Optional[Dict[str, List[Tuple[int, int, List[str], List[str]]]]] = None
        self.chrom_style = None

    def _open(self):
        return gzip.open(self.path, "rt") if self.path.suffix == ".gz" else open(self.path)

    def _load(self, wanted: Optional[set] = None) -> None:
        recs: Dict[Tuple[str, int], List[List[str]]] = {}
        blocks: Dict[str, List[Tuple[int, int, List[str], List[str]]]] = {}
        with self._open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                c, pos = _norm_chrom(f[0]), int(f[1])
                if self.chrom_style is None:
                    self.chrom_style = "chr" if f[0].lower().startswith("chr") else ""
                end = re.search(r"(?:^|;)END=(\d+)", f[7])
                if end and f[4] in ("<NON_REF>", "<*>", "."):
                    blocks.setdefault(c, []).append((pos, int(end.group(1)), f[8].split(":"), f[9].split(":")))
                    continue
                if wanted is None or (c, pos) in wanted:
                    recs.setdefault((c, pos), []).append(f)
        self._records, self._blocks = recs, blocks

    def _tabix_query(self, chrom: str, pos: int) -> List[List[str]]:
        out = []
        for c in (f"chr{chrom}", chrom):
            try:
                r = subprocess.run(["tabix", str(self.path), f"{c}:{pos}-{pos}"], capture_output=True, text=True, timeout=30)
            except (OSError, subprocess.TimeoutExpired):
                return []
            for line in r.stdout.splitlines():
                out.append(line.split("\t"))
            if out:
                break
        return out

    def _quality_ok(self, fmt: List[str], val: List[str]) -> Tuple[bool, str]:
        d = dict(zip(fmt, val))
        gq, dp = d.get("GQ"), d.get("DP")
        try:
            if gq not in (None, ".") and float(gq) < self.rules["min_gq"]:
                return False, "low_gq"
            if dp not in (None, ".") and float(dp) < self.rules["min_dp"]:
                return False, "low_dp"
        except ValueError:
            return False, "bad_quality_field"
        return True, "ok" if (gq not in (None, ".") and dp not in (None, ".")) else "quality_unverified"

    def alt_dosage(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """{'dosage': 0|1|2|None, 'status': called|ref_block|not_called|filtered|low_gq|...}."""
        chrom = _norm_chrom(chrom)
        rows = None
        if self._tabix:
            rows = self._tabix_query(chrom, pos)
        else:
            if self._records is None:
                self._load()
            rows = self._records.get((chrom, pos), [])
        for f in rows:
            if int(f[1]) != pos:
                continue
            end = re.search(r"(?:^|;)END=(\d+)", f[7])
            if end and f[4] in ("<NON_REF>", "<*>", "."):
                continue
            if f[3].upper() != ref.upper():
                continue
            alts = [a.upper() for a in f[4].split(",")]
            if f[6] not in ("PASS", ".", ""):
                return {"dosage": None, "status": f"filtered:{f[6]}"}
            fmt, val = f[8].split(":"), f[9].split(":")
            ok, why = self._quality_ok(fmt, val)
            if not ok:
                return {"dosage": None, "status": why}
            gt = re.split(r"[/|]", dict(zip(fmt, val)).get("GT", "./."))
            if "." in gt or len(gt) != 2:
                return {"dosage": None, "status": "no_call"}
            if alt.upper() in alts:
                idx = str(alts.index(alt.upper()) + 1)
                return {"dosage": sum(1 for g in gt if g == idx), "status": why if why != "ok" else "called"}
            return {"dosage": 0, "status": "called_other_alt" if any(g != "0" for g in gt) else "called"}
        # gVCF reference block
        blocks = []
        if self._tabix:
            for f in self._tabix_query(chrom, pos):
                end = re.search(r"(?:^|;)END=(\d+)", f[7])
                if end and f[4] in ("<NON_REF>", "<*>", "."):
                    blocks.append((int(f[1]), int(end.group(1)), f[8].split(":"), f[9].split(":")))
        else:
            blocks = self._blocks.get(chrom, [])
        for s, e, fmt, val in blocks:
            if s <= pos <= e:
                ok, why = self._quality_ok(fmt, val)
                return {"dosage": 0 if ok else None, "status": "ref_block" if ok else why}
        return {"dosage": None, "status": "not_called"}

    def records_in(self, chrom: str, start: int, end: int) -> List[List[str]]:
        chrom = _norm_chrom(chrom)
        if self._tabix:
            out = []
            for c in (f"chr{chrom}", chrom):
                try:
                    r = subprocess.run(["tabix", str(self.path), f"{c}:{start}-{end}"], capture_output=True, text=True, timeout=60)
                except (OSError, subprocess.TimeoutExpired):
                    return []
                out = [l.split("\t") for l in r.stdout.splitlines()]
                if out:
                    break
            return out
        if self._records is None:
            self._load()
        return [f for (c, p), fs in self._records.items() if c == chrom and start <= p <= end for f in fs]


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


# ------------------------------------------------------------------ score
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


def explain(st: Dict[str, Any], ws: Path) -> Dict[str, Any]:
    if (st["member"].get("answers") or {}).get("genetic_disclosure") != "yes":
        raise LAError("genetic results need genetic_disclosure=yes (`la.py member <ws> genetic_disclosure=yes|no`)", EXIT_INPUT)
    from .methods import pick
    variants = pick(st, "variants")
    if not variants:
        raise LAError("no usable VCF for this member (intake, identity, primary)", EXIT_INPUT)
    tm = data("trait_map.json")
    rules = tm["score_rules"]
    asm = _assembly(st)
    vcf = Vcf(Path(variants["path"]), rules)
    labs = member_labs(st)
    results: Dict[str, Any] = {}
    readouts: List[Dict[str, Any]] = []
    failed = {}
    for key, lab in labs.items():
        a = tm["analytes"][key]
        try:
            g = pubdata.gwas_trait(ws, a["efo"], pmax=rules["gwas_pmax"])
        except LAError as e:                       # one source down never blocks the other analytes
            failed[key] = str(e)[:200]
            results[key] = {"analyte": key, "label_zh": a["label_zh"], "member_lab": lab, "efo": a["efo"], "not_retrieved": str(e)[:200],
                            "loci_tested": 0, "loci_called": 0, "coverage": 0.0, "loci": [],
                            "monogenic_scan": {"genes": a.get("monogenic_genes", []), "pathogenic_or_likely": [], "not_retrieved": True}}
            continue
        best: Dict[str, Dict[str, Any]] = {}
        for r in g["records"]:
            if r["rsid"] not in best or r["p_value"] < best[r["rsid"]]["p_value"]:
                best[r["rsid"]] = r
        mv = pubdata.myvariant_rsids(ws, list(best)[:400], assembly=asm)
        loci = []
        for rs, r in best.items():
            for h in mv.get(rs, []):
                if not h["pos"] or len(h["ref_allele"]) != 1 or len(h["alt_allele"]) != 1:
                    continue
                if r["effect_allele"] not in (h["ref_allele"], h["alt_allele"]):
                    continue
                if h["af_eas"] is None:
                    continue
                p_alt = h["af_eas"]
                p_eff = p_alt if r["effect_allele"] == h["alt_allele"] else 1 - p_alt
                loci.append({**r, "chrom": _norm_chrom(h["chrom"]), "pos": h["pos"], "ref": h["ref_allele"], "alt": h["alt_allele"],
                             "p_eff_eas": p_eff, "variant_ref": h["ref"]})
                break
        loci = _clump(loci, rules["clump_kb"])[:rules["max_loci"]]
        mean = var = score = 0.0
        covered = 0
        detail = []
        for l in loci:
            gt = vcf.alt_dosage(l["chrom"], l["pos"], l["ref"], l["alt"])
            raising = l["direction"] > 0                  # effect allele raises the trait
            risk_allele = l["effect_allele"] if raising else (l["alt"] if l["effect_allele"] == l["ref"] else l["ref"])
            p_risk = l["p_eff_eas"] if raising else 1 - l["p_eff_eas"]
            dos = None
            if gt["dosage"] is not None:
                dos = gt["dosage"] if risk_allele == l["alt"] else 2 - gt["dosage"]
                covered += 1
                score += dos
                mean += 2 * p_risk
                var += 2 * p_risk * (1 - p_risk)
            detail.append({"rsid": l["rsid"], "genes": l["mapped_genes"], "risk_allele": risk_allele, "risk_allele_freq_eas": round(p_risk, 4),
                           "member_risk_alleles": dos, "call_status": gt["status"], "p_value": l["p_value"], "effect_text": l["beta_text"],
                           "gwas_ref": l["ref"], "pubmed_id": l["pubmed_id"]})
        coverage = covered / len(loci) if loci else 0.0
        res = {"analyte": key, "label_zh": a["label_zh"], "member_lab": lab, "efo": a["efo"], "loci_tested": len(loci),
               "loci_called": covered, "coverage": round(coverage, 3), "loci": detail}
        if loci and coverage >= rules["min_coverage"] and var > 0:
            z = (score - mean) / math.sqrt(var)
            pct = 100 * 0.5 * (1 + math.erf(z / math.sqrt(2)))
            res.update(score=score, expected=round(mean, 2), z=round(z, 2), pct_eas=round(pct, 1))
            readouts.append({"id": f"gen.{key}.grs_pct", "label_zh": f"{a['label_zh']}遗传倾向（东亚人群百分位）", "value": round(pct, 1),
                             "unit": "%", "kind": "genetic_score", "method": "native.genotype_phenotype", "analyte": key,
                             "coverage_pct": round(100 * coverage, 1), "note_zh": f"{covered} 个全基因组显著位点的风险等位基因计数，按 gnomAD 东亚频率换算",
                             "systems": []})
        else:
            res["score_not_computed"] = f"coverage {coverage:.0%} of {len(loci)} lead loci is below {rules['min_coverage']:.0%}"
        # monogenic scan: the member's own variants inside each gene, annotated with ClinVar
        plp = []
        for gene in a.get("monogenic_genes", []):
            reg = pubdata.gene_region(ws, gene, "GRCh38" if asm == "hg38" else "GRCh37")
            if not reg:
                continue
            recs = vcf.records_in(reg["chrom"], reg["start"], reg["end"])
            ids = []
            for f in recs:
                if re.search(r"(?:^|;)END=", f[7]):
                    continue
                gt = f[9].split(":")[0] if len(f) > 9 else "./."
                if not re.search(r"[1-9]", gt):
                    continue
                for alt in f[4].split(","):
                    if len(f[3]) == 1 and len(alt) == 1:
                        ids.append(f"chr{_norm_chrom(f[0])}:g.{f[1]}{f[3]}>{alt}")
            ann = pubdata.myvariant_hgvs(ws, ids, assembly=asm) if ids else {}
            for h, x in ann.items():
                sig = " ".join(x["clinvar_significance"]).lower()
                if "pathogenic" in sig and "conflicting" not in sig and "benign" not in sig:
                    plp.append({"gene": gene, "hgvs": h, "significance": x["clinvar_significance"], "conditions": x["conditions"],
                                "clinvar_ref": x["ref"], "af_eas": x["af_eas"]})
        res["monogenic_scan"] = {"genes": a.get("monogenic_genes", []), "pathogenic_or_likely": plp}
        if a.get("monogenic_genes"):
            readouts.append({"id": f"gen.{key}.clinvar_plp", "label_zh": f"{a['label_zh']}相关单基因 ClinVar 致病/可能致病变异", "value": len(plp),
                             "unit": "个", "kind": "genetic_finding", "method": "native.genotype_phenotype", "analyte": key, "systems": []})
        results[key] = res
    out = {"generated_at": now_iso(), "assembly": asm, "vcf": variants["path"], "analytes": results, "not_retrieved": failed}
    p = ws / "work" / "insights" / "genotype_phenotype.json"
    write_json(p, out)
    return {"path": str(p), "analytes": {k: {"flag": v["member_lab"]["flag"], "coverage": v["coverage"], "pct_eas": v.get("pct_eas"),
                                             "plp": len(v["monogenic_scan"]["pathogenic_or_likely"])} for k, v in results.items()},
            "not_retrieved": failed, "readouts": readouts}
