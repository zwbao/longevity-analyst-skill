"""Genotype ↔ phenotype: the member's genotype at the loci that published GWAS tie to a lab analyte, a risk-allele
count placed in the East Asian distribution, and a ClinVar scan of the member's own variants in monogenic genes.

The harness gathers and computes; whether an abnormal value is "mainly inherited" is the agent's judgment
(workflows/04c-insights.md). Rules that keep a statement from being falsely reassuring:
- a site the VCF does not cover is 'not called', never reference, unless a gVCF reference block covers it or the
  agent recorded that absent sites are reference for this delivery (`--absent-as-ref "<reason>"`);
- the monogenic scan reports 'not scanned' for a gene whose region or annotation could not be retrieved (never 0),
  and counts only alleles the member actually carries, from calls that pass the quality rules; indels are annotated;
- ClinVar P/LP means every submission says Pathogenic/Likely pathogenic with review criteria provided;
- a heterozygous P/LP variant in a recessive gene is a carrier finding, not a cause.
"""
from __future__ import annotations

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


def _norm_chrom(c: str) -> str:
    m = CHR.match(str(c))
    return m.group(1).upper().replace("MT", "M") if m else str(c)


def _is_block(f: List[str]) -> bool:
    return bool(END.search(f[7])) and f[4] in ("<NON_REF>", "<*>", ".")


# ------------------------------------------------------------------ VCF access
class Vcf:
    """Genotype lookups on one sample: tabix when an index and the tabix binary both exist, otherwise one streaming
    pass that keeps only the positions and regions asked for (a whole-genome VCF is never held in memory)."""

    def __init__(self, path: Path, rules: Dict[str, Any], absent_is_ref: bool = False):
        self.path = Path(path)
        self.rules = rules
        self.absent_is_ref = absent_is_ref
        indexed = self.path.suffix == ".gz" and (Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists())
        self._tabix = indexed and shutil.which("tabix") is not None
        self._rows: Dict[Tuple[str, int], List[List[str]]] = {}
        self._blocks: Dict[str, List[Tuple[int, int, List[str], List[str], str]]] = {}
        self._loaded = False

    def _open(self):
        return gzip.open(self.path, "rt") if self.path.suffix == ".gz" else open(self.path)

    def prefetch(self, positions: Set[Tuple[str, int]], regions: List[Tuple[str, int, int]]) -> None:
        """One streaming pass (no tabix): keep rows at the wanted positions or inside the wanted regions, plus gVCF blocks."""
        if self._tabix or self._loaded:
            return
        reg: Dict[str, List[Tuple[int, int]]] = {}
        for c, s, e in regions:
            reg.setdefault(_norm_chrom(c), []).append((int(s), int(e)))
        with self._open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                c, pos = _norm_chrom(f[0]), int(f[1])
                if _is_block(f):
                    self._blocks.setdefault(c, []).append((pos, int(END.search(f[7]).group(1)), f[8].split(":"), f[9].split(":"), f[6]))
                elif (c, pos) in positions or any(s <= pos <= e for s, e in reg.get(c, ())):
                    self._rows.setdefault((c, pos), []).append(f)
        self._loaded = True

    def _tabix_rows(self, chrom: str, start: int, end: int) -> List[List[str]]:
        for c in (f"chr{chrom}", chrom):
            try:
                r = subprocess.run(["tabix", str(self.path), f"{c}:{start}-{end}"], capture_output=True, text=True, timeout=60)
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

    def alt_dosage(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """Dosage of this exact alt allele. Every record at the position is read (split multi-allelic rows too)."""
        chrom = _norm_chrom(chrom)
        if self._tabix:
            rows = self._tabix_rows(chrom, pos, pos)
            blocks = [(int(f[1]), int(END.search(f[7]).group(1)), f[8].split(":"), f[9].split(":"), f[6]) for f in rows if _is_block(f)]
            rows = [f for f in rows if not _is_block(f)]
        else:
            self._need_load()
            rows, blocks = self._rows.get((chrom, pos), []), self._blocks.get(chrom, [])
        other = None
        for f in rows:
            if int(f[1]) != pos or f[3].upper() != ref.upper():
                continue
            alts = [a.upper() for a in f[4].split(",")]
            fmt, val = f[8].split(":"), f[9].split(":")
            ok, why = self.quality(fmt, val, f[6])
            gt = re.split(r"[/|]", dict(zip(fmt, val)).get("GT", "./."))
            if alt.upper() in alts:
                if not ok:
                    return {"dosage": None, "status": why}
                if "." in gt or len(gt) != 2:
                    return {"dosage": None, "status": "no_call"}
                idx = str(alts.index(alt.upper()) + 1)
                return {"dosage": sum(1 for g in gt if g == idx), "status": why}
            if ok and "." not in gt and len(gt) == 2:
                other = {"dosage": 0, "status": "called_other_alt"}     # this row is another allele; keep looking
        if other:
            return other
        for s, e, fmt, val, filt in blocks:
            if s <= pos <= e:
                ok, why = self.quality(fmt, val, filt)
                return {"dosage": 0 if ok else None, "status": "ref_block" if ok else why}
        if self.absent_is_ref:
            return {"dosage": 0, "status": "assumed_ref"}
        return {"dosage": None, "status": "not_called"}

    def ref_dosage(self, chrom: str, pos: int, ref: str) -> Dict[str, Any]:
        """Copies of the reference base: 2 minus every alt allele called at the position, over all rows (joined or
        split multi-allelic records give the same count)."""
        chrom = _norm_chrom(chrom)
        if self._tabix:
            rows = self._tabix_rows(chrom, pos, pos)
            blocks = [(int(f[1]), int(END.search(f[7]).group(1)), f[8].split(":"), f[9].split(":"), f[6]) for f in rows if _is_block(f)]
            rows = [f for f in rows if not _is_block(f)]
        else:
            self._need_load()
            rows, blocks = self._rows.get((chrom, pos), []), self._blocks.get(chrom, [])
        rows = [f for f in rows if int(f[1]) == pos and f[3].upper() == ref.upper()]
        if rows:
            alt_copies, status = 0, "called"
            for f in rows:
                fmt, val = f[8].split(":"), f[9].split(":")
                ok, why = self.quality(fmt, val, f[6])
                gt = re.split(r"[/|]", dict(zip(fmt, val)).get("GT", "./."))
                if not ok:
                    return {"dosage": None, "status": why}
                if "." in gt or len(gt) != 2:
                    return {"dosage": None, "status": "no_call"}
                alt_copies += sum(1 for g in gt if g != "0")
                status = why if why != "called" else status
            return {"dosage": max(0, 2 - alt_copies), "status": status}
        for s_, e, fmt, val, filt in blocks:
            if s_ <= pos <= e:
                ok, why = self.quality(fmt, val, filt)
                return {"dosage": 2 if ok else None, "status": "ref_block" if ok else why}
        if self.absent_is_ref:
            return {"dosage": 2, "status": "assumed_ref"}
        return {"dosage": None, "status": "not_called"}

    def carried(self, chrom: str, start: int, end: int) -> List[Dict[str, Any]]:
        """Alt alleles the member carries in a region, each with zygosity and a quality verdict."""
        chrom = _norm_chrom(chrom)
        if self._tabix:
            rows = self._tabix_rows(chrom, start, end)
        else:
            self._need_load()
            rows = [f for (c, p), fs in self._rows.items() if c == chrom and start <= p <= end for f in fs]
        out = []
        for f in rows:
            if len(f) < 10 or _is_block(f):
                continue
            fmt, val = f[8].split(":"), f[9].split(":")
            d = dict(zip(fmt, val))
            gt = re.split(r"[/|]", d.get("GT", "./."))
            present = sorted({g for g in gt if g not in (".", "0")})
            if not present:
                continue
            ok, why = self.quality(fmt, val, f[6])
            alts = f[4].split(",")
            ad = d.get("AD", "")
            for g in present:
                i = int(g)
                if i > len(alts) or alts[i - 1] in ("*", "<NON_REF>", "<*>") or alts[i - 1].startswith("<"):
                    continue
                reads = None
                if "," in ad:
                    try:
                        reads = int(ad.split(",")[i])
                    except (ValueError, IndexError):
                        reads = None
                verdict = why
                if ok and reads is not None and reads < self.rules["min_alt_reads"]:
                    ok, verdict = False, "few_alt_reads"
                out.append({"chrom": chrom, "pos": int(f[1]), "ref": f[3].upper(), "alt": alts[i - 1].upper(),
                            "zygosity": "hom" if gt.count(g) == 2 else "het", "quality_ok": ok, "quality": verdict, "alt_reads": reads})
        return out


def hgvs_id(chrom: str, pos: int, ref: str, alt: str) -> Optional[str]:
    """myvariant.info genomic HGVS id for a VCF allele: SNV, deletion, insertion or delins. None for a duplication-
    style insertion that myvariant keys differently (counted as 'not annotated')."""
    ref, alt = ref.upper(), alt.upper()
    c = f"chr{_norm_chrom(chrom)}"
    if not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
        return None
    if len(ref) == 1 and len(alt) == 1:
        return f"{c}:g.{pos}{ref}>{alt}"
    if ref[0] == alt[0]:
        if len(alt) == 1:                                           # deletion after the anchor base
            s, e = pos + 1, pos + len(ref) - 1
            return f"{c}:g.{s}del" if s == e else f"{c}:g.{s}_{e}del"
        if len(ref) == 1:                                           # insertion after the anchor base
            return f"{c}:g.{pos}_{pos + 1}ins{alt[1:]}"
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return f"{c}:g.{pos}delins{alt}" if len(ref) == 1 else f"{c}:g.{pos}_{pos + len(ref) - 1}delins{alt}"


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


def pick_hit(hits: List[Dict[str, Any]], effect: str) -> Tuple[Optional[Dict[str, Any]], str]:
    """The allele the GWAS effect refers to, with its East Asian frequency.
    - effect allele is one alt: that alt (its frequency);
    - effect allele is the reference base: the reference, whose frequency is 1 minus every alt's (a multi-allelic site
      is fine: the member's reference copies are 2 minus all alt copies)."""
    snv = [h for h in hits if h.get("pos") and len(h["ref_allele"]) == 1 and len(h["alt_allele"]) == 1]
    if not snv:
        return None, "no_snv_record"
    if len({(h["chrom"], h["pos"], h["ref_allele"]) for h in snv}) > 1:
        return None, "rsid_maps_to_several_positions"
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


def plp_unanimous(rec: Dict[str, Any]) -> bool:
    """Every submission P/LP, and at least one with review criteria (ClinVar ≥1 star); any conflict disqualifies."""
    sig = {s.lower() for s in rec.get("clinvar_significance", [])}
    rev = [r.lower().strip() for r in rec.get("review_status", [])]
    starred = any(r.startswith("criteria provided") or r.startswith("reviewed by expert panel") or r == "practice guideline" for r in rev)
    return bool(sig) and sig <= PLP_OK and starred and not any("conflicting" in r for r in rev)


# ------------------------------------------------------------------ explain
def explain(st: Dict[str, Any], ws: Path, absent_as_ref: Optional[str] = None) -> Dict[str, Any]:
    if (st["member"].get("answers") or {}).get("genetic_disclosure") != "yes":
        raise LAError("genetic results need genetic_disclosure=yes (`la.py member <ws> genetic_disclosure=yes|no`)", EXIT_INPUT)
    if absent_as_ref is not None and len(absent_as_ref.strip()) < 8:
        raise LAError("--absent-as-ref needs the reason (what says this VCF lists every non-reference site)", EXIT_INPUT)
    from .methods import pick
    variants = pick(st, "variants")
    if not variants:
        raise LAError("no usable VCF for this member (intake, identity, primary)", EXIT_INPUT)
    tm = data("trait_map.json")
    rules, modes = tm["score_rules"], tm.get("gene_modes", {})
    asm = _assembly(st)
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
                it["regions"][gene] = reg
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
        # monogenic scan
        plp: List[Dict[str, Any]] = []
        lowq: List[Dict[str, Any]] = []
        not_annotated = 0
        for gene, reg in it["regions"].items():
            try:
                carried = vcf.carried(reg["chrom"], reg["start"], reg["end"])
            except LAError as e:
                it["unscanned"][gene] = str(e)[:120]
                continue
            ids: Dict[str, Dict[str, Any]] = {}
            for c in carried:
                h = hgvs_id(c["chrom"], c["pos"], c["ref"], c["alt"])
                if h:
                    ids[h] = c
                else:
                    not_annotated += 1
            try:
                ann = pubdata.myvariant_hgvs(ws, list(ids), assembly=asm) if ids else {}
            except LAError as e:
                it["unscanned"][gene] = f"ClinVar annotation not retrieved: {str(e)[:100]}"
                continue
            not_annotated += sum(1 for h, c in ids.items() if h not in ann and (len(c["ref"]) > 1 or len(c["alt"]) > 1))
            gene_hits = []
            for h, x in ann.items():
                c = ids.get(h)
                if not c or not plp_unanimous(x):
                    continue
                gene_hits.append({"gene": gene, "hgvs": h, "significance": x["clinvar_significance"], "review_status": x.get("review_status", []),
                                  "conditions": x["conditions"], "clinvar_ref": x["ref"], "af_eas": x["af_eas"], "zygosity": c["zygosity"],
                                  "inheritance": modes.get(gene, "unknown"), "quality": c["quality"], "alt_reads": c["alt_reads"],
                                  "_ok": c["quality_ok"]})
            good = [x for x in gene_hits if x["_ok"]]
            for x in gene_hits:
                recessive = x["inheritance"] in ("AR", "XLR")
                # a single heterozygous hit in a recessive gene is a carrier finding; two good hits may be compound heterozygous
                x["carrier_only"] = recessive and x["zygosity"] == "het" and len(good) < 2
                x["possible_compound_het"] = recessive and x["zygosity"] == "het" and len(good) >= 2
                ok = x.pop("_ok")
                (plp if ok else lowq).append(x)
        res["monogenic_scan"] = {"genes": a.get("monogenic_genes", []), "not_scanned": it["unscanned"],
                                 "pathogenic_or_likely": plp, "low_quality_not_counted": lowq, "variants_not_annotated": not_annotated,
                                 "rule": "ClinVar: every submission P/LP with review criteria; call: PASS, GQ/DP and alt-read thresholds"}
        if a.get("monogenic_genes") and not it["unscanned"]:
            readouts.append({"id": f"gen.{key}.clinvar_plp", "label_zh": f"{a['label_zh']}相关单基因致病/可能致病变异（显性、纯合或复合杂合）",
                             "value": sum(1 for x in plp if not x["carrier_only"]), "unit": "个", "kind": "genetic_finding",
                             "method": "native.genotype_phenotype", "analyte": key, "systems": []})
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
                             "not_scanned": sorted(v["monogenic_scan"]["not_scanned"])} for k, v in results.items()},
            "readouts": readouts}
