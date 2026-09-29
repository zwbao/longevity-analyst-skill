"""Methods implemented here because no longevity-skills skill covers them.

Each function applies a published, fixed rule to one person's data and
returns readouts. None of them estimates an age or a risk the source did not
define.
"""
from __future__ import annotations

import gzip
import io
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import LAError, EXIT_INPUT, data

# APOE: rs429358 (C112R) and rs7412 (R158C). Positions on GRCh38 / GRCh37.
APOE_SITES = {
    "GRCh38": {"rs429358": ("19", 44908684, "T", "C"), "rs7412": ("19", 44908822, "C", "T")},
    "GRCh37": {"rs429358": ("19", 45411941, "T", "C"), "rs7412": ("19", 45412079, "C", "T")},
}
# haplotype by (rs429358 allele, rs7412 allele)
APOE_HAP = {("T", "T"): "ε2", ("T", "C"): "ε3", ("C", "C"): "ε4", ("C", "T"): "ε1(罕见)"}


def _open(path: Path):
    return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def _scan_vcf(path: Path, wanted: Dict[Tuple[str, int], str]) -> Dict[str, Any]:
    """All records at the wanted positions, plus gVCF reference blocks (INFO END=) that span them."""
    found: Dict[str, List[Dict[str, Any]]] = {rs: [] for rs in wanted.values()}
    blocks: Dict[str, Dict[str, Any]] = {}
    samples: List[str] = []
    last = max(p for (_, p) in wanted)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0][3:] if parts[0].startswith("chr") else parts[0]
            if chrom != "19":
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos > last:
                break
            fmt = parts[8].split(":")
            sample = parts[9].split(":")
            def fval(key: str) -> Optional[str]:
                if key in fmt and len(sample) > fmt.index(key):
                    v = sample[fmt.index(key)]
                    return None if v in (".", "") else v
                return None
            gt = fval("GT")
            gq = fval("GQ")
            dp = fval("MIN_DP") or fval("DP")
            ad_raw = fval("AD")
            try:
                ad = [int(x) for x in ad_raw.split(",")] if ad_raw else None
            except ValueError:
                ad = None
            rec = {"ref": parts[3], "alt": parts[4].split(","), "gt": gt, "pos": pos, "filter": parts[6], "ad": ad,
                   "gq": _num_or_none(gq), "gq_raw": gq,
                   "dp": _num_or_none(dp)}
            if ("19", pos) in wanted:
                found[wanted[("19", pos)]].append(rec)
            m = re.search(r"(?:^|;)END=(\d+)", parts[7])
            if m and all(a in ("<NON_REF>", "<*>", ".") for a in rec["alt"]):
                stop = int(m.group(1))
                for (c, p), rs in wanted.items():
                    if pos <= p <= stop:
                        blocks[rs] = rec
    return {"samples": samples, "sites": found, "blocks": blocks}


def _alleles(site: Dict[str, Any]) -> Optional[List[str]]:
    gt = site.get("gt")
    if not gt or "." in gt:
        return None
    idx = [int(x) for x in re.split(r"[|/]", gt)]
    al = [site["ref"], *site["alt"]]
    return [al[i] for i in idx]


MIN_GQ, MIN_DP = 20, 10


def _num_or_none(v: Optional[str]) -> Optional[float]:
    try:
        return float(v) if v not in (None, "", ".") else None
    except ValueError:
        return None


def _quality(r: Dict[str, Any]) -> Optional[str]:
    """Why a call is not trustworthy, or None. Missing GQ/DP counts as not trustworthy for reference calls."""
    if r.get("filter") not in ("PASS", ".", None):
        return f"FILTER {r['filter']}"
    if r.get("gq") is not None and r["gq"] < MIN_GQ:
        return f"GQ {r['gq']:.0f} < {MIN_GQ}"
    if r.get("dp") is not None and r["dp"] < MIN_DP:
        return f"depth {r['dp']:.0f} < {MIN_DP}"
    return None


def _call_site(recs: List[Dict[str, Any]], block: Optional[Dict[str, Any]], ref: str, alt: str) -> Tuple[Optional[List[str]], Optional[str]]:
    """Genotype at one site from possibly several records (split multi-allelics) or a gVCF reference block.

    A reference call (block or 0/0 record) needs GQ and depth fields that pass; a zero-coverage block is not a call.
    """
    real = [r for r in recs if not all(a in ("<NON_REF>", "<*>") for a in r["alt"])]
    if not real:
        src = block or (recs[0] if recs else None)
        if src and src.get("gt") and set(re.split(r"[|/]", src["gt"])) == {"0"}:
            if src.get("gq") is None or src.get("dp") is None:
                return None, "reference block without GQ/depth fields; coverage unknown"
            q = _quality(src)
            return (None, f"reference block not trustworthy: {q}") if q else ([ref, ref], None)
        return None, "absent"
    for r in real:
        q = _quality(r)
        if q:
            return None, f"call not trustworthy: {q}"
        if r.get("gt") and len(re.split(r"[|/]", r["gt"])) != 2:
            return None, f"genotype {r['gt']} is not diploid"
        ad = r.get("ad")
        if ad and r.get("gt") and "." not in r["gt"]:
            idx = [int(x) for x in re.split(r"[|/]", r["gt"])]
            if any(i < len(ad) and ad[i] == 0 for i in set(idx)) and sum(ad) >= 10:
                return None, f"allele depths {ad} contradict genotype {r['gt']}"
    calls = []
    for r in real:
        if r["ref"] != ref:
            return None, f"REF {r['ref']} != expected {ref}"
        a = _alleles(r)
        if a is None:
            return None, "no genotype"
        calls.append([x for x in a if x != "<NON_REF>"])
    carrying = [a for a in calls if alt in a]
    other = [a for a in calls if any(x not in (ref, alt) for x in a)]
    if other:
        return None, "a different ALT allele is called at this position"
    if len({tuple(sorted(a)) for a in carrying}) > 1:
        return None, "conflicting records at the same position"
    return (carrying[0] if carrying else [ref, ref]), None


def apoe(vcf: Path, assembly: Optional[str]) -> Dict[str, Any]:
    asm = assembly or "GRCh38"
    sites = APOE_SITES[asm]
    wanted = {(c, p): rs for rs, (c, p, _, _) in sites.items()}
    scan = _scan_vcf(vcf, wanted)
    if len(scan["samples"]) != 1:
        raise LAError(f"{vcf.name} has {len(scan['samples'])} samples; APOE genotyping needs a single-person VCF", EXIT_INPUT)
    calls: Dict[str, Optional[List[str]]] = {}
    notes = []
    why_zh = {"absent": "该位点不在 VCF 里；只含变异位点的 VCF 无法区分纯合参考和没测到，因此不下结论"}
    for rs, (c, p, ref, alt) in sites.items():
        a, why = _call_site(scan["sites"].get(rs, []), scan["blocks"].get(rs), ref, alt)
        if why and why.startswith("REF"):
            raise LAError(f"{rs}: VCF {why} on {asm}; wrong assembly?", EXIT_INPUT)
        if a is None:
            notes.append(f"{rs}（{asm} chr{c}:{p}）无法判定：{why_zh.get(why, why)}")
        calls[rs] = a
    if any(v is None for v in calls.values()):
        return {"status": "not_determined", "readouts": [], "notes": notes, "assembly": asm}
    a1, a2 = calls["rs429358"], calls["rs7412"]
    c = a1.count("C")   # rs429358-C marks ε4 (and rare ε1)
    t = a2.count("T")   # rs7412-T marks ε2 (and rare ε1)
    table = {(0, 0): "ε3/ε3", (0, 1): "ε2/ε3", (0, 2): "ε2/ε2", (1, 0): "ε3/ε4", (2, 0): "ε4/ε4",
             (1, 1): None, (2, 1): "ε1/ε4（罕见）", (1, 2): "ε1/ε2（罕见）", (2, 2): "ε1/ε1（罕见）"}
    geno = table[(c, t)]
    n_e4: Optional[int]
    if geno is None:
        geno = "ε2/ε4（也可能是极罕见的 ε1/ε3：两个位点都杂合，未定相时不能区分，需要定相或 Sanger 确认）"
        n_e4 = None
        notes.append("两个位点都杂合，ε4 个数未定")
    else:
        n_e4 = geno.split("（")[0].count("ε4")
    unverified = any(r.get("gq") is None and r.get("dp") is None for rs in scan["sites"].values() for r in rs)
    if unverified:
        notes.append("VCF 记录没有 GQ/DP 等质量字段（例如参考面板或推断数据），基因型质量无法核对")
    return {
        "status": "ok", "assembly": asm, "notes": notes,
        "readouts": [
            {"id": "native.apoe.genotype", "label_zh": "APOE 基因型", "value": geno, "unit": "", "kind": "computed_quality_unverified" if unverified else "computed",
             "basis": "rs429358 + rs7412 两位点组合（标准 APOE 定义）", "calls": {k: "/".join(v) for k, v in calls.items()}},
            {"id": "native.apoe.e4_count", "label_zh": "APOE ε4 等位基因个数", "value": n_e4, "unit": "个", "kind": "computed_quality_unverified" if unverified else "computed",
             "basis": "由上面的基因型计数"},
        ],
    }


def vcf_summary(vcf: Path) -> Dict[str, Any]:
    n = snv = indel = het = hom = 0
    contigs = set()
    with _open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t", 10)
            if len(p) > 4 and all(a in ("<NON_REF>", "<*>") for a in p[4].split(",")):
                continue                     # gVCF reference block, not a variant
            n += 1
            contigs.add(p[0])
            ref, alts = p[3], [a for a in p[4].split(",") if a not in ("<NON_REF>", "<*>")]
            if all(len(a) == len(ref) == 1 for a in alts):
                snv += 1
            else:
                indel += 1
            if len(p) > 9:
                gt = p[9].split(":")[0]
                al = re.split(r"[|/]", gt)
                if len(al) == 2 and "." not in al:
                    het += al[0] != al[1]
                    hom += al[0] == al[1] and al[0] != "0"
    whole_genome = n > 3_000_000
    return {"status": "ok", "readouts": [
        {"id": "native.vcf.records", "label_zh": "VCF 记录数", "value": n, "unit": "条", "kind": "descriptive"},
        {"id": "native.vcf.contigs", "label_zh": "覆盖的染色体/contig 数", "value": len(contigs), "unit": "个", "kind": "descriptive"},
        {"id": "native.vcf.snv", "label_zh": "SNV 记录数", "value": snv, "unit": "条", "kind": "descriptive"},
        {"id": "native.vcf.indel_or_other", "label_zh": "InDel/其他记录数", "value": indel, "unit": "条", "kind": "descriptive"},
    ], "notes": [] if whole_genome else [f"只有 {n} 条记录、{len(contigs)} 个 contig：这不是全基因组 VCF，PGx/PRS/ACMG 结果会不完整，报告必须说明"]}


# ---------------------------------------------------------------- gut
def read_metaphlan(path: Path) -> Dict[str, float]:
    """Species-level relative abundances (percent) from a MetaPhlAn-style table."""
    sp: Dict[str, float] = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            clade = parts[0]
            last = clade.split("|")[-1]
            if not last.startswith("s__") or "|t__" in clade:
                continue
            val_s = parts[-1] if len(parts) == 2 else (parts[2] if len(parts) >= 3 and re.match(r"^[\d.eE+-]+$", parts[2] or "") else parts[1])
            try:
                v = float(val_s)
            except ValueError:
                continue
            sp[last] = sp.get(last, 0.0) + v
    if not sp:
        raise LAError(f"{path.name}: no species-level (s__) rows found", EXIT_INPUT)
    return sp


def gut_diversity(path: Path) -> Dict[str, Any]:
    sp = read_metaphlan(path)
    total = sum(sp.values())
    p = [v / total for v in sp.values() if v > 0]
    shannon = -sum(x * math.log(x) for x in p)
    simpson = 1 - sum(x * x for x in p)
    return {"status": "ok", "notes": [f"物种层 {len(p)} 个；丰度总和 {total:.2f}%（已按物种层重新归一化）"], "readouts": [
        {"id": "native.gut.richness", "label_zh": "物种丰富度（检出物种数）", "value": len(p), "unit": "个", "kind": "computed", "basis": "MetaPhlAn 物种层检出数；受测序深度影响"},
        {"id": "native.gut.shannon", "label_zh": "Shannon 多样性指数", "value": round(shannon, 3), "unit": "", "kind": "computed", "basis": "−Σ p·ln p（自然对数）"},
        {"id": "native.gut.simpson", "label_zh": "Gini-Simpson 多样性指数", "value": round(simpson, 3), "unit": "", "kind": "computed", "basis": "1 − Σ p²"},
    ]}


def gmhi(path: Path) -> Dict[str, Any]:
    """GMHI exactly as GMHI.R (Gupta et al. 2020)."""
    ref = data("gmhi.json")
    sp = read_metaphlan(path)
    keep = {k: v for k, v in sp.items() if "unclassified" not in k and "virus" not in k.lower()}
    total = sum(keep.values())
    if total <= 0:
        raise LAError("no classified species for GMHI", EXIT_INPUT)
    prof = {k: v / total for k, v in keep.items()}
    prof = {k: (0.0 if v < ref["zero_threshold"] else v) for k, v in prof.items()}
    mh = [prof[s] for s in ref["mh_species"] if prof.get(s, 0) > 0]
    mn = [prof[s] for s in ref["mn_species"] if prof.get(s, 0) > 0]
    alpha = lambda xs: -sum(math.log(x) * x for x in xs if x > 0)
    psi_mh = (len(mh) / ref["mh_prime"]) * alpha(mh)
    psi_mn = (len(mn) / ref["mn_prime"]) * alpha(mn)
    val = math.log10((psi_mh + ref["pseudocount"]) / (psi_mn + ref["pseudocount"]))
    matched = sum(1 for s in ref["mh_species"] + ref["mn_species"] if s in sp)
    return {"status": "ok", "notes": [
        f"GMHI 名单 50 个物种里本样本检出 {len(mh)}/7 个健康富集种、{len(mn)}/43 个健康稀缺种",
        "GMHI 的 50 个物种名来自 MetaPhlAn2；本方法只在 MetaPhlAn2 谱上运行（其它版本在计划阶段被拦下）",
    ], "readouts": [
        {"id": "native.gut.gmhi", "label_zh": "肠道菌群健康指数 GMHI", "value": round(val, 3), "unit": "", "kind": "computed",
         "basis": f"Gupta 2020 Nat Commun GMHI.R 原式；>0 偏向健康人群谱，<0 偏向疾病人群谱（原文切点 0）", "citation": ref["doi"]},
        {"id": "native.gut.gmhi_mh_detected", "label_zh": "GMHI 健康富集种检出数", "value": len(mh), "unit": "/7", "kind": "computed"},
        {"id": "native.gut.gmhi_mn_detected", "label_zh": "GMHI 健康稀缺种检出数", "value": len(mn), "unit": "/43", "kind": "computed"},
    ]}


# --------------------------------------------------------- proteomics
def proteomics_describe(path: Path, platform: str) -> Dict[str, Any]:
    """QC description of the member's derived protein column (protein<TAB>value). No protein age."""
    n = missing = 0
    with _open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            n += 1
            try:
                float(parts[1])
            except ValueError:
                missing += 1
    return {"status": "ok", "notes": [f"平台：{platform}。本期只做质控描述，不计算任何蛋白年龄或器官年龄"], "readouts": [
        {"id": "native.prot.proteins", "label_zh": "蛋白条目数", "value": n, "unit": "个", "kind": "descriptive"},
        {"id": "native.prot.quantified", "label_zh": "有定量值的蛋白数", "value": n - missing, "unit": "个", "kind": "descriptive"},
        {"id": "native.prot.missing_pct", "label_zh": "缺失值比例", "value": round(100 * missing / n, 1) if n else None, "unit": "%", "kind": "descriptive"},
    ]}
