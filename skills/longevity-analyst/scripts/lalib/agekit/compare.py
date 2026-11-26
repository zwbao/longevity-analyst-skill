"""Which of two omics profiles looks older, from published per-feature age directions (a method of this library).

For each protein (Olink NPX) or CpG (methylation beta or M value) that both profiles have, the profile with the higher
value gets a vote towards "older" when the feature rises with age in the reference cohort, towards "younger" when it
falls; each vote is weighted by the strength of that feature's age association (|z|). No reference means are needed,
so values from any study or lab can be compared as long as both profiles were measured the same way.

References (both CC BY, neither is a LongevityBench source):
  proteins  Iona et al., Eur J Epidemiol 2025 (doi 10.1007/s10654-025-01311-z), eTable 5: China Kadoorie Biobank,
            n=2,006, Olink Explore 3072, SD change per year of age (sex, area, fasting time, temperature, plate adjusted)
  CpGs      McCartney et al., Genome Med 2020 (doi 10.1186/s13073-019-0693-z; data doi 10.7488/ds/2709): Generation
            Scotland EPIC age EWAS; kept: direction replicated in the second wave and discovery P < 1e-5

Standard library only.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from typing import Any, Dict, List, Optional, Tuple

from . import store
from .library import Library, LKError

REFS = {
    "proteomics": {"file": "data/reference/olink_age_direction_ckb.csv.gz",
                   "source": "Iona et al., Eur J Epidemiol 2025, doi 10.1007/s10654-025-01311-z (China Kadoorie Biobank, n=2,006, CC BY)",
                   "check": "LongevityBench mini Olink pairwise (public young-vs-old cohort, not UK Biobank; 45 pairs): 39/45 = 87%"},
    "methylation": {"file": "data/reference/dnam_age_direction_gs.csv.gz",
                    "source": "McCartney et al., Genome Med 2020, doi 10.1186/s13073-019-0693-z; data doi 10.7488/ds/2709 (Generation Scotland, CC BY)",
                    "check": "LongevityBench mini DNAm pairwise (GEO cohorts, 45 pairs): 40/45 = 89%"},
}
_CACHE: Dict[Tuple[str, str], Dict[str, Tuple[float, float, float]]] = {}
_PROTEIN_SIG_P = 2.305e-6  # the paper's Bonferroni threshold over principal components
_NUM = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"


def _load(lib: Library, modality: str) -> Dict[str, Tuple[float, float, float]]:
    """feature -> (direction weight: signed |z|, raw slope, p)."""
    key = (str(lib.root), modality)
    if key in _CACHE:
        return _CACHE[key]
    path = lib.root / REFS[modality]["file"]
    if not path.is_file():
        raise LKError(f"reference {REFS[modality]['file']} is missing from this library", "library_broken")
    table: Dict[str, Tuple[float, float]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if modality == "proteomics":
                b, se, p = float(r["sd_per_year"]), float(r["se"]), float(r["p"])
                if p < 0.05 and se > 0:
                    table[r["gene"].upper()] = (b / se, b, p)
            else:
                table[r["cpg"]] = (float(r["z"]), float(r["m_per_year"]), 0.0)
    _CACHE[key] = table
    return table


def parse_pairs(text: str) -> Dict[str, Tuple[float, float]]:
    """Feature -> (value A, value B) from text such as 'ELN (A: 0.31, B: -0.12)' or 'cg16867657 (...) ... A: 0.73, B: 0.54',
    or from two blocks 'Participant A: ... GENE (0.3); ...' / 'Participant B: ...'."""
    pairs: Dict[str, Tuple[float, float]] = {}
    for m in re.finditer(rf"\b(cg\d{{8}})\b.*?\bA:\s*({_NUM}),\s*B:\s*({_NUM})", text, re.S):
        pairs.setdefault(m.group(1), (float(m.group(2)), float(m.group(3))))
    for m in re.finditer(rf"\b([A-Za-z][A-Za-z0-9_.\-]*)\s*\(\s*A:\s*({_NUM}),\s*B:\s*({_NUM})", text):
        if not m.group(1).startswith("cg"):
            pairs.setdefault(m.group(1).upper(), (float(m.group(2)), float(m.group(3))))
    if pairs:
        return pairs
    blocks = re.split(r"(?im)^\s*(?:participant|individual|person|sample|subject)\s+B\s*[:：]", text)
    if len(blocks) == 2:
        a, b = (_single(x) for x in blocks)
        return {k: (a[k], b[k]) for k in a if k in b}
    return {}


def _single(text: str) -> Dict[str, float]:
    out = {m.group(1).upper(): float(m.group(2)) for m in re.finditer(rf"\b([A-Za-z][A-Za-z0-9_.\-]*)\s*\(\s*({_NUM})\s*\)", text)}
    for m in re.finditer(rf"\b(cg\d{{8}})\b[^;\n]*?(?:beta|β|value)?\s*[:=]?\s*({_NUM})\s*(?:[;\n]|$)", text, re.I):
        out.setdefault(m.group(1), float(m.group(2)))
    return out


def compare_age(lib: Library, pairs: Dict[str, Tuple[float, float]], modality: Optional[str] = None) -> Dict[str, Any]:
    if not pairs:
        return {"status": "input_problem", "message_zh": "没有读到两份可以逐项比较的数值（例如 ELN (A: 0.3, B: -0.1) 或两段 Participant A/B）。"}
    if modality is None:
        modality = "methylation" if sum(k.startswith("cg") for k in pairs) > len(pairs) / 2 else "proteomics"
    if modality not in REFS:
        raise LKError("modality must be proteomics or methylation", "bad_input", 3)
    ref = _load(lib, modality)
    present = [(k, a, b) for k, (a, b) in pairs.items() if k in ref and a == a and b == b]
    strict = _PROTEIN_SIG_P if modality == "proteomics" else 1.0
    used = [x for x in present if ref[x[0]][2] < strict]
    low_confidence = False
    if len(used) < 5 and len(present) >= 5:  # too few strongly age-related features: fall back to nominal ones, say so
        used, low_confidence = present, True
    head = {"method": f"age-direction vote ({modality})",
            "label_zh": "本库的方法：按已发表队列里每项随年龄升降的方向投票（不是已发表时钟，只比较谁更“老”，不估年龄）",
            "reference": REFS[modality]["source"], "independent_check": REFS[modality]["check"], **lib.stamp()}
    if len(used) < 5:
        return {**head, "status": "input_problem",
                "message_zh": f"两份数据里能在参考表中找到的只有 {len(used)} 项，至少要 5 项。", "features_given": len(pairs)}
    older_a = sum(abs(ref[k][0]) for k, a, b in used if (a > b) == (ref[k][0] > 0) and a != b)
    older_b = sum(abs(ref[k][0]) for k, a, b in used if (a < b) == (ref[k][0] > 0) and a != b)
    total = older_a + older_b
    share_a = older_a / total if total else 0.5
    top = sorted(used, key=lambda x: -abs(ref[x[0]][0]))[:8]
    return {**head, "status": "ok",
            "outputs": {"looks_older": "A" if share_a > 0.5 else "B" if share_a < 0.5 else "tie",
                        "weighted_share_for_A_older": round(share_a, 3),
                        "features_used": len(used), "features_given": len(pairs), "low_confidence": low_confidence},
            "strongest_features": [{"feature": k, "A": a, "B": b, "rises_with_age": ref[k][0] > 0,
                                    "points_to": ("A" if (a > b) == (ref[k][0] > 0) else "B") if a != b else "tie"}
                                   for k, a, b in top],
            "caveats_zh": (["强年龄相关的项目不到 5 个，改用名义显著（P<0.05）的项目投票，把握较低。"] if low_confidence else []) + ["两份数据必须用同一平台、同样的标准化方法测得，否则差别可能来自批次而不是年龄。",
                           "只判断谁“看起来”更老；得分接近 0.5 时说明区分不开。不给具体年龄，也不是生物年龄。",
                           f"独立核对（参考队列与题源不重叠）：{REFS[modality]['check']}。"]}


def compare_from_payload(lib: Library, payload: Dict[str, Any]) -> Dict[str, Any]:
    if payload.get("profile_a") and payload.get("profile_b"):
        a, b = payload["profile_a"], payload["profile_b"]
        if not isinstance(a, dict) or not isinstance(b, dict):
            raise LKError("profile_a and profile_b must be objects {feature: value}", "bad_input", 3)
        norm = lambda k: k if str(k).startswith("cg") else str(k).upper()  # noqa: E731
        a = {norm(k): float(v) for k, v in a.items()}
        b = {norm(k): float(v) for k, v in b.items()}
        pairs = {k: (a[k], b[k]) for k in a if k in b}
    elif isinstance(store.resolve(payload), str):
        pairs = parse_pairs(store.resolve(payload))
    else:
        raise LKError("give panel_ref, panel_text, or profile_a and profile_b", "bad_input", 3)
    return compare_age(lib, pairs, payload.get("modality"))


def looks_pairwise_omics(text: str) -> bool:
    return len(re.findall(rf"\(\s*A:\s*{_NUM},\s*B:\s*{_NUM}", text)) >= 10 or \
        len(re.findall(rf"\bcg\d{{8}}\b.*?\bA:\s*{_NUM},\s*B:\s*{_NUM}", text[:200000], re.S)) >= 10
