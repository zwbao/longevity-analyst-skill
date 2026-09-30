"""Lab rows: which analyte a row could be (lookup), and whether the agent confirmed it (judgment).

Name matching alone cannot tell 血肌酐 from 尿肌酐, a fasting from a 2-hour glucose, or triglycerides from
thyroglobulin written "TG": that is a judgment about meaning. So the harness only proposes candidates by name
(whole name, the part outside the parentheses, the parenthetical), and a row reaches a method or an organ index only
after the agent confirmed it with ``la.py labs confirm`` ("yes: this row is <candidate>, serum/plasma, fasting where
it matters") or mapped it explicitly with ``la.py labs map``. A confirmation is bound to the row's content; editing
the row voids it.
"""
from __future__ import annotations

import hashlib
import re
from typing import Any, Dict, List, Optional, Set

from .common import data, load_json, longevity_skills_home, skillkit

_INDEX: Optional[Dict[str, Set[str]]] = None
_LABEL: Dict[str, str] = {}
_GROUPS: List[Set[str]] = []


def _index() -> Dict[str, Set[str]]:
    """Folded name -> analyte codes (LOINC, or the input key when a method gives none)."""
    global _INDEX
    if _INDEX is None:
        kit = skillkit()
        idx: Dict[str, Set[str]] = {}
        home = longevity_skills_home()
        for p in sorted((home / "skills").glob("*/skill.json")):
            try:
                man = load_json(p)
            except Exception:
                continue
            for spec in kit.input_specs(man, "measurements"):
                if spec.get("key", "").startswith("cg") and spec["key"][2:].isdigit():
                    continue
                lo = spec.get("loinc") or [spec["key"]]
                keys = {str(k).casefold() for k in (lo if isinstance(lo, list) else [lo])}
                _GROUPS.append(keys)
                for k in keys:
                    _LABEL.setdefault(k, spec.get("label_zh") or spec["key"])
                for n in [spec["key"], spec.get("label_zh", ""), *spec.get("aliases", [])]:
                    f = kit.fold_name(n)
                    if f:
                        idx.setdefault(f, set()).update(keys)
        for analyte, a in data("organs.json")["lab_aliases"].items():
            keys = {str(k).casefold() for k in a.get("loinc", [analyte])}
            for k in keys:
                _LABEL[k] = a.get("confirm_zh") or _LABEL.get(k) or analyte
            for n in a["names"] + a.get("abbrev", []):
                idx.setdefault(kit.fold_name(n), set()).update(keys)
        _INDEX = idx
    return _INDEX


def _content_key(r: Dict[str, Any]) -> str:
    raw = "\x1f".join(str(r.get(k, "")) for k in ("marker", "value", "unit", "ref_range", "page", "source_file", "maps_to"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:12]


def row_key(r: Dict[str, Any], labs: Optional[List[Dict[str, Any]]] = None) -> str:
    """Content hash plus the occurrence number among identical rows (a report may print the same row twice)."""
    k = _content_key(r)
    if labs is None:
        return r.get("_row_key") or k
    n = sum(1 for x in labs[:next(i for i, x in enumerate(labs) if x is r)] if _content_key(x) == k)
    return f"{k}-{n + 1}"


def keyed(st: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    out = {}
    for r in st["labs"]:
        k = row_key(r, st["labs"])
        r["_row_key"] = k
        out[k] = r
    return out


def numeric(r: Dict[str, Any]) -> bool:
    return any(ch.isdigit() for ch in str(r.get("value", "")))


def name_for_kit(name: str) -> str:
    """'白蛋白 ALB' -> '白蛋白(ALB)': the exact string the method matcher receives."""
    name = str(name).replace("\u3000", " ").replace("\xa0", " ").strip()
    if "(" in name or "（" in name:
        return name
    parts = name.split()
    return f"{parts[0]}({' '.join(parts[1:])})" if len(parts) > 1 else name


def name_keys(marker: str) -> Set[str]:
    """Every folded form a lab name can match on: whole, the name_for_kit rewrite, and each space-separated part."""
    kit = skillkit()
    raw = str(marker or "")
    names = {raw, name_for_kit(raw)} | set(raw.replace("\u3000", " ").replace("\xa0", " ").split())
    out: Set[str] = set()
    for n in names:
        out |= set(kit.name_variants(n))
    return out


def analytes(r: Dict[str, Any]) -> Set[str]:
    """Analyte codes a row could be: every variant the method matcher or the organ calculators could use."""
    kit = skillkit()
    idx = _index()
    if r.get("maps_to"):
        return set(idx.get(kit.fold_name(r["maps_to"]), set()))
    raw = str(r["marker"])
    names = {raw, name_for_kit(raw)} | set(raw.replace("\u3000", " ").replace("\xa0", " ").split())
    out: Set[str] = set()
    for n in names:
        for v in kit.name_variants(n):
            out |= idx.get(v, set())
    return out


def candidates_zh(codes: Set[str]) -> List[str]:
    """One label per analyte group; an organ-calculator label (which says serum / fasting) wins over a method label."""
    _index()
    organ = set()
    for a in data("organs.json")["lab_aliases"].values():
        lo = {str(c).casefold() for c in a.get("loinc", [])}
        if lo & codes:
            organ.add(a.get("confirm_zh"))
            for g in _GROUPS:                          # the method's wider group (e.g. any glucose) is the same question
                if g & lo:
                    lo = lo | g
            codes = codes - lo
    return sorted(organ | {_LABEL.get(c, c) for c in codes})


def _answers(st: Dict[str, Any]) -> List[tuple]:
    """(key, row, current answer or None) for every row, keys computed once."""
    out = []
    for k, r in keyed(st).items():
        c = r.get("confirm") or {}
        out.append((k, r, c.get("answer") if c.get("row_key") == k and c.get("answer") in ("yes", "no") else None))
    return out


def usable_rows(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Numeric rows the agent confirmed 'yes' for exactly this content and position (a map is a yes, unless refused)."""
    return [r for _, r, ans in _answers(st) if numeric(r) and (ans == "yes" or (r.get("maps_to") and ans != "no"))]


def refused_rows(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    ok = {id(r) for r in usable_rows(st)}
    return [r for r in st["labs"] if numeric(r) and id(r) not in ok]


def pending(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Every numeric row without a current answer; a mapped row counts as answered."""
    out = []
    for k, r, ans in _answers(st):
        if not numeric(r) or r.get("maps_to") or ans:
            continue
        out.append({"row_key": k, "position": st["labs"].index(r) + 1, "page": r.get("page", ""),
                    "marker": r["marker"], "value": r.get("value"), "unit": r.get("unit", ""),
                    "ref_range": r.get("ref_range", ""), "source_file": r.get("source_file"),
                    "candidates": candidates_zh(analytes(r))})
    return out


FLAG = re.compile(r"^\s*[↑↓▲▼*]+|[↑↓▲▼*!]+\s*$|(?<=\d)\s*(?:HH|LL|H|L)\s*$|\s*[（(]\s*(?:[↑↓▲▼]|高|低|偏高|偏低|HH|LL|H|L)\s*[)）]\s*$|"
                  r"(?<=\d)\s*(?:偏高|偏低|升高|降低|高|低|危急值?)\s*$")


def clean_value(raw: Any) -> str:
    """'8.4↑', '8.4 H', '250*', '6.8（偏高）' -> the number: the lab's out-of-range flag is not part of the value."""
    t = str(raw).strip()
    for _ in range(3):
        t = FLAG.sub("", t).strip()
        t = re.sub(r"(?<=\d)\s*[%％]$", "", t)             # '7.2%' - the unit belongs in the unit column
    return t
