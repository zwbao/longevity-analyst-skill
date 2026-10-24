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


def usable_rows(st: Dict[str, Any], history: bool = False) -> List[Dict[str, Any]]:
    """Numeric rows the agent confirmed 'yes' for exactly this content and position (a map is a yes, unless refused).

    Only rows of this visit (see `visit_date`) unless `history`: methods, indices and estimators read one visit;
    the twin keeps every confirmed row with its own date for `twin series`."""
    rows = [r for _, r, ans in _answers(st) if numeric(r) and (ans == "yes" or (r.get("maps_to") and ans != "no"))]
    if history:
        return rows
    v = visit_date(st)
    return [r for r in rows if not v or collected_at(r, st)[0] in (None, v)]


# ------------------------------------------------------------ sampling dates (collected_at)
DATE_COLUMNS = ("date", "collectedat", "collectiondate", "sampledate", "sampletime", "采样日期", "采样时间", "检验日期",
                "采集日期", "采集时间", "送检日期", "检测日期", "体检日期", "日期")
_DATE = re.compile(r"^\s*(\d{4})\s*(?:[-/.年]\s*(\d{1,2})\s*[-/.月]\s*(\d{1,2})\s*日?|(\d{2})(\d{2}))(?:[\sT].*)?$")


def parse_date(raw: Any) -> Optional[str]:
    """'2024-09-03', '2024/9/3', '2024.09.03 08:10', '20240903', '2024年9月3日' -> ISO date; anything else None."""
    import datetime as _d
    m = _DATE.match(str(raw or ""))
    if not m:
        return None
    y, mo, d = (m.group(1), m.group(2) or m.group(4), m.group(3) or m.group(5))
    try:
        return _d.date(int(y), int(mo), int(d)).isoformat()
    except ValueError:
        return None


def collected_at(r: Dict[str, Any], st: Dict[str, Any]) -> tuple:
    """(date, source) for a lab row: its own date column (row), else the date the agent recorded for its file with a
    quote (file), else the workspace sample date (workspace), else (None, None)."""
    d = parse_date(r.get("date_raw")) if r.get("date_raw") else None
    if d:
        return d, "row"
    f = next((x for x in st.get("files", []) if x["id"] == r.get("source_file")), None) or {}
    if (f.get("collected_at") or {}).get("date"):
        return f["collected_at"]["date"], "file"
    if st["member"].get("sample_date"):
        return st["member"]["sample_date"], "workspace"
    return None, None


def visit_date(st: Dict[str, Any]) -> Optional[str]:
    """The visit this workspace reports on: the sample date, else the latest row date (None: rows carry no dates)."""
    if st["member"].get("sample_date"):
        return st["member"]["sample_date"]
    ds = [collected_at(r, st)[0] for r in st.get("labs", [])]
    ds = [d for d in ds if d]
    return max(ds) if ds else None


def stamp_dates(st: Dict[str, Any]) -> List[str]:
    """Write collected_at / collected_at_source on every row; ask (pending judgment) where a file's date is not clear:
    a date cell that is not a date, rows of one file with and without dates, or an undated file in a workspace whose
    other rows carry dates and that has no sample date. Never defaulted."""
    from .common import add_judgment
    asked = []
    by_file: Dict[str, List[Dict[str, Any]]] = {}
    for r in st.get("labs", []):
        d, src = collected_at(r, st)
        r["collected_at"], r["collected_at_source"] = d, src
        by_file.setdefault(r.get("source_file") or "", []).append(r)
    dated_elsewhere = any(r.get("collected_at_source") in ("row", "file") for r in st.get("labs", []))
    for fid, rows in by_file.items():
        f = next((x for x in st.get("files", []) if x["id"] == fid), None)
        if not f or (f.get("collected_at") or {}).get("date"):
            continue
        bad = [r.get("date_raw") for r in rows if r.get("date_raw") and not parse_date(r["date_raw"])]
        has = [r for r in rows if parse_date(r.get("date_raw"))]
        why = None
        if bad:
            why = f"date cells that are not dates ({', '.join(map(str, bad[:3]))})"
        elif has and len(has) < len(rows):
            why = f"{len(rows) - len(has)} of {len(rows)} rows have no date while the others do"
        elif not has and not st["member"].get("sample_date") and dated_elsewhere:
            why = "no date in the file, other files carry dates, and the workspace has no sample date"
        if why:
            add_judgment(st, "collected_at", fid,
                         f"{f['name']}: when were these samples collected? {why}. Read the report header (采样/检验日期) or ask "
                         "the user, then `la.py labs date <ws> --file " + fid + " --date YYYY-MM-DD --source \"<quote>\" --reason \"...\"`.")
            asked.append(fid)
    return asked


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
                    "collected_at": collected_at(r, st)[0], "candidates": candidates_zh(analytes(r))})
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
