"""Read a pasted lab panel (any language, free text) into measurement rows, per person.

Structural parsing only: it finds "name (unit) value", "name: value unit", "name value unit" and "name=value unit"
pieces, the age ("82-year-old", "age 55", "55岁") and sex of each person, and splits text that describes several
people ("Participant A:", "Participant B:", "患者甲："). Which method to run stays with the caller; `detect` lists
the methods whose required inputs the panel covers, using the same name matching as `run`.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

from . import aliases
from .library import Library

_NUM = r"[-+]?\d+(?:[.,]\d+)?"
_UNIT = r"(?:[A-Za-zµμ%/^*×·0-9.\-]+(?:\s?(?:cells|million|thousand)[A-Za-zµμ/^0-9]*)?(?:/[A-Za-zµμ0-9.^]+)?|%|‰)"
# "Albumin (g/dL) 4.90"  /  "White blood cell count (1000 cells/uL) 6.10"
_PAREN = re.compile(rf"(?P<item>[^;:\n()=]{{2,80}}?(?:\s*\([A-Za-z][^()]{{0,14}}\))?)\s*\((?P<unit>[^()]{{1,30}})\)\s*[:=]?\s*(?P<value>{_NUM})(?!\d)(?![.,]\d)")
# NHANES codebook style with the unit after the name: "Systolic: Blood pres (1st rdg) mm Hg 136"
_UNIT_AFTER = re.compile(rf"(?P<item>[A-Za-z][^;\n=]{{2,60}}?)\s+(?P<unit>mm ?Hg)\s+(?P<value>{_NUM})(?!\d)(?![.,]\d)", re.I)
# "Creatinine: 0.9 mg/dL", "肌酐=78 umol/L", "BMI 24.2 kg/m2", "肌酐 78 μmol/L"
_PLAIN = re.compile(rf"(?P<item>[^\s;:,，；：=()（）\d][^;:,，；：=()（）\n]{{0,60}}?)\s*[:=：]?\s*(?P<value>{_NUM})\s*(?P<unit>{_UNIT})?(?=\s*(?:[;,，；\n]|$))")
_PERSON = re.compile(r"^\s*(?:participant|patient|person|subject|individual|case|患者|受试者|参与者|对象)\s*([A-Z0-9甲乙丙丁一二三四])\s*[:：]",
                     re.I | re.M)
_AGE = [re.compile(r"(\d{1,3}(?:\.\d+)?)[-\s]*(?:year|yr)s?[-\s]*old", re.I),
        re.compile(r"\bage[d]?\s*[:=]?\s*(\d{1,3}(?:\.\d+)?)\b", re.I),
        re.compile(r"(?:年龄|年纪)\s*[:：=]?\s*(\d{1,3}(?:\.\d+)?)"),
        re.compile(r"(\d{1,3})\s*岁")]
_SEX = [(re.compile(r"\b(?:female|woman|women|girl|mother|mom|wife|daughter|sister)\b", re.I), "female"),
        (re.compile(r"\b(?:male|man|men|boy|father|dad|husband|son|brother)\b", re.I), "male"),
        (re.compile(r"女|我妈|妈妈|母亲|老妈|婆婆|奶奶|外婆|姥姥|妻子|老婆|姐姐|妹妹"), "female"),
        (re.compile(r"男|我爸|爸爸|父亲|老爸|公公|爷爷|外公|姥爷|丈夫|老公|哥哥|弟弟"), "male")]
_SKIP_ITEMS = {"age", "aged", "年龄", "year", "years", "option", "options"}


def _split_people(text: str) -> List[Dict[str, str]]:
    marks = list(_PERSON.finditer(text))
    if len(marks) < 2:
        return [{"label": None, "text": text}]
    people = []
    for i, m in enumerate(marks):
        end = marks[i + 1].start() if i + 1 < len(marks) else len(text)
        people.append({"label": m.group(1), "text": text[m.end():end]})
    return people


def _age(text: str) -> Optional[float]:
    for pattern in _AGE:
        m = pattern.search(text)
        if m:
            value = float(m.group(1))
            if 0 < value < 125:
                return value
    return None


def _sex(text: str) -> Optional[str]:
    found = {sex for pattern, sex in _SEX if pattern.search(text)}
    return found.pop() if len(found) == 1 else None


def _looks_like_unit(text: str) -> bool:
    words = text.split()
    return len(words) <= 2 and not re.search(r"[*=]|\b(?:sec|rdg|reading|average|mean of)\b", text, re.I)


def _clean_item(item: str) -> str:
    item = re.sub(r"^.*?(?:blood panel|biomarkers|labs?|panel|化验|检验|生化|血常规)\s*[:：]\s*", "", item, flags=re.I)
    return item.strip(" \t-–—•*·.;:,，；：")


def parse_rows(text: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    taken = []
    for m in _UNIT_AFTER.finditer(text):
        rows.append({"item": _clean_item(m.group("item")), "value": m.group("value").replace(",", "."), "unit": m.group("unit")})
        taken.append((m.start(), m.end()))
    for m in _PAREN.finditer(text):
        if any(a <= m.start() < b for a, b in taken):
            continue
        item = _clean_item(m.group("item"))
        unit = m.group("unit").strip()
        if not _looks_like_unit(unit):  # "60 sec. pulse (30 sec. pulse * 2) 68": the parenthesis is part of the name
            item, unit = f"{item} ({unit})", ""
        if item and item.casefold() not in _SKIP_ITEMS:
            rows.append({"item": item, "value": m.group("value").replace(",", "."), "unit": unit})
            taken.append((m.start(), m.end()))
    rest = list(text)
    for a, b in taken:
        rest[a:b] = [" "] * (b - a)
    rest_text = re.sub(r",\s*((?:frozen |refrigerated )?(?:serum|plasma|whole blood|blood|fasting))\b", r" \1", "".join(rest), flags=re.I)
    for piece in re.split(r"[;\n；。，]|(?<!\d)\.(?=\s)", rest_text):
        for m in _PLAIN.finditer(piece):
            item = _clean_item(m.group("item"))
            if not item or item.casefold() in _SKIP_ITEMS or re.search(r"(?:year|岁|old|taken for|months?|days?|\bat age\b|\bage at\b)", piece[m.start():m.end()], re.I):
                continue
            if len(item.split()) > 7:
                continue
            rows.append({"item": item, "value": m.group("value").replace(",", "."), "unit": (m.group("unit") or "").strip().rstrip(".")})
    return rows


def parse_panel(text: str) -> List[Dict[str, Any]]:
    """[{label, age, sex, rows}] — one entry per person described in the text."""
    out = []
    for person in _split_people(text):
        out.append({"label": person["label"], "age": _age(person["text"]), "sex": _sex(person["text"]),
                    "rows": parse_rows(person["text"])})
    return out


def detect(lib: Library, text: str) -> Dict[str, Any]:
    """Which runnable methods a pasted panel covers, per person, with the rows the library recognised."""
    from .methods import methods  # local import: methods imports nothing from here
    people = []
    for person in parse_panel(text):
        names = [r["item"] for r in person["rows"]]
        english = {n: aliases.analyte_of(n) for n in names}
        have = names + [a for a in english.values() if a] + [lib_name for a in english.values() if a
                                                               for lib_name in aliases.ANALYTES[a]["lib"][:2]]
        if person["age"] is not None:
            have.append("age")
        if person["sex"]:
            have.append("sex")
        listed = methods(lib, have=have)["methods"] if names else []
        complete = []
        for m in listed:
            if not (m.get("complete") and m["status"] == "runnable"):
                continue
            paper = (lib.skills[m["name"]].get("paper") or {})
            outputs = [o["key"] for o in (lib.manifest(m["name"]).get("outputs") or []) if isinstance(o, dict) and o.get("key")]
            complete.append({"name": m["name"], "label_zh": m.get("label_zh"), "uses": m.get("uses_from_have", []),
                             "outputs": outputs,
                             "paper": " ".join(str(paper[k]) for k in ("journal", "year") if paper.get(k)) or None})
        from .reference import analyte_for  # local import: reference imports this module
        ref_items = sorted({a for a in (analyte_for(r["item"]) for r in person["rows"]) if a})
        people.append({"label": person["label"], "age": person["age"], "sex": person["sex"],
                       "rows": len(person["rows"]), "complete_methods": complete,
                       "age_reference_items": len(ref_items) if person["sex"] and len(ref_items) >= 3 else 0})
    from .compare import looks_pairwise_omics
    return {"people": people, "pairwise_omics": looks_pairwise_omics(text), **lib.stamp()}


def hook_context(result: Dict[str, Any]) -> Optional[str]:
    """One short paragraph for an agent host's prompt hook, or None when the text has no usable panel."""
    lines = []
    for p in result["people"]:
        who = f"Person {p['label']}" if p["label"] else "The person"
        if p.get("age_reference_items"):
            lines.append(f"- {who}: {p['age_reference_items']} routine values fit the age reference model -> MCP tool "
                         f"estimate_age gives the age these values look like, a 90% interval and age-band "
                         f"probabilities (fitted on NHANES 2011-2016, held-out mean error ~10 years; not a published method).")
        if not p["complete_methods"]:
            continue
        methods_text = "; ".join(
            f"{m['name']} [{m['paper'] or 'paper in read_method'}] -> {', '.join(m['outputs'][:8]) or 'see read_method'}"
            for m in p["complete_methods"][:6])
        lines.append(f"- {who} ({p['rows']} values read): {methods_text}.")
    if result.get("pairwise_omics"):
        lines.append("- The prompt compares two proteomic or methylation profiles -> MCP tool compare_age votes over every "
                     "protein/CpG with its direction of change with age in a published cohort (gives which looks older, not an "
                     "age; independent check about 87-89% correct on public benchmark pairs).")
    if not lines:
        return None
    return ("longevity plugin: the prompt contains measurements that methods in the longevity library can compute from "
            "directly.\n" + "\n".join(lines) + "\n"
            "Before answering a question about these people's age, which is older, biological age, mortality or survival, "
            "kidney function or body composition, run the relevant tool (run_method with name=<method>, estimate_age, "
            "compare_age). A fitted model is more reliable than judging the values by eye. Then weigh its result together "
            "with the rest of the profile (history, medicines), cite its source and caveats, and say if you chose not to use it.")
