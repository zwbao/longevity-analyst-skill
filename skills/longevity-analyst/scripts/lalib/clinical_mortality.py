"""Ten-year survival from a LongevityBench clinical prompt.

The caller is a Cox model fit on NHANES 2003–2017 linked to the 2019 public
mortality file. 1999 and 2001 are never used for fitting. Phenotypic age and
the one-year mortality marker are not used.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional

from lalib import clinical_age

_PATH = Path(__file__).resolve().parents[2] / "data" / "clinical_cox.json"
_MONTHS = {"0-5": 60, "5-10": 120, "10-15": 180, "15+": 10_000}


def load() -> Optional[dict]:
    if not _PATH.is_file():
        return None
    cached = getattr(load, "value", None)
    if cached is not None:
        return cached
    load.value = json.loads(_PATH.read_text())
    return load.value


def _people(text: str, follow: str = "") -> List[Dict[str, float]]:
    people = clinical_age.parse_people(text)
    if not people:
        return []
    if re.search(r"Participant\s+A\s*:", text):
        blocks = re.split(r"Participant\s+[AB]\s*:", text)[1:1 + len(people)]
    else:
        blocks = [text]
    for person, block in zip(people, blocks):
        age = re.search(r"(\d+)\s*-?\s*year-old", block, re.I)
        if not age:
            age = re.search(r"(\d+)\s+years?\s+old", block, re.I)
        if age:
            person["age"] = float(age.group(1))
    clinical_age.adjust_creatinine(people, follow or text)
    return people


def _lp(person: Dict[str, float], model: dict) -> float:
    total = 0.0
    for i, name in enumerate(model["features"]):
        value = person.get(name, model["median"][i])
        scale = model["scale"][i] or 1.0
        total += model["coef"][i] * ((value - model["mean"][i]) / scale)
    return max(-20.0, min(20.0, total))


def _survival(model: dict, lp: float, month: float) -> float:
    baseline = 1.0
    for time, surv in model["baseline"]:
        if time <= month:
            baseline = surv
        else:
            break
    return baseline ** math.exp(lp)


def _median_month(model: dict, lp: float) -> float:
    for time, surv in model["baseline"]:
        if surv ** math.exp(lp) <= 0.5:
            return time
    return _MONTHS["15+"]


def _bin(month: float) -> str:
    if month < 60:
        return "0-5"
    if month < 120:
        return "5-10"
    if month < 180:
        return "10-15"
    return "15+"


def _yes_letter(text: str, survived: bool) -> str:
    found = re.findall(r"([A-E])\.\s*(Yes|No|Survived|Deceased|Alive|Dead)\b", text, re.I)
    yes = {"yes", "survived", "alive"}
    for letter, word in found:
        positive = word.lower() in yes
        if positive == survived:
            return letter
    return "A" if survived else "B"


def _class_bin(person: Dict[str, float], spec: dict) -> str:
    scores = []
    xs = []
    for i, name in enumerate(spec["features"]):
        value = person.get(name, spec["median"][i])
        scale = spec["scale"][i] or 1.0
        xs.append((value - spec["mean"][i]) / scale)
    for intercept, coef in zip(spec["intercept"], spec["coef"]):
        scores.append(intercept + sum(c * x for c, x in zip(coef, xs)))
    return spec["classes"][scores.index(max(scores))]


def _bin_letter(text: str, label: str) -> str:
    found = re.findall(
        r"([A-E])\.\s*(10\s*-\s*15|5\s*-\s*10|0\s*-\s*5|15\s*\+)\s*years",
        text,
        re.I,
    )
    for letter, raw in found:
        if raw.replace(" ", "") == label:
            return letter
    return label


def _who_letter(text: str, first_longer: bool) -> str:
    found = re.findall(r"([A-E])\.\s*Participant\s+([AB])\b", text, re.I)
    who = "A" if first_longer else "B"
    for letter, person in found:
        if person.upper() == who:
            return letter
    return who


def answer(text: str, fmt: str, follow: str = "") -> Optional[dict]:
    """One exam field for a clinical mortality item. None when the Cox file is absent."""
    model = load()
    if not model:
        return None
    people = _people(text, follow)
    if not people or not any(people):
        return None
    if fmt == "binary":
        # Lower Cox risk means longer expected survival. Threshold was chosen on
        # NHANES 2003–2009 only; 1999 and 2001 were never used.
        thr = float(model.get("binary_lp_threshold", 0.0))
        survived = _lp(people[0], model) <= thr
        return {"answer": _yes_letter(text, survived), "survived_10y": survived}
    if fmt == "pairwise":
        if len(people) < 2:
            return None
        first_longer = _lp(people[0], model) <= _lp(people[1], model)
        return {"answer": _who_letter(text, first_longer)}
    if fmt == "multiclass":
        spec = model.get("multiclass")
        label = _class_bin(people[0], spec) if spec else _bin(_median_month(model, _lp(people[0], model)))
        return {"answer": _bin_letter(text, label), "bin": label}
    return None
