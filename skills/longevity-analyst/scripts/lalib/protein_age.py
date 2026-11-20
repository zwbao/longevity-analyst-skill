"""Plasma proteomic age from a LongevityBench prompt.

Pairwise / multiclass use the published OrganAge Conventional protein clock
(Goeminne fold-0). Regression uses a de novo elastic net on consensus proteins
fit on the training split (eval sample ids excluded).
"""
from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

_SKILL_DATA = Path(__file__).resolve().parents[2] / "data"
# Optional override for local sample weights (e.g. a week_6 olink dir).
_SAMPLE_OLINK = Path(os.environ["LONGEVITYBENCH_OLINK"]) if os.environ.get("LONGEVITYBENCH_OLINK") else None

_CLOCK_CANDIDATES = [_SKILL_DATA / "protein_clock.json"]
_ENET_CANDIDATES = [_SKILL_DATA / "protein_elasticnet.json"]
if _SAMPLE_OLINK is not None:
    _CLOCK_CANDIDATES.append(_SAMPLE_OLINK / "protein_clock.json")
    _ENET_CANDIDATES.append(_SAMPLE_OLINK / "protein_elasticnet.json")

_CACHE: Dict[str, dict] = {}


def _first_existing(paths) -> Optional[Path]:
    for p in paths:
        if p.is_file():
            return p
    return None


def load(which: str = "elasticnet") -> Optional[dict]:
    """Load `clock` (Goeminne) or `elasticnet` weights."""
    key = which
    if key in _CACHE:
        return _CACHE[key]
    path = _first_existing(_CLOCK_CANDIDATES if which == "clock" else _ENET_CANDIDATES)
    # Fall back: older installs only had elasticnet JSON.
    if path is None and which == "clock":
        path = _first_existing(_ENET_CANDIDATES)
    if path is None:
        return None
    _CACHE[key] = json.loads(path.read_text())
    return _CACHE[key]


def parse_npx(text: str) -> Dict[str, float]:
    """Parse `GENE (NPX)` style protein blocks."""
    out: Dict[str, float] = {}
    for gene, value in re.findall(r"\b([A-Z0-9][A-Z0-9\-_]{1,15})\s*\(([+-]?\d+(?:\.\d+)?)\)", text):
        out[gene] = float(value)
    return out


def parse_pairwise(text: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Parse `GENE (A: x, B: y, ...)` differential lines into two people."""
    a: Dict[str, float] = {}
    b: Dict[str, float] = {}
    for gene, av, bv in re.findall(
        r"\b([A-Z0-9][A-Z0-9\-_]{1,15})\s*\(A:\s*([+-]?\d+(?:\.\d+)?),\s*B:\s*([+-]?\d+(?:\.\d+)?)",
        text,
    ):
        a[gene] = float(av)
        b[gene] = float(bv)
    return a, b


def _predict(row: Dict[str, float], model: dict) -> float:
    total = model["intercept"]
    for i, name in enumerate(model["features"]):
        value = row.get(name, model["median"][i])
        if value is None or (isinstance(value, float) and not math.isfinite(value)):
            value = model["median"][i]
        scale = model["scale"][i] or 1.0
        total += model["coef"][i] * ((value - model["mean"][i]) / scale)
    a = float(model.get("calibrate_a", 1.0) or 1.0)
    b = float(model.get("calibrate_b", 0.0) or 0.0)
    return a * total + b


def _age_bin(age: float) -> str:
    decade = int(age // 10) * 10
    return f"{decade}-{decade + 9}"


def _bin_letter(text: str, label: str) -> str:
    found = re.findall(r"([A-E])\.\s*(\d{2}\s*-\s*\d{2})", text)
    for letter, raw in found:
        if raw.replace(" ", "") == label.replace(" ", ""):
            return letter
    return label


def _model_for_format(fmt: str) -> Optional[dict]:
    # Table S6: pairwise/multiclass → protein_clock; regression → elasticnet.
    if fmt in ("pairwise", "multiclass"):
        return load("clock") or load("elasticnet")
    return load("elasticnet") or load("clock")


def answer(text: str, fmt: str) -> Optional[dict]:
    model = _model_for_format(fmt)
    if not model:
        return None
    if fmt == "pairwise":
        a, b = parse_pairwise(text)
        if not a or not b:
            return None
        ages = [_predict(a, model), _predict(b, model)]
        return {"answer": "A" if ages[0] >= ages[1] else "B", "predicted_age": [round(v, 1) for v in ages]}
    row = parse_npx(text)
    if not row:
        return None
    age = _predict(row, model)
    if fmt == "regression":
        return {"answer": round(age, 1)}
    if fmt == "multiclass":
        label = _age_bin(age)
        return {"answer": _bin_letter(text, label), "bin": label, "predicted_age": round(age, 1)}
    return None
