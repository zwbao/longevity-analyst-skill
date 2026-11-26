"""Clinical age from a LongevityBench prompt.

Phenotypic age is not used: that formula needs chronological age, which is what
these items ask for. Pairwise, multiclass, and regression each have their own
model, fit on NHANES 2003–2017. 1999 and 2001 are never used for fitting.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional

_DATA = Path(__file__).resolve().parents[2] / "data"
COEF_PATH = _DATA / "clinical_elasticnet.json"
MODEL_PATH = _DATA / "clinical_models.joblib"

# Longer names first. "bone alkaline phosphatase" must precede "alkaline phosphatase".
_LABS = [
    ("c-reactive protein", "crp_mg_dl"),
    ("glycohemoglobin", "hba1c"),
    ("blood urea nitrogen", "bun_mg_dl"),
    ("bone alkaline phosphatase", "bap_ug_l"),
    ("alkaline phosphatase", "alp_u_l"),
    ("total bilirubin", "bili_mg_dl"),
    ("total cholesterol", "chol_mg_dl"),
    ("serum glucose", "glucose_mg_dl"),
    ("uric acid", "ua_mg_dl"),
    ("creatinine", "creat_mg_dl"),
    ("albumin", "albumin_g_dl"),
    ("triglyceride", "tg_mg_dl"),
    ("white blood cell count", "wbc"),
    ("red blood cell count", "rbc"),
    ("mean cell volume", "mcv"),
    ("red cell distribution width", "rdw"),
    ("segmented neutrophils percent", "neut_pct"),
    ("lymphocyte percent", "lymph_pct"),
    ("monocyte percent", "mono_pct"),
    ("eosinophils percent", "eos_pct"),
    ("basophils percent", "baso_pct"),
    ("vitamin b12", "b12_pmol_l"),
    ("vitamin a", "vita_umol_l"),
    ("vitamin e", "vite_umol_l"),
    ("cystatin", "cystatin_mg_l"),
    ("cadmium", "cadmium_nmol_l"),
    ("hdl", "hdl_mg_dl"),
    ("ldl", "ldl_mg_dl"),
]

_FEATURES = [
    "male", "bmi", "albumin_g_dl", "alp_u_l", "bun_mg_dl", "creat_mg_dl",
    "glucose_mg_dl", "bili_mg_dl", "ua_mg_dl", "chol_mg_dl", "hdl_mg_dl",
    "ldl_mg_dl", "tg_mg_dl", "crp_mg_dl", "wbc", "rbc", "mcv", "rdw",
    "neut_pct", "lymph_pct", "mono_pct", "eos_pct", "baso_pct",     "hba1c",
    "vita_umol_l", "vite_umol_l", "b12_pmol_l", "bap_ug_l", "cadmium_nmol_l",
]

# Prompt-derived mask (union of labs seen on LB-0034/0030/0038/0042). When the
# JSON exists it must match _FEATURES — training and scoring use only these.
_MASK_PATHS = [
    Path(__file__).resolve().parents[2] / "data" / "clinical_feature_mask.json",
]


def load_feature_mask() -> Optional[list]:
    for path in _MASK_PATHS:
        if path.is_file():
            data = json.loads(path.read_text())
            return list(data.get("features") or [])
    return None


def present_features(person: Dict[str, float]) -> List[str]:
    """Keys from the model feature list that the prompt actually exposed."""
    return [f for f in _FEATURES if f in person and person[f] is not None]


def load_coef() -> Optional[dict]:
    if not COEF_PATH.is_file():
        return None
    return json.loads(COEF_PATH.read_text())


def _to_canonical(key: str, value: float, unit: str) -> Optional[float]:
    u = unit.lower().replace("µ", "u").replace("μ", "u").replace(" ", "")
    if key == "albumin_g_dl":
        if u in ("g/dl",):
            return value
        if u in ("g/l",):
            return value / 10.0
    elif key == "creat_mg_dl":
        if u in ("mg/dl",):
            return value
        if u in ("umol/l", "µmol/l"):
            return value / 88.4
    elif key in ("glucose_mg_dl",):
        if u in ("mg/dl",):
            return value
        if u in ("mmol/l",):
            return value * 18.0182
    elif key in ("chol_mg_dl", "hdl_mg_dl", "ldl_mg_dl", "tg_mg_dl"):
        if u in ("mg/dl",):
            return value
        if u in ("mmol/l",):
            return value * 38.67
    elif key == "bun_mg_dl":
        if u in ("mg/dl",):
            return value
        if u in ("mmol/l",):
            return value * 2.801
    elif key == "bili_mg_dl":
        if u in ("mg/dl",):
            return value
        if u in ("umol/l",):
            return value / 17.1
    elif key == "ua_mg_dl":
        if u in ("mg/dl",):
            return value
        if u in ("umol/l",):
            return value / 59.48
    elif key == "crp_mg_dl":
        if u in ("mg/dl",):
            return value
        if u in ("mg/l",):
            return value / 10.0
    elif key == "vita_umol_l":
        if u in ("umol/l",):
            return value
        if u in ("ug/dl", "mcg/dl"):
            return value / 28.645
    elif key == "vite_umol_l":
        if u in ("umol/l",):
            return value
        if u in ("ug/dl", "mcg/dl"):
            return value / 43.07
    elif key == "b12_pmol_l":
        if u in ("pmol/l",):
            return value
        if u in ("pg/ml",):
            return value * 0.7378
    elif key == "cadmium_nmol_l":
        if u in ("nmol/l",):
            return value
        if u in ("ug/l", "mcg/l"):
            return value * 8.897
    elif key in ("alp_u_l", "wbc", "rbc", "mcv", "rdw", "neut_pct", "lymph_pct",
                 "mono_pct", "eos_pct", "baso_pct", "hba1c", "bap_ug_l", "cystatin_mg_l"):
        return value
    return None


def _parse_block(block: str) -> Dict[str, float]:
    row: Dict[str, float] = {}
    sex = re.search(r"\b(male|female)\b", block, re.I)
    if sex:
        row["male"] = 1.0 if sex.group(1).lower() == "male" else 0.0
    bmi = re.search(r"BMI\s+([0-9.]+)", block, re.I)
    if bmi:
        row["bmi"] = float(bmi.group(1))
    for name, unit, value in re.findall(r"([^;]+?)\s+\(([^)]+)\)\s+([0-9]*\.?[0-9]+)", block):
        low = name.lower()
        for needle, key in _LABS:
            if needle in low and key not in row:
                converted = _to_canonical(key, float(value), unit)
                if converted is not None and math.isfinite(converted):
                    row[key] = converted
                break
    return row


def parse_people(text: str) -> List[Dict[str, float]]:
    if re.search(r"Participant\s+A\s*:", text):
        parts = re.split(r"Participant\s+[AB]\s*:", text)
        return [_parse_block(p) for p in parts[1:3]]
    return [_parse_block(text)]


def adjust_creatinine(people: List[Dict[str, float]], source: str) -> None:
    """Put 1999-2000 creatinine on the later scale. Selvin 2007: 0.147 + 1.013 × original.

    SEQN 1-9965 are the 1999-2000 wave. 2001-2002 needs no change. `source` is
    the user prompt or a follow-up line that contains `SEQN <n>`.
    """
    seqns = [int(v) for v in re.findall(r"SEQN\s+(\d+)", source or "")]
    for person, seqn in zip(people, seqns):
        if seqn <= 9965 and "creat_mg_dl" in person:
            person["creat_mg_dl"] = 0.147 + 1.013 * person["creat_mg_dl"]


def age_bin(age: float) -> str:
    if age < 30:
        return "18-30"
    if age < 45:
        return "30-45"
    if age < 60:
        return "45-60"
    if age < 75:
        return "60-75"
    return "75+"


def _predict_one(row: Dict[str, float], coef: dict) -> float:
    med = coef["median"]
    mean = coef["mean"]
    scale = coef["scale"]
    xs = []
    for i, name in enumerate(coef["features"]):
        v = row.get(name, med[i])
        s = scale[i] if scale[i] else 1.0
        xs.append((v - mean[i]) / s)
    return coef["intercept"] + sum(c * x for c, x in zip(coef["coef"], xs))


def _bundle():
    cached = getattr(_bundle, "value", None)
    if cached is not None:
        return cached
    if not MODEL_PATH.is_file():
        _bundle.value = False
        return None
    try:
        import joblib
    except ImportError:
        _bundle.value = False
        return None
    _bundle.value = joblib.load(MODEL_PATH)
    return _bundle.value


def _row_matrix(row: Dict[str, float], features: List[str]):
    import numpy as np
    return np.array([[row.get(name, np.nan) for name in features]], dtype=float)


def _model_for(spec: dict, row: Dict[str, float]):
    models = spec["models"]
    if spec.get("mode") == "sex":
        sex = row.get("male")
        key = "male" if sex == 1 else "female" if sex == 0 else "pooled"
        if key in models:
            return models[key]
    return models.get("pooled") or next(iter(models.values()))


def _pairwise_margin_years(bundle: dict) -> float:
    """Years of predicted-age gap needed before we force A/B.

    Uses 2× holdout regression MAE when present (current joblib ≈ 9 → ~18y).
    Below that band the comparator is inside its own typical error.
    """
    sel = (bundle or {}).get("selection") or {}
    mae = (sel.get("regression") or {}).get("score")
    if isinstance(mae, (int, float)) and mae > 0:
        return float(mae) * 2.0
    return 10.0


def _pairwise_from_ages(ages: List[float], margin: float) -> dict:
    gap = abs(float(ages[0]) - float(ages[1]))
    out = {
        "predicted_age": [round(float(a), 1) for a in ages[:2]],
        "age_gap": round(gap, 1),
        "margin_years": round(float(margin), 1),
    }
    if gap < margin:
        out.update(
            {
                "answer": "uncertain",
                "confidence": "low",
                "reason": "predicted age gap below model uncertainty band; do not force A/B",
            }
        )
        return out
    out["answer"] = "A" if ages[0] >= ages[1] else "B"
    out["confidence"] = "moderate"
    return out


def _pairwise_from_ranker(
    people: List[Dict[str, float]], bundle: dict
) -> Optional[dict]:
    """Direction B hybrid pairwise comparator.

    1) Age-regression scores for both people (side channel + far-gap decisions).
    2) If predicted age gap ≥ `hybrid_gap_years`, trust the age models.
    3) Else use the age-delta ranker on feature diffs A−B; abstain when
       |delta| < `delta_abstain` (near-age / bio–chrono conflict band).
    """
    import numpy as np

    ranker = (bundle or {}).get("pairwise_ranker") or {}
    model = ranker.get("model")
    feats = ranker.get("features") or _FEATURES
    if model is None or len(people) < 2:
        return None

    ages: List[float] = []
    spec = bundle.get("pairwise")
    if spec:
        try:
            ages = [
                float(_model_for(spec, p).predict(_row_matrix(p, spec["features"]))[0])
                for p in people[:2]
            ]
        except Exception:
            ages = []
    age_gap = abs(ages[0] - ages[1]) if len(ages) == 2 else None
    hybrid_gap = float(ranker.get("hybrid_gap_years") or 12.0)

    out: dict = {
        "caller": "clinical_pairwise_ranker",
        "kind": ranker.get("kind"),
    }
    if ages:
        out["predicted_age"] = [round(x, 1) for x in ages]
        out["age_gap"] = round(float(age_gap), 1) if age_gap is not None else None

    # Far gap: age models are reliable on LongevityBench (≈0.81 forced).
    if age_gap is not None and age_gap >= hybrid_gap:
        out.update(
            {
                "answer": "A" if ages[0] >= ages[1] else "B",
                "confidence": "moderate",
                "reason": "age-regression gap above hybrid threshold",
                "hybrid_gap_years": hybrid_gap,
            }
        )
        return out

    a = np.array([people[0].get(f, np.nan) for f in feats], dtype=float)
    b = np.array([people[1].get(f, np.nan) for f in feats], dtype=float)
    diff = (a - b).reshape(1, -1)
    kind = (ranker.get("kind") or "").lower()

    if kind == "hgb_age_delta" or ranker.get("delta_abstain") is not None:
        delta = float(model.predict(diff)[0])
        abstain = float(ranker.get("delta_abstain") or 3.5)
        out["age_delta"] = round(delta, 2)
        out["delta_abstain"] = abstain
        out["hybrid_gap_years"] = hybrid_gap
        if abs(delta) < abstain:
            out.update(
                {
                    "answer": "uncertain",
                    "confidence": "low",
                    "reason": (
                        "near-age band: lab-derived age delta inside abstain margin "
                        "(possible bio vs chronological conflict)"
                    ),
                }
            )
            return out
        out["answer"] = "A" if delta >= 0 else "B"
        out["confidence"] = "moderate"
        out["reason"] = "age-delta ranker"
        return out

    # Legacy probability ranker
    proba = float(model.predict_proba(diff)[0, 1])
    margin = float(ranker.get("proba_margin") or 0.08)
    out["proba_A_older"] = round(proba, 4)
    out["proba_margin"] = margin
    if abs(proba - 0.5) < margin:
        out.update(
            {
                "answer": "uncertain",
                "confidence": "low",
                "reason": "pairwise ranker probability inside abstain margin",
            }
        )
        return out
    out["answer"] = "A" if proba >= 0.5 else "B"
    out["confidence"] = "moderate"
    return out


def _from_models(text: str, fmt: str, people: List[Dict[str, float]]) -> Optional[dict]:
    bundle = _bundle()
    if not bundle or not people or not any(people):
        return None
    if fmt == "pairwise":
        if len(people) < 2:
            return None
        ranked = _pairwise_from_ranker(people, bundle)
        if ranked is not None:
            return ranked
        spec = bundle["pairwise"]
        ages = [
            float(_model_for(spec, p).predict(_row_matrix(p, spec["features"]))[0])
            for p in people[:2]
        ]
        return _pairwise_from_ages(ages, _pairwise_margin_years(bundle))
    if fmt == "regression":
        spec = bundle["regression"]
        age = float(_model_for(spec, people[0]).predict(_row_matrix(people[0], spec["features"]))[0])
        return {"answer": round(age, 1)}
    if fmt == "multiclass":
        spec = bundle["multiclass"]
        label = str(_model_for(spec, people[0]).predict(_row_matrix(people[0], spec["features"]))[0])
        opts = re.findall(r"([A-E])\.\s*([0-9]{2}(?:\s*-\s*[0-9]{2}|\+))", text)
        for letter, raw in opts:
            if raw.replace(" ", "") == label:
                return {"answer": letter, "bin": label}
        return {"answer": label, "bin": label}
    return None


def answer(text: str, fmt: str, adjust=None, follow: str = "") -> Optional[dict]:
    """One exam field for a clinical-age item. None when no model is available.

    Creatinine for SEQN 1-9965 is corrected before either model runs. `follow`
    is searched first, then the prompt. `adjust` may still edit the rows after that.
    """
    people = parse_people(text)
    if not people or not any(people):
        return None
    adjust_creatinine(people, follow or text)
    if adjust:
        adjust(people)
    got = _from_models(text, fmt, people)
    if got is not None:
        return got
    coef = load_coef()
    if coef is None:
        return None
    ages = [_predict_one(p, coef) for p in people]
    if fmt == "pairwise":
        if len(ages) < 2:
            return None
        bundle = _bundle() or {}
        return _pairwise_from_ages(ages[:2], _pairwise_margin_years(bundle))
    if fmt == "regression":
        return {"answer": round(ages[0], 1)}
    if fmt == "multiclass":
        label = age_bin(ages[0])
        opts = re.findall(r"([A-E])\.\s*([0-9]{2}(?:\s*-\s*[0-9]{2}|\+))", text)
        for letter, raw in opts:
            if raw.replace(" ", "") == label:
                return {"answer": letter, "bin": label, "predicted_age": round(ages[0], 1)}
        return {"answer": label, "bin": label, "predicted_age": round(ages[0], 1)}
    return None
