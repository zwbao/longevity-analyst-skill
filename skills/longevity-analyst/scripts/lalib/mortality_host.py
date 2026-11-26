"""Host-facing mortality evidence: Cox scores are evidence, not the final product answer.

Product path for clinical mortality pairwise/multiclass remains Host. This module
exposes structured Cox comparisons plus history-keyword flags for the Host to weigh.
"""
from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional

from . import clinical_mortality as CM

_HISTORY_HEAVY = (
    "alzheimer", "cancer", "malignant", "chemotherapy", "dialysis", "cirrhosis",
    "heart failure", "copd", "oxygen", "hospice", "metastat", "leukemia", "lymphoma",
)


def _history_flags(text: str) -> List[str]:
    low = text.lower()
    return [k for k in _HISTORY_HEAVY if k in low]


def evidence(text: str, fmt: str, follow: str = "") -> Dict[str, Any]:
    """Return Cox evidence + optional suggested field for Host.

    `fmt`: binary | pairwise | multiclass
    """
    fmt = (fmt or "").strip().lower()
    if fmt not in {"binary", "pairwise", "multiclass"}:
        return {"status": "bad_input", "message": "fmt must be binary|pairwise|multiclass"}
    model = CM.load()
    if not model:
        return {"status": "unavailable", "message": "clinical_cox.json missing"}
    people = CM._people(text, follow)
    if not people or not any(people):
        return {"status": "input_problem", "message": "could not parse clinical panel"}

    flags = _history_flags(text + "\n" + (follow or ""))
    got = CM.answer(text, fmt, follow=follow)
    out: Dict[str, Any] = {
        "kind": "mortality-host",
        "status": "ok",
        "provided_by": "longevity-analyst mortality-evidence",
        "format": fmt,
        "note": (
            "Cox linear predictor / survival curves are Host evidence only. "
            "Weigh time-band wording, medicines and disease history yourself; "
            "do not treat survival_probability alone as the exam answer."
        ),
        "history_flags": flags,
        "cox_answer_field": (got or {}).get("answer"),
        "cox_extra": {k: v for k, v in (got or {}).items() if k != "answer"},
    }

    # Pairwise strength from LP gap
    if fmt == "pairwise" and len(people) >= 2:
        lp0 = CM._lp(people[0], model)
        lp1 = CM._lp(people[1], model)
        # P(A outlives B) ≈ logistic(-(lpA-lpB)) with lower LP = lower hazard
        p_a = 1 / (1 + math.exp(lp0 - lp1))
        gap = abs(lp0 - lp1)
        out["cox_pairwise"] = {
            "lp_A": round(lp0, 4),
            "lp_B": round(lp1, 4),
            "p_A_outlives_B": round(p_a, 3),
            "lp_gap": round(gap, 4),
        }
        if flags:
            out["suggested_letter"] = None
            out["confidence"] = "low"
            out["host_hint"] = (
                "history keywords present — do not trust Cox alone; weigh disease course"
            )
        elif gap < 0.15 or 0.4 <= p_a <= 0.6:
            out["suggested_letter"] = None
            out["confidence"] = "low"
            out["host_hint"] = "emit uncertain — Cox cannot separate the pair"
        else:
            out["suggested_letter"] = got["answer"] if got else ("A" if p_a >= 0.5 else "B")
            out["confidence"] = "medium" if gap < 0.4 else "high"
            out["host_hint"] = f"prefer {out['suggested_letter']} from Cox risk ordering"
        return out

    if fmt == "binary":
        if flags:
            out["suggested_letter"] = None
            out["confidence"] = "low"
            out["host_hint"] = "history keywords present — Host must weigh beyond Cox 10y threshold"
        else:
            out["suggested_letter"] = (got or {}).get("answer")
            out["confidence"] = "medium"
            out["host_hint"] = "Cox 10y threshold suggestion; still check age and comorbidities in prose"
        return out

    # multiclass: Cox band is weak per age_evidence caveat — always low confidence suggestion
    out["suggested_letter"] = (got or {}).get("answer")
    out["suggested_bin"] = (got or {}).get("bin")
    out["confidence"] = "low"
    out["host_hint"] = (
        "Cox does not map well onto survival-time bins; treat band as weak prior only"
    )
    return out
