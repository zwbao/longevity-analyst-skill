"""Validate Host LLM outputs for capability A (LongevityBench) and B (Jarchow).

These checks are format gates only. They do not score correctness against labels.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

JARCHOW_KEYS = (
    "comprehensiveness",
    "correctness",
    "usefulness",
    "explainability",
    "safety",
    "recommendation",
    "rationale",
)

_LETTER = re.compile(r"^[A-Da-d]$")
_INT = re.compile(r"^-?\d+$")
_FLOAT = re.compile(r"^-?\d+(\.\d+)?$")
# Common decade / bin labels seen on LongevityBench multiclass prompts.
_BIN = re.compile(
    r"^(0-5|5-10|10-15|15-20|20-25|25-30|30-35|35-40|40-45|45-50|"
    r"50-55|55-60|60-65|65-70|70-75|75-80|80-85|85-90|90\+|90-\+|90\+|"
    r"[A-Da-d]|decade_\d+|age_bin_\d+)$"
)


# Hard cap so Host prose cannot sneak through as a "short" token.
_MAX_FIELD_CHARS = {
    "pairwise": 24,
    "binary": 24,
    "multiclass": 32,
    "regression": 32,
}
_PROSE_HINT = re.compile(
    r"\b(because|therefore|however|based on|in conclusion|the answer is)\b",
    re.I,
)


def check_longevitybench(text: str, fmt: str) -> Dict[str, Any]:
    """Return {ok, format, value, errors[]} for a single exam field."""
    raw = (text or "").strip()
    errors: List[str] = []
    if not raw:
        return {"ok": False, "kind": "longevitybench", "format": fmt, "value": None,
                "errors": ["empty answer"]}
    if "\n" in raw.strip():
        # Allow a single trailing newline; reject multi-line prose.
        lines = [ln for ln in raw.splitlines() if ln.strip()]
        if len(lines) > 1:
            return {"ok": False, "kind": "longevitybench", "format": fmt, "value": None,
                    "errors": ["answer must be a single field, not multi-line prose"]}
        raw = lines[0].strip()

    fmt = (fmt or "").strip().lower()
    value: Any = raw
    max_chars = _MAX_FIELD_CHARS.get(fmt, 32)
    if len(raw) > max_chars:
        return {
            "ok": False,
            "kind": "longevitybench",
            "format": fmt,
            "value": None,
            "errors": [
                f"answer longer than {max_chars} chars looks like prose, not a single field "
                f"(got {len(raw)} chars)"
            ],
        }
    if _PROSE_HINT.search(raw):
        return {
            "ok": False,
            "kind": "longevitybench",
            "format": fmt,
            "value": None,
            "errors": ["answer contains explanatory prose; emit only the exam field"],
        }

    if fmt in ("pairwise", "binary"):
        # Letter choice, or yes/no / increase/decrease style short tokens.
        # `uncertain` / `cannot_tell` = honest abstain (product short-QA; not a forced guess).
        if _LETTER.match(raw):
            value = raw.upper()
        elif raw.lower() in (
            "yes",
            "no",
            "increase",
            "decrease",
            "up",
            "down",
            "a",
            "b",
            "uncertain",
            "cannot_tell",
        ):
            value = raw.lower()
        elif len(raw) <= 4 and raw.isalpha():
            value = raw.upper() if len(raw) == 1 else raw.lower()
        else:
            errors.append(
                f"pairwise/binary expects a single letter, short token, or uncertain, got {raw!r}"
            )
    elif fmt == "multiclass":
        if _LETTER.match(raw):
            value = raw.upper()
        elif _BIN.match(raw.replace(" ", "")):
            value = raw.replace(" ", "")
        else:
            errors.append(
                f"multiclass expects a letter or age/survival bin label, got {raw!r}"
            )
    elif fmt == "regression":
        cleaned = raw.replace(",", "").rstrip("%")
        if _FLOAT.match(cleaned) or _INT.match(cleaned):
            value = float(cleaned)
        else:
            errors.append(f"regression expects a single number, got {raw!r}")
    else:
        errors.append(
            f"unknown longevitybench format {fmt!r}; "
            "use pairwise|binary|multiclass|regression"
        )

    return {
        "ok": not errors,
        "kind": "longevitybench",
        "format": fmt,
        "value": value if not errors else None,
        "errors": errors,
    }


def check_jarchow(obj: Union[Dict[str, Any], str, Path]) -> Dict[str, Any]:
    """Return {ok, missing[], errors[]} for a Jarchow five-axis JSON object."""
    errors: List[str] = []
    if isinstance(obj, (str, Path)):
        path = Path(obj)
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as e:
            return {"ok": False, "kind": "jarchow", "missing": list(JARCHOW_KEYS),
                    "errors": [f"cannot read JSON: {e}"]}
    else:
        data = obj

    if not isinstance(data, dict):
        return {"ok": False, "kind": "jarchow", "missing": list(JARCHOW_KEYS),
                "errors": ["payload must be a JSON object"]}

    missing = [k for k in JARCHOW_KEYS if k not in data]
    if missing:
        errors.append(f"missing keys: {', '.join(missing)}")

    rec = data.get("recommendation")
    if "recommendation" in data and rec not in (True, False, "yes", "no", "Yes", "No", 0, 1):
        errors.append("recommendation must be binary (yes/no or true/false)")

    for axis in ("comprehensiveness", "correctness", "usefulness", "explainability", "safety"):
        if axis in data and data[axis] is None:
            errors.append(f"{axis} must not be null")

    if "rationale" in data and not str(data.get("rationale") or "").strip():
        errors.append("rationale must be non-empty text")

    return {
        "ok": not errors,
        "kind": "jarchow",
        "missing": missing,
        "errors": errors,
        "keys_present": [k for k in JARCHOW_KEYS if k in data],
    }
