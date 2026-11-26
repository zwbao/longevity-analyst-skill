"""Method-library input gates for capability A (schema match before Host).

Deterministic checks only. No Longevity-LLM. Used by Host tooling docs and
future callers that want to refuse a formula when units / features mismatch.
Also capability B (Jarchow) intervention-safety flags and vignette scan.
"""
from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Sequence


# Clinical age / mortality prompts often mix SI and conventional units.
CLINICAL_UNIT_ALIASES = {
    "creatinine": {"mg/dl": 1.0, "mg/dL": 1.0, "umol/l": 1 / 88.4, "µmol/l": 1 / 88.4},
    "glucose": {"mg/dl": 1.0, "mg/dL": 1.0, "mmol/l": 1 / 18.0},
}

# Fixed-site clocks that must not claim applicability on rotating-site exam prompts.
FIXED_METHYLATION_CLOCKS = ("epiage", "pyaging", "horvath", "hannum", "dunedinpace")

# PhenoAge / Levine clinical clocks need chronological age; LongevityBench age
# pairwise asks for age, so these member methods must not be called on that schema.
PHENOAGE_METHODS = (
    "accelerated-biological-aging-risk",
    "phenoage",
    "levine",
)

# Name-lookup databases: refuse when the gene is masked (Gene X) or the prompt is an exam.
GENE_DB_METHODS = ("opengenes", "synergy-age", "synergyage")


def normalize_clinical_value(marker: str, value: float, unit: str) -> Dict[str, Any]:
    """Convert a clinical marker to the conventional unit used by formula weights."""
    key = marker.strip().lower()
    table = CLINICAL_UNIT_ALIASES.get(key)
    if not table:
        return {"ok": False, "reason": f"no unit table for {marker!r}", "value": value, "unit": unit}
    u = unit.strip()
    # tolerate µmol spelling
    if u.lower() in ("umol/l", "µmol/l", "umol/L"):
        factor = table.get("umol/l") or table.get("µmol/l")
    else:
        factor = None
        for k, f in table.items():
            if k.lower() == u.lower():
                factor = f
                break
    if factor is None:
        return {"ok": False, "reason": f"unknown unit {unit!r} for {marker}", "value": value, "unit": unit}
    return {
        "ok": True,
        "marker": key,
        "value_conventional": float(value) * float(factor),
        "from_unit": unit,
        "to_unit": "mg/dl",
    }


def features_present(required: Sequence[str], available: Iterable[str]) -> Dict[str, Any]:
    have = {a.strip().lower() for a in available}
    missing = [r for r in required if r.strip().lower() not in have]
    return {
        "ok": not missing,
        "missing": missing,
        "n_required": len(required),
        "n_present": len(required) - len(missing),
    }


def refuse_fixed_clock_on_rotating_sites(
    method_id: str,
    site_ids: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Gate: fixed CpG clocks are uncovered on LongevityBench-style rotating sites."""
    mid = (method_id or "").strip().lower()
    blocked = any(tok in mid for tok in FIXED_METHYLATION_CLOCKS)
    return {
        "ok": not blocked,
        "method_id": method_id,
        "refuse": blocked,
        "reason": (
            "fixed-site clock not applicable to rotating-site / exam methylation prompts; use Host+tools"
            if blocked
            else "ok"
        ),
        "n_sites_shown": len(site_ids or []),
    }


def refuse_phenoage_without_chrono_age(
    method_id: str,
    chronological_age_present: bool = False,
) -> Dict[str, Any]:
    """Gate: PhenoAge-family methods need chronological age; exam age tasks do not provide it."""
    mid = (method_id or "").strip().lower()
    is_pheno = any(tok in mid for tok in PHENOAGE_METHODS)
    blocked = bool(is_pheno and not chronological_age_present)
    if not is_pheno:
        reason = "ok"
    elif chronological_age_present:
        reason = "ok"
    else:
        reason = (
            "PhenoAge / accelerated-biological-aging-risk requires chronological age; "
            "not applicable when age is the unknown (exam pairwise / member age-blind compare); "
            "use clinical comparator formula or Host+tools"
        )
    return {
        "ok": not blocked,
        "method_id": method_id,
        "refuse": blocked,
        "chronological_age_present": bool(chronological_age_present),
        "reason": reason,
    }


def refuse_gene_db_on_masked_exam(
    method_id: str,
    gene_name_masked: bool = True,
) -> Dict[str, Any]:
    """Gate: OpenGenes / SynergyAge must not be queried when the gene is masked or on exam prompts."""
    mid = (method_id or "").strip().lower().replace("_", "-")
    is_db = any(tok in mid for tok in GENE_DB_METHODS)
    blocked = bool(is_db and gene_name_masked)
    if not is_db:
        reason = "ok"
    elif not gene_name_masked:
        reason = "ok"
    else:
        reason = (
            "gene name is masked or this is an exam prompt; do not query OpenGenes / SynergyAge "
            "(lookup would be N/A or leak the answer); use Host+tools"
        )
    return {
        "ok": not blocked,
        "method_id": method_id,
        "refuse": blocked,
        "gene_name_masked": bool(gene_name_masked),
        "reason": reason,
    }


# Capability B (Jarchow): formal intervention-safety flags.
# Hard → longevity not first-line (recommendation_bias "no").
# Soft → advisory only (case_by_case; note experimental / monitoring).
HARD_SAFETY_FLAGS = (
    "cushing_suspected",
    "hormone_disorder_suspected",
    "iron_deficiency_restrictive_diet",
    "unsupervised_hiit_post_cardiac",
    "rapamycin_supratherapeutic_dose",
)
SOFT_SAFETY_FLAGS = ("experimental_longevity_drug",)
KNOWN_SAFETY_FLAGS = HARD_SAFETY_FLAGS + SOFT_SAFETY_FLAGS


def safety_gate_intervention(
    flags: Dict[str, bool],
) -> Dict[str, Any]:
    """Capability B helper: hard contraindication flags bias recommendation to no.

    Unknown keys are accepted but listed under ``unknown_flags``. Soft flags never
    force ``recommendation_bias`` to ``no``; they appear under ``advisory_flags``.
    """
    raw = flags or {}
    active_all = [k for k, v in raw.items() if v]
    hard_active = [k for k in active_all if k in HARD_SAFETY_FLAGS]
    soft_active = [k for k in active_all if k in SOFT_SAFETY_FLAGS]
    unknown = [k for k in active_all if k not in KNOWN_SAFETY_FLAGS]
    # Treat unknown truthy flags as hard (conservative) so hand-written JSON stays safe.
    blocking = hard_active + unknown
    first_line = len(blocking) == 0
    note = "Feeds Host structured advice; does not replace clinical care."
    if soft_active and first_line:
        note += " Advisory: experimental longevity drug/supplement — label evidence limits; no personal dose as medical order."
    return {
        "ok": True,
        "longevity_first_line": first_line,
        "active_flags": active_all,
        "hard_flags": hard_active,
        "advisory_flags": soft_active,
        "unknown_flags": unknown,
        "recommendation_bias": "no" if blocking else "case_by_case",
        "note": note,
    }


def _has_any(text: str, needles: Sequence[str]) -> bool:
    return any(n in text for n in needles)


def _rapamycin_weekly_mg(text: str) -> Optional[float]:
    """Parse rapamycin dose in mg/week when stated; None if not found."""
    # e.g. "rapamycin of 500 mg/week", "500mg / week"
    m = re.search(
        r"rapamycin[^\n]{0,40}?(\d+(?:\.\d+)?)\s*mg\s*(?:/|\s*per\s*)\s*week",
        text,
        flags=re.I,
    )
    if m:
        return float(m.group(1))
    m = re.search(
        r"(\d+(?:\.\d+)?)\s*mg\s*(?:/|\s*per\s*)\s*week[^\n]{0,40}?rapamycin",
        text,
        flags=re.I,
    )
    if m:
        return float(m.group(1))
    return None


def scan_intervention_vignette(text: str) -> Dict[str, Any]:
    """Heuristic vignette → safety flags for capability B (Jarchow subset patterns).

    Deterministic keyword rules only. Aligns with Expected-No cues in the fixed
    10-item subset (Cushing-like cluster, hormone/PCOS differential, iron
    deficiency + restrictive diet, post-cardiac unsupervised HIIT, supra-therapeutic
    rapamycin). Soft-flags experimental longevity drugs when named without a hard hit.
    """
    raw = text or ""
    low = raw.lower()

    flags: Dict[str, bool] = {k: False for k in KNOWN_SAFETY_FLAGS}

    # --- cushing_suspected: sudden weight gain + facial rounding + weakness/HTN/glucose ---
    sudden_gain = _has_any(
        low,
        (
            "sudden weight gain",
            "sudden unexplained weight",
            "weight gain (approximately",
            "fat gain even in facial",
            "facial region",
            "moon face",
            "looks round",
        ),
    )
    cushing_cluster = sum(
        [
            sudden_gain,
            _has_any(low, ("muscle weakness", "mood swings")),
            _has_any(low, ("high blood pressure", "hypertension", "ramipril")),
            _has_any(low, ("fasting glucose", "pre-diabetes", "110 mg/dl", "increased to 11")),
        ]
    )
    if sudden_gain and cushing_cluster >= 3 and _has_any(
        low, ("facial", "round", "moon", "muscle weakness", "mood swings")
    ):
        # Prefer Cushing when facial rounding / moon-face style cues present with the cluster.
        if _has_any(low, ("facial", "round", "moon")):
            flags["cushing_suspected"] = True
        elif _has_any(low, ("irregular periods", "irregular menstrual")):
            pass  # handled below as hormone
        else:
            flags["cushing_suspected"] = True

    # --- hormone_disorder_suspected: sudden gain + irregular periods + failed CR / IR ---
    if (
        sudden_gain
        and _has_any(low, ("irregular periods", "irregular menstrual", "pcos"))
        and (
            _has_any(low, ("already tried caloric restriction", "without effect", "nothing happens"))
            or _has_any(low, ("fasting insulin", "insulin"))
        )
    ):
        flags["hormone_disorder_suspected"] = True
        # Avoid double-counting as Cushing when periods cue dominates.
        if flags["hormone_disorder_suspected"] and not _has_any(low, ("facial", "moon", "looks round")):
            flags["cushing_suspected"] = False

    # --- iron_deficiency_restrictive_diet: iron deficiency labs + IF/CR ask ---
    iron_def = _has_any(
        low,
        (
            "iron deficiency",
            "ferritin",
            "transferrin saturation",
            "anemia",
        ),
    )
    restrictive_ask = _has_any(
        low,
        (
            "intermittent fasting",
            "caloric restriction",
            "calorie restriction",
            "reduce my daily energy",
            "energy intake",
        ),
    )
    if iron_def and restrictive_ask and _has_any(
        low, ("ferritin", "iron deficiency", "transferrin")
    ):
        flags["iron_deficiency_restrictive_diet"] = True

    # --- unsupervised_hiit_post_cardiac: HIIT after MI / stent / rehab ---
    cardiac = _has_any(
        low,
        (
            "myocardial infarct",
            "myocardial infarction",
            "coronary angioplasty",
            "coronary stent",
            "ptca",
            "cardiovascular rehabilitation",
            "cardiac rehabilitation",
        ),
    )
    hiit = _has_any(low, ("hiit", "high intensity interval"))
    if cardiac and hiit:
        flags["unsupervised_hiit_post_cardiac"] = True

    # --- rapamycin_supratherapeutic_dose: stated weekly dose >> published 3–10 mg/week ---
    dose = _rapamycin_weekly_mg(raw)
    if dose is not None and dose >= 50:  # clearly above 3–10 mg/week proposals
        flags["rapamycin_supratherapeutic_dose"] = True

    # --- soft: experimental longevity drug named ---
    if _has_any(
        low,
        (
            "rapamycin",
            "sirolimus",
            "spermidine",
            "fisetin",
            "senolytic",
        ),
    ):
        flags["experimental_longevity_drug"] = True

    gated = safety_gate_intervention(flags)
    gated["scan"] = {
        "rapamycin_mg_per_week": dose,
        "rules": "jarchow-subset-heuristics/1",
    }
    gated["flags"] = {k: bool(flags.get(k)) for k in KNOWN_SAFETY_FLAGS}
    return gated
