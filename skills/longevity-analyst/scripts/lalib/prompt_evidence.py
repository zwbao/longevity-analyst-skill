"""Host tools that extract structured evidence from exam prompt text only.

No OpenGenes / SynergyAge database lookups. Used when method-gate refuses gene-db
on masked names, or when SynergyAge effects are already printed in the prompt.
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

KINDS = ("opengenes-masked", "synergyage-pairwise")


# ------------------------------------------------------------------ OpenGenes v2
# Section-aware cues + pathway→direction priors (no gene-symbol reverse lookup).

_PATHWAY_INC = (
    "inflam", "cytokine", "chemokine", "interleukin", "immune response",
    "nf-kappa", "nfkb", "toll-like", "complement", "senescence",
    "extracellular matrix", "collagen", "angiogenesis",
)
_PATHWAY_DEC = (
    "mitochondri", "oxidative phosphoryl", "electron transport", "respiration",
    "dna repair", "base excision", "nucleotide excision", "homologous recombin",
    "cytoskeletal motor", "microtubule motor", "axon transport", "synaptic vesicle",
    "muscle contraction", "sarcomere",
)

_DEC_CUES = (
    "decreas", "downregulat", "down-regulat", "reduced expression", "reduction",
    "apoptosis", "dna repair", "mitochond", "muscle atrophy", "cytoskeletal motor",
    "microtubule", "motor activity", "axonal transport",
)
_INC_CUES = (
    "increas", "upregulat", "up-regulat", "elevat",
    "inflammatory", "inflammation", "cytokine", "chemokine", "senescence",
    "immune activation",
)


def _section(text: str, title: str) -> str:
    m = re.search(
        rf"{re.escape(title)}\s*:\s*(.*?)(?=\n[A-Z][^\n]{{0,40}}:\s|\nOptions:|\Z)",
        text,
        re.I | re.S,
    )
    return m.group(1).strip() if m else ""


def _option_letters_opengenes(text: str) -> Tuple[str, str]:
    if re.search(r"A\.\s*increases", text, re.I) and re.search(r"B\.\s*decreases", text, re.I):
        return "A", "B"
    if re.search(r"A\.\s*decreases", text, re.I) and re.search(r"B\.\s*increases", text, re.I):
        return "B", "A"
    return "A", "B"


def _pathway_scores(blob: str) -> Tuple[float, float, List[str], List[str]]:
    low = blob.lower()
    inc_hits = [k for k in _PATHWAY_INC if k in low]
    dec_hits = [k for k in _PATHWAY_DEC if k in low]
    return float(len(inc_hits)), float(len(dec_hits)), inc_hits, dec_hits


def opengenes_masked(text: str) -> Dict[str, Any]:
    """v2: segment scores + pathway priors; multi-segment agreement required."""
    animal = _section(text, "Animal model evidence")
    partners = _section(text, "Interaction partners")
    funcs = _section(text, "Functional annotations")
    tissue = _section(text, "Tissue expression data")
    # Remainder often has the age question
    rest = text

    has_animal = bool(animal) and not re.search(
        r"no model organism data available", animal, re.I
    )
    explicit_dec = bool(
        re.search(r"(decreased|decreases|reduction|down).{0,60}(with age|between ages)", rest, re.I)
    )
    explicit_inc = bool(
        re.search(r"(increased|increases|elevation|up).{0,60}(with age|between ages)", rest, re.I)
    )

    segments: Dict[str, Dict[str, Any]] = {}
    score_inc = score_dec = 0.0

    def add_seg(name: str, blob: str, weight: float = 1.0) -> None:
        nonlocal score_inc, score_dec
        if not blob:
            segments[name] = {"empty": True}
            return
        low = blob.lower()
        inc = [c for c in _INC_CUES if c in low]
        dec = [c for c in _DEC_CUES if c in low]
        p_inc, p_dec, ph_inc, ph_dec = _pathway_scores(blob)
        si = weight * (len(inc) + 1.5 * p_inc)
        sd = weight * (len(dec) + 1.5 * p_dec)
        score_inc += si
        score_dec += sd
        segments[name] = {
            "increase_cues": inc,
            "decrease_cues": dec,
            "pathway_increase": ph_inc,
            "pathway_decrease": ph_dec,
            "score_increase": round(si, 2),
            "score_decrease": round(sd, 2),
        }

    add_seg("animal_model", animal, 1.5)
    add_seg("functional_annotations", funcs, 2.0)
    add_seg("interaction_partners", partners, 0.8)
    add_seg("tissue_expression", tissue, 0.5)

    if has_animal and re.search(
        r"(knockout|knock-out|loss.of.function|null|deficient|knockdown)", animal, re.I
    ):
        score_dec += 1.5
        segments["animal_model"]["lof_bonus"] = 1.5

    if explicit_dec:
        score_dec += 5.0
    if explicit_inc:
        score_inc += 5.0

    # Agreement: count segments that lean one way
    lean_inc = lean_dec = 0
    for seg in segments.values():
        if seg.get("empty"):
            continue
        si, sd = seg.get("score_increase", 0), seg.get("score_decrease", 0)
        if si > sd + 0.5:
            lean_inc += 1
        elif sd > si + 0.5:
            lean_dec += 1

    inc_letter, dec_letter = _option_letters_opengenes(text)
    gap = abs(score_dec - score_inc)
    suggested = None
    confidence = "low"
    agreeing = 0
    if explicit_dec and not explicit_inc:
        suggested, confidence, agreeing = dec_letter, "high", max(lean_dec, 1)
    elif explicit_inc and not explicit_dec:
        suggested, confidence, agreeing = inc_letter, "high", max(lean_inc, 1)
    elif score_dec > score_inc and gap >= 2.0 and lean_dec >= 2:
        suggested = dec_letter
        confidence = "high" if gap >= 4 and lean_dec >= 2 else "medium"
        agreeing = lean_dec
    elif score_inc > score_dec and gap >= 2.0 and lean_inc >= 2:
        suggested = inc_letter
        confidence = "high" if gap >= 4 and lean_inc >= 2 else "medium"
        agreeing = lean_inc
    elif gap >= 3.0 and (lean_dec >= 1 or lean_inc >= 1):
        # single strong segment — medium at best
        suggested = dec_letter if score_dec > score_inc else inc_letter
        confidence = "medium"
        agreeing = max(lean_dec, lean_inc)
    else:
        suggested, confidence, agreeing = None, "low", 0

    return {
        "kind": "opengenes-masked",
        "status": "ok",
        "version": 2,
        "provided_by": "longevity-analyst prompt-evidence",
        "note": "Prompt text only (segmented). Do not query OpenGenes. Masked Gene X stays masked.",
        "has_animal_model": has_animal,
        "segments": segments,
        "totals": {
            "score_increase": round(score_inc, 2),
            "score_decrease": round(score_dec, 2),
            "segments_lean_increase": lean_inc,
            "segments_lean_decrease": lean_dec,
            "explicit_decrease_with_age": explicit_dec,
            "explicit_increase_with_age": explicit_inc,
            "agreeing_segments": agreeing,
        },
        "option_map": {"increase_letter": inc_letter, "decrease_letter": dec_letter},
        "suggested_letter": suggested,
        "confidence": confidence,
        "host_hint": (
            "emit uncertain"
            if confidence == "low" or suggested is None
            else f"prefer {suggested}; still weigh remaining prose"
        ),
    }


# ------------------------------------------------------------------ SynergyAge
def _effect_table(text: str) -> Dict[str, float]:
    effects: Dict[str, float] = {}
    for name, pct in re.findall(
        r"([A-Za-z0-9_();\+\-]+):\s*([+\-]?\d+(?:\.\d+)?)%\s*lifespan",
        text,
    ):
        effects[name] = float(pct)
    for name, pct in re.findall(
        r"([A-Za-z0-9_();\+\-]+):\s*([+\-]?\d+(?:\.\d+)?)%\s*(?:,|\.|Mean)",
        text,
    ):
        effects.setdefault(name, float(pct))
    for name, pct in re.findall(
        r"([A-Za-z0-9_();\+\-]+):\s*[+\-]?\d+(?:\.\d+)?%(?:,\s*[+\-]?\d+(?:\.\d+)?%)*\.\s*Mean:\s*([+\-]?\d+(?:\.\d+)?)%",
        text,
    ):
        effects[name] = float(pct)
    return effects


def _score_label(label: str, effects: Dict[str, float]) -> Optional[float]:
    lab = label.replace(" ", "")
    if lab in effects:
        return effects[lab]
    parts = [p for p in re.split(r"[;+]", lab) if p]
    if not parts:
        return None
    vals = []
    for p in parts:
        if p in effects:
            vals.append(effects[p])
        else:
            hit = next((effects[k] for k in effects if k.lower() == p.lower()), None)
            if hit is None:
                return None
            vals.append(hit)
    return float(sum(vals))


def synergyage_pairwise(text: str) -> Dict[str, Any]:
    effects = _effect_table(text)
    a_m = re.search(r"Options:\s*A\.\s*(\S+)", text)
    b_m = re.search(r"B\.\s*(\S+)", text)
    if not a_m or not b_m:
        return {
            "kind": "synergyage-pairwise",
            "status": "input_problem",
            "message": "could not parse Options A/B",
            "effects": effects,
        }
    label_a, label_b = a_m.group(1), b_m.group(1)
    sa = effects.get(label_a.replace(" ", ""))
    sb = effects.get(label_b.replace(" ", ""))
    additive_a = additive_b = False
    if sa is None:
        sa = _score_label(label_a, effects)
        additive_a = sa is not None and (";" in label_a or "+" in label_a)
    if sb is None:
        sb = _score_label(label_b, effects)
        additive_b = sb is not None and (";" in label_b or "+" in label_b)
    out: Dict[str, Any] = {
        "kind": "synergyage-pairwise",
        "status": "ok",
        "provided_by": "longevity-analyst prompt-evidence",
        "note": "Printed % lifespan only (exact compound if present; else additive). Do not query SynergyAge DB.",
        "effects": effects,
        "options": {"A": label_a, "B": label_b},
        "scores": {"A": sa, "B": sb},
        "additive_prior": {"A": additive_a, "B": additive_b},
    }
    if sa is None or sb is None:
        out["suggested_letter"] = None
        out["confidence"] = "low"
        out["host_hint"] = "emit uncertain — missing effect for an allele"
        return out
    gap = abs(sa - sb)
    if additive_a or additive_b:
        if gap < 20.0:
            out["suggested_letter"] = None
            out["confidence"] = "low"
            out["host_hint"] = (
                "emit uncertain — compound allele scored only by additivity; "
                "true synergy may flip the winner"
            )
            return out
    if gap < 1.0:
        out["suggested_letter"] = None
        out["confidence"] = "low"
        out["host_hint"] = "emit uncertain — lifespan delta near zero"
        return out
    out["suggested_letter"] = "A" if sa >= sb else "B"
    out["confidence"] = "high" if gap >= 5.0 and not (additive_a or additive_b) else "medium"
    out["host_hint"] = f"prefer {out['suggested_letter']} (score A={sa}, B={sb})"
    return out


def run(kind: str, text: str) -> Dict[str, Any]:
    if kind not in KINDS:
        return {"status": "bad_input", "message": f"--kind must be one of {', '.join(KINDS)}"}
    if not text.strip():
        return {"status": "bad_input", "message": "empty input"}
    if kind == "opengenes-masked":
        return opengenes_masked(text)
    return synergyage_pairwise(text)
