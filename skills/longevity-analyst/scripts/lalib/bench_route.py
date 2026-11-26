"""Route one LongevityBench prompt to a formula or the host LLM.

Member reports stay on workflows 01–06. This module does not read a workspace,
does not call a clock, and does not query OpenGenes or SynergyAge.

Default paths are formula | host | uncovered. Paper Longevity-LLM names appear
only under `control` metadata for optional对照 runs (--control-model).
lb_id / detect_task are eval-adapter hints for scoring tables; the product path
is schema (family/format + whether a formula can run), not "task id → download model".
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Dict, Optional

# Paper names and the repositories on the Longevity-LLM collection (control only).
REPOS = {
    "L-Qwen3.5-9B": "insilicomedicine/longevity-llm",
    "L-Qwen3-0.6B": "insilicomedicine/Qwen3-0.6B-Longevity",
    "L-Qwen3-1.7B": "insilicomedicine/Qwen3-1.7B-Longevity",
    "L-LFM2-1.2B": "LiquidAI/LFM2-1.2B-Longevity",
    "L-LFM2-2.6B": "LiquidAI/LFM2-2.6B-Longevity",
}

_DATA = Path(__file__).resolve().parents[2] / "data"
_MODEL = _DATA / "clinical_models.joblib"
_COX = _DATA / "clinical_cox.json"
_PROTEIN = _DATA / "protein_elasticnet.json"
_PROTEIN_CLOCK = _DATA / "protein_clock.json"
_PROTEIN_ALT = (
    Path(os.environ["LONGEVITYBENCH_OLINK"]) / "protein_elasticnet.json"
    if os.environ.get("LONGEVITYBENCH_OLINK")
    else None
)
# Optional local checkpoints via LONGEVITY_LLM_HOME/<repo-name> (control only).
_LLM_HOME = Path(os.environ["LONGEVITY_LLM_HOME"]) if os.environ.get("LONGEVITY_LLM_HOME") else None
_LOCAL_MODELS = {
    "L-Qwen3.5-9B": (_LLM_HOME / "longevity-llm") if _LLM_HOME else None,
    "L-LFM2-1.2B": (_LLM_HOME / "LFM2-1.2B-Longevity") if _LLM_HOME else None,
}


def _control(name: str) -> Dict[str, Optional[str]]:
    local = _LOCAL_MODELS.get(name)
    return {
        "name": name,
        "repo": REPOS[name],
        "local_path": str(local) if local else None,
    }


# paper_best is a published reference number for对照 tables (not a Skill pass bar).
_TASKS = {
    "LB-0034": {
        "family": "clinical_age", "format": "pairwise", "kind": "formula",
        "target": "clinical_hgb",
        "do_not_call": ["accelerated-biological-aging-risk", "phenoage"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.81, "higher_is_better": True,
            "source": "Table S6 ML:hgb",
        },
        "control": None,
        "note": "Gradient boosting fit on NHANES 2003-2017. Eval years are 1999 and 2001.",
    },
    "LB-0030": {
        "family": "clinical_age", "format": "multiclass", "kind": "formula",
        "target": "clinical_elasticnet",
        "do_not_call": ["accelerated-biological-aging-risk", "phenoage"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.43, "higher_is_better": True,
            "source": "Table S6 ML:logreg",
        },
        "control": None,
        "note": "Pooled elastic-net logistic regression fit on NHANES 2003-2017.",
    },
    "LB-0038": {
        "family": "clinical_age", "format": "regression", "kind": "formula",
        "target": "clinical_elasticnet",
        "do_not_call": ["accelerated-biological-aging-risk", "phenoage"],
        "wired": True,
        "paper_best": {
            "metric": "mae_years", "value": 11.55, "higher_is_better": False,
            "source": "Table S6 ML:elasticnet",
        },
        "control": None,
        "note": "Pooled elastic net fit on NHANES 2003-2017.",
    },
    "LB-0042": {
        "family": "clinical_mortality", "format": "binary", "kind": "formula",
        "target": "clinical_cox",
        "do_not_call": ["one-year-mortality-aging-marker"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.81, "higher_is_better": True,
            "source": "Table S6 ML:cox (best LLM Sonnet-4.6 0.808)",
        },
        "control": None,
        "note": "Cox fit on NHANES 2003-2017. Score with balanced accuracy, not raw accuracy.",
    },
    "LB-0050": {
        "family": "clinical_mortality", "format": "pairwise", "kind": "host",
        "target": "host",
        "do_not_call": ["one-year-mortality-aging-marker"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.684, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B (ML:cox 0.65)",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Default path is Host+tools. Longevity-LLM is control only.",
    },
    "LB-0046": {
        "family": "clinical_mortality", "format": "multiclass", "kind": "host",
        "target": "host",
        "do_not_call": ["one-year-mortality-aging-marker"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.389, "higher_is_better": True,
            "source": "Table S6 best LLM gemini-3-pro (ML:cox_rmst_bin 0.30)",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Default path is Host+tools. Longevity-LLM is control only.",
    },
    "LB-0006": {
        "family": "methylation", "format": "pairwise", "kind": "host",
        "target": "host",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.868, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Rotating sites; fixed clocks stay in do_not_call. Default is Host+tools.",
    },
    "LB-0010": {
        "family": "methylation", "format": "regression", "kind": "host",
        "target": "host",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "mae_years", "value": 12.96, "higher_is_better": False,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Rotating sites; fixed clocks stay in do_not_call. Default is Host+tools.",
    },
    "LB-0002": {
        "family": "methylation", "format": "multiclass", "kind": "host",
        "target": "host",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": None, "higher_is_better": True,
            "source": "parity is 9B's own Table S3 score; task-best LLM grok-4 is 0.203",
            "task_best_llm": {"model": "grok-4", "value": 0.203},
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Default path is Host+tools. Do not treat grok-4 0.203 as a Skill pass bar.",
    },
    "LB-0018": {
        "family": "transcriptome", "format": "pairwise", "kind": "host",
        "target": "host",
        "do_not_call": [],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.676, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Default path is Host+tools. Longevity-LLM is control only.",
    },
    "LB-0014": {
        "family": "transcriptome", "format": "multiclass", "kind": "host",
        "target": "host",
        "do_not_call": [],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.223, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Default path is Host+tools. Longevity-LLM is control only.",
    },
    "LB-0299": {
        "family": "protein", "format": "pairwise", "kind": "formula",
        "target": "protein_clock",
        "do_not_call": ["organ-aging-plasma-proteome"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.75, "higher_is_better": True,
            "source": "Table S6 ML:protein_clock",
        },
        "control": None,
        "note": "Goeminne OrganAge Conventional fold-0 on prompt NPX.",
    },
    "LB-0074": {
        "family": "protein", "format": "multiclass", "kind": "formula",
        "target": "protein_clock",
        "do_not_call": ["organ-aging-plasma-proteome"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.56, "higher_is_better": True,
            "source": "Table S6 ML:protein_clock",
        },
        "control": None,
        "note": "Goeminne OrganAge Conventional fold-0; age decade from predicted age.",
    },
    "LB-0082": {
        "family": "protein", "format": "regression", "kind": "formula",
        "target": "protein_elasticnet",
        "do_not_call": ["organ-aging-plasma-proteome"],
        "wired": True,
        "paper_best": {
            "metric": "mae_years", "value": 3.94, "higher_is_better": False,
            "source": "Table S6 ML:elasticnet",
        },
        "control": None,
        "note": "De novo elastic net on consensus proteins; eval kitGuids excluded from fit.",
    },
    "LB-0086": {
        "family": "opengenes", "format": "binary", "kind": "host",
        "target": "host",
        "do_not_call": ["opengenes"],
        "wired": True,
        "paper_best": {
            "metric": "accuracy", "value": 0.780, "higher_is_better": True,
            "source": "main text L-LFM2-1.2B (9B 0.710)",
        },
        "control": _control("L-LFM2-1.2B"),
        "note": "Masked gene name is not looked up in OpenGenes. Default is Host+tools.",
    },
    "LB-0122": {
        "family": "synergyage", "format": "pairwise", "kind": "host",
        "target": "host",
        "do_not_call": ["synergy-age"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.846, "higher_is_better": True,
            "source": "Table S3 Claude-Opus-4.5 (L-Qwen3.5-9B 0.731)",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Do not query SynergyAge on exam prompts. Default is Host+tools.",
    },
    "LB-0126": {
        "family": "synergyage", "format": "regression", "kind": "host",
        "target": "host",
        "do_not_call": ["synergy-age"],
        "wired": True,
        "paper_best": {
            "metric": "mae_percent_lifespan", "value": 11.96, "higher_is_better": False,
            "source": "Table S3 Claude-Opus-4.5 (L-Qwen3.5-9B 12.38)",
        },
        "control": _control("L-Qwen3.5-9B"),
        "note": "Do not query SynergyAge on exam prompts. Default is Host+tools.",
    },
}

_MINI = {
    "LB-0301": "LB-0002",
    "LB-0302": "LB-0006",
    "LB-0303": "LB-0010",
    "LB-0304": "LB-0014",
    "LB-0305": "LB-0018",
    "LB-0306": "LB-0030",
    "LB-0307": "LB-0034",
    "LB-0308": "LB-0038",
    "LB-0309": "LB-0042",
    "LB-0310": "LB-0046",
    "LB-0311": "LB-0050",
    "LB-0312": "LB-0074",
    "LB-0313": "LB-0082",
    "LB-0314": "LB-0086",
    "LB-0315": "LB-0122",
    "LB-0316": "LB-0126",
    "LB-0317": "LB-0299",
    # Eval parquet uses LB-0220 for transcriptome pairwise (paper id LB-0018).
    "LB-0220": "LB-0018",
}


def canonical_id(lb_id: Optional[str]) -> Optional[str]:
    if not lb_id:
        return None
    key = lb_id.strip().upper()
    return _MINI.get(key, key if key in _TASKS else None)


def detect_task(text: str) -> Optional[str]:
    """Identify a paper task from the user prompt alone. Unknown text returns None."""
    t = text.lower()
    if re.search(r"\bgene\s*x\b", t) and "with age" in t:
        return "LB-0086"
    if "lifespan" in t and ("perturb" in t or "synergy" in t or "gene-pair" in t or "gene pair" in t):
        if "percent" in t or "%" in t:
            return "LB-0126"
        return "LB-0122"
    if "methylation" in t or re.search(r"\bcg\d{6,}\b", t):
        if "decade" in t or "age group" in t or "which age bin" in t:
            return "LB-0002"
        if "how old" in t or "chronological age" in t or "predict" in t and "age" in t and "older" not in t:
            return "LB-0010"
        return "LB-0006"
    if "olink" in t or "proteom" in t or "npx" in t:
        if "decade" in t or "age group" in t:
            return "LB-0074"
        if "older" in t:
            return "LB-0299"
        return "LB-0082"
    if "gtex" in t or "transcript" in t:
        if "decade" in t or "age group" in t:
            return "LB-0014"
        return "LB-0018"
    mortality = any(w in t for w in ("deceased", "survived longer", "still alive", "years of follow", "mortality"))
    clinical = any(w in t for w in ("clinical", "blood panel", "blood test", "nhanes"))
    if mortality:
        if "which" in t and ("survived" in t or "longer" in t):
            return "LB-0050"
        if "0-5" in t or "5-10" in t or "bin" in t or "category" in t:
            return "LB-0046"
        return "LB-0042"
    if clinical or "which participant is older" in t:
        if "age group" in t or "decade" in t:
            return "LB-0030"
        if "how old" in t or "chronological age" in t or ("predict" in t and "older" not in t):
            return "LB-0038"
        return "LB-0034"
    return None


def _formula_wired(target: str) -> tuple[bool, str]:
    if target in ("clinical_hgb", "clinical_elasticnet"):
        if _MODEL.is_file():
            return True, ""
        return False, f"missing formula weights: {_MODEL}"
    if target == "clinical_cox":
        if _COX.is_file():
            return True, ""
        return False, f"missing formula weights: {_COX}"
    if target == "protein_elasticnet":
        alt_ok = bool(_PROTEIN_ALT and _PROTEIN_ALT.is_file())
        if _PROTEIN.is_file() or alt_ok:
            return True, ""
        return False, f"missing formula weights: {_PROTEIN}"
    if target == "protein_clock":
        if _PROTEIN_CLOCK.is_file():
            return True, ""
        return False, f"missing formula weights: {_PROTEIN_CLOCK}"
    return False, f"unknown formula target: {target}"


def _families() -> Dict[str, Dict[str, object]]:
    """Collapse task table into schema keys (family, format) → route policy.

    Product path uses this map. LongevityBench ids are eval metadata only.
    """
    out: Dict[str, Dict[str, object]] = {}
    for lb_id, spec in _TASKS.items():
        key = f"{spec['family']}|{spec['format']}"
        # First canonical task wins; mini aliases share the same family/format.
        if key not in out:
            out[key] = {
                "family": spec["family"],
                "format": spec["format"],
                "kind": spec["kind"],
                "target": spec["target"],
                "do_not_call": list(spec["do_not_call"]),
                "wired": bool(spec["wired"]),
                "note": spec["note"],
                "example_lb_id": lb_id,  # eval hint only
                "control": spec.get("control"),
                "paper_best": spec.get("paper_best"),
            }
    return out


_SCHEMA = _families()


def route_by_schema(family: str, fmt: str) -> Dict[str, object]:
    """Route by input schema alone — the product-facing entry point.

    Does not require an LB-xxxx id. Does not download Longevity-LLM.
    """
    family = (family or "").strip().lower()
    fmt = (fmt or "").strip().lower()
    key = f"{family}|{fmt}"
    spec = _SCHEMA.get(key)
    if not spec:
        known = sorted(_SCHEMA.keys())
        return {
            "status": "unrouted",
            "capability": "A",
            "family": family,
            "format": fmt,
            "path": "uncovered",
            "callable": False,
            "reason": f"unknown schema {key!r}; known: {known}",
        }
    kind = spec["kind"]
    target = spec["target"]
    missing = ""
    if kind == "formula":
        wired, missing = _formula_wired(target)
        path = "formula" if wired else "uncovered"
    else:
        wired = bool(spec["wired"])
        path = "host"
    out: Dict[str, object] = {
        "status": "routed",
        "capability": "A",
        "family": family,
        "format": fmt,
        "kind": kind,
        "target": target,
        "path": path,
        "do_not_call": spec["do_not_call"],
        "wired": wired,
        "callable": False,
        "writes_report": False,
        "host_protocol": fmt if path == "host" else None,
        "note": spec["note"],
        "schema_key": key,
        # Eval-only fields kept optional for adapters that still tabulate by id.
        "example_lb_id": spec.get("example_lb_id"),
        "paper_best": spec.get("paper_best"),
        "control": spec.get("control"),
    }
    if missing:
        out["reason"] = missing
    if path == "uncovered" and kind == "formula":
        out["note"] = (spec["note"] + " " if spec["note"] else "") + (
            "Formula weights missing; treat as uncovered until weights are installed."
        )
    return out


def list_schemas() -> list:
    return sorted(_SCHEMA.keys())


def route(text: str = "", lb_id: Optional[str] = None) -> Dict[str, object]:
    """Return the route for one exam prompt (eval adapter).

    Prefer `route_by_schema` for product paths. `path` is formula | host |
    uncovered. `callable` stays false until a formula caller fills an exam
    field. Host answers are produced by the host LLM outside this module.
    Longevity-LLM is never the default path.
    """
    paper_id = canonical_id(lb_id) or detect_task(text)
    if paper_id not in _TASKS:
        return {
            "status": "unrouted",
            "capability": "A",
            "lb_id": paper_id,
            "path": "uncovered",
            "callable": False,
            "reason": "prompt is not one of the 17 LongevityBench tasks",
        }
    spec = _TASKS[paper_id]
    # Delegate policy to schema route so product and eval stay consistent.
    base = route_by_schema(spec["family"], spec["format"])
    base["lb_id"] = paper_id
    base["paper_best"] = spec["paper_best"]
    base["control"] = spec.get("control")
    base["note"] = spec["note"]
    if base.get("path") == "uncovered" and spec["kind"] == "formula":
        wired, missing = _formula_wired(spec["target"])
        if not wired and missing:
            base["reason"] = missing
    if spec.get("candidates"):
        base["candidates"] = spec["candidates"]
    return base


def route_json(text: str = "", lb_id: Optional[str] = None) -> str:
    return json.dumps(route(text, lb_id), ensure_ascii=False, indent=2)
