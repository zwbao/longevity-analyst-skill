"""Route one LongevityBench prompt to a formula or an open model.

Member reports stay on workflows 01–06. This module does not read a workspace,
does not call a clock, and does not query OpenGenes or SynergyAge.
Open models are named here. Inference runs only after a runner is wired.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Dict, Optional

# Paper names and the repositories on the Longevity-LLM collection.
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
_PROTEIN_ALT = (
    Path(os.environ["LONGEVITYBENCH_OLINK"]) / "protein_elasticnet.json"
    if os.environ.get("LONGEVITYBENCH_OLINK")
    else None
)
# Optional local checkpoints via LONGEVITY_LLM_HOME/<repo-name>.
_LLM_HOME = Path(os.environ["LONGEVITY_LLM_HOME"]) if os.environ.get("LONGEVITY_LLM_HOME") else None
_LOCAL_MODELS = {
    "L-Qwen3.5-9B": (_LLM_HOME / "longevity-llm") if _LLM_HOME else None,
    "L-LFM2-1.2B": (_LLM_HOME / "LFM2-1.2B-Longevity") if _LLM_HOME else None,
}

# paper_best is the better of Table S6 ML and Best LLM (lower MAE for regression).
# SynergyAge has no ML in S6; absolute scores come from Table S3 (mmc3.xlsx).
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
        "repo": None,
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
        "repo": None,
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
        "repo": None,
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
        "repo": None,
        "note": "Cox fit on NHANES 2003-2017. Score with balanced accuracy, not raw accuracy.",
    },
    "LB-0050": {
        "family": "clinical_mortality", "format": "pairwise", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["one-year-mortality-aging-marker"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.684, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B (ML:cox 0.65)",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "Published best is the 9B. Local Cox stays a control only and is not the parity target.",
    },
    "LB-0046": {
        "family": "clinical_mortality", "format": "multiclass", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["one-year-mortality-aging-marker"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.389, "higher_is_better": True,
            "source": "Table S6 best LLM gemini-3-pro (ML:cox_rmst_bin 0.30)",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "Published best is an LLM. Local bin classifier is a control only.",
    },
    "LB-0006": {
        "family": "methylation", "format": "pairwise", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.868, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "One user prompt in, one exam field out.",
    },
    "LB-0010": {
        "family": "methylation", "format": "regression", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "mae_years", "value": 12.96, "higher_is_better": False,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "One user prompt in, one exam field out.",
    },
    "LB-0002": {
        "family": "methylation", "format": "multiclass", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["epiage", "pyaging", "horvath"],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": None, "higher_is_better": True,
            "source": "parity is 9B's own Table S3 score; task-best LLM grok-4 is 0.203",
            "task_best_llm": {"model": "grok-4", "value": 0.203},
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "Route stays on 9B. Do not use grok-4's 0.203 as the 9B parity number.",
    },
    "LB-0018": {
        "family": "transcriptome", "format": "pairwise", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": [],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.676, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "One user prompt in, one exam field out.",
    },
    "LB-0014": {
        "family": "transcriptome", "format": "multiclass", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": [],
        "wired": True,
        "paper_best": {
            "metric": "balanced_accuracy", "value": 0.223, "higher_is_better": True,
            "source": "Table S6 best LLM L-Qwen3.5-9B",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "One user prompt in, one exam field out.",
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
        "repo": None,
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
        "repo": None,
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
        "repo": None,
        "note": "De novo elastic net on consensus proteins; eval kitGuids excluded from fit.",
    },
    "LB-0086": {
        "family": "opengenes", "format": "binary", "kind": "model",
        "target": "L-LFM2-1.2B",
        "do_not_call": ["opengenes"],
        "wired": True,
        "paper_best": {
            "metric": "accuracy", "value": 0.780, "higher_is_better": True,
            "source": "main text L-LFM2-1.2B (9B 0.710)",
        },
        "repo": REPOS["L-LFM2-1.2B"],
        "note": "A masked gene name is not looked up in OpenGenes.",
    },
    "LB-0122": {
        "family": "synergyage", "format": "pairwise", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["synergy-age"],
        "wired": True,
        "paper_best": {
            "metric": "concordance", "value": 0.846, "higher_is_better": True,
            "source": "Table S3 Claude-Opus-4.5 (L-Qwen3.5-9B 0.731)",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "Exam prompts are not answered by querying SynergyAge. Score via MLX 9B when run.",
    },
    "LB-0126": {
        "family": "synergyage", "format": "regression", "kind": "model",
        "target": "L-Qwen3.5-9B",
        "do_not_call": ["synergy-age"],
        "wired": True,
        "paper_best": {
            "metric": "mae_percent_lifespan", "value": 11.96, "higher_is_better": False,
            "source": "Table S3 Claude-Opus-4.5 (L-Qwen3.5-9B 12.38)",
        },
        "repo": REPOS["L-Qwen3.5-9B"],
        "note": "Exam prompts are not answered by querying SynergyAge. Score via MLX 9B when run.",
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


def route(text: str = "", lb_id: Optional[str] = None) -> Dict[str, object]:
    """Return the route for one exam prompt.

    `callable` stays false until a connected caller returns the exam field.
    """
    paper_id = canonical_id(lb_id) or detect_task(text)
    if paper_id not in _TASKS:
        return {
            "status": "unrouted",
            "lb_id": paper_id,
            "callable": False,
            "reason": "prompt is not one of the 17 LongevityBench tasks",
        }
    spec = _TASKS[paper_id]
    if spec["target"] in ("clinical_hgb", "clinical_elasticnet"):
        wired = bool(spec["wired"] and _MODEL.is_file())
    elif spec["target"] == "clinical_cox":
        wired = bool(spec["wired"] and _COX.is_file())
    elif spec["target"] == "protein_elasticnet":
        alt_ok = bool(_PROTEIN_ALT and _PROTEIN_ALT.is_file())
        wired = bool(spec["wired"] and (_PROTEIN.is_file() or alt_ok))
    elif spec["target"] == "protein_clock":
        clock = _DATA / "protein_clock.json"
        wired = bool(spec["wired"] and clock.is_file())
    elif spec["kind"] == "model":
        local = _LOCAL_MODELS.get(spec["target"])
        wired = bool(spec["wired"] and local and (local / "config.json").is_file())
    else:
        wired = bool(spec["wired"])
    local_path = _LOCAL_MODELS.get(spec["target"])
    out = {
        "status": "routed",
        "lb_id": paper_id,
        "family": spec["family"],
        "format": spec["format"],
        "kind": spec["kind"],
        "target": spec["target"],
        "do_not_call": spec["do_not_call"],
        "wired": wired,
        "callable": False,
        "writes_report": False,
        "paper_best": spec["paper_best"],
        "repo": spec["repo"],
        "local_path": str(local_path) if local_path else None,
        "note": spec["note"],
    }
    if spec.get("candidates"):
        out["candidates"] = spec["candidates"]
    return out


def route_json(text: str = "", lb_id: Optional[str] = None) -> str:
    return json.dumps(route(text, lb_id), ensure_ascii=False, indent=2)
