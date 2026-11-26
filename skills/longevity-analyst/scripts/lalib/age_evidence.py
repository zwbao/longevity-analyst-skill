"""Evidence for "how old does this person look", "which of two people is older" and "how long will they survive",
computed from measured values the user pasted (blood tests, DNA methylation, plasma proteins, gene expression).

Calculators only. The host reads the evidence together with everything else in the message and makes the call.
No benchmark ids, no answer letters: every method here also runs on a member's own report.

Methods (each result names its source and how it was checked):
  labs, age        NHANES age reference model (fitted on NHANES 2011-2016, held out 2017-2018): the age a set of
                   routine values looks like, a 90% interval, band probabilities, P(A older) for two people.
  labs, survival   Cox model on NHANES 2003-2017 with linked mortality (from PR #2's data/clinical_cox.json): survival
                   probability at 5/10/15 years and P(A outlives B). Needs the chronological age.
  methylation/proteome, two people   age-direction vote over every CpG/protein with a published age association
                   (Generation Scotland EWAS; China Kadoorie Biobank Olink), neither a LongevityBench source.
  proteome, one person   Goeminne et al. OrganAge conventional plasma clock (data/protein_clock.json).
  methylation, one person   no calibrated single-profile estimator: strongest age-associated CpGs listed for the host.
  expression, two top-lists   literature gene-set direction vote on Sample A/B Top-N genes (not GTEx-fitted);
                   low_confidence when few hits or near-tie → Host should prefer uncertain.

Standard library only.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from .agekit import compare, expression_vote, panel, reference
from .agekit.library import Library

SKILL_DIR = Path(__file__).resolve().parents[2]
LIB = Library(SKILL_DIR)
MODALITIES = ("labs", "methylation", "proteome", "expression")


# ------------------------------------------------------------------ input
def _strip_wrappers(text: str) -> str:
    m = re.search(r"<command-args>(.*)</command-args>", text, re.S)
    return m.group(1) if m else text


def _claude_transcripts(cwd: str) -> Optional[str]:
    root = Path.home() / ".claude" / "projects"
    now = time.time()
    files = [f for f in glob.glob(str(root / "*" / "*.jsonl")) if now - os.path.getmtime(f) < 6 * 3600]
    files.sort(key=os.path.getmtime, reverse=True)
    for f in files[:80]:
        last = None
        with open(f, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if '"user"' not in line:
                    continue
                try:
                    e = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if e.get("type") != "user" or not e.get("cwd") or not (cwd == e["cwd"] or cwd.startswith(e["cwd"].rstrip("/") + "/")):
                    continue
                c = (e.get("message") or {}).get("content")
                text = c if isinstance(c, str) else "\n".join(x.get("text", "") for x in (c or []) if isinstance(x, dict) and x.get("type") == "text")
                if text and not text.startswith("Base directory for this skill") and len(text) > 40:
                    last = text
        if last:
            return _strip_wrappers(last)
    return None


def _codex_transcripts(cwd: str) -> Optional[str]:
    root = Path.home() / ".codex" / "sessions"
    now = time.time()
    files = [f for f in glob.glob(str(root / "**" / "rollout-*.jsonl"), recursive=True) if now - os.path.getmtime(f) < 6 * 3600]
    files.sort(key=os.path.getmtime, reverse=True)
    for f in files[:40]:
        same, last = False, None
        with open(f, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                try:
                    e = json.loads(line)
                except json.JSONDecodeError:
                    continue
                p = e.get("payload") or {}
                if e.get("type") == "session_meta" and p.get("cwd"):
                    same = cwd == p["cwd"] or cwd.startswith(p["cwd"].rstrip("/") + "/")
                if same and p.get("type") == "message" and p.get("role") == "user":
                    text = "\n".join(x.get("text", "") for x in p.get("content") or [] if isinstance(x, dict))
                    if len(text) > 40 and not text.lstrip().startswith("<"):
                        last = text
        if same and last:
            return last
    return None


def read_input(file: Optional[str] = None, text: Optional[str] = None, from_transcript: bool = False) -> Dict[str, Any]:
    if text:
        return {"text": text, "source": "--text"}
    if file:
        return {"text": Path(file).read_text(encoding="utf-8", errors="ignore"), "source": f"file {file}"}
    if from_transcript:
        cwd = os.getcwd()
        got = _claude_transcripts(cwd) or _codex_transcripts(cwd)
        if got:
            return {"text": got, "source": "latest user message in this session's transcript"}
        return {"text": "", "source": "transcript", "problem": "no transcript for this working directory; write the data to a file and pass --file"}
    return {"text": "", "source": "none", "problem": "pass --file, --text or --from-transcript"}


# ------------------------------------------------------------------ labs
def _p_older(pa: Dict[str, Any], pb: Dict[str, Any]) -> float:
    ages, a = pa["ages"], pa["p"]
    b = pb["p"]
    cum_b, below = [], 0.0
    for v in b:
        cum_b.append(below)
        below += v
    return sum(a[i] * (cum_b[i] + 0.5 * b[i]) for i in range(len(ages)))


def _cox(block: str, age: Optional[float]) -> Optional[Dict[str, Any]]:
    from . import clinical_age, clinical_mortality
    model = clinical_mortality.load()
    if not model or age is None:
        return None
    person = clinical_age._parse_block(block)
    person["age"] = float(age)
    missing = [f for f in model["features"] if f not in person]
    lp = clinical_mortality._lp(person, model)
    surv = {y: clinical_mortality._survival(model, lp, 12 * y) for y in (5, 10, 15)}
    median = None
    for t, s in model["baseline"]:
        if s ** math.exp(lp) <= 0.5:
            median = round(t / 12, 1)
            break
    return {"linear_predictor": round(lp, 3),
            "survival_probability": {f"{y}y": round(v, 3) for y, v in surv.items()},
            "median_survival_years": median if median is not None else ">15 (beyond the follow-up curve)",
            "features_used": len(model["features"]) - len(missing),
            "features_missing_filled_with_population_median": missing}


def labs(text: str, bands: Optional[str]) -> Dict[str, Any]:
    blocks = panel._split_people(text)
    people = panel.parse_panel(text)
    out: Dict[str, Any] = {"modality": "labs", "people": []}
    posts = []
    for blk, person in zip(blocks, people):
        entry: Dict[str, Any] = {"label": person["label"], "sex": person["sex"], "chronological_age_in_text": person["age"],
                                 "values_read": len(person["rows"])}
        if person["age"] is None:
            est = reference.estimate_age(LIB, person["rows"], person["sex"], bands, with_posterior=True)
            posts.append(est.pop("_posterior", None))
            entry["lab_age_reference"] = est
        else:
            posts.append(None)
            entry["lab_age_reference"] = {"status": "skipped", "reason": "the age is given; how old the labs look is not a "
                                          "validated predictor of survival, so it is not computed here"}
        cox = _cox(blk["text"], person["age"])
        if cox:
            entry["survival_cox"] = {"method": "Cox proportional hazards, NHANES 2003-2017 with public linked mortality "
                                               "(fitted by the PR #2 author; age, sex, BMI and 27 labs)",
                                     **cox,
                                     "caveat": "US population; uses age plus labs only. Medical history, medicines, "
                                               "smoking and frailty in the text are NOT in this model: weigh them yourself. "
                                               "It does not predict how long someone survives (time bands) well: do not map "
                                               "these probabilities onto survival-time options."}
        out["people"].append(entry)
    if len(people) == 2 and all(posts):
        out["p_A_older_by_lab_age"] = round(_p_older(posts[0], posts[1]), 3)
    if len(people) == 2 and all("survival_cox" in p for p in out["people"]):
        la, lb = (p["survival_cox"]["linear_predictor"] for p in out["people"])
        pa = 1 / (1 + math.exp(la - lb))
        out["p_A_outlives_B_by_cox"] = round(pa, 3)
        out["cox_comparison_strength"] = ("cannot separate them (0.40-0.60): decide from history, medicines and the rest"
                                          if 0.4 <= pa <= 0.6 else "moderate" if 0.3 <= pa <= 0.7 else "clear")
    if any("survival_cox" in p for p in out["people"]):
        out["cox_discrimination"] = ("modest: age-and-lab Cox models rank two people correctly about 2 times in 3 "
                                     "(LongevityBench paper's Cox baseline, concordance 0.65); history, medicines, "
                                     "smoking and frailty often add more")
    out["how_to_use"] = ("lab_age_reference answers 'what age do these values look like' (held-out mean error ~10 years). "
                         "When the age is given in the text, use survival_cox instead and combine it with history and medicines.")
    return out


# ------------------------------------------------------------------ omics
def _single_methylation(text: str) -> Dict[str, Any]:
    vals = compare._single(text)
    ref = compare._load(LIB, "methylation")
    present = [(k, v, ref[k][0]) for k, v in vals.items() if k in ref]
    present.sort(key=lambda x: -abs(x[2]))
    return {"status": "no_calibrated_estimator",
            "message": "No single-profile methylation age estimator here: the listed sites are not a fixed clock's sites. "
                       "Below are the most strongly age-associated CpGs in this profile (Generation Scotland EWAS, n~18k); "
                       "reason from them and from the rest of the text.",
            "cpgs_read": len(vals), "cpgs_with_age_association": len(present),
            "strongest": [{"cpg": k, "beta": v, "rises_with_age": z > 0, "z": round(z, 1)} for k, v, z in present[:15]]}


def _single_proteome(text: str, bands: Optional[str]) -> Dict[str, Any]:
    from . import protein_age
    model = protein_age._model_for_format("regression")
    row = protein_age.parse_npx(text)
    if not model or not row:
        return {"status": "input_problem", "message": "no 'PROTEIN (value)' values read or clock weights missing"}
    used = [f for f in model["features"] if f in row]
    age = protein_age._predict(row, model)
    out = {"status": "ok", "method": "Goeminne et al. OrganAge conventional plasma protein clock (data/protein_clock.json), "
                                     "linear calibration fitted by the PR #2 author on a training split",
           "predicted_age": round(age, 1), "proteins_read": len(row), "clock_proteins_present": len(used),
           "validation": "author-reported mean absolute error 4.4 years on 93 profiles of a 396-protein Olink panel held "
                         "out from the calibration (PR #2 data/scores_formula.json); not re-checked here",
           "caveat": "Clock proteins absent from the input are filled with the reference median, which pulls the "
                     "estimate towards the reference mean age when far fewer than ~300 are present."}
    if bands:
        try:
            custom = reference.parse_bands(bands)
            out["band_of_point_estimate"] = reference.band_of(int(round(age)), custom)
        except Exception as exc:  # noqa: BLE001
            out["band_problem"] = str(exc)
    return out


def omics(text: str, modality: str, bands: Optional[str]) -> Dict[str, Any]:
    pairs = compare.parse_pairs(text)
    kind = "methylation" if modality == "methylation" else "proteomics"
    if len(pairs) >= 5:
        res = compare.compare_age(LIB, pairs, kind)
        return {"modality": modality, "two_profiles": True, "age_direction_vote": res}
    if modality == "methylation":
        return {"modality": modality, "two_profiles": False, "single_profile": _single_methylation(text)}
    return {"modality": modality, "two_profiles": False, "single_profile": _single_proteome(text, bands)}


def expression(text: str) -> Dict[str, Any]:
    """Pairwise top-gene lists → literature direction vote; else honest single-profile stub."""
    genes_a, genes_b = expression_vote.parse_top_lists(text)
    if len(genes_a) >= 10 and len(genes_b) >= 10:
        vote = expression_vote.compare_top_lists(SKILL_DIR, text)
        return {"modality": "expression", "two_profiles": True, "age_direction_vote": vote}
    return {
        "modality": "expression",
        "two_profiles": False,
        "status": "no_calibrated_estimator",
        "message": "No single-profile transcriptomic age estimator independent of GTEx. "
                   "For two Top-N gene lists, paste both Sample A/B blocks so the literature "
                   "gene-set vote can run. Otherwise reason from tissue biology and say how sure you are.",
        "tissue": expression_vote.parse_tissue(text),
        "n_genes_parsed": max(len(genes_a), len(genes_b)),
    }


def run(modality: str, text: str, bands: Optional[str] = None) -> Dict[str, Any]:
    if modality not in MODALITIES:
        return {"status": "bad_input", "message": f"--modality must be one of {', '.join(MODALITIES)}"}
    if not text.strip():
        return {"status": "bad_input", "message": "empty input"}
    if modality == "labs":
        return labs(text, bands)
    if modality == "expression":
        return expression(text)
    return omics(text, modality, bands)
