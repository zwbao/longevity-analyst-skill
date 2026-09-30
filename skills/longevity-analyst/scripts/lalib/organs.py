"""Organ checkup table: deterministic organ indices, AI estimates, and the per-organ view.

Three layers, never mixed:
  measured / model  - organ ages from published methods (readouts from `methods run`)
  calibrated index  - published formulas computed here from lab rows (eGFR, FIB-4, TyG, AIP) and China-PAR
  AI estimate       - an estimator subagent's ranges for organ age and disease probability, registered through
                      `la.py organ register`, kind `llm_estimate`, always shown in its own column

The harness checks the estimates' form (ranges, widths, basis ids, retrieved PMIDs, confidence); whether an
estimate is sensible is the reviewer's judgment.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import EXIT_INPUT, LAError, data, load_json, now_iso, sha256_file, skillkit, write_json


# ------------------------------------------------------------ deterministic calculators
def _convert(v: float, f: Any) -> float:
    return v * f[0] + f[1] if isinstance(f, list) else v * f


def lab_values(st: Dict[str, Any], analyte: str, blank_unit: Optional[str] = None) -> List[float]:
    """Every usable (confirmed, numeric, plausible) value for an analyte. `blank_unit` reads a row without a unit in
    that unit - used only where reading it errs on the safe side (a diagnostic-threshold override)."""
    from .labnames import analytes, clean_value, usable_rows
    spec = data("organs.json")["lab_aliases"][analyte]
    kit = skillkit()
    codes = {str(c).casefold() for c in spec.get("loinc", [analyte])}
    units = {kit.normalize_unit(u): f for u, f in spec["units"].items()}
    lo, hi = spec.get("plausible", [0, float("inf")])
    out = []
    for r in usable_rows(st):
        if not (analytes(r) & codes):
            continue
        try:
            v = kit.parse_number(clean_value(r["value"]))
        except ValueError:
            continue
        u = kit.normalize_unit(r.get("unit", "")) or (kit.normalize_unit(blank_unit) if blank_unit else "")
        if u not in units:
            continue
        v = _convert(v, units[u])
        if math.isfinite(v) and lo <= v <= hi:
            out.append(v)
    return out


def _lab_value(st: Dict[str, Any], analyte: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """(value in the canonical unit, problem). A row counts only when its name makes it a candidate for this analyte
    AND the agent confirmed it (`la.py labs confirm`) or mapped it (`la.py labs map`). Several confirmed rows with
    different values, a censored, non-finite or implausible value, or an unknown unit give None plus the reason;
    nothing is guessed. Reasons are member-facing Chinese."""
    from .labnames import analytes, clean_value, usable_rows
    ok = {id(r) for r in usable_rows(st)}
    spec = data("organs.json")["lab_aliases"][analyte]
    kit = skillkit()
    codes = {str(c).casefold() for c in spec.get("loinc", [analyte])}
    units = {kit.normalize_unit(u): f for u, f in spec["units"].items()}
    lo, hi = spec.get("plausible", [0, float("inf")])
    used: List[Dict[str, Any]] = []
    why: List[str] = []
    for r in st["labs"]:
        if not (analytes(r) & codes):
            continue
        if id(r) not in ok:
            c = (r.get("confirm") or {}).get("answer")
            if c != "no":
                why.append(f"「{r['marker']}」尚未确认是否为{spec.get('confirm_zh', analyte)}")
            continue
        raw = clean_value(r["value"])
        try:
            v = kit.parse_number(raw)
        except ValueError:
            why.append(f"「{r['marker']}」的结果「{raw}」不是可计算的数（删失值或非数字）")
            continue
        u = kit.normalize_unit(r.get("unit", ""))
        if u not in units:
            why.append(f"「{r['marker']}」的单位「{r.get('unit', '')}」不在可换算单位里")
            continue
        v = _convert(v, units[u])
        if not (math.isfinite(v) and lo <= v <= hi):
            why.append(f"「{r['marker']}」换算后超出生理可能范围，可能单位或数值有误")
            continue
        used.append({"value": v, "marker": r["marker"], "unit_in": r.get("unit", "")})
    vals = {round(x["value"], 6) for x in used}
    if len(vals) > 1:
        return None, f"{'、'.join(x['marker'] for x in used)} 有多行且数值不同，未自动取其一"
    if used:
        return used[0], ("未采用的行：" + "；".join(why)) if why else None
    return None, "；".join(why) or None


def egfr_ckd_epi_2021(scr_umol: float, age: int, sex: str) -> float:
    scr = scr_umol / 88.4
    if sex not in ("female", "male"):
        raise ValueError(f"sex must be female or male, got {sex!r}")   # never fall back to the male equation silently
    female = sex == "female"
    k, a = (0.7, -0.241) if female else (0.9, -0.302)
    v = 142 * min(scr / k, 1) ** a * max(scr / k, 1) ** -1.200 * 0.9938 ** age
    return v * 1.012 if female else v


def kdigo_g(egfr: float) -> str:
    if not math.isfinite(egfr):
        raise ValueError("eGFR must be finite")
    for lim, g in ((90, "G1"), (60, "G2"), (45, "G3a"), (30, "G3b"), (15, "G4")):
        if egfr >= lim:
            return g
    return "G5"


def fib4(age: int, ast: float, alt: float, plt: float) -> float:
    return age * ast / (plt * math.sqrt(alt))


def fib4_tier(v: float, age: int) -> Optional[str]:
    if age < 35:
        return None                                        # not validated under 35
    low = 2.0 if age >= 65 else 1.30
    if v < low:
        return "低于排除界值（在脂肪肝人群中进展期纤维化可能性低）"
    if v > 2.67:
        return "高于界值（建议就医进一步评估，如弹性成像）"
    return "界值之间（建议就医进一步评估）"


def compute_indices(st: Dict[str, Any]) -> Dict[str, Any]:
    """Readouts from published organ formulas; each missing input is named, nothing is imputed."""
    m = st["member"]
    age, sex = m.get("age"), m.get("sex")
    calc = data("organs.json")["calculators"]
    ro: List[Dict[str, Any]] = []
    notes: List[str] = []
    missing: Dict[str, List[str]] = {}

    def add(rid, key, value, unit="", kind="computed", **extra):
        c = calc[key]
        system = {"egfr": ["liver_kidney"], "fib4": ["liver_kidney"], "tyg": ["cardiometabolic"], "aip": ["cardiometabolic"]}[key]
        ro.append({"id": f"native.organ.{rid}", "label_zh": extra.pop("label_zh", c["label_zh"]), "value": value, "unit": unit,
                   "kind": kind, "basis": c["citation"], "doi": c["doi"], "systems": system, **extra})

    probs: Dict[str, str] = {}

    def lab(a):
        v, p = _lab_value(st, a)
        if p:
            probs[a] = p
        return v

    adult = bool(age) and int(age) >= 18
    cr = lab("creatinine")
    if cr and adult and sex in ("female", "male"):
        e = round(egfr_ckd_epi_2021(cr["value"], int(age), sex), 1)
        add("egfr", "egfr", e, calc["egfr"]["unit"], inputs=[cr["marker"]])
        add("egfr_stage", "egfr", kdigo_g(e), "", label_zh="肾功能分期（KDIGO GFR 分类，单次测量）", inputs=[cr["marker"]])
        notes.append("eGFR 分期来自单次肌酐：诊断慢性肾脏病需要间隔 3 个月以上的重复检查和尿白蛋白")
    else:
        missing["egfr"] = [x for x, ok in (("肌酐", cr), ("年龄（成年人）", adult), ("性别", sex)) if not ok]
    ast, alt, plt = lab("ast"), lab("alt"), lab("platelets")
    if ast and alt and plt and adult:
        raw = fib4(int(age), ast["value"], alt["value"], plt["value"])
        v = round(raw, 2)
        add("fib4", "fib4", v, "", inputs=[ast["marker"], alt["marker"], plt["marker"]], value_raw=raw)
        tier = fib4_tier(raw, int(age))
        if tier:
            add("fib4_tier", "fib4", tier, "", label_zh="FIB-4 界值提示（筛查用，非诊断）", inputs=[ast["marker"], alt["marker"], plt["marker"]])
        else:
            notes.append("FIB-4 未在 35 岁以下人群验证，不给分层")
    else:
        missing["fib4"] = [x for x, ok in (("谷草转氨酶", ast), ("谷丙转氨酶", alt), ("血小板", plt), ("年龄（成年人）", adult)) if not ok]
    tg, glu, hdl = lab("tg"), lab("glucose"), lab("hdl")
    if tg and glu:
        v = math.log((tg["value"] / 0.01129) * (glu["value"] * 18.0) / 2)
        add("tyg", "tyg", round(v, 2), "", kind="descriptive", inputs=[tg["marker"], glu["marker"]])
    else:
        missing["tyg"] = [x for x, ok in (("甘油三酯", tg), ("空腹血糖", glu)) if not ok]
    if tg and hdl:
        add("aip", "aip", round(math.log10(tg["value"] / hdl["value"]), 3), "", kind="descriptive", inputs=[tg["marker"], hdl["marker"]])
    else:
        missing["aip"] = [x for x, ok in (("甘油三酯", tg), ("高密度脂蛋白胆固醇", hdl)) if not ok]
    notes += [f"{a}：{p}" for a, p in probs.items()]
    notes.append("FIB-4 " + calc["fib4"]["note_zh"])
    return {"status": "ok", "readouts": ro, "notes": notes, "missing_inputs": {k: v for k, v in missing.items() if v}}


# ------------------------------------------------------------ AI estimate layer
def _organ_readouts(ws: Path) -> List[Dict[str, Any]]:
    p = ws / "work" / "organs" / "organ_readouts.json"
    return load_json(p)["readouts"] if p.exists() else []


def _unusable_basis(st: Dict[str, Any], ro: List[Dict[str, Any]]) -> List[str]:
    """Readout ids and lab markers an estimate may not rest on: low-quality or not confirmed as the member's."""
    bad_kinds = {"implausible", "computed_low_coverage", "computed_quality_unverified", "computed_coverage_unknown"}
    out = [r["id"] for r in ro if r.get("kind") in bad_kinds or r.get("provenance_uncertain")]
    from .labnames import refused_rows, usable_rows
    unc = {f["id"] for f in st["files"] if f.get("provenance_uncertain") or f.get("excluded")}
    out += [l["marker"] for l in st["labs"] if l.get("source_file") in unc]
    ok_markers = {l["marker"] for l in usable_rows(st)}         # a name shared with a confirmed row means that row
    out += [l["marker"] for l in refused_rows(st) if l["marker"] not in ok_markers]
    return sorted(set(out))


def bundle(st: Dict[str, Any], ws: Path) -> Dict[str, Any]:
    from .labnames import usable_rows
    spec = data("organs.json")
    ro = load_json(ws / "work" / "readouts.json")["readouts"]
    out = ws / "work" / "organs" / "bundles"
    old = {b["organ"]: b.get("bundle_sha256") for b in (st.get("organs") or {}).get("bundles") or []}
    made = []
    for organ, o in spec["organs"].items():
        measured = [r for r in ro if r.get("method") in o["age_readout_methods"] and r.get("unit") == "a"
                    and any(t in r["id"].lower() for t in o["age_readout_tokens"])]
        indices = [r for r in ro if r["id"] in o["index_readouts"]]
        b = {"organ": organ, "label_zh": o["label_zh"], "diseases": o["diseases"],
             "member": {k: st["member"].get(k) for k in ("id", "age", "sex", "sample_date", "answers")},
             "identity": st.get("identity"),
             "measured_organ_age": measured, "calibrated_or_published_indices": indices,
             "all_readouts": ro, "labs": [{k: v for k, v in l.items() if k not in ("confirm", "_row_key")} for l in usable_rows(st)],
             "rules": spec["ai_estimate_rules"], "contract": "references/organ-estimator.md",
             "unusable_basis": _unusable_basis(st, ro),
             "write_to": str(ws / "work" / "organs" / "estimates" / f"{organ}.json")}
        write_json(out / f"{organ}.json", b)
        made.append({"organ": organ, "measured_age": bool(measured), "indices": len(indices), "bundle": str(out / f"{organ}.json"),
                     "bundle_sha256": sha256_file(out / f"{organ}.json", limit=None)})
    st.setdefault("organs", {})["bundles"] = made
    for b in made:                      # a registration or skip made from a different bundle no longer stands
        if old.get(b["organ"]) != b["bundle_sha256"]:
            _void(st, ws, b["organ"])
            (st["organs"].get("skipped") or {}).pop(b["organ"], None)
    return {"bundles": made, "note": "give each estimator its bundle path and bundle_sha256; the estimate must carry that bundle_sha256"}


def _num(x: Any, name: str, problems: List[str]) -> Optional[float]:
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        problems.append(f"{name} must be a finite number")
        return None
    return float(x)


EST_KEYS = {"organ", "method", "confidence", "bundle_sha256", "age_estimate", "age_estimate_null_reason", "disease_risks"}
AGE_KEYS = {"low", "point", "high", "basis"}
RISK_KEYS = {"disease", "horizon_years", "low", "point", "high", "basis", "evidence", "rationale_zh"}


def active_overrides(ws: Path, organ: str, st: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Calibrated results that make an AI probability meaningless (a diagnostic threshold is already met)."""
    if st is None:
        st = load_json(ws / "state.json") if (ws / "state.json").exists() else {"labs": []}
    ro = {r["id"]: r for r in load_json(ws / "work" / "readouts.json")["readouts"]} if (ws / "work" / "readouts.json").exists() else {}
    seen, out = set(), []
    for o in data("organs.json").get("calibrated_overrides", []):
        if o["organ"] != organ or o["disease"] in seen:
            continue
        if "readout" in o:
            r = ro.get(o["readout"]) or {}
            vals = [r.get("value_raw", r.get("value"))]
            if o["readout"] == "native.organ.egfr" and vals == [None] and st.get("labs") and st.get("member", {}).get("age") \
                    and int(st["member"]["age"]) >= 18 and st["member"].get("sex") in ("female", "male"):
                vals += [egfr_ckd_epi_2021(c, int(st["member"]["age"]), st["member"]["sex"])   # every confirmed creatinine row
                         for c in lab_values(st, "creatinine", blank_unit="umol/L")]
        else:                                        # any confirmed value past the threshold counts (duplicates too)
            vals = lab_values(st, o["lab"], blank_unit=o.get("blank_unit")) if st.get("labs") else []
        for v in vals:
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                continue
            if ("below" in o and v < o["below"]) or ("at_least" in o and v >= o["at_least"] - 1e-9) or ("above" in o and v > o["above"]):
                out.append(o)
                seen.add(o["disease"])
                break
    return out


def unreadable_check_rows(st: Dict[str, Any], organ: str) -> List[str]:
    """Confirmed rows of an organ's threshold analytes whose value or unit the harness cannot read."""
    from .labnames import analytes, clean_value, usable_rows
    kit = skillkit()
    out = []
    spec = data("organs.json")
    for a in spec.get("check_labs", {}).get(organ, []):
        al = spec["lab_aliases"][a]
        codes = {str(c).casefold() for c in al.get("loinc", [a])}
        units = {kit.normalize_unit(u) for u in al["units"]}
        for r in usable_rows(st):
            if not (analytes(r) & codes):
                continue
            try:
                v = kit.parse_number(clean_value(r["value"]))
                lo, hi = al.get("plausible", [0, float("inf")])
                # a blank unit is readable only when the value is plausible in the canonical unit (151 glucose is mg/dL)
                ok = kit.normalize_unit(r.get("unit", "")) in units if r.get("unit") else lo <= v <= hi
            except ValueError:
                ok = False
            if not ok:
                out.append(r["marker"])
    return out


def check_labs(st: Dict[str, Any], organ: str) -> List[str]:
    """Confirmed lab values relevant to an organ, for the reviewer's side-by-side check."""
    out = []
    for a in data("organs.json").get("check_labs", {}).get(organ, []):
        v = _lab_value(st, a)[0]
        if v:
            out.append(f"{v['marker']}: {round(v['value'], 3)}")
    return out


def _void(st: Dict[str, Any], ws: Path, organ: str) -> None:
    o = st.setdefault("organs", {})
    if organ in (o.get("estimates") or {}):
        o["estimates"].pop(organ)
        _write_organ_readouts(st, ws, organ, [])


def register(st: Dict[str, Any], ws: Path, organ: str, pmid_check=None) -> Dict[str, Any]:
    try:
        return _register(st, ws, organ, pmid_check)
    except LAError:
        _void(st, ws, organ)            # a rejected re-registration never leaves the previous estimate in force
        raise


def _register(st: Dict[str, Any], ws: Path, organ: str, pmid_check=None) -> Dict[str, Any]:
    spec = data("organs.json")
    if organ not in spec["organs"]:
        raise LAError(f"unknown organ {organ!r}; one of {sorted(spec['organs'])}", EXIT_INPUT)
    bundles = {b["organ"]: b for b in (st.get("organs") or {}).get("bundles") or []}
    if organ not in bundles:
        raise LAError("run `la.py organ bundle` first; estimates are made from the current bundle", EXIT_INPUT)
    bpath = Path(bundles[organ]["bundle"])
    if not bpath.exists() or sha256_file(bpath, limit=None) != bundles[organ]["bundle_sha256"]:
        raise LAError(f"{bpath} changed after `organ bundle`; run `la.py organ bundle` again", EXIT_INPUT)
    rules = spec["ai_estimate_rules"]
    path = ws / "work" / "organs" / "estimates" / f"{organ}.json"
    if not path.exists():
        raise LAError(f"{path} not written", EXIT_INPUT)
    try:
        est = load_json(path)
    except ValueError as e:
        raise LAError(f"{path} is not JSON ({e})", EXIT_INPUT)
    if not isinstance(est, dict):
        raise LAError(f"{path} must be a JSON object", EXIT_INPUT)
    from .report import member_readouts
    readout_ids = {r["id"] for r in load_json(ws / "work" / "readouts.json")["readouts"]} | set(member_readouts(ws))
    lab_names = {l["marker"] for l in st["labs"]}
    unusable = set(load_json(bpath).get("unusable_basis") or [])
    age = st["member"].get("age")
    problems: List[str] = []
    extra = set(est) - EST_KEYS
    if extra:
        problems.append(f"unknown fields {sorted(extra)}; the contract allows {sorted(EST_KEYS)}")
    if est.get("bundle_sha256") != bundles[organ]["bundle_sha256"]:
        problems.append("`bundle_sha256` must equal the current bundle's (copy it from the bundle's `bundle_sha256` "
                        "in `organ bundle` output); an estimate made from an older bundle is not registered")
    if est.get("organ") != organ:
        problems.append("`organ` does not match")
    if est.get("method") != "llm_estimate":
        problems.append("`method` must be llm_estimate")
    if est.get("confidence") not in rules["confidence_allowed"]:
        problems.append(f"`confidence` must be one of {rules['confidence_allowed']}")

    def check_basis(basis: Any, where: str, need_member_data: bool = False) -> None:
        if not isinstance(basis, list) or not basis or not all(isinstance(b, str) for b in basis):
            problems.append(f"{where}: `basis` must be a list of readout ids or lab markers")
            return
        answers = {f"member.answers.{k}" for k in (st["member"].get("answers") or {})}
        bad = [b for b in basis if b not in readout_ids and b not in lab_names and b not in ("member.age", "member.sex") and b not in answers]
        if bad:
            problems.append(f"{where}: basis {bad} are not readouts or lab markers of this member")
        low = [b for b in basis if b in unusable]
        if low:
            problems.append(f"{where}: basis {low} are low-quality or not confirmed as the member's (bundle `unusable_basis`)")
        if need_member_data and all(b.startswith("member.") for b in basis):
            problems.append(f"{where}: an organ age needs organ-relevant member data in `basis`, not only age/sex; "
                            "with none, set age_estimate to null with a reason")

    pmids: List[str] = []
    out: List[Dict[str, Any]] = []
    ae = est.get("age_estimate")
    if ae is not None:
        if not isinstance(ae, dict) or set(ae) - AGE_KEYS:
            problems.append(f"age_estimate must be an object with {sorted(AGE_KEYS)}")
        elif not age:
            problems.append("the member's age is unknown: age_estimate must be null")
        else:
            lo, pt, hi = (_num(ae.get(k), f"age_estimate.{k}", problems) for k in ("low", "point", "high"))
            if None not in (lo, pt, hi):
                if not lo <= pt <= hi:
                    problems.append("age_estimate: need low ≤ point ≤ high")
                if hi - lo < rules["min_age_width_years"]:
                    problems.append(f"age_estimate: range narrower than {rules['min_age_width_years']} years claims precision an estimate does not have")
                if max(abs(lo - age), abs(hi - age)) > rules["max_age_offset_years"]:
                    problems.append(f"age_estimate: range reaches more than {rules['max_age_offset_years']} years from chronological age")
            check_basis(ae.get("basis"), "age_estimate", need_member_data=True)
            if not problems:
                out.append({"id": f"organ.{organ}.age_est", "label_zh": f"{spec['organs'][organ]['label_zh']}年龄（AI 估计）",
                            "value": round(pt, 1), "low": round(lo, 1), "high": round(hi, 1), "unit": "a", "kind": "llm_estimate",
                            "confidence": est.get("confidence"), "organ": organ, "basis": ae.get("basis"), "method": "llm_estimate",
                            "systems": [], "source": str(path)})
    elif not isinstance(est.get("age_estimate_null_reason"), str) or not est["age_estimate_null_reason"].strip():
        problems.append("age_estimate is null: give `age_estimate_null_reason`")
    risks = est.get("disease_risks")
    overridden = {o["disease"]: o for o in active_overrides(ws, organ, st)}
    unread = unreadable_check_rows(st, organ)
    if unread:
        problems.append(f"confirmed rows {unread} cannot be read (value or unit); a diagnostic threshold may be hidden in them. "
                        "Correct a mis-transcribed value or unit with `la.py labs fix`; a censored value (<5) must not be rewritten as a "
                        "number - skip this organ with `la.py organ skip` and say why")
    all_overridden = set(spec["organs"][organ]["diseases"]) <= set(overridden)
    if not isinstance(risks, list) or (not risks and not all_overridden):
        problems.append("`disease_risks` must be a non-empty list")
        risks = []
    seen: set = set()
    for i, r in enumerate(risks):
        w = f"disease_risks[{i}]"
        if not isinstance(r, dict):
            problems.append(f"{w} must be an object")
            continue
        if set(r) - RISK_KEYS:
            problems.append(f"{w}: unknown fields {sorted(set(r) - RISK_KEYS)}")
        dis = r.get("disease") if isinstance(r.get("disease"), str) else ""
        dis = dis.strip()
        if dis not in spec["organs"][organ]["diseases"]:
            problems.append(f"{w}: disease {dis!r} is not one of this organ's diseases {spec['organs'][organ]['diseases']}")
        elif dis in seen:
            problems.append(f"{w}: disease {dis!r} listed twice")
        elif dis in overridden:
            problems.append(f"{w}: {dis} is not estimated: {overridden[dis]['message_zh']} (leave it out; the table states this)")
        seen.add(dis)
        hz = r.get("horizon_years")
        if isinstance(hz, bool) or not isinstance(hz, int) or not rules["horizon_years"][0] <= hz <= rules["horizon_years"][1]:
            problems.append(f"{w}: horizon_years must be a whole number in {rules['horizon_years']}")
        lo, pt, hi = (_num(r.get(k), f"{w}.{k}", problems) for k in ("low", "point", "high"))
        if None not in (lo, pt, hi):
            if not 0 <= lo <= pt <= hi <= 1:
                problems.append(f"{w}: need 0 ≤ low ≤ point ≤ high ≤ 1 (probabilities, not percent)")
            if hi - lo < max(rules["min_prob_width_abs"], rules["min_prob_width_rel"] * pt):
                problems.append(f"{w}: range too narrow for an estimate (≥ {rules['min_prob_width_abs']} absolute and ≥ half the point value)")
        check_basis(r.get("basis"), w)
        if not isinstance(r.get("rationale_zh"), str):
            problems.append(f"{w}: `rationale_zh` must be a string")
        else:                                   # the same number rules the trace applies, caught here instead of later
            from .report import PLACEHOLDER, trace_text
            bad_ids = [m.group(1) for m in PLACEHOLDER.finditer(r["rationale_zh"]) if m.group(1) not in readout_ids]
            if bad_ids:
                problems.append(f"{w}: rationale_zh cites {bad_ids}, which are not readout ids ({{{{r:…}}}} takes an id, never a value)")
            if age and re.search(rf"\{{\{{n:[^{{}}]*(?<!\d){age}(?!\d)", r["rationale_zh"]):
                problems.append(f"{w}: the member's own age is written as a literal; use {{{{r:…}}}} or leave it out")
            tt = trace_text(r["rationale_zh"])
            if tt:
                problems.append(f"{w}: rationale_zh breaks the number rules: " + "; ".join(tt[:3]))
        ev = r.get("evidence")
        if not isinstance(ev, list):
            problems.append(f"{w}: `evidence` must be a list of {{\"type\": \"pubmed\", \"ref\": \"<PMID>\"}}")
            ev = []
        for j, e in enumerate(ev):
            if not (isinstance(e, dict) and set(e) == {"type", "ref"} and e["type"] == "pubmed"
                    and isinstance(e["ref"], (str, int)) and not isinstance(e["ref"], bool) and re.fullmatch(r"[1-9][0-9]{0,8}", str(e["ref"]))):
                problems.append(f"{w}.evidence[{j}]: only {{\"type\": \"pubmed\", \"ref\": \"<PMID digits>\"}} is accepted")
            else:
                pmids.append(str(e["ref"]))
        if not problems:
            out.append({"id": f"organ.{organ}.risk.{i + 1}", "label_zh": f"{dis}（{hz} 年，AI 估计）", "disease": dis,
                        "value": round(pt, 4), "low": lo, "high": hi, "unit": "概率", "horizon_years": hz,
                        "kind": "llm_estimate", "confidence": est.get("confidence"), "organ": organ, "basis": r.get("basis"),
                        "rationale_zh": r["rationale_zh"], "evidence": [{"type": "pubmed", "ref": str(e["ref"])} for e in ev],
                        "method": "llm_estimate", "systems": [], "source": str(path)})
    uniq = sorted(set(pmids))
    if risks and len(uniq) < rules["min_pubmed_per_organ"]:
        problems.append(f"cite at least {rules['min_pubmed_per_organ']} PubMed article retrieved in this run (baseline incidence, cohort data)")
    if uniq and pmid_check is not None:
        live = pmid_check(uniq)
        dead = [p for p in uniq if p not in live]
        if dead:
            problems.append(f"PMIDs {dead} do not exist or are retracted")
        from .evidence import _cached_pmids
        unretrieved = [p for p in uniq if p not in _cached_pmids(ws)]
        if unretrieved:
            problems.append(f"PMIDs {unretrieved} were not retrieved with `la.py evidence pubmed` in this run")
    if problems:
        raise LAError(f"{organ} estimate rejected:\n- " + "\n- ".join(problems), EXIT_INPUT)
    reg = st.setdefault("organs", {}).setdefault("estimates", {})
    reg[organ] = {"path": str(path), "sha256": sha256_file(path, limit=None), "readouts": [r["id"] for r in out], "at": now_iso()}
    (st["organs"].get("skipped") or {}).pop(organ, None)
    _write_organ_readouts(st, ws, organ, out)
    return reg[organ]


def _write_organ_readouts(st: Dict[str, Any], ws: Path, organ: str, new: List[Dict[str, Any]]) -> None:
    """Replace one organ's rows; rows of organs no longer registered are dropped."""
    reg = (st.get("organs") or {}).get("estimates") or {}
    keep = [r for r in _organ_readouts(ws) if r.get("organ") != organ and r.get("organ") in reg]
    dst = ws / "work" / "organs" / "organ_readouts.json"
    write_json(dst, {"generated_at": now_iso(), "readouts": keep + new})
    st["organs"]["organ_readouts_sha256"] = sha256_file(dst, limit=None)


def skip(st: Dict[str, Any], ws: Path, organ: str, reason: str) -> Dict[str, Any]:
    if organ not in data("organs.json")["organs"]:
        raise LAError(f"unknown organ {organ!r}", EXIT_INPUT)
    o = st.setdefault("organs", {})
    o.setdefault("skipped", {})[organ] = {"reason": reason, "at": now_iso()}
    (o.get("estimates") or {}).pop(organ, None)
    _write_organ_readouts(st, ws, organ, [])
    return {"skipped": organ}


def table(ws: Path, readouts: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One row per organ for the report and the twin."""
    spec = data("organs.json")
    rows = []
    for organ, o in spec["organs"].items():
        measured = [r for r in readouts.values() if r.get("method") in o["age_readout_methods"] and r.get("unit") == "a"
                    and any(t in r["id"].lower() for t in o["age_readout_tokens"]) and r.get("kind") not in ("llm_estimate", "implausible")]
        idx = [readouts[i] for i in o["index_readouts"] if i in readouts]
        est_age = readouts.get(f"organ.{organ}.age_est")
        est_risks = [r for r in readouts.values() if r.get("organ") == organ and r["id"].startswith(f"organ.{organ}.risk.")]
        rows.append({"organ": organ, "label_zh": o["label_zh"], "measured": measured, "indices": idx,
                     "ai_age": est_age, "ai_risks": est_risks, "overrides": active_overrides(ws, organ)})
    return rows
