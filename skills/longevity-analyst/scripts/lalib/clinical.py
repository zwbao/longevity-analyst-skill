"""Diagnostic conclusions and medication suggestions (PLAN WS-L.3, WS-L.4), each a labelled item.

`clinical dx` evaluates data/diagnostic_criteria.json (la-dx-criteria/1) for one guideline set (cn by default, or
ada) against the member's confirmed lab rows (all dates: the current snapshot and work/twin/series.json) and method
readouts (KDIGO GFR category, FIB-4). A hit is written in the fixed wording
"按《…》的标准，符合……的诊断标准。" with the guideline quote as its source. AI judgements from the clinical-estimator
subagent (references/clinical-estimator.md) are registered with --ai and kept in their own column:
"AI 判断：可能存在……（概率区间……，依据……）。"

`clinical rx` evaluates data/medication_rules.json (la-rx-rules/1), turns every CPIC recommendation in
work/pgx/pgx.json into a `cpic` item with CPIC's text verbatim (dose-adjustment text included), and registers the
estimator's AI suggestions (always labelled AI 判断 with a confidence). Thresholds are compared in the guideline's
own unit; the member's value is converted here.
"""
from __future__ import annotations

import datetime as _dt
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import labels as L
from .common import DATA_DIR, EXIT_INPUT, EXIT_USAGE, LAError, data, load_json, now_iso, sha256_file, skillkit, write_json

ACTIONS = ("consider_start", "consider_stop", "consider_switch", "adjust", "avoid", "monitor")
ACTION_ZH = {"consider_start": "考虑开始", "consider_stop": "考虑停用", "consider_switch": "考虑换药", "adjust": "调整剂量",
             "avoid": "避免使用", "monitor": "监测"}
DX_WORDING = "按《{g}》的标准，符合{c}的诊断标准。"
AI_DX_WORDING = "AI 判断：可能存在{c}（概率区间 {lo}–{hi}，依据{b}）。"
EXECUTOR = ["physician", "member"]
DRUG_ZH = {"simvastatin": "辛伐他汀", "atorvastatin": "阿托伐他汀", "rosuvastatin": "瑞舒伐他汀", "pravastatin": "普伐他汀",
           "fluvastatin": "氟伐他汀", "lovastatin": "洛伐他汀", "pitavastatin": "匹伐他汀", "clopidogrel": "氯吡格雷",
           "warfarin": "华法林", "metformin": "二甲双胍"}       # names only, to match what the member wrote


# ------------------------------------------------------------ member values with dates
def _spec(analyte: str, extra: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    al = data("organs.json")["lab_aliases"]
    return al.get(analyte) or extra.get(analyte)


def _convert(v: float, f: Any) -> float:
    return v * f[0] + f[1] if isinstance(f, list) else v * f


def _row_is(r_marker: str, codes: set, spec: Dict[str, Any]) -> bool:
    from .labnames import analytes
    kit = skillkit()
    lo = {str(c).casefold() for c in spec.get("loinc", [])}
    if lo and analytes({"marker": r_marker}) & lo:
        return True
    names = {kit.fold_name(n) for n in spec.get("names", [])}
    return bool(names & set(kit.name_variants(r_marker)))


def marker_values(st: Dict[str, Any], ws: Path, analyte: str, extra: Dict[str, Any]) -> List[Dict[str, Any]]:
    """[{t, value (canonical unit), marker, from}] sorted by date, one per date; rows whose unit cannot be read are left out."""
    from .labnames import clean_value, usable_rows
    spec = _spec(analyte, extra)
    if not spec:
        raise LAError(f"marker {analyte!r} has no unit table (organs.json lab_aliases or the rule file's markers)", EXIT_INPUT)
    if spec.get("derived_from"):
        return _derived(st, ws, analyte, spec, extra)
    kit = skillkit()
    units = {kit.normalize_unit(u): f for u, f in spec["units"].items()}
    lo, hi = spec.get("plausible", [0, float("inf")])
    snap = st["member"].get("sample_date") or now_iso()[:10]
    by_t: Dict[str, Dict[str, Any]] = {}

    def add(t: str, raw: Any, unit: str, marker: str, src: str) -> None:
        try:
            v = kit.parse_number(clean_value(str(raw)))
        except ValueError:
            return
        u = kit.normalize_unit(unit or "")
        if u not in units:
            return
        v = _convert(v, units[u])
        if math.isfinite(v) and lo <= v <= hi and t not in by_t:
            by_t[t] = {"t": t, "value": v, "marker": marker, "unit_in": unit, "from": src}
    for r in usable_rows(st):
        if _row_is(r["marker"], set(), spec):
            add(str(r.get("collected_at") or snap)[:10], r["value"], r.get("unit", ""), r["marker"], "lab")
    from .predict import _series
    for s in _series(ws / "work" / "twin" / "series.json"):
        if _row_is(str(s.get("marker", "")), set(), spec):
            for p in s.get("points", []):
                add(str(p.get("t", ""))[:10], p.get("value"), s.get("unit", ""), s["marker"], "series")
    return [by_t[t] for t in sorted(by_t)]


def _derived(st: Dict[str, Any], ws: Path, analyte: str, spec: Dict[str, Any], extra: Dict[str, Any]) -> List[Dict[str, Any]]:
    """eGFR at every date a creatinine was measured (CKD-EPI 2021, the formula the organ table uses); age on that date."""
    if analyte != "egfr" or spec["derived_from"] != "creatinine":
        raise LAError(f"no derivation for {analyte}", EXIT_INPUT)
    from .organs import egfr_ckd_epi_2021
    m = st["member"]
    if not m.get("age") or m.get("sex") not in ("male", "female"):
        return []
    snap = _dt.date.fromisoformat((m.get("sample_date") or now_iso())[:10])
    out = []
    for x in marker_values(st, ws, "creatinine", extra):
        age = int(m["age"]) - round((snap - _dt.date.fromisoformat(x["t"])).days / 365.25)
        if age >= 18:
            out.append(dict(x, value=egfr_ckd_epi_2021(x["value"], age, m["sex"]), marker=f"eGFR（{x['marker']}）"))
    return out


def _readout_values(ws: Path) -> Dict[str, Dict[str, Any]]:
    out = {}
    for p in (ws / "work" / "readouts.json", ws / "work" / "insights" / "insight_readouts.json"):
        if p.exists():
            out.update({r["id"]: r for r in load_json(p).get("readouts", [])})
    return out


def _test(v: Any, c: Dict[str, Any]) -> bool:
    if "in" in c:
        return str(v) in [str(x) for x in c["in"]]
    if "eq" in c:
        return str(v) == str(c["eq"])
    if not isinstance(v, (int, float)):
        return False
    ok = True
    for op, f in (("gte", lambda a, b: a >= b), ("gt", lambda a, b: a > b), ("lte", lambda a, b: a <= b), ("lt", lambda a, b: a < b)):
        if op in c:
            ok = ok and f(v, c[op])
    return ok


def _thr(c: Dict[str, Any], spec: Dict[str, Any]) -> Dict[str, Any]:
    """A condition with its thresholds converted from the guideline's unit into the canonical unit."""
    kit = skillkit()
    units = {kit.normalize_unit(u): f for u, f in spec["units"].items()}
    f = units.get(kit.normalize_unit(c.get("unit", "")))
    if f is None:
        raise LAError(f"condition unit {c.get('unit')!r} is not convertible for {c.get('marker')}", EXIT_INPUT)
    return {**c, **{op: _convert(float(c[op]), f) for op in ("gte", "gt", "lte", "lt") if op in c}}


def evaluate(cond: Dict[str, Any], ctx: Dict[str, Any]) -> Tuple[Optional[bool], List[Dict[str, Any]]]:
    """(True | False | None when the data are not there, the values used)."""
    if "all_of" in cond or "any_of" in cond:
        parts = [evaluate(c, ctx) for c in cond.get("all_of") or cond.get("any_of")]
        ev = [e for _, es in parts for e in es]
        res = [r for r, _ in parts]
        if "all_of" in cond:
            return (False if False in res else None if None in res else True), ev
        return (True if True in res else None if None in res else False), ev
    if "dx" in cond:
        hit = cond["dx"] in ctx.get("dx_hits", set())
        return hit, ([{"type": "rule", "ref": cond["dx"]}] if hit else [])
    if "member" in cond:
        v = (ctx["st"]["member"] or {}).get(cond["member"])
        if v is None:
            return None, []
        return _test(v, cond), [{"type": "member", "ref": cond["member"], "value": v}]
    if "readout" in cond:
        r = ctx["readouts"].get(cond["readout"])
        if r is None or r.get("kind") in ("implausible",) or r.get("provenance_uncertain"):
            return None, []
        return _test(r.get("value"), cond), [{"type": "readout", "ref": r["id"], "value": r.get("value"), "unit": r.get("unit", "")}]
    if "marker" in cond:
        spec = _spec(cond["marker"], ctx["extra"])
        vals = ctx["values"].setdefault(cond["marker"], marker_values(ctx["st"], ctx["ws"], cond["marker"], ctx["extra"]))
        if not vals:
            return None, []
        c = _thr(cond, spec)
        canon = next(iter(spec["units"]))
        latest = vals[-1]
        ev = [{"type": "lab", "ref": latest["marker"], "value": round(latest["value"], 3), "unit": canon, "t": latest["t"]}]
        if not _test(latest["value"], c):
            return False, ev
        need = int(cond.get("distinct_days", 1))
        hits = [x for x in vals if _test(x["value"], c)]
        if len(hits) < need:
            return False, ev
        if cond.get("min_span_days"):
            span = (_dt.date.fromisoformat(hits[-1]["t"]) - _dt.date.fromisoformat(hits[0]["t"])).days
            if span < int(cond["min_span_days"]):
                return False, ev
        return True, [{"type": "lab", "ref": x["marker"], "value": round(x["value"], 3), "unit": canon, "t": x["t"]} for x in hits[-need:]]
    raise LAError(f"unknown condition {cond}", EXIT_INPUT)


def _ctx(st: Dict[str, Any], ws: Path, extra: Dict[str, Any]) -> Dict[str, Any]:
    return {"st": st, "ws": ws, "extra": extra, "values": {}, "readouts": _readout_values(ws)}


def _src_label(rule: Dict[str, Any], basis: str, evidence: List[Dict[str, Any]]) -> Dict[str, Any]:
    s = rule["source"]
    sources = [{"type": "guideline", "ref": s["title"], "quote": s["quote"],
                "url": f"https://doi.org/{s['doi']}" if s.get("doi") else s["url"]},
               {"type": "rule", "ref": rule["id"]}]
    sources += [{"type": e["type"], "ref": str(e["ref"])} for e in evidence if e["type"] in ("lab", "readout")]
    return L.make_label(basis, sources, guideline_zh=rule["guideline_zh"])


def _sources_ok(d: Dict[str, Any], kind: str) -> List[str]:
    bad = []
    for r in d.get("rules", []):
        s = r.get("source") or {}
        if not (s.get("title") and s.get("quote") and (s.get("doi") or s.get("url"))):
            bad.append(f"{r.get('id')}: source needs title, quote and doi or url")
        if not r.get("guideline_zh"):
            bad.append(f"{r.get('id')}: guideline_zh missing")
        if kind == "rx" and r.get("action") not in ACTIONS:
            bad.append(f"{r.get('id')}: action {r.get('action')!r}")
    return bad


# ------------------------------------------------------------ clinical dx
def inputs_digest(st: Dict[str, Any], ws: Path, kind: str) -> str:
    from .labnames import usable_rows
    sha = lambda p: sha256_file(p, limit=None) if p.exists() else None   # noqa: E731
    parts = {"labs": L.digest([[r["marker"], r["value"], r.get("unit"), r.get("collected_at")] for r in usable_rows(st)]),
             "member": L.digest({k: st["member"].get(k) for k in ("age", "sex", "sample_date", "answers")}),
             "readouts": sha(ws / "work" / "readouts.json"), "insight": sha(ws / "work" / "insights" / "insight_readouts.json"),
             "series": sha(ws / "work" / "twin" / "series.json")}
    if kind == "dx":
        parts.update(rules=sha(DATA_DIR / "diagnostic_criteria.json"), ai=sha(ws / "work" / "clinical" / "ai_dx.json"),
                     set=((st.get("clinical") or {}).get("dx") or {}).get("guideline_set"))
    else:
        parts.update(rules=sha(DATA_DIR / "medication_rules.json"), ai=sha(ws / "work" / "clinical" / "ai_rx.json"),
                     dx=sha(L.item_path(ws, "dx")), pgx=sha(ws / "work" / "pgx" / "pgx.json"))
    return L.digest(parts)


def _dx_caveats(rule: Dict[str, Any], ev: List[Dict[str, Any]], d: Dict[str, Any], ctx: Dict[str, Any]) -> List[str]:
    """What the reader must know next to a met criterion: it rests on one measurement; the same condition in another
    guideline set of the file is (not) met (data-driven: the rules with the same label in other sets)."""
    out = []
    labs = [e for e in ev if e.get("type") == "lab"]
    if len(labs) == 1:
        out.append(f"这是 {labs[0].get('t')} 一次测量得出的结论，需要复查确认。")
    for other in d["rules"]:
        if other is rule or other["label_zh"] != rule["label_zh"] or other["set"] in (rule["set"], "any"):
            continue
        res, _ = evaluate(other["when"], ctx)
        if res is False:
            out.append(f"按《{other['guideline_zh']}》的界值，这次的结果不符合{other['label_zh']}。")
    return out


def dx(st: Dict[str, Any], ws: Path, guideline_set: str = "cn", ai_file: Optional[Path] = None) -> Dict[str, Any]:
    d = data("diagnostic_criteria.json")
    if guideline_set not in d["guideline_sets"]:
        raise LAError(f"--guideline-set must be one of {d['guideline_sets']}", EXIT_USAGE)
    bad = _sources_ok(d, "dx")
    if bad:
        raise LAError("diagnostic_criteria.json: " + "; ".join(bad), EXIT_INPUT)
    ctx = _ctx(st, ws, d.get("markers", {}))
    items, evaluated = [], []
    for rule in d["rules"]:
        if rule["set"] not in (guideline_set, "any"):
            continue
        res, ev = evaluate(rule["when"], ctx)
        evaluated.append({"rule_id": rule["id"], "label_zh": rule["label_zh"], "guideline_zh": rule["guideline_zh"],
                          "result": {True: "met", False: "not_met", None: "no_data"}[res], "values": ev})
        if res:
            items.append({"id": f"DX{len(items) + 1}", "kind": "guideline_criteria", "rule_id": rule["id"], "label_zh": rule["label_zh"],
                          "conclusion_zh": DX_WORDING.format(g=rule["guideline_zh"], c=rule["label_zh"]),
                          "note_zh": rule.get("note_zh"), "caveats_zh": _dx_caveats(rule, ev, d, ctx), "evidence": ev,
                          "executor": ["physician"], "label": _src_label(rule, "guideline_criteria", ev)})
    if ai_file:
        items += _ai_dx(st, ws, ai_file, len(items))
    elif (ws / "work" / "clinical" / "ai_dx.json").exists():
        (ws / "work" / "clinical" / "ai_dx.json").unlink()
    used = sorted({r["guideline_zh"] for r in d["rules"] if r["set"] in (guideline_set, "any")})
    out = {"schema": "la-diagnoses/1", "guideline_set": guideline_set, "guidelines_used_zh": used,
           "generated_at": now_iso(), "items": items, "evaluated": evaluated}
    write_json(L.item_path(ws, "dx"), out)
    st.setdefault("clinical", {})["dx"] = {"guideline_set": guideline_set}
    L.record(st, ws, "dx", inputs_digest(st, ws, "dx"))
    return {"path": str(L.item_path(ws, "dx")), "guideline_set": guideline_set, "guidelines_used_zh": used,
            "met": [x["rule_id"] for x in items if x["kind"] == "guideline_criteria"],
            "not_met": [e["rule_id"] for e in evaluated if e["result"] == "not_met"],
            "no_data": [e["rule_id"] for e in evaluated if e["result"] == "no_data"],
            "ai": sum(1 for x in items if x["kind"] == "ai_estimate"),
            "next": "dispatch references/clinical-estimator.md for the AI column, then re-run with --ai <file>" if not ai_file else None}


def _pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def _ai_common(ws: Path, ai_file: Path) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    from .evidence import _cached_pmids, verify_pmids
    doc = load_json(ai_file)
    rows = doc.get("items") if isinstance(doc, dict) else None
    if not isinstance(rows, list):
        raise LAError("the estimator's file needs {\"items\": [...]} (references/clinical-estimator.md)", EXIT_INPUT)
    pm = sorted({str(p) for r in rows if isinstance(r, dict) for p in r.get("pmids", [])})
    return rows, {"cached": _cached_pmids(ws), "live": verify_pmids(pm) if pm else {}, "readouts": set(_readout_values(ws))}


def _ai_row_problems(r: Any, n: int, ctx: Dict[str, Any], st: Dict[str, Any], prose: Tuple[str, ...]) -> List[str]:
    from .labnames import usable_rows
    from .report import trace_text
    if not isinstance(r, dict):
        return [f"item {n}: not an object"]
    bad = []
    if r.get("confidence") not in L.AI_LEVELS:
        bad.append(f"item {n}: confidence must be low or very_low")
    pm = [str(p) for p in r.get("pmids", [])]
    basis = r.get("basis") or []
    if not pm and not basis:
        bad.append(f"item {n}: give its basis (readout ids / lab markers) and PMIDs retrieved in this run")
    for p in pm:
        if p not in ctx["cached"] or p not in ctx["live"]:
            bad.append(f"item {n}: PMID {p} was not retrieved in this run or PubMed does not confirm it")
    labs = {x["marker"] for x in usable_rows(st)}
    for b in basis:
        if b not in ctx["readouts"] and b not in labs and not str(b).startswith(("member.", "rule:", "DX")):
            bad.append(f"item {n}: basis {b!r} is not a readout id or a confirmed lab marker of this member")
    for k in prose:
        if not str(r.get(k) or "").strip():
            bad.append(f"item {n}: missing {k}")
    prose_text = "\n".join(str(r.get(k) or "") for k in prose)
    tb = trace_text(prose_text) + (["{{…}}"] if "{{" in prose_text else [])
    if tb:
        bad.append(f"item {n}: {'/'.join(prose)} carry numbers, doses or links ({tb[0]}); put a dose in dose_text_zh")
    return bad


def _ai_sources(r: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [{"type": "pmid", "ref": str(p)} for p in r.get("pmids", [])] + \
        [{"type": "readout" if "." in str(b) else "lab", "ref": str(b)} for b in r.get("basis", [])]


def _ai_dx(st: Dict[str, Any], ws: Path, ai_file: Path, start: int) -> List[Dict[str, Any]]:
    import json as _j
    rows, ctx = _ai_common(ws, ai_file)
    out, bad = [], []
    for n, r in enumerate(rows, 1):
        b = _ai_row_problems(r, n, ctx, st, ("condition_zh", "basis_zh"))
        pr = r.get("probability_range") if isinstance(r, dict) else None
        if not (isinstance(pr, list) and len(pr) == 2 and all(isinstance(x, (int, float)) for x in pr) and 0 <= pr[0] <= pr[1] <= 1):
            b.append(f"item {n}: probability_range must be [low, high] between 0 and 1")
        if b:
            bad += b
            continue
        lab = L.make_label("ai_estimate", _ai_sources(r), confidence={"level": r["confidence"], "range": pr})
        out.append({"id": f"DX{start + len(out) + 1}", "kind": "ai_estimate", "label_zh": r["condition_zh"],
                    "condition_zh": r["condition_zh"], "basis_zh": r["basis_zh"], "probability_range": pr,
                    "conclusion_zh": AI_DX_WORDING.format(c=r["condition_zh"], lo=_pct(pr[0]), hi=_pct(pr[1]), b=r["basis_zh"]),
                    "executor": ["physician"], "label": lab})
    if bad:
        raise LAError("AI diagnostic judgements rejected:\n- " + "\n- ".join(bad), EXIT_INPUT)
    write_json(ws / "work" / "clinical" / "ai_dx.json", _j.loads(Path(ai_file).read_text(encoding="utf-8")))
    return out


# ------------------------------------------------------------ clinical rx
def current_medications(st: Dict[str, Any]) -> List[str]:
    a = st["member"].get("answers") or {}
    v = a.get("medications") or a.get("current_medications") or ""      # `la.py member <ws> medications=…` (the key pgx reads)
    return [x.strip() for x in re.split(r"[,，、;；/]+", str(v)) if x.strip() and x.strip() not in ("无", "none", "no")]


def _taking(meds: List[str], *names: str) -> bool:
    f = lambda s: re.sub(r"\s+", "", str(s)).casefold()   # noqa: E731
    return any(n and (f(n) in f(m) or f(m) in f(n)) for m in meds for n in names)


def _cpic_action(rec: str, taking: bool, category: Optional[str] = None) -> str:
    """The action CPIC's own words call for; the words themselves are shown verbatim. The skill's category of the
    row (standard / caution / avoid) decides when it is there; the text only refines caution into a dose change."""
    t = rec.casefold().strip()
    dose = re.search(r"[≤<]\s*\d|\bdose\b|\d+\s*%|\btimes\b", t)
    if category == "avoid":
        return "consider_switch" if taking else "avoid"
    if category == "standard":
        return "monitor"
    if category == "caution":
        return "adjust" if dose and "desired starting dose" not in t else "monitor"
    if re.search(r"\bavoid\b|not recommended|contraindicated|prescribe an alternative|alternative (?:statin|drug|agent)", t):
        return "consider_switch" if taking else "avoid"
    if t.startswith(("prescribe desired starting dose", "initiate therapy with recommended starting dose")):
        return "monitor"
    if re.search(r"[≤<]\s*\d+\s*mg|\b(?:reduce|lower|limit|decrease)\w*\b.{0,30}\bdose", t):
        return "adjust"
    return "monitor"


def _cpic_items(st: Dict[str, Any], ws: Path, start: int) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    p = ws / "work" / "pgx" / "pgx.json"
    if not p.exists():
        return [], "work/pgx/pgx.json not found (pharmacogenomics not run)"
    d = load_json(p)
    if d.get("schema") != "la-pgx/1":
        raise LAError("work/pgx/pgx.json is not la-pgx/1", EXIT_INPUT)
    if d.get("status") != "ran":
        return [], d.get("reason_zh") or "pharmacogenomics not run"
    meds = current_medications(st)
    out = []
    for res in d.get("results", []):
        for x in res.get("drugs", []):
            rec = str(x.get("recommendation_text") or "").strip()
            if not rec:
                continue
            dz = x.get("drug_zh") or DRUG_ZH.get(str(x.get("drug", "")).casefold())
            taking = _taking(meds, x.get("drug"), dz)
            src = {"type": "cpic", "ref": x.get("guideline") or x.get("guideline_url") or "CPIC", "quote": rec}
            if x.get("guideline_url"):
                src["url"] = x["guideline_url"]
            lk = x.get("lookup") or {}
            basis = "；".join(f"{g} {v}" for g, v in lk.items()) if len(lk) > 1 else f"{res.get('gene')} {res.get('phenotype')}"
            lab = L.make_label("cpic", [src, {"type": "lab", "ref": f"基因检测 {basis}"}])
            gene, pheno = (("、".join(lk), "；".join(f"{g} {v}" for g, v in lk.items())) if len(lk) > 1
                           else (res.get("gene"), res.get("phenotype")))
            out.append({"id": f"RX{start + len(out) + 1}", "kind": "cpic", "gene": gene, "phenotype": pheno,
                        "drug_or_class_zh": dz or x.get("drug"), "drug": x.get("drug"),
                        "action": _cpic_action(rec, taking, x.get("category")), "cpic_category": x.get("category"),
                        "currently_taking": taking,
                        **({"caveats_zh": [f"{'、'.join(x['untested_genes'])} 没有检测结果，这条推荐没有考虑它；"
                                           f"用这个药之前，请医生确认{'、'.join(x['untested_genes'])}。"]} if x.get("untested_genes") else {}),
                        "label_zh": f"{gene} 与 {dz or x.get('drug')}",
                        "text_zh": rec, **({"dose_text_zh": x["dose_text"]} if x.get("dose_text") else {}),
                        "strength": x.get("strength"), "executor": list(EXECUTOR), "label": lab})
    return out, None


def rx(st: Dict[str, Any], ws: Path, ai_file: Optional[Path] = None) -> Dict[str, Any]:
    d = data("medication_rules.json")
    bad = _sources_ok(d, "rx")
    if bad:
        raise LAError("medication_rules.json: " + "; ".join(bad), EXIT_INPUT)
    dxp = L.item_path(ws, "dx")
    dxd = load_json(dxp) if dxp.exists() else {"items": [], "guideline_set": "cn"}
    gset = dxd.get("guideline_set", "cn")
    ctx = _ctx(st, ws, d.get("markers", {}))
    ctx["dx_hits"] = {x["rule_id"] for x in dxd.get("items", []) if x.get("kind") == "guideline_criteria"}
    meds = current_medications(st)
    items, evaluated = [], []
    for rule in d["rules"]:
        if rule["set"] not in (gset, "any"):
            continue
        res, ev = evaluate(rule["when"], ctx)
        evaluated.append({"rule_id": rule["id"], "result": {True: "met", False: "not_met", None: "no_data"}[res], "values": ev})
        if res:
            taking = _taking(meds, *rule.get("drug_names_zh", []))
            items.append({"id": f"RX{len(items) + 1}", "kind": "guideline_recommendation", "rule_id": rule["id"],
                          "drug_or_class_zh": rule["drug_class_zh"], "action": rule["action"], "currently_taking": taking,
                          "label_zh": f"{rule['drug_class_zh']}：{ACTION_ZH[rule['action']]}", "text_zh": rule["text_zh"],
                          **({"dose_text_zh": rule["dose_text_zh"]} if rule.get("dose_text_zh") else {}),
                          "note_zh": rule.get("note_zh"), "evidence": ev, "executor": list(rule.get("executor") or EXECUTOR),
                          "label": _src_label(rule, "guideline_recommendation", ev)})
    cp, why = _cpic_items(st, ws, len(items))
    items += cp
    if ai_file:
        items += _ai_rx(st, ws, ai_file, len(items))
    elif (ws / "work" / "clinical" / "ai_rx.json").exists():
        (ws / "work" / "clinical" / "ai_rx.json").unlink()
    out = {"schema": "la-medications/1", "guideline_set": gset, "generated_at": now_iso(), "current_medications": meds,
           "items": items, "evaluated": evaluated, **({"cpic_not_used": why} if why else {})}
    write_json(L.item_path(ws, "rx"), out)
    L.record(st, ws, "rx", inputs_digest(st, ws, "rx"))
    return {"path": str(L.item_path(ws, "rx")), "guideline_set": gset,
            "guideline": [x["rule_id"] for x in items if x["kind"] == "guideline_recommendation"],
            "cpic": [x["drug"] for x in items if x["kind"] == "cpic"], "cpic_not_used": why,
            "ai": sum(1 for x in items if x["kind"] == "ai_estimate"),
            "next": "dispatch references/clinical-estimator.md (medications part), then re-run with --ai <file>" if not ai_file else None}


def _ai_rx(st: Dict[str, Any], ws: Path, ai_file: Path, start: int) -> List[Dict[str, Any]]:
    import json as _j
    rows, ctx = _ai_common(ws, ai_file)
    meds = current_medications(st)
    out, bad = [], []
    for n, r in enumerate(rows, 1):
        b = _ai_row_problems(r, n, ctx, st, ("drug_or_class_zh", "text_zh", "basis_zh"))
        if isinstance(r, dict) and r.get("action") not in ACTIONS:
            b.append(f"item {n}: action must be one of {ACTIONS}")
        if isinstance(r, dict) and r.get("dose_text_zh") is not None and not str(r["dose_text_zh"]).strip():
            b.append(f"item {n}: dose_text_zh is empty")
        if b:
            bad += b
            continue
        lab = L.make_label("ai_estimate", _ai_sources(r), confidence={"level": r["confidence"]})
        out.append({"id": f"RX{start + len(out) + 1}", "kind": "ai_estimate", "drug_or_class_zh": r["drug_or_class_zh"],
                    "action": r["action"], "currently_taking": _taking(meds, r["drug_or_class_zh"]),
                    "label_zh": f"{r['drug_or_class_zh']}：{ACTION_ZH[r['action']]}", "text_zh": r["text_zh"], "basis_zh": r["basis_zh"],
                    **({"dose_text_zh": str(r["dose_text_zh"])} if r.get("dose_text_zh") else {}),
                    "executor": list(EXECUTOR), "label": lab})
    if bad:
        raise LAError("AI medication suggestions rejected:\n- " + "\n- ".join(bad), EXIT_INPUT)
    write_json(ws / "work" / "clinical" / "ai_rx.json", _j.loads(Path(ai_file).read_text(encoding="utf-8")))
    return out
