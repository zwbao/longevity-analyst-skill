"""Personal response prediction (PLAN WS-L.2, appendix B.7): will this plan item work for this member?

Statistical model (label `statistical_model`): prior from the trial mean effect in the library's effects table
(scaled by adherence) and the person-to-person spread (tau_pct, or |mu0| by default); updated with the member's own
earlier experiments of the same kind (work/twin/series.json); noise from the library's biological variation table.
Three outputs: P(the next experiment shows a favourable change beyond noise), P(the true effect is favourable),
and the 80% interval of the expected change. Every assumption is written into the label.

Items the model cannot cover (no effects row, an effect not in the marker's units, no baseline) are listed in
work/predict/ai_requests.json for the response-estimator subagent (references/response-estimator.md); its answers are
registered with `predict response <ws> --ai <file>`: intervals only, confidence low or very_low, labelled AI 判断.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Any, Dict, List, Optional, Tuple

from . import labels as L
from .common import EXIT_INPUT, LAError, load_json, longevity_skills_home, now_iso, sha256_file, skillkit, write_json

SCHEMA = "la-predictions/1"
N = NormalDist()
Z80 = 1.2816


# ------------------------------------------------------------ the model (appendix B.7)
def sd_pct(cvi: float, cva: Optional[float]) -> float:
    cva = 0.5 * cvi if cva is None else cva
    return math.sqrt(cvi ** 2 + cva ** 2)


def model(effect: float, ci: Tuple[float, float], baseline: float, cvi: float, cva: Optional[float], k1: int, k2: int, *,
          lower_better: bool = True, personal: List[Tuple[float, int, int]] = (), a: float = 1.0, a_next: float = 1.0,
          tau: Optional[float] = None, z: float = 1.96) -> Dict[str, Any]:
    """Effects in percent of baseline. personal: (observed change %, k before, k after) per earlier experiment."""
    mu0 = 100 * effect / baseline * a
    se = 100 * (ci[1] - ci[0]) / (2 * 1.96) / baseline
    tau = abs(mu0) if tau is None else tau
    s0sq = se ** 2 + tau ** 2

    def v(n1: int, n2: int) -> float:
        return sd_pct(cvi, cva) ** 2 * (1 / n1 + 1 / n2)
    prec, num = 1 / s0sq, mu0 / s0sq
    for d, n1, n2 in personal:
        vi = v(n1, n2)
        prec += 1 / vi
        num += d / vi
    ssq = 1 / prec
    mu = ssq * num
    vn = v(k1, k2)
    band = z * math.sqrt(vn)
    sd_d = math.sqrt(ssq + vn)
    m = mu * a_next
    if lower_better:
        p_detect, p_true = N.cdf((-band - m) / sd_d), N.cdf(-mu / math.sqrt(ssq))
    else:
        p_detect, p_true = 1 - N.cdf((band - m) / sd_d), 1 - N.cdf(-mu / math.sqrt(ssq))
    return {"mu0_pct": mu0, "se_pct": se, "tau_pct": tau, "s0_pct": math.sqrt(s0sq), "band_pct": band,
            "posterior_mu_pct": mu, "posterior_s_pct": math.sqrt(ssq), "p_detect_favorable": p_detect,
            "p_true_favorable": p_true, "expected_change_pct": m, "interval_80_pct": [m - Z80 * sd_d, m + Z80 * sd_d]}


# ------------------------------------------------------------ library and member inputs
def _library() -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    home = longevity_skills_home()
    effects = {}
    with open(home / "data" / "effects.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                r = json.loads(line)
                effects[r["id"]] = r
    return effects, load_json(home / "data" / "biological_variation.json")


def _bv_match(bv: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
    kit = skillkit()
    idx = {}
    for m in bv["markers"]:
        for n in [m["key"], m.get("label_zh", ""), *m.get("aliases", [])]:
            if n:
                idx[kit.fold_name(n)] = m
    for v in kit.name_variants(str(name)):
        if v in idx:
            return idx[v]
    return None


def _num(v: Any) -> Optional[float]:
    s = str(v).strip().replace(",", ".")
    if s[:1] in "<>≤≥":
        return None
    try:
        x = float(s)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def _baseline(st: Dict[str, Any], ws: Path, bv: Dict[str, Any], marker: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """The member's latest value of this marker in the marker's own unit: confirmed lab rows, else the series."""
    from .labnames import usable_rows
    kit = skillkit()
    rows = [r for r in usable_rows(st) if _bv_match(bv, r["marker"]) is marker
            and kit.normalize_unit(r.get("unit", "")) == kit.normalize_unit(marker["unit"]) and _num(r["value"])]
    vals = {_num(r["value"]) for r in rows}
    if len(vals) == 1:
        r = rows[0]
        return {"value": _num(r["value"]), "marker": r["marker"], "from": "lab"}
    sp = ws / "work" / "twin" / "series.json"
    for s in _series(sp):
        if _bv_match(bv, s.get("marker", "")) is marker and s.get("points") \
                and kit.normalize_unit(s.get("unit", "")) == kit.normalize_unit(marker["unit"]):
            v = _num(s["points"][-1].get("value"))
            if v:
                return {"value": v, "marker": s["marker"], "from": "series"}
    return None


def _series(sp: Path) -> List[Dict[str, Any]]:
    if not sp.exists():
        return []
    try:
        d = load_json(sp)
    except (OSError, ValueError):
        return []
    ms = d.get("markers", d.get("series", []))
    if isinstance(ms, dict):
        ms = [dict(v, marker=v.get("marker", k)) for k, v in ms.items()]
    return [m for m in ms if isinstance(m, dict)]


def _mentions(lst: Any, item: Dict[str, Any], row: Dict[str, Any]) -> bool:
    for e in lst or []:
        if not isinstance(e, dict):
            continue
        if e.get("plan_item") == item["id"] or (e.get("label_zh") and e["label_zh"] in [row.get("intervention_zh"), *row.get("keywords", [])]):
            return True
    return False


def personal_experiments(ws: Path, bv: Dict[str, Any], marker: Dict[str, Any], item: Dict[str, Any],
                         row: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Earlier before/after pairs of this marker during which this intervention ran (twin series)."""
    out = []
    for s in _series(ws / "work" / "twin" / "series.json"):
        if _bv_match(bv, s.get("marker", "")) is not marker:
            continue
        k_at = {p.get("t"): int(p.get("k") or 1) for p in s.get("points", [])}
        pairs = [p for p in s.get("pairs", []) if isinstance(p, dict) and p.get("change_pct") is not None]
        if any("concurrent_exposures" in p for p in pairs):
            hits = [p for p in pairs if _mentions(p.get("concurrent_exposures"), item, row)]
        else:                              # marker-level list: the exposure window is the latest comparison
            hits = pairs[-1:] if pairs and _mentions(s.get("concurrent_exposures"), item, row) else []
        for p in hits:
            out.append({"ref": f"{s.get('marker')}:{p.get('from')}→{p.get('to')}", "change_pct": float(p["change_pct"]),
                        "k1": k_at.get(p.get("from"), 1), "k2": k_at.get(p.get("to"), 1)})
    return out


def _exposures(ws: Path) -> List[Dict[str, Any]]:
    for p in (ws / "work" / "predict" / "exposures.json", ws / "work" / "twin" / "twin.json"):
        if p.exists():
            try:
                d = load_json(p)
            except (OSError, ValueError):
                continue
            items = d.get("items") if d.get("schema") == "la-exposures/1" else d.get("exposures")
            if items:
                return [x for x in items if isinstance(x, dict)]
    return []


def adherence(exp: Optional[Dict[str, Any]], planned_per_week: float, today: str) -> Tuple[float, str]:
    """Check-ins / planned (cap 1), with the assumption in plain Chinese."""
    if not exp or exp.get("status") in (None, "proposed") or not exp.get("adopted_at"):
        return 1.0, "尚未开始执行，依从比例按 1 计"
    if exp.get("cumulative") is None:
        return 1.0, "宿主只记录是否采纳、不记录打卡次数，依从比例按 1 计"
    try:
        weeks = max((_dt.date.fromisoformat(today[:10]) - _dt.date.fromisoformat(exp["adopted_at"])).days / 7, 1 / 7)
    except ValueError:
        return 1.0, "采纳日期无法读取，依从比例按 1 计"
    planned = planned_per_week * weeks
    a = min(1.0, exp["cumulative"] / planned) if planned > 0 else 1.0
    return a, f"依从比例按打卡次数线性缩放：累计打卡 {exp['cumulative']} 次，计划约 {planned:.0f} 次（每周 {planned_per_week:g} 次），比例 {a:.2f}"


def _plan(ws: Path) -> List[Dict[str, Any]]:
    p = ws / "work" / "intervene" / "plan.json"
    return load_json(p).get("items", []) if p.exists() else []


def inputs_digest(st: Dict[str, Any], ws: Path) -> str:
    from .labnames import usable_rows
    home = longevity_skills_home(required=False)
    parts = {"plan": _sha(ws / "work" / "intervene" / "plan.json"), "series": _sha(ws / "work" / "twin" / "series.json"),
             "exposures": L.digest(_exposures(ws)), "labs": L.digest([[r["marker"], r["value"], r.get("unit")] for r in usable_rows(st)]),
             "ai": _sha(ws / "work" / "predict" / "ai_responses.json"),
             "effects": _sha(home / "data" / "effects.jsonl") if home else None,
             "biovar": _sha(home / "data" / "biological_variation.json") if home else None,
             "opts": (st.get("predict") or {}).get("opts")}
    return L.digest(parts)


def _sha(p: Path) -> Optional[str]:
    return sha256_file(p, limit=None) if p.exists() else None


# ------------------------------------------------------------ predict response
def predict_item(item: Dict[str, Any], row: Dict[str, Any], marker: Dict[str, Any], base: Dict[str, Any], k: int, z: float,
                 exp: Optional[Dict[str, Any]], personal: List[Dict[str, Any]], planned_per_week: float, today: str,
                 pid: str) -> Dict[str, Any]:
    eff = row["effect"]
    a, a_note = adherence(exp, planned_per_week, today)
    tau = row.get("tau_pct")
    res = model(eff["value"], tuple(eff["ci"]), base["value"], marker["cvi_pct"], marker.get("cva_pct"), k, k,
                lower_better=marker.get("better") == "lower", personal=[(p["change_pct"], p["k1"], p["k2"]) for p in personal],
                a=a, a_next=1.0, tau=tau, z=z)
    assumptions = [a_note, "下次实验按完全依从计算",
                   "个人差异按效应表给出的值" if tau is not None else "个人差异默认与平均效应同量级（效应表未给出）",
                   "效应和测量误差按正态分布近似", "效应换算为相对基线的百分比，基线取你最近一次的测量值",
                   "测量误差取个体内变异和分析变异" + ("（分析变异按 0.5 倍个体内变异缺省）" if marker.get("cva_pct") is None else "")]
    if personal:
        assumptions += ["你本人做过的各次实验相互独立",
                        "你本人的前后对比没有对照组，同一段时间里其他因素带来的变化也被算作这项干预的效果"]
    iv = [round(res["interval_80_pct"][0], 1), round(res["interval_80_pct"][1], 1)]
    sources = [{"type": "effects", "ref": row["id"]}, {"type": "biovar", "ref": marker["key"]}]
    if row.get("doi"):
        sources.append({"type": "doi", "ref": row["doi"]})
    sources += [{"type": "personal_experiment", "ref": p["ref"]} for p in personal]
    label = L.make_label("statistical_model", sources, assumptions=assumptions, confidence={"interval_80": iv, "unit": "%"})
    return {"schema": "la-prediction/1", "id": pid, "plan_item": item["id"], "intervention_zh": row.get("intervention_zh") or item["action_zh"],
            "marker": marker.get("label_zh") or marker["key"],
            "design": {"repeats_each_side": k, "weeks": (item.get("retest") or {}).get("after_weeks"), "adherence": round(a, 3)},
            "p_detect_favorable": round(res["p_detect_favorable"], 3), "p_true_favorable": round(res["p_true_favorable"], 3),
            "expected_change_pct": round(res["expected_change_pct"], 1), "interval_80_pct": iv, "label": label}, \
        {"marker_key": marker["key"], "baseline": {"value": base["value"], "unit": marker["unit"], "from": base["from"], "marker": base["marker"]},
         "model": {k_: round(v, 4) for k_, v in res.items() if isinstance(v, float)}}


def run(st: Dict[str, Any], ws: Path, *, repeats: Optional[int] = None, planned_per_week: float = 7.0,
        ai_file: Optional[Path] = None, exposures_file: Optional[Path] = None) -> Dict[str, Any]:
    """Predict every plan item that has an effects row and a judgeable target; list the rest for the AI estimator."""
    effects, bv = _library()
    z = float(bv.get("z", 1.96))
    today = str(st["member"].get("sample_date") or now_iso())[:10]     # adherence up to this visit, not up to the run date
    if exposures_file:
        d = load_json(exposures_file)
        if d.get("schema") != "la-exposures/1":
            raise LAError(f"{exposures_file} is not la-exposures/1", EXIT_INPUT)
        write_json(ws / "work" / "predict" / "exposures.json", d)
    exps = {e.get("plan_item"): e for e in _exposures(ws)}
    items, requests, skipped, details = [], [], [], {}
    for it in _plan(ws):
        if it.get("category") in ("test", "referral", "medication"):
            continue
        rows = [effects[e["ref"]] for e in it.get("evidence", []) if e.get("type") == "effects" and e.get("ref") in effects]
        targets = list(it.get("targets", [])) + ([str((it.get("retest") or {}).get("what"))] if (it.get("retest") or {}).get("what") else [])
        tmarkers = [m for m in (_bv_match(bv, t) for t in targets) if m]
        done = False
        why = "效应表没有对应条目" if not rows else ""
        for row in rows:
            m = bv_m = None
            if row.get("marker_key"):
                bv_m = next((x for x in bv["markers"] if x["key"] == row["marker_key"]), None)
            m = bv_m or _bv_match(bv, row.get("marker_zh", "")) or _bv_match(bv, row.get("marker", ""))
            if not m or m not in tmarkers:
                why = "效应表条目的指标不是这一项的目标，或该指标没有个体内变异数据（不可判定）"
                continue
            if (row.get("effect") or {}).get("kind") != "mean_difference" or kit_unit(row["effect"].get("unit")) != kit_unit(m["unit"]):
                why = "效应表给的是标准化效应或单位与指标不同，不能换算为相对变化"
                continue
            base = _baseline(st, ws, bv, m)
            if not base:
                why = "没有你本人这项指标的基线值"
                continue
            k = repeats or int(m.get("average_days") or 1)
            pers = personal_experiments(ws, bv, m, it, row)
            pi, det = predict_item(it, row, m, base, k, z, exps.get(it["id"]), pers, planned_per_week, today, f"P{len(items) + 1}")
            items.append(pi)
            details[pi["id"]] = det
            done = True
            break
        if not done:
            requests.append({"plan_item": it["id"], "action_zh": it.get("action_zh"), "targets": it.get("targets"),
                             "why_not_model_zh": why})
    if not ai_file and (ws / "work" / "predict" / "ai_responses.json").exists():
        (ws / "work" / "predict" / "ai_responses.json").unlink()      # AI answers belong to the plan they were asked about
    if ai_file:
        for r in _ai_items(ws, ai_file, {x["plan_item"] for x in requests}):
            items.append(dict(r, id=f"P{len(items) + 1}"))
    out = {"schema": SCHEMA, "generated_at": now_iso(), "z": z, "items": items, "details": details}
    write_json(ws / "work" / "predict" / "predictions.json", out)
    write_json(ws / "work" / "predict" / "ai_requests.json", {"schema": "la-predict-requests/1", "items": requests,
                                                               "estimator": "references/response-estimator.md"})
    st["predict"] = dict(st.get("predict") or {}, opts={"repeats": repeats, "planned_per_week": planned_per_week})
    L.record(st, ws, "p", inputs_digest(st, ws))
    return {"path": str(ws / "work" / "predict" / "predictions.json"), "statistical": sum(1 for x in items if x["label"]["basis_kind"] == "statistical_model"),
            "ai": sum(1 for x in items if x["label"]["basis_kind"] == "ai_estimate"), "ai_requests": requests, "skipped": skipped}


def kit_unit(u: Any) -> str:
    return skillkit().normalize_unit(str(u or ""))


def _ai_items(ws: Path, ai_file: Path, wanted: set) -> List[Dict[str, Any]]:
    """Register the response-estimator's answers: interval only, confidence low/very_low, PMIDs retrieved this run."""
    from .evidence import _cached_pmids, verify_pmids
    from .report import trace_text
    doc = load_json(ai_file)
    rows = doc.get("items") if isinstance(doc, dict) else None
    if not isinstance(rows, list):
        raise LAError("the estimator's file needs {\"items\": [...]} (references/response-estimator.md)", EXIT_INPUT)
    cached = _cached_pmids(ws)
    want_pm = sorted({str(p) for r in rows if isinstance(r, dict) for p in r.get("pmids", [])})
    live = verify_pmids(want_pm) if want_pm else {}
    out, bad = [], []
    for r in rows:
        if not isinstance(r, dict):
            bad.append("each item is an object")
            continue
        pi = r.get("plan_item")
        if pi not in wanted:
            bad.append(f"{pi}: not an item the model could not cover (see work/predict/ai_requests.json)")
            continue
        for k in ("p_detect_favorable", "p_true_favorable"):
            v = r.get(k)
            if v is not None and not (isinstance(v, list) and len(v) == 2 and all(isinstance(x, (int, float)) for x in v) and 0 <= v[0] <= v[1] <= 1):
                bad.append(f"{pi}: {k} must be an interval [low, high] between 0 and 1, never a point")
        if r.get("p_detect_favorable") is None and r.get("p_true_favorable") is None:
            bad.append(f"{pi}: give p_detect_favorable or p_true_favorable as an interval")
        if r.get("confidence") not in L.AI_LEVELS:
            bad.append(f"{pi}: confidence must be low or very_low")
        pm = [str(p) for p in r.get("pmids", [])]
        if not pm:
            bad.append(f"{pi}: anchor the estimate on at least one PMID retrieved with `la.py evidence pubmed`")
        for p in pm:
            if p not in cached or p not in live:
                bad.append(f"{pi}: PMID {p} was not retrieved in this run or PubMed does not confirm it")
        for k in ("basis_zh", "intervention_zh", "marker"):
            if not str(r.get(k) or "").strip():
                bad.append(f"{pi}: missing {k}")
        prose = str(r.get("basis_zh") or "") + "\n" + str(r.get("intervention_zh") or "")
        tb = trace_text(prose) + (["{{…}}"] if "{{" in prose else [])
        if tb:
            bad.append(f"{pi}: basis_zh/intervention_zh carry numbers or links ({tb[0]})")
        if bad:
            continue
        label = L.make_label("ai_estimate", [{"type": "pmid", "ref": p} for p in pm], confidence={"level": r["confidence"]},
                             assumptions=[r["basis_zh"]])          # the estimator's basis, in words
        out.append({"schema": "la-prediction/1", "plan_item": pi, "intervention_zh": r["intervention_zh"], "marker": r.get("marker_zh") or r["marker"],
                    **{k: r[k] for k in ("p_detect_favorable", "p_true_favorable") if r.get(k) is not None}, "label": label})
    if bad:
        raise LAError("AI predictions rejected:\n- " + "\n- ".join(bad), EXIT_INPUT)
    write_json(ws / "work" / "predict" / "ai_responses.json", doc)
    return out
