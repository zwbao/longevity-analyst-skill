"""Several twin snapshots of one member: adjacent pairs, first vs latest, and a trend test (PLAN appendix B).

Noise comes from published within-person biological variation (longevity-skills data/biological_variation.json)
for lab markers, and from technical-replicate noise (data/readout_noise.json, when the library has it) for readouts
such as clocks. A marker without a noise model is shown side by side and never called better or worse.

The trend test asks whether the points are consistent with a stable set point given the published noise: a
weighted least-squares slope with known per-point error, z = b / SE_b, |z| > z_trend (default 2.576), with at least
3 points over at least 180 days. `twin compare` is the two-point case: it uses `pair()` below.
"""
from __future__ import annotations

import datetime as _dt
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import EXIT_INPUT, LAError, load_json, longevity_skills_home, now_iso, skillkit

DEFAULT_MIN_DAYS = 28
TREND_MIN_POINTS, TREND_MIN_SPAN_DAYS = 3, 180
Z_TREND = 2.576
IDENTITY_MAX_MISMATCH, IDENTITY_MAX_FRACTION = 2, 0.01

CAVEAT_ZH = {
    "censored_or_non_numeric": "删失值（<、>）、非数值或 0 值不参与判断",
    "unconfirmed_source": "来源未确认属于本人的值不参与判断",
    "same_day_duplicate_first_kept": "同一天有多个值，只取文件中的第一个",
    "unit_not_convertible": "单位无法换算到噪声表的单位，这些值不参与判断",
    "unit_unverified": "单位与噪声表不同但各次一致，按相对变化判断",
    "returned_within_noise": "曾有一次变化超出波动，之后又回到了波动范围内",
    "seasonal": "25(OH)D 有季节性波动，不同季节的测量差异可能来自日照而非干预",
    "snapshot_excluded_identity": "有一次检测的胚系基因型与其他几次不一致，那一次的数据不参与对比",
}
IDENTITY_ZH = "两次检测的胚系基因型不一致，胚系 DNA 不会改变，可能拿错了样本，这次的对比暂不采用。"
READOUT_BEYOND_ZH = "超出测量误差"
READOUT_NOTE_ZH = "超出测量误差不等于生物学上的真实变化。"


# ------------------------------------------------------------ noise tables
def readout_noise_path() -> Path:
    return longevity_skills_home() / "data" / "readout_noise.json"


def noise() -> Dict[str, Any]:
    """Biological-variation index (folded name -> marker) plus the readout-noise table when the library has one."""
    home = longevity_skills_home()
    bv = load_json(home / "data" / "biological_variation.json")
    kit = skillkit()
    idx: Dict[str, Dict[str, Any]] = {}
    for m in bv["markers"]:
        for n in [m["key"], m.get("label_zh", ""), *m.get("aliases", [])]:
            if n:
                idx[kit.fold_name(n)] = m
    rp = readout_noise_path()
    rn = load_json(rp) if rp.exists() else None
    return {"z": bv["z"], "z_trend": bv.get("z_trend") or Z_TREND, "idx": idx, "kit": kit, "readouts": rn}


def match(N: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
    for v in N["kit"].name_variants(str(name or "")):
        if v in N["idx"]:
            return N["idx"][v]
    return None


def match_readout(N: Dict[str, Any], rid: str, label: str = "") -> Optional[Dict[str, Any]]:
    rn = N.get("readouts")
    if not rn:
        return None
    kit = N["kit"]
    names = {kit.fold_name(x) for x in (rid, label) if x}
    for r in rn.get("readouts", []):
        keys = {kit.fold_name(x) for x in [r.get("key", ""), *r.get("match", [])] if x}
        if names & keys or any(str(rid).endswith("." + k) or str(rid).endswith("_" + k) for k in keys if k):
            return r
    return None


def resolve(N: Dict[str, Any], what: str) -> Optional[Dict[str, Any]]:
    """What a retest target resolves to: {table, key, min_retest_days} or None (judgeable = found)."""
    m = match(N, what)
    if m:
        return {"table": "biological_variation", "key": m["key"], "min_retest_days": m.get("min_retest_days")}
    r = match_readout(N, what, what)
    if r and _sd_abs(r) is not None:
        return {"table": "readout_noise", "key": r["key"], "min_retest_days": r.get("min_retest_days")}
    return None


def _sd_abs(r: Dict[str, Any]) -> Optional[float]:
    try:
        return float(r.get("sd_abs"))
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------ preprocessing
def num(v: Any) -> Optional[float]:
    """Plain numbers only: censored values such as '<0.5' or '>90' are not judged."""
    s = str(v).strip().replace(",", ".")
    if s[:1] in "<>≤≥":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def days(a: str, b: str) -> Optional[int]:
    try:
        return (_dt.date.fromisoformat(str(b)[:10]) - _dt.date.fromisoformat(str(a)[:10])).days
    except ValueError:
        return None


def to_unit(m: Dict[str, Any], unit: str, kit) -> Optional[float]:
    """Factor from `unit` to the table unit: 1 for the table unit, the declared `convert` factor, else None."""
    u = kit.normalize_unit(unit or "")
    if u == kit.normalize_unit(m.get("unit", "")):
        return 1.0
    for k, f in (m.get("convert") or {}).items():
        if kit.normalize_unit(k) == u:
            return float(f)
    return None


def _caveat(out: List[Dict[str, Any]], code: str, **kw: Any) -> None:
    if not any(c["code"] == code for c in out):
        out.append({"code": code, "zh": CAVEAT_ZH.get(code, code), **kw})


def prepare(obs: List[Dict[str, Any]], m: Optional[Dict[str, Any]], kit) -> Tuple[List[Tuple[str, float]], List[Dict[str, Any]]]:
    """Observations -> (date, value in the table unit) in input order, same-day duplicates dropped (first kept)."""
    caveats: List[Dict[str, Any]] = []
    rows = []
    for o in obs:
        if o.get("provenance_uncertain"):
            _caveat(caveats, "unconfirmed_source")
            continue
        v = num(o.get("value"))
        if v is None or v == 0:
            _caveat(caveats, "censored_or_non_numeric")
            continue
        rows.append((str(o["t"])[:10], v, o.get("unit", "")))
    if m:
        units = {kit.normalize_unit(u) for _, _, u in rows}
        conv = []
        for t, v, u in rows:
            f = to_unit(m, u, kit)
            if f is None and len(units) == 1:
                f = 1.0                                # one consistent unit throughout: relative change still holds
                _caveat(caveats, "unit_unverified")
            if f is None:
                _caveat(caveats, "unit_not_convertible")
                continue
            conv.append((t, v * f))
        rows2 = conv
    else:
        rows2 = [(t, v) for t, v, _ in rows]
    seen, out = set(), []
    for t, v in rows2:
        if t in seen:
            _caveat(caveats, "same_day_duplicate_first_kept", t=t)
            continue
        seen.add(t)
        out.append((t, v))
    return out, caveats


def cluster(series: List[Tuple[str, float]], min_days: int, log_normal: bool) -> List[Dict[str, Any]]:
    """Windows of `min_days` from the earliest observation; one point per window (median date, mean of values, k)."""
    pts = sorted(((_dt.date.fromisoformat(t), v) for t, v in series), key=lambda x: x[0])
    if not pts:
        return []
    start = pts[0][0]
    win: Dict[int, List[Tuple[_dt.date, float]]] = {}
    for d, v in pts:
        win.setdefault((d - start).days // min_days, []).append((d, v))
    out = []
    for w in sorted(win):
        g = win[w]
        ds = sorted(d for d, _ in g)
        k = len(g)
        med = ds[(k - 1) // 2] if k % 2 else ds[k // 2 - 1] + (ds[k // 2] - ds[k // 2 - 1]) / 2
        vals = [v for _, v in g]
        val = math.exp(sum(math.log(x) for x in vals) / k) if log_normal else sum(vals) / k
        out.append({"t": med.isoformat(), "value": round(val, 6), "k": k})
    return out


# ------------------------------------------------------------ judgments
def _cva(m: Dict[str, Any]) -> float:
    return (m["cva_pct"] if m.get("cva_pct") is not None else 0.5 * m["cvi_pct"]) / 100


def band(m: Dict[str, Any], z: float, k1: int = 1, k2: int = 1) -> Tuple[float, float]:
    """(down_pct, up_pct): the RCV band for means of k1 and k2 separate-day results."""
    f = math.sqrt(1 / k1 + 1 / k2)
    cvi = m["cvi_pct"] / 100
    if m.get("log_normal"):
        s = math.sqrt(math.log(cvi ** 2 + 1) + math.log(_cva(m) ** 2 + 1))
        return 100 * (math.exp(-z * f * s) - 1), 100 * (math.exp(z * f * s) - 1)
    b = 100 * z * f * math.sqrt(cvi ** 2 + _cva(m) ** 2)
    return -b, b


def pair(a: Dict[str, Any], b: Dict[str, Any], m: Dict[str, Any], z: float) -> Dict[str, Any]:
    """One comparison of two points (dicts with t, value, k)."""
    out = {"from": a["t"], "to": b["t"]}
    gap = days(a["t"], b["t"])
    min_days = m.get("min_retest_days") or DEFAULT_MIN_DAYS
    if gap is not None and gap < 1:
        return {**out, "verdict": "not_judged", "why": "same_day"}
    if gap is not None and gap < min_days:
        return {**out, "verdict": "not_judged", "why": "too_close", "days": gap, "min_retest_days": min_days}
    lo, hi = band(m, z, a.get("k", 1), b.get("k", 1))
    ch = 100 * (b["value"] - a["value"]) / a["value"]
    return {**out, "change_pct": round(ch, 1), "band_pct": [round(lo, 1), round(hi, 1)],
            "verdict": "increase_beyond_noise" if ch > hi else "decrease_beyond_noise" if ch < lo else "within_noise"}


def trend(points: List[Dict[str, Any]], m: Dict[str, Any], z_trend: float) -> Dict[str, Any]:
    n = len(points)
    if n < TREND_MIN_POINTS:
        return {"n": n, "verdict": "not_judged", "why": "fewer_than_3_points"}
    span = days(points[0]["t"], points[-1]["t"]) or 0
    if span < TREND_MIN_SPAN_DAYS:
        return {"n": n, "span_days": span, "verdict": "not_judged", "why": "span_under_180_days"}
    t0 = _dt.date.fromisoformat(points[0]["t"])
    t = [(_dt.date.fromisoformat(p["t"]) - t0).days / 365.25 for p in points]
    cvi = (m.get("long_term_cvi_pct") or m["cvi_pct"]) / 100
    mean = sum(p["value"] for p in points) / n
    if m.get("log_normal"):
        sl = math.sqrt(math.log(cvi ** 2 + 1) + math.log(_cva(m) ** 2 + 1))
        y = [math.log(p["value"]) for p in points]
        s = [sl / math.sqrt(p["k"]) for p in points]
    else:
        y = [p["value"] for p in points]
        s = [math.sqrt(cvi ** 2 + _cva(m) ** 2) * mean / math.sqrt(p["k"]) for p in points]
    w = [1 / x ** 2 for x in s]
    W = sum(w)
    tw = sum(a * b for a, b in zip(w, t)) / W
    yw = sum(a * b for a, b in zip(w, y)) / W
    sxx = sum(wi * (ti - tw) ** 2 for wi, ti in zip(w, t))
    b = sum(wi * (ti - tw) * (yi - yw) for wi, ti, yi in zip(w, t, y)) / sxx
    z = b / (1 / math.sqrt(sxx))
    annual = (math.exp(b) - 1) * 100 if m.get("log_normal") else b / mean * 100
    verdict = ("rising_beyond_noise" if z > 0 else "falling_beyond_noise") if abs(z) > z_trend else "no_trend_detected"
    return {"n": n, "span_days": span, "slope_per_year": round(b, 4), "annual_change_pct": round(annual, 1),
            "z": round(z, 2), "z_crit": z_trend, "scale": "ln" if m.get("log_normal") else "linear", "verdict": verdict}


BEYOND = ("increase_beyond_noise", "decrease_beyond_noise")


def summary(pairs: List[Dict[str, Any]], bvl: Dict[str, Any], tr: Dict[str, Any], caveats: List[Dict[str, Any]]) -> str:
    if tr.get("verdict") in ("rising_beyond_noise", "falling_beyond_noise"):
        return "trend_beyond_noise"
    if bvl.get("verdict") in BEYOND:
        return "change_beyond_noise"
    if any(p.get("verdict") in BEYOND for p in pairs):
        if bvl.get("verdict") == "within_noise":
            _caveat(caveats, "returned_within_noise")
        return "change_beyond_noise"
    if any(p.get("verdict") == "within_noise" for p in pairs) or bvl.get("verdict") == "within_noise" \
            or tr.get("verdict") == "no_trend_detected":
        return "within_noise"
    return "not_judged"


def marker_series(name: str, obs: List[Dict[str, Any]], m: Optional[Dict[str, Any]], N: Dict[str, Any]) -> Dict[str, Any]:
    """One lab marker across snapshots: A.4 of the plan."""
    kit = N["kit"]
    pts, caveats = prepare(obs, m, kit)
    unit = m["unit"] if m else next((o.get("unit", "") for o in reversed(obs) if o.get("unit")), "")
    out: Dict[str, Any] = {"marker": name, "marker_key": m["key"] if m else None, "unit": unit}
    if not m:
        out.update(points=[{"t": t, "value": v, "k": 1} for t, v in sorted(pts)], pairs=[],
                   baseline_vs_latest={"verdict": "not_judged", "why": "no_noise_model"},
                   trend={"verdict": "not_judged", "why": "no_noise_model"}, summary_verdict="not_judged",
                   why="no within-person biological variation data for this marker; shown side by side only",
                   caveats=caveats)
        out["observed_at"] = out["points"][-1]["t"] if out["points"] else None
        return out
    points = cluster(pts, m.get("min_retest_days") or DEFAULT_MIN_DAYS, bool(m.get("log_normal")))
    z = N["z"]
    pairs = [pair(points[i], points[i + 1], m, z) for i in range(len(points) - 1)]
    if len(points) >= 2:
        bvl = pair(points[0], points[-1], m, z)
    else:
        bvl = {"verdict": "not_judged", "why": "fewer_than_2_points"}
    tr = trend(points, m, N["z_trend"])
    if m["key"] == "vitd" and len({p["t"][5:7] for p in points}) > 1:
        _caveat(caveats, "seasonal")
    out.update(points=points, pairs=pairs, baseline_vs_latest=bvl, trend=tr,
               summary_verdict=summary(pairs, bvl, tr, caveats), bv_source=(m.get("cvi_source") or {}).get("doi"),
               caveats=caveats, observed_at=points[-1]["t"] if points else None)
    if m.get("better"):
        out["better"] = m["better"]
    return out


def readout_pair(a: Dict[str, Any], b: Dict[str, Any], r: Optional[Dict[str, Any]], z: float) -> Dict[str, Any]:
    """Readouts (clocks …) are judged only against a technical-replicate noise table: |Δ| > z·√2·sd_abs."""
    va, vb = num(a.get("value")), num(b.get("value"))
    sd = _sd_abs(r) if r else None
    if sd is None or va is None or vb is None:
        why = ("no published within-person noise model is loaded for this readout; show the two values, do not call it better or worse"
               if sd is None else "censored or non-numeric value")
        return {"verdict": "not_judged", "why": why}
    d = vb - va
    lim = z * math.sqrt(2) * sd
    v = "increase_beyond_noise" if d > lim else "decrease_beyond_noise" if d < -lim else "within_noise"
    res = {"change_abs": round(d, 2), "band_abs": [round(-lim, 2), round(lim, 2)], "verdict": v, "noise_key": r["key"],
           "noise_kind": r.get("noise_kind"), "noise_source": (r.get("source") or {}).get("doi") or (r.get("source") or {}).get("title")}
    if v != "within_noise":
        res.update(wording_zh=READOUT_BEYOND_ZH, note_zh=READOUT_NOTE_ZH)
    return res


# ------------------------------------------------------------ identity
def _calls(snap: Dict[str, Any]) -> Dict[str, str]:
    fp = snap.get("genotype_fingerprint") or {}
    return {k: str(v) for k, v in (fp.get("calls") or {}).items()}


def _apoe(snap: Dict[str, Any]) -> Optional[str]:
    for r in snap.get("readouts", []):
        if str(r.get("id", "")).startswith("native.apoe"):
            return str(r.get("value"))
    return None


def fingerprint_mismatch(a: Dict[str, str], b: Dict[str, str], max_mismatch: int = IDENTITY_MAX_MISMATCH,
                         max_fraction: float = IDENTITY_MAX_FRACTION) -> Optional[Dict[str, Any]]:
    """None when the overlap is empty; else {overlap, mismatches, alert}."""
    common = sorted(set(a) & set(b))
    if not common:
        return None
    mm = sum(1 for k in common if "".join(sorted(a[k])) != "".join(sorted(b[k])))
    return {"overlap": len(common), "mismatches": mm, "alert": mm >= max_mismatch or mm / len(common) > max_fraction}


def identity(snaps: List[Dict[str, Any]], paths: List[str], max_mismatch: int = IDENTITY_MAX_MISMATCH,
             max_fraction: float = IDENTITY_MAX_FRACTION) -> Tuple[List[Dict[str, Any]], set]:
    """Which snapshots disagree with the rest on germline genotype. The reference is the snapshot that agrees with
    the most others (ties: the earliest); a snapshot that disagrees with it leaves every comparison."""
    n = len(snaps)

    def agree(i: int, j: int) -> Optional[Tuple[bool, Dict[str, Any]]]:
        fm = fingerprint_mismatch(_calls(snaps[i]), _calls(snaps[j]), max_mismatch, max_fraction)
        if fm is not None:
            return (not fm["alert"], {"basis": "fingerprint", **fm})
        a, b = _apoe(snaps[i]), _apoe(snaps[j])
        if a is not None and b is not None:
            return (a == b, {"basis": "apoe", "previous": a, "current": b})
        return None

    if n < 2:
        return [], set()
    score = []
    for i in range(n):
        score.append(sum(1 for j in range(n) if j != i and (agree(i, j) or (True,))[0]))
    ref = max(range(n), key=lambda i: (score[i], -i))
    alerts, excluded = [], set()
    for j in range(n):
        if j == ref:
            continue
        res = agree(ref, j)
        if res and not res[0]:
            excluded.add(j)
            alerts.append({"snapshot": paths[j], "t": (snaps[j].get("snapshot") or {}).get("t"), "against": paths[ref],
                           **res[1], "verdict": "identity_alert", "wording_zh": IDENTITY_ZH})
    return alerts, excluded


# ------------------------------------------------------------ exposures
def _target_keys(N: Dict[str, Any], targets: List[str]) -> set:
    kit = N["kit"]
    out = set()
    for t in targets or []:
        m = match(N, t)
        out.add(("bv", m["key"]) if m else ("name", kit.fold_name(t)))
    return out


def exposures_of(snaps: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Exposure items across snapshots; a later snapshot's record of the same plan item wins."""
    by: Dict[str, Dict[str, Any]] = {}
    for s in snaps:
        for e in s.get("exposures") or []:
            by[str(e.get("plan_item"))] = e
    return list(by.values())


def concurrent(N: Dict[str, Any], exposures: List[Dict[str, Any]], marker: str, key: Optional[str], t0: str, t1: str) -> List[Dict[str, Any]]:
    """Adopted or graduated items targeting this marker that were under way between t0 and t1."""
    kit = N["kit"]
    me = ("bv", key) if key else ("name", kit.fold_name(marker))
    out = []
    for e in exposures:
        if e.get("status") not in ("adopted", "graduated") or me not in _target_keys(N, e.get("targets") or []):
            continue
        if e.get("adopted_at") and str(e["adopted_at"]) > t1:
            continue
        weeks = [w for w in e.get("weekly_checkins") or [] if t0 <= str(w.get("week_start")) <= t1]
        out.append({"plan_item": e.get("plan_item"), "label_zh": e.get("label_zh"), "status": e.get("status"),
                    "adopted_at": e.get("adopted_at"),
                    "cumulative_in_window": sum(int(w["count"]) for w in weeks) if weeks else None,
                    "cumulative_total": e.get("cumulative")})
    return out


def _window(s: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    if s["summary_verdict"] == "trend_beyond_noise" or s["baseline_vs_latest"].get("verdict") in BEYOND:
        return s["points"][0]["t"], s["points"][-1]["t"]
    for p in s["pairs"]:
        if p.get("verdict") in BEYOND:
            return p["from"], p["to"]
    return None


# ------------------------------------------------------------ concordance transitions
TRANSITION_ZH = {
    "expressed": "有先天倾向的这一项，这次已经表现出来了，建议列入问医生的单子。",
    "back_real": "这一项回到了参考范围内，而且变化超出了你的正常波动，是真实的变化。",
    "back_within": "这一项回到了参考范围内，但变化还在你的正常波动之内。",
    "back_no_model": "这一项回到了参考范围内；这一项没有噪声模型，不能判断变化是不是超出了你的正常波动。",
}
PUBLIC_DATA_ZH = "检索的公共数据更新了，遗传分位随之变化；这不是你身上的变化。"


def transitions(snaps: List[Dict[str, Any]], kept: List[int], by_key: Dict[str, Dict[str, Any]], N: Dict[str, Any],
                exposures: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    out, alerts = [], []
    prev: Dict[str, Tuple[int, Dict[str, Any]]] = {}
    for i in kept:
        s = snaps[i]
        for c in s.get("concordance") or []:
            an = c.get("analyte")
            if not an:
                continue
            if an in prev:
                j, p = prev[an]
                t0 = p.get("t") or (snaps[j].get("snapshot") or {}).get("t")
                t1 = c.get("t") or (s.get("snapshot") or {}).get("t")
                base = {"analyte": an, "from_t": t0, "to_t": t1}
                if p.get("genetic_band") != c.get("genetic_band") and "not_computed" not in (p.get("genetic_band"), c.get("genetic_band")):
                    public = (p.get("gwas_retrieved_on") or p.get("retrieved_on")) != (c.get("gwas_retrieved_on") or c.get("retrieved_on")) \
                        or p.get("trait_map_version") != c.get("trait_map_version")
                    if public:
                        out.append({**base, "kind": "public_data_updated", "from": p.get("class"), "to": c.get("class"),
                                    "genetic_band": [p.get("genetic_band"), c.get("genetic_band")], "member_change": False,
                                    "wording_zh": PUBLIC_DATA_ZH})
                    else:
                        alerts.append({"basis": "concordance", "analyte": an, "t": t1, "verdict": "identity_alert",
                                       "wording_zh": IDENTITY_ZH})
                elif p.get("class") != c.get("class"):
                    out.append(_transition(base, p.get("class"), c.get("class"), an, by_key, N, exposures))
            prev[an] = (i, c)
    return out, alerts


def _transition(base, fr, to, analyte, by_key, N, exposures) -> Dict[str, Any]:
    m = match(N, analyte)
    s = by_key.get(m["key"]) if m else None
    tr = {**base, "marker": (s or {}).get("marker") or analyte, "from": fr, "to": to, "member_change": True}
    if (fr, to) == ("latent", "concordant_adverse"):
        return {**tr, "kind": "tendency_expressed", "ask_doctor": True, "wording_zh": TRANSITION_ZH["expressed"]}
    if (fr, to) in (("concordant_adverse", "latent"), ("acquired_adverse", "concordant_typical")):
        sv = (s or {}).get("summary_verdict")
        real = True if sv in ("trend_beyond_noise", "change_beyond_noise") else False if sv == "within_noise" else None
        key = {True: "back_real", False: "back_within", None: "back_no_model"}[real]
        tr.update(kind="back_in_range", real_change=real, wording_zh=TRANSITION_ZH[key])
        if fr == "acquired_adverse":
            tr["concurrent_exposures"] = concurrent(N, exposures, analyte, m["key"] if m else None, base["from_t"] or "", base["to_t"] or "9999")
        return tr
    return {**tr, "kind": "class_changed"}


# ------------------------------------------------------------ the command
def load_snapshot(p: Path) -> Dict[str, Any]:
    try:
        s = load_json(p)
    except (OSError, ValueError) as e:
        raise LAError(f"{p} is not a twin snapshot ({e})", EXIT_INPUT)
    if not isinstance(s, dict) or "observations" not in s:
        raise LAError(f"{p} is not a twin snapshot (no observations)", EXIT_INPUT)
    return s


def _snap_t(s: Dict[str, Any]) -> str:
    return str((s.get("snapshot") or {}).get("t") or max((str(o.get("t", "")) for o in s.get("observations", [])), default=""))


def run(snapshots: List[Dict[str, Any]], paths: List[str], max_mismatch: int = IDENTITY_MAX_MISMATCH) -> Dict[str, Any]:
    if not snapshots:
        raise LAError("twin series needs at least one snapshot", EXIT_INPUT)
    ids = {(s.get("member") or {}).get("id") for s in snapshots}
    if len(ids) != 1:
        raise LAError(f"the snapshots belong to different members ({sorted(map(str, ids))})", EXIT_INPUT)
    if not os.environ.get("LONGEVITY_SKILLS_HOME"):
        lib = next((s.get("longevity_skills_home") for s in reversed(snapshots) if s.get("longevity_skills_home")), None)
        if lib:
            os.environ["LONGEVITY_SKILLS_HOME"] = lib
    order = sorted(range(len(snapshots)), key=lambda i: (_snap_t(snapshots[i]), i))
    snaps = [snapshots[i] for i in order]
    paths = [paths[i] for i in order]
    N = noise()
    warnings = []
    for p, s in zip(paths, snaps):
        if s.get("schema") != "la-twin/2":
            warnings.append(f"{p}: {s.get('schema') or 'no schema'} snapshot: observation dates are the snapshot date, "
                            "and it has no exposures, concordance or genotype fingerprint")
    alerts, excluded = identity(snaps, paths, max_mismatch)
    kept = [i for i in range(len(snaps)) if i not in excluded]
    obs_by: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for i in kept:
        s = snaps[i]
        for o in s.get("observations", []):
            o = dict(o)
            if s.get("schema") != "la-twin/2" or not o.get("t"):
                o["t"] = _snap_t(s)
            m = match(N, o["marker"])
            k = ("bv", m["key"]) if m else ("name", N["kit"].fold_name(o["marker"]))
            e = obs_by.setdefault(k, {"m": m, "name": o["marker"], "obs": []})
            e["name"] = o["marker"]                       # the member's own name for it, latest snapshot
            e["obs"].append(o)
    exposures = exposures_of([snaps[i] for i in kept])
    markers, by_key = [], {}
    for k, e in obs_by.items():
        s = marker_series(e["name"], e["obs"], e["m"], N)
        last = max((o for o in e["obs"] if not o.get("provenance_uncertain") and num(o.get("value")) is not None),
                   key=lambda o: str(o["t"]), default=None)
        if last is not None:
            s["latest_value"], s["latest_unit"], s["latest_ref_range"] = num(last["value"]), last.get("unit", ""), last.get("ref_range", "")
        if excluded:
            _caveat(s["caveats"], "snapshot_excluded_identity")
        w = _window(s) if s["summary_verdict"] in ("trend_beyond_noise", "change_beyond_noise") else None
        if w:
            ce = concurrent(N, exposures, s["marker"], s["marker_key"], *w)
            if ce:
                s["concurrent_exposures"] = ce
        markers.append(s)
        if s["marker_key"]:
            by_key[s["marker_key"]] = s
    readouts = []
    rdict: Dict[str, List[Dict[str, Any]]] = {}
    for i in kept:
        for r in snaps[i].get("readouts", []):
            rdict.setdefault(r["id"], []).append({**r, "t": r.get("t") or _snap_t(snaps[i])})
    for rid, rs in rdict.items():
        if len(rs) < 2 or rid.startswith("native.apoe"):
            continue
        nm = match_readout(N, rid, rs[-1].get("label_zh", ""))
        pr = [{"from": a["t"], "to": b["t"], **readout_pair(a, b, nm, (N.get("readouts") or {}).get("z") or N["z"])}
              for a, b in zip(rs, rs[1:])]
        bvl = readout_pair(rs[0], rs[-1], nm, (N.get("readouts") or {}).get("z") or N["z"])
        readouts.append({"readout": rid, "label_zh": rs[-1].get("label_zh"), "unit": rs[-1].get("unit", ""),
                         "points": [{"t": r["t"], "value": r.get("value")} for r in rs], "pairs": pr,
                         "baseline_vs_latest": bvl, "observed_at": rs[-1]["t"]})
    trans, conc_alerts = transitions(snaps, kept, by_key, N, exposures)
    obs_dates = [p["t"] for s in markers for p in s["points"]] + [r["observed_at"] for r in readouts]
    return {"schema": "la-series/1", "generated_at": now_iso(), "member": next(iter(ids)),
            "snapshots": [{"path": p, "t": _snap_t(s), "schema": s.get("schema"), "excluded": i in excluded}
                          for i, (p, s) in enumerate(zip(paths, snaps))],
            "observed_at": max(obs_dates) if obs_dates else None,
            "z": N["z"], "z_trend": N["z_trend"], "readout_noise": bool(N.get("readouts")),
            "identity_alerts": alerts + conc_alerts, "warnings": warnings,
            "markers": sorted(markers, key=lambda s: (s["summary_verdict"] == "not_judged", str(s["marker"]))),
            "readouts": readouts, "transitions": trans,
            "cva_note": "CVA 缺省按 EFLM 期望分析性能取 0.5×CVI（Sandberg 2024）；实验室给出实测 CVA 时应以实测为准"}
