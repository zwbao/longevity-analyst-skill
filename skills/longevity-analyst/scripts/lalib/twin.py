"""Digital twin v0: a provenance-carrying longitudinal state of one member.

twin.json holds every observation (lab rows, method readouts) with its time
and source, the intervention register and the retest plan. `compare` decides
whether a change between two snapshots exceeds normal within-person
variation, using the reference change value (RCV) from published
within-subject biological variation. Where no variation data exists the
change is reported as "not judged", never as improved or worsened.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional

from .labnames import usable_rows
from .common import EXIT_INPUT, LAError, load_json, longevity_skills_home, now_iso, skillkit, write_json

SCHEMA = "la-twin/1"


def _kit():
    return skillkit()


def build(st: Dict[str, Any], ws_root: Path) -> Dict[str, Any]:
    t = st["member"].get("sample_date") or now_iso()[:10]
    ro_path = ws_root / "work" / "readouts.json"
    readouts = load_json(ro_path)["readouts"] if ro_path.exists() else []
    ip = ws_root / "work" / "insights" / "insight_readouts.json"
    if ip.exists():
        readouts = readouts + [r for r in load_json(ip)["readouts"] if r.get("kind") in ("genetic_score", "population_position", "descriptive")]
    plan_path = ws_root / "work" / "intervene" / "plan.json"
    plan = load_json(plan_path) if plan_path.exists() else {"items": []}
    op = ws_root / "work" / "organs" / "organ_readouts.json"
    organ_ro = load_json(op)["readouts"] if op.exists() else []
    unc = {f["id"] for f in st["files"] if f.get("provenance_uncertain")}
    twin = {
        "schema": SCHEMA,
        "organ_estimates": [{"t": t, **{k: r.get(k) for k in ("id", "organ", "label_zh", "value", "low", "high", "unit",
                                                             "horizon_years", "confidence", "kind")}} for r in organ_ro],
        "longevity_skills_home": st.get("longevity_skills_home"),
        "member": {k: st["member"].get(k) for k in ("id", "age", "sex")},
        "identity": st.get("identity") or {"answer": "not_checked"},
        "snapshot": {"t": t, "built_at": now_iso(), "workspace": str(ws_root)},
        "observations": [{"t": t, "kind": "lab", "marker": l["marker"], "value": l["value"], "unit": l.get("unit", ""),
                          "source_file": l.get("source_file"), "page": l.get("page", ""),
                          "provenance_uncertain": l.get("source_file") in unc} for l in usable_rows(st)],
        "readouts": [{"t": t, "id": r["id"], "label_zh": r["label_zh"], "value": r["value"], "unit": r.get("unit", ""),
                      "kind": r.get("kind"), "method": r.get("method"), "systems": r.get("systems", []),
                      "source": r.get("source"), "provenance_uncertain": bool(r.get("provenance_uncertain"))} for r in readouts],
        "data_files": [{"id": f["id"], "name": f["name"], "kind": f["kind"], "sha256": f["sha256"], "status": f.get("status"),
                        "excluded": bool(f.get("excluded")), "provenance_uncertain": bool(f.get("provenance_uncertain"))}
                       for f in st["files"]],
        "pipelines": [{"run_id": p["run_id"], "pipeline": p["pipeline"], "revision": p["revision"], "status": p["status"]} for p in st["pipelines"]],
        "interventions": [{"id": i["id"], "category": i["category"], "action_zh": i["action_zh"], "targets": i["targets"],
                           "executor": i["executor"], "status": "proposed", "since": None} for i in plan.get("items", [])],
        "retest_plan": [{"intervention": i["id"], "what": i["retest"].get("what"), "after_weeks": i["retest"].get("after_weeks")}
                        for i in plan.get("items", []) if i.get("retest")],
    }
    write_json(ws_root / "work" / "twin" / "twin.json", twin)
    st["twin"] = {"path": str(ws_root / "work" / "twin" / "twin.json"), "built_at": now_iso()}
    return {"path": st["twin"]["path"], "observations": len(twin["observations"]), "readouts": len(twin["readouts"]),
            "interventions": len(twin["interventions"])}


def _bv_index() -> Dict[str, Dict[str, Any]]:
    home = longevity_skills_home()
    bv = load_json(home / "data" / "biological_variation.json")
    kit = _kit()
    idx = {}
    for m in bv["markers"]:
        for n in [m["key"], m.get("label_zh", ""), *m.get("aliases", [])]:
            if n:
                idx[kit.fold_name(n)] = m
    return {"z": bv["z"], "idx": idx, "kit": kit}


def _num(v: Any) -> Optional[float]:
    """Plain numbers only: censored values such as '<0.5' or '>90' are not judged."""
    s = str(v).strip().replace(",", ".")
    if s[:1] in "<>≤≥":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _days(a: str, b: str) -> Optional[int]:
    import datetime as _d
    try:
        return (_d.date.fromisoformat(b[:10]) - _d.date.fromisoformat(a[:10])).days
    except ValueError:
        return None


def rcv(m: Dict[str, Any], z: float) -> Dict[str, float]:
    cvi = m["cvi_pct"] / 100
    cva = (m["cva_pct"] if m.get("cva_pct") is not None else 0.5 * m["cvi_pct"]) / 100
    if m.get("log_normal"):
        sigma = math.sqrt(math.log(cvi ** 2 + 1) + math.log(cva ** 2 + 1))
        return {"up_pct": 100 * (math.exp(z * math.sqrt(2) * sigma) - 1), "down_pct": 100 * (math.exp(-z * math.sqrt(2) * sigma) - 1)}
    r = 100 * z * math.sqrt(2) * math.sqrt(cvi ** 2 + cva ** 2)
    return {"up_pct": r, "down_pct": -r}


def compare(prev_path: Path, cur_path: Path) -> Dict[str, Any]:
    prev, cur = load_json(prev_path), load_json(cur_path)
    import os as _os
    if not _os.environ.get("LONGEVITY_SKILLS_HOME"):
        lib = cur.get("longevity_skills_home") or prev.get("longevity_skills_home")
        if lib:
            _os.environ["LONGEVITY_SKILLS_HOME"] = lib
    if prev.get("member", {}).get("id") != cur.get("member", {}).get("id"):
        raise LAError("the two twin snapshots belong to different members", EXIT_INPUT)
    B = _bv_index()
    kit = B["kit"]
    rows: List[Dict[str, Any]] = []

    def match(name: str) -> Optional[Dict[str, Any]]:
        for v in kit.name_variants(name):
            if v in B["idx"]:
                return B["idx"][v]
        return None

    from collections import Counter
    pcount = Counter(o["marker"] for o in prev["observations"])
    ccount = Counter(o["marker"] for o in cur["observations"])
    pmap = {o["marker"]: o for o in prev["observations"]}
    for o in cur["observations"]:
        p = pmap.get(o["marker"])
        if not p:
            continue
        if o.get("provenance_uncertain") or p.get("provenance_uncertain"):
            rows.append({"marker": o["marker"], "verdict": "not_judged", "why": "a value whose source is not confirmed as the member's"})
            continue
        if pcount[o["marker"]] > 1 or ccount[o["marker"]] > 1:
            rows.append({"marker": o["marker"], "verdict": "not_judged", "why": "the marker appears more than once in a snapshot"})
            continue
        a, b = _num(p["value"]), _num(o["value"])
        row = {"marker": o["marker"], "prev": p["value"], "cur": o["value"], "unit": o.get("unit", ""),
               "prev_t": p["t"], "cur_t": o["t"]}
        m = match(o["marker"])
        gap = _days(p["t"], o["t"])
        if gap is not None and gap < 1:
            row.update(verdict="not_judged", why="both snapshots are from the same day")
            rows.append(row)
            continue
        if a is None or b is None or a == 0:
            row["verdict"] = "not_judged"
            row["why"] = "censored (<, >), non-numeric or zero value"
        elif not m:
            row["verdict"] = "not_judged"
            row["why"] = "no within-person biological variation data for this marker"
        elif kit.normalize_unit(o.get("unit", "")) != kit.normalize_unit(p.get("unit", "")):
            row["verdict"] = "not_judged"
            row["why"] = f"units differ between visits ({p.get('unit')} vs {o.get('unit')})"
        elif m.get("min_retest_days") and gap is not None and gap < m["min_retest_days"]:
            row["verdict"] = "not_judged"
            row["why"] = f"only {gap} days apart; this marker needs at least {m['min_retest_days']} days between measurements"
        else:
            r = rcv(m, B["z"])
            ch = 100 * (b - a) / a
            if m["key"] == "vitd" and o["t"][5:7] != p["t"][5:7]:
                row["caveat"] = "25(OH)D 有季节性波动，不同季节的两次测量差异可能来自日照而非干预"
            row.update(change_pct=round(ch, 1), rcv_up_pct=round(r["up_pct"], 1), rcv_down_pct=round(r["down_pct"], 1),
                       bv_source=m["cvi_source"]["doi"],
                       verdict="increase_beyond_noise" if ch > r["up_pct"] else "decrease_beyond_noise" if ch < r["down_pct"] else "within_noise")
        rows.append(row)
    rprev = {r["id"]: r for r in prev["readouts"]}
    alerts = []
    for r in cur["readouts"]:
        p = rprev.get(r["id"])
        if not p:
            continue
        if r["id"].startswith("native.apoe") and str(p["value"]) != str(r["value"]):
            alerts.append(f"germline genotype changed between visits ({r['id']}: {p['value']} → {r['value']}): germline DNA does not change, so one of the samples is probably someone else's")
            rows.append({"readout": r["id"], "label_zh": r["label_zh"], "prev": p["value"], "cur": r["value"], "verdict": "identity_alert"})
            continue
        rows.append({"readout": r["id"], "label_zh": r["label_zh"], "prev": p["value"], "cur": r["value"], "unit": r.get("unit", ""),
                     "verdict": "not_judged", "why": "no published within-person noise model is loaded for this readout; show the two values, do not call it better or worse"})
    cva_note = "CVA 缺省按 EFLM 期望分析性能取 0.5×CVI（Sandberg 2024）；实验室给出实测 CVA 时应以实测为准"
    return {"prev": str(prev_path), "cur": str(cur_path), "z": B["z"], "alerts": alerts, "cva_note": cva_note, "rows": rows}
