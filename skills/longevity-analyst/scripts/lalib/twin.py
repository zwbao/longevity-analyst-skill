"""Digital twin: a provenance-carrying longitudinal state of one member (la-twin/2).

twin.json holds every observation (lab rows with their own sampling date, method readouts) with its source, the
intervention register with the exposures a host recorded (la-exposures/1), the gene-lab concordance, the germline
fingerprint used to check that two visits are the same person, the surveillance plan and the retest plan with
whether each retest can be judged. `series` (lalib/series.py) decides whether changes across snapshots exceed
normal within-person variation; `compare` is its two-point case. Where no noise model exists the change is
reported as "not judged", never as improved or worsened.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import series as S
from .labnames import usable_rows
from .common import (EXIT_INPUT, LAError, SKILL_DIR, load_json, now_iso, sha256_file, skillkit, write_json)

SCHEMA = "la-twin/2"
EXPOSURE_SCHEMA = SKILL_DIR / "schema" / "la-exposures-1.schema.json"
JUDGEABLE_FALSE_ZH = "这一项复测时只能并排看两次结果，不能判断是不是真的变化。"


def _kit():
    return skillkit()


# ------------------------------------------------------------ exposures (la-exposures/1)
def _schema_errors(doc: Any, sch: Dict[str, Any], path: str = "$") -> List[str]:
    """The subset of JSON Schema the shared schemas use: type, const, enum, pattern, min/maxLength, minimum,
    required, properties, additionalProperties: false, items."""
    errs: List[str] = []
    types = sch.get("type")
    if types:
        tl = types if isinstance(types, list) else [types]
        ok = {"object": isinstance(doc, dict), "array": isinstance(doc, list), "string": isinstance(doc, str),
              "integer": isinstance(doc, int) and not isinstance(doc, bool), "null": doc is None,
              "number": isinstance(doc, (int, float)) and not isinstance(doc, bool), "boolean": isinstance(doc, bool)}
        if not any(ok.get(t) for t in tl):
            return [f"{path}: expected {'/'.join(tl)}"]
    if "const" in sch and doc != sch["const"]:
        errs.append(f"{path}: must be {sch['const']!r}")
    if "enum" in sch and doc not in sch["enum"]:
        errs.append(f"{path}: {doc!r} is not one of {sch['enum']}")
    if isinstance(doc, str):
        if "pattern" in sch and not re.search(sch["pattern"], doc):
            errs.append(f"{path}: {doc!r} does not match {sch['pattern']}")
        if len(doc) < sch.get("minLength", 0) or len(doc) > sch.get("maxLength", 10 ** 9):
            errs.append(f"{path}: length {len(doc)} outside [{sch.get('minLength', 0)}, {sch.get('maxLength', '∞')}]")
    if isinstance(doc, (int, float)) and not isinstance(doc, bool) and "minimum" in sch and doc < sch["minimum"]:
        errs.append(f"{path}: below {sch['minimum']}")
    if isinstance(doc, dict):
        for k in sch.get("required", []):
            if k not in doc:
                errs.append(f"{path}: missing {k}")
        props = sch.get("properties", {})
        if sch.get("additionalProperties") is False:
            errs += [f"{path}: unknown field {k!r}" for k in doc if k not in props]
        for k, v in doc.items():
            if k in props:
                errs += _schema_errors(v, props[k], f"{path}.{k}")
    if isinstance(doc, list) and "items" in sch:
        for i, v in enumerate(doc):
            errs += _schema_errors(v, sch["items"], f"{path}[{i}]")
    return errs


def load_exposures(path: Path, member_id: Optional[str]) -> Dict[str, Any]:
    try:
        doc = load_json(path)
    except (OSError, ValueError) as e:
        raise LAError(f"{path} is not JSON ({e})", EXIT_INPUT)
    errs = _schema_errors(doc, load_json(EXPOSURE_SCHEMA))
    if errs:
        raise LAError("exposures file rejected (schema la-exposures/1):\n- " + "\n- ".join(errs[:20]), EXIT_INPUT)
    warnings = []
    if member_id and doc["member"] != member_id:
        warnings.append(f"exposures are for member {doc['member']!r}, the workspace is {member_id!r}; check it is the same person")
    return {"doc": doc, "warnings": warnings}


# ------------------------------------------------------------ retests
def retest_check(N: Dict[str, Any], what: str, targets: List[str], after_weeks: Optional[int],
                 labels: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Can this retest be judged next time? `what` resolves against the biological-variation table or the readout
    noise table; failing that, every target must resolve. Warnings, never a refusal."""
    labels = labels or {}
    text = re.sub(r"\{\{r:([^|}]+)(?:\|[^}]*)?\}\}", lambda m: labels.get(m.group(1), m.group(1)), str(what or ""))
    text = re.sub(r"\{\{n:([^}]*)\}\}", r"\1", text)
    hit = S.resolve(N, text)
    if not hit and targets:
        hits = [S.resolve(N, labels.get(t, t)) or S.resolve(N, t) for t in targets]
        hit = hits[0] if all(hits) and len({h["key"] for h in hits}) == 1 else None
    out: Dict[str, Any] = {"judgeable": bool(hit), "warnings": []}
    if hit:
        out["noise"] = hit
        mind = hit.get("min_retest_days") or (S.DEFAULT_MIN_DAYS if hit["table"] == "biological_variation" else None)
        if after_weeks and mind and after_weeks * 7 < mind:
            out["warnings"].append(f"retest after {after_weeks} weeks is sooner than the {mind} days this marker needs between measurements")
    else:
        out["warnings"].append(f"{text or '(empty)'!r} has no within-person noise model: the retest can only be shown side by side, not judged")
    return out


def plan_retests(ws_root: Path, plan: Dict[str, Any]) -> List[Dict[str, Any]]:
    """judgeable / warnings for every plan item with a retest (intervene register and twin build share this)."""
    labels = {}
    for f in (ws_root / "work" / "readouts.json", ws_root / "work" / "insights" / "insight_readouts.json",
              ws_root / "work" / "organs" / "organ_readouts.json"):
        if f.exists():
            labels.update({r["id"]: r.get("label_zh", r["id"]) for r in load_json(f)["readouts"]})
    N = S.noise()
    out = []
    for i in plan.get("items", []):
        rt = i.get("retest") or {}
        if not rt.get("what"):
            continue
        chk = retest_check(N, rt.get("what"), i.get("targets") or [], rt.get("after_weeks"), labels)
        out.append({"item": i["id"], "what": rt.get("what"), "after_weeks": rt.get("after_weeks"), **chk})
    return out


def _surveillance(ws_root: Path) -> List[Dict[str, Any]]:
    p = ws_root / "work" / "surveillance" / "plan.json"
    if not p.exists():
        return []
    try:
        return list(load_json(p).get("items") or [])
    except (OSError, ValueError, AttributeError):
        return []


def _fingerprint(ws_root: Path) -> Optional[Dict[str, Any]]:
    for p, reused in ((ws_root / "work" / "insights" / "fingerprint.json", False),
                      (ws_root / "work" / "reuse" / "germline" / "fingerprint.json", True)):
        if p.exists():
            fp = load_json(p)
            out = {"sites": fp.get("sites", len(fp.get("calls") or {})), "sha256": fp.get("sha256"), "calls": fp.get("calls") or {}}
            if reused:
                out["reused_from"] = (load_json(ws_root / "work" / "reuse" / "germline.json") or {}).get("from") \
                    if (ws_root / "work" / "reuse" / "germline.json").exists() else None
            return out
    return None


def _concordance(ws_root: Path, t: str) -> List[Dict[str, Any]]:
    p = ws_root / "work" / "insights" / "concordance.json"
    if not p.exists():
        return []
    doc = load_json(p)
    return [{"t": t, **{k: v for k, v in it.items() if k != "basis"}, "gwas_retrieved_on": doc.get("retrieved_on"),
             "trait_map_version": doc.get("trait_map_version")} for it in doc.get("items") or []]


def build(st: Dict[str, Any], ws_root: Path, exposures: Optional[Path] = None) -> Dict[str, Any]:
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
    warnings: List[str] = []
    exp_items: List[Dict[str, Any]] = []
    exp_dst = ws_root / "work" / "twin" / "exposures.json"
    if exposures is not None:
        ex = load_exposures(Path(exposures), st["member"].get("id"))
        warnings += ex["warnings"]
        exp_items = ex["doc"]["items"]
        write_json(exp_dst, ex["doc"])               # the bound copy: review trace hashes it
        st["exposures"] = {"source_path": str(Path(exposures).resolve()), "path": str(exp_dst),
                           "sha256": sha256_file(exp_dst, limit=None), "source": ex["doc"]["source"],
                           "generated_at": ex["doc"]["generated_at"], "items": len(exp_items)}
    elif st.get("exposures") and exp_dst.exists() and sha256_file(exp_dst, limit=None) == st["exposures"].get("sha256"):
        exp_items = load_json(exp_dst)["items"]       # a rebuild keeps the exposures registered earlier
        warnings.append("kept the exposures registered with the previous build; pass --exposures to replace them")
    else:
        exp_dst.unlink(missing_ok=True)
        st.pop("exposures", None)
    by_item = {e["plan_item"]: e for e in exp_items}
    from .labnames import collected_at
    obs = []
    for l in usable_rows(st, history=True):
        ct, src = collected_at(l, st)
        obs.append({"t": ct or t, "t_source": src if ct else "workspace", "kind": "lab", "marker": l["marker"],
                    "value": l["value"], "unit": l.get("unit", ""), "ref_range": l.get("ref_range", ""),
                    "source_file": l.get("source_file"), "page": l.get("page", ""),
                    "provenance_uncertain": l.get("source_file") in unc})
    judge = {r["item"]: r for r in plan_retests(ws_root, plan)}
    twin = {
        "schema": SCHEMA,
        "organ_estimates": [{"t": t, **{k: r.get(k) for k in ("id", "organ", "label_zh", "value", "low", "high", "unit",
                                                             "horizon_years", "confidence", "kind")}} for r in organ_ro],
        "longevity_skills_home": st.get("longevity_skills_home"),
        "member": {k: st["member"].get(k) for k in ("id", "age", "sex")},
        "identity": st.get("identity") or {"answer": "not_checked"},
        "snapshot": {"t": t, "built_at": now_iso(), "workspace": str(ws_root)},
        "observations": obs,
        "readouts": [{"t": t, "id": r["id"], "label_zh": r["label_zh"], "value": r["value"], "unit": r.get("unit", ""),
                      "kind": r.get("kind"), "method": r.get("method"), "systems": r.get("systems", []),
                      "source": r.get("source"), "provenance_uncertain": bool(r.get("provenance_uncertain"))} for r in readouts],
        "data_files": [{"id": f["id"], "name": f["name"], "kind": f["kind"], "sha256": f["sha256"], "status": f.get("status"),
                        "excluded": bool(f.get("excluded")), "provenance_uncertain": bool(f.get("provenance_uncertain"))}
                       for f in st["files"]],
        "pipelines": [{"run_id": p["run_id"], "pipeline": p["pipeline"], "revision": p["revision"], "status": p["status"]} for p in st["pipelines"]],
        "interventions": [{"id": i["id"], "category": i["category"], "action_zh": i["action_zh"], "targets": i["targets"],
                           "executor": i["executor"],
                           "status": by_item[i["id"]]["status"] if i["id"] in by_item else "proposed",
                           "since": (by_item[i["id"]].get("status_since") or by_item[i["id"]].get("adopted_at")) if i["id"] in by_item else None}
                          for i in plan.get("items", [])],
        "retest_plan": [{"intervention": i["id"], "what": i["retest"].get("what"), "after_weeks": i["retest"].get("after_weeks"),
                         **({"judgeable": judge[i["id"]]["judgeable"]} if i["id"] in judge else {})}
                        for i in plan.get("items", []) if i.get("retest")],
        "exposures": exp_items,
        "concordance": _concordance(ws_root, t),
        "genotype_fingerprint": _fingerprint(ws_root),
        "surveillance": _surveillance(ws_root),
    }
    write_json(ws_root / "work" / "twin" / "twin.json", twin)
    st["twin"] = {"path": str(ws_root / "work" / "twin" / "twin.json"), "built_at": now_iso()}
    res = {"path": st["twin"]["path"], "schema": SCHEMA, "observations": len(twin["observations"]),
           "observation_dates": sorted({o["t"] for o in obs}), "readouts": len(twin["readouts"]),
           "interventions": len(twin["interventions"]), "exposures": len(exp_items),
           "concordance": len(twin["concordance"]), "fingerprint_sites": (twin["genotype_fingerprint"] or {}).get("sites"),
           "surveillance": len(twin["surveillance"])}
    if warnings:
        res["warnings"] = warnings
    return res


def rcv(m: Dict[str, Any], z: float) -> Dict[str, float]:
    lo, hi = S.band(m, z)
    return {"up_pct": hi, "down_pct": lo}


def _current(snap: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Each marker's observations at its latest date in this snapshot (a la-twin/2 snapshot may carry history rows)."""
    last: Dict[str, str] = {}
    for o in snap["observations"]:
        t = str(o.get("t") or "")
        if t >= last.get(o["marker"], ""):
            last[o["marker"]] = t
    return [o for o in snap["observations"] if str(o.get("t") or "") == last[o["marker"]]]


def compare(prev_path: Path, cur_path: Path) -> Dict[str, Any]:
    """Two snapshots: the two-point case of `series` (same pair judgment), in the original row format."""
    prev, cur = load_json(prev_path), load_json(cur_path)
    import os as _os
    if not _os.environ.get("LONGEVITY_SKILLS_HOME"):
        lib = cur.get("longevity_skills_home") or prev.get("longevity_skills_home")
        if lib:
            _os.environ["LONGEVITY_SKILLS_HOME"] = lib
    if prev.get("member", {}).get("id") != cur.get("member", {}).get("id"):
        raise LAError("the two twin snapshots belong to different members", EXIT_INPUT)
    N = S.noise()
    kit = N["kit"]
    rows: List[Dict[str, Any]] = []
    from collections import Counter
    pobs, cobs = _current(prev), _current(cur)
    pcount = Counter(o["marker"] for o in pobs)
    ccount = Counter(o["marker"] for o in cobs)
    pmap = {o["marker"]: o for o in pobs}
    seen = set()
    for o in cobs:
        if o["marker"] in seen:
            continue
        seen.add(o["marker"])
        m = S.match(N, o["marker"])
        p = pmap.get(o["marker"]) or (next((x for x in pobs if S.match(N, x["marker"]) is m), None) if m else None)
        if not p:
            continue
        if o.get("provenance_uncertain") or p.get("provenance_uncertain"):
            rows.append({"marker": o["marker"], "verdict": "not_judged", "why": "a value whose source is not confirmed as the member's"})
            continue
        if pcount[p["marker"]] > 1 or ccount[o["marker"]] > 1:
            rows.append({"marker": o["marker"], "verdict": "not_judged", "why": "the marker appears more than once in a snapshot"})
            continue
        a, b = S.num(p["value"]), S.num(o["value"])
        row = {"marker": o["marker"], "prev": p["value"], "cur": o["value"], "unit": o.get("unit", ""),
               "prev_t": p["t"], "cur_t": o["t"]}
        gap = S.days(p["t"], o["t"])
        if gap is not None and gap < 1:
            row.update(verdict="not_judged", why="both snapshots are from the same day")
            rows.append(row)
            continue
        fa = fb = 1.0
        if m and kit.normalize_unit(o.get("unit", "")) != kit.normalize_unit(p.get("unit", "")):
            fa, fb = S.to_unit(m, p.get("unit", ""), kit), S.to_unit(m, o.get("unit", ""), kit)
        if a is None or b is None or a == 0:
            row["verdict"] = "not_judged"
            row["why"] = "censored (<, >), non-numeric or zero value"
        elif not m:
            row["verdict"] = "not_judged"
            row["why"] = "no within-person biological variation data for this marker"
        elif fa is None or fb is None:
            row["verdict"] = "not_judged"
            row["why"] = f"units differ between visits ({p.get('unit')} vs {o.get('unit')})"
        else:
            pr = S.pair({"t": p["t"], "value": a * fa, "k": 1}, {"t": o["t"], "value": b * fb, "k": 1}, m, N["z"])
            if pr["verdict"] == "not_judged":
                row.update(verdict="not_judged",
                           why=f"only {gap} days apart; this marker needs at least {pr.get('min_retest_days')} days between measurements")
                rows.append(row)
                continue
            if m["key"] == "vitd" and o["t"][5:7] != p["t"][5:7]:
                row["caveat"] = "25(OH)D 有季节性波动，不同季节的两次测量差异可能来自日照而非干预"
            row.update(change_pct=pr["change_pct"], rcv_up_pct=pr["band_pct"][1], rcv_down_pct=pr["band_pct"][0],
                       bv_source=m["cvi_source"]["doi"], verdict=pr["verdict"])
        rows.append(row)
    rprev = {r["id"]: r for r in prev["readouts"]}
    alerts = []
    rz = (N.get("readouts") or {}).get("z") or N["z"]
    for r in cur["readouts"]:
        p = rprev.get(r["id"])
        if not p:
            continue
        if r["id"].startswith("native.apoe") and str(p["value"]) != str(r["value"]):
            alerts.append(f"germline genotype changed between visits ({r['id']}: {p['value']} → {r['value']}): germline DNA does not change, so one of the samples is probably someone else's")
            rows.append({"readout": r["id"], "label_zh": r["label_zh"], "prev": p["value"], "cur": r["value"], "verdict": "identity_alert"})
            continue
        row = {"readout": r["id"], "label_zh": r["label_zh"], "prev": p["value"], "cur": r["value"], "unit": r.get("unit", "")}
        row.update(S.readout_pair(p, r, S.match_readout(N, r["id"], r.get("label_zh", "")), rz))
        rows.append(row)
    fm = S.fingerprint_mismatch(S._calls(prev), S._calls(cur))
    if fm and fm["alert"]:
        alerts.append(S.IDENTITY_ZH)
        rows.append({"readout": "genotype_fingerprint", "label_zh": "胚系基因指纹", "overlap": fm["overlap"],
                     "mismatches": fm["mismatches"], "verdict": "identity_alert"})
    cva_note = "CVA 缺省按 EFLM 期望分析性能取 0.5×CVI（Sandberg 2024）；实验室给出实测 CVA 时应以实测为准"
    return {"prev": str(prev_path), "cur": str(cur_path), "z": N["z"], "alerts": alerts, "cva_note": cva_note, "rows": rows}


# ------------------------------------------------------------ series in a workspace
def provisional(st: Dict[str, Any], ws_root: Path) -> Dict[str, Any]:
    """This workspace's confirmed lab rows as an in-memory snapshot, before `twin build` (the labs-only fast route
    uses it to see which systems changed). Never written to disk."""
    from .labnames import collected_at
    t = st["member"].get("sample_date") or now_iso()[:10]
    unc = {f["id"] for f in st["files"] if f.get("provenance_uncertain")}
    obs = []
    for l in usable_rows(st, history=True):
        ct, src = collected_at(l, st)
        obs.append({"t": ct or t, "t_source": src or "workspace", "kind": "lab", "marker": l["marker"], "value": l["value"],
                    "unit": l.get("unit", ""), "provenance_uncertain": l.get("source_file") in unc})
    return {"schema": SCHEMA, "provisional": True, "member": {k: st["member"].get(k) for k in ("id", "age", "sex")},
            "longevity_skills_home": st.get("longevity_skills_home"), "snapshot": {"t": t, "workspace": str(ws_root)},
            "observations": obs, "readouts": [], "exposures": []}


def series_problems(ws_root: Path, st: Dict[str, Any]) -> List[str]:
    """work/twin/series.json is rendered into the report: it must be what `twin series --ws` wrote."""
    p = ws_root / "work" / "twin" / "series.json"
    rec = (st.get("twin") or {}).get("series") or {}
    if p.exists() and (not rec or sha256_file(p, limit=None) != rec.get("sha256")):
        return ["work/twin/series.json does not match what `la.py twin series --ws` wrote; run it again"]
    if rec and not p.exists():
        return ["work/twin/series.json is missing; run `la.py twin series --ws` again"]
    return []
