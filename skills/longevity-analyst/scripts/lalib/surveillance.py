"""Surveillance plan (WS-E): guideline rules (data/surveillance_rules.json, la-surveillance/1) evaluated on this member,
merged with the retests of the registered intervention plan into work/surveillance/plan.json (la-surveillance-plan/1).

A rule is applied only when its source page was fetched in this workspace (`la.py evidence fetch`). "once" means the
marker needs measuring once in adulthood: once it has been measured (now or in an earlier snapshot) it leaves the
retest list. Where a rule and an intervention retest the same marker the shorter interval wins and both reasons stay.
No rules: the plan holds the intervention retests only, as before.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import EXIT_INPUT, LAError, data, load_json, now_iso, skillkit, write_json

SCHEMA = "la-surveillance-plan/1"
WHEN_KEYS = {"age_min", "lab_measured", "readout", "concordance_class", "clinvar_gene"}


def rules() -> List[Dict[str, Any]]:
    try:
        doc = data("surveillance_rules.json")
    except FileNotFoundError:
        return []
    if doc.get("schema") != "la-surveillance/1":
        raise LAError("data/surveillance_rules.json is not la-surveillance/1", EXIT_INPUT)
    return doc.get("rules") or []


def check_rules(rs: List[Dict[str, Any]]) -> List[str]:
    """Structural problems of the rules file: unknown conditions, no source quote, no fetchable url."""
    out = []
    ids = set()
    for r in rs:
        rid = r.get("id", "?")
        if rid in ids:
            out.append(f"{rid}: duplicate rule id")
        ids.add(rid)
        if set(r.get("when") or {}) - WHEN_KEYS:
            out.append(f"{rid}: unknown when keys {sorted(set(r['when']) - WHEN_KEYS)}")
        iv = (r.get("action") or {}).get("interval")
        if not (iv == "once" or (isinstance(iv, int) and not isinstance(iv, bool) and 1 <= iv <= 520)):
            out.append(f"{rid}: action.interval must be weeks (1-520) or \"once\"")
        src = r.get("source") or {}
        if not (src.get("title") and src.get("quote") and re.match(r"^https?://[^\s?#]+$", str(src.get("url", "")))):
            out.append(f"{rid}: source needs title, quote and the fetched url")
        if r.get("reason") not in ("genetic", "guideline") or r.get("executor") not in ("member", "physician"):
            out.append(f"{rid}: reason genetic|guideline and executor member|physician")
    return out


def unfetched(ws: Path, rs: List[Dict[str, Any]]) -> List[str]:
    from .evidence import _cached_urls
    have = _cached_urls(ws) if (ws / "work" / "evidence").exists() else set()
    return [r["id"] for r in rs if (r.get("source") or {}).get("url") not in have]


def _names(analyte: str) -> List[str]:
    a = data("trait_map.json")["analytes"].get(analyte) or {}
    return [a.get("label_zh", analyte), *a.get("names", [])]


def _mentions(text: str, analyte: str) -> bool:
    kit = skillkit()
    t = kit.fold_name(re.sub(r"\{\{r:([^}|]*)[^}]*\}\}", r" \1 ", str(text or "")))
    return any(kit.fold_name(n) and kit.fold_name(n) in t for n in _names(analyte))


def measured_before(analyte: str, history: List[Path]) -> bool:
    kit = skillkit()
    want = {kit.fold_name(n) for n in _names(analyte)}
    for p in history:
        for o in load_json(p).get("observations") or []:
            if o.get("kind", "lab") == "lab" and kit.fold_name(str(o.get("marker", ""))) in want:
                return True
    return False


def _holds(when: Dict[str, Any], st: Dict[str, Any], ws: Path, labs: Dict[str, Any], history: List[Path]) -> bool:
    for k, v in (when or {}).items():
        if k == "age_min":
            if not isinstance(st["member"].get("age"), (int, float)) or st["member"]["age"] < v:
                return False
        elif k == "lab_measured":
            if v not in labs and not measured_before(v, history):
                return False
        elif k == "readout":
            ro = {}
            for f in (ws / "work" / "readouts.json", ws / "work" / "insights" / "insight_readouts.json"):
                if f.exists():
                    ro.update({r["id"]: r for r in load_json(f)["readouts"]})
            if v.get("id") not in ro or ("value" in v and ro[v["id"]].get("value") != v["value"]):
                return False
        elif k == "concordance_class":
            cp = ws / "work" / "insights" / "concordance.json"
            items = load_json(cp)["items"] if cp.exists() else []
            if not any(i["analyte"] == v.get("analyte") and i["class"] == v.get("class") for i in items):
                return False
        elif k == "clinvar_gene":
            gp = ws / "work" / "insights" / "genotype_phenotype.json"
            an = load_json(gp)["analytes"] if gp.exists() else {}
            if not any(x["gene"] == v and not x.get("carrier_only") for a in an.values()
                       for x in a["monogenic_scan"]["pathogenic_or_likely"]):
                return False
    return True


def plan(st: Dict[str, Any], ws: Path, history: Optional[List[Path]] = None) -> Dict[str, Any]:
    from .genomics import member_labs
    history = history or []
    rs = rules()
    bad = check_rules(rs)
    if bad:
        raise LAError("surveillance rules invalid:\n- " + "\n- ".join(bad), EXIT_INPUT)
    labs = member_labs(st)
    applied = [r for r in rs if _holds(r.get("when") or {}, st, ws, labs, history)]
    missing = unfetched(ws, applied)
    if missing:
        urls = sorted({r["source"]["url"] for r in applied if r["id"] in missing})
        raise LAError("fetch each rule's source first (`la.py evidence fetch <ws> --url <url>`): " + ", ".join(urls), EXIT_INPUT)
    items: List[Dict[str, Any]] = []
    for r in applied:
        a = r["action"]
        once = a["interval"] == "once"
        done = once and bool(a.get("analyte")) and (a["analyte"] in labs or measured_before(a["analyte"], history))
        items.append({"id": r["id"], "marker": "、".join(a.get("markers") or []), "analyte": a.get("analyte"),
                      "interval_weeks": None if once else int(a["interval"]), "once": once, "done": done,
                      "reason": r["reason"], "reasons": [r["reason"]], "rule_ids": [r["id"]], "interventions": [],
                      "executor": r["executor"], "note_zh": r.get("note_zh"),
                      "source": {k: r["source"].get(k) for k in ("title", "doi", "url", "quote", "locator") if r["source"].get(k)}})
    if (st.get("processed") or {}).get("variants"):
        items.append({"id": "S-germline", "marker": "胚系基因型", "analyte": None, "interval_weeks": None, "once": True, "done": True,
                      "reason": "genetic", "reasons": ["genetic"], "rule_ids": [], "interventions": [], "executor": "member",
                      "note_zh": "胚系基因型一生不变，测一次即可，以后复测不用再测。", "source": None})
    pp = ws / "work" / "intervene" / "plan.json"
    for it in (load_json(pp)["items"] if pp.exists() else []):
        rt = it.get("retest") or {}
        if not rt.get("what") or not rt.get("after_weeks"):
            continue
        wk = int(rt["after_weeks"])
        hit = next((x for x in items if x.get("analyte") and _mentions(rt["what"], x["analyte"])), None)
        if hit:                                            # same marker: shorter interval wins, both reasons kept
            hit["interventions"].append(it["id"])
            if "intervention" not in hit["reasons"]:
                hit["reasons"].append("intervention")
            if hit["once"] or hit["interval_weeks"] is None or wk < hit["interval_weeks"]:
                hit.update(interval_weeks=wk, once=False, done=False, reason="intervention")
            continue
        items.append({"id": f"R-{it['id']}", "marker": rt["what"], "analyte": None, "interval_weeks": wk, "once": False, "done": False,
                      "reason": "intervention", "reasons": ["intervention"], "rule_ids": [], "interventions": [it["id"]],
                      "executor": it.get("executor"), "note_zh": None, "source": None})
    doc = {"schema": SCHEMA, "generated_at": now_iso(), "history": [str(p) for p in history], "items": items}
    write_json(ws / "work" / "surveillance" / "plan.json", doc)
    return doc


def retests(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Items that still need a test: everything except a once-marker already measured."""
    return [x for x in doc.get("items", []) if not (x["once"] and x["done"])]
