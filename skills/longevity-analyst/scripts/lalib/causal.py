"""Would lowering X help *this* member? Published Mendelian-randomisation estimates (EpiGraphDB) projected onto the
member's own value, a target value, the population SD and the member's baseline risk.

The projection is arithmetic on retrieved numbers: OR per SD from the MR record, the change in SD units from the
member's value to the target (SD from the NHANES same-sex, same-age-band stratum), and the odds conversion applied
to a registered baseline risk. Which exposure, outcome and target make sense is the agent's choice; that the MR
estimate transfers to this member is an assumption the report states.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Optional

from . import pubdata
from .common import EXIT_INPUT, LAError, data, load_json, now_iso, skillkit, write_json
from .reference import _band


def _member_value(st: Dict[str, Any], analyte: str) -> Dict[str, Any]:
    from .labnames import clean_value, name_keys, usable_rows
    kit = skillkit()
    a = data("ref_population.json")["labs_nhanes"]["analytes"].get(analyte)
    if not a:
        raise LAError(f"no population reference for {analyte!r}; one of {sorted(data('ref_population.json')['labs_nhanes']['analytes'])}", EXIT_INPUT)
    names = {kit.fold_name(n) for n in a["names"]}
    for r in usable_rows(st):
        if (name_keys(r["marker"]) | (name_keys(r["maps_to"]) if r.get("maps_to") else set())) & names \
                and kit.normalize_unit(r.get("unit", "")) == kit.normalize_unit(a["unit"]):
            return {"value": kit.parse_number(clean_value(r["value"])), "marker": r["marker"], "unit": a["unit"]}
    raise LAError(f"the member has no confirmed {a['label_zh']} row in {a['unit']}", EXIT_INPUT)


def project(st: Dict[str, Any], ws: Path, mr_ref: str, analyte: str, target: float, baseline_readout: str,
            exposure_match: Optional[str] = None) -> Dict[str, Any]:
    refs = pubdata.known_refs(ws)
    if mr_ref not in refs or not mr_ref.startswith("mr:"):
        raise LAError(f"{mr_ref!r} is not an MR record retrieved in this workspace (`la.py insights mr` first)", EXIT_INPUT)
    mr = refs[mr_ref]
    cur = _member_value(st, analyte)                       # also rejects an analyte with no population reference
    kit = skillkit()
    accepted = {kit.fold_name(t) for t in (data("trait_map.json").get("mr_exposure_traits") or {}).get(analyte, [])}
    if kit.fold_name(str(mr.get("exposure", ""))) not in accepted:
        if not exposure_match or len(exposure_match.strip()) < 8:
            raise LAError(f"MR exposure {mr.get('exposure')!r} is not a recognised name for {analyte}; if it measures the same quantity, "
                          "pass --exposure-match \"<reason>\"", EXIT_INPUT)
    m = st["member"]
    band = _band(int(m.get("age") or 0))
    stratum = data("ref_population.json")["labs_nhanes"]["analytes"][analyte]["strata"].get(f"{m.get('sex')}:{band}") if band else None
    if not stratum or not stratum.get("sd"):
        raise LAError(f"no population SD for {analyte} in {m.get('sex')} {band}", EXIT_INPUT)
    ro = {r["id"]: r for r in load_json(ws / "work" / "readouts.json")["readouts"]}
    op = ws / "work" / "organs" / "organ_readouts.json"
    if op.exists():
        ro.update({r["id"]: r for r in load_json(op)["readouts"]})
    if baseline_readout not in ro:
        raise LAError(f"baseline {baseline_readout!r} is not a readout (e.g. china-par-ascvd-risk.risk_10y_pct or organ.<o>.risk.<n>)", EXIT_INPUT)
    b = ro[baseline_readout]
    if b.get("unit") != "%" or not isinstance(b.get("value"), (int, float)) or not 0 < float(b["value"]) < 100:
        raise LAError(f"baseline {baseline_readout!r} must be a risk in percent (unit '%', 0-100); it is {b.get('value')!r} {b.get('unit')!r}", EXIT_INPUT)
    p0 = float(b["value"]) / 100.0
    d_sd = (target - cur["value"]) / stratum["sd"]

    def after(beta: float) -> float:
        orr = math.exp(beta * d_sd)
        return p0 * orr / (1 - p0 + p0 * orr)
    lo_b, hi_b = pubdata.se_ci(mr["b"], mr["se"])
    p1, pa, pb = after(mr["b"]), after(lo_b), after(hi_b)
    caveats = ["孟德尔随机化估计的是终生暴露差异的效应；成年后短期降低同样幅度，获益通常更小，推算值偏乐观",
               f"研究的结局是「{mr.get('outcome')}」，基线风险来自「{b.get('label_zh') or baseline_readout}」，两者定义不完全相同时只能作量级参考"]
    rec = {"ref": f"proj:{mr_ref}:{analyte}:{target}:{baseline_readout}", "mr_ref": mr_ref, "exposure": mr["exposure"], "outcome": mr["outcome"],
           "analyte": analyte, "member_value": cur["value"], "target": target, "unit": cur["unit"], "population_sd": stratum["sd"],
           "delta_sd": round(d_sd, 3), "baseline_readout": baseline_readout, "baseline_label": b.get("label_zh"), "baseline_risk": round(p0, 4),
           "risk_after": round(p1, 4), "risk_after_ci": [round(min(pa, pb), 4), round(max(pa, pb), 4)],
           "absolute_change": round(p1 - p0, 4), "exposure_match": exposure_match, "caveats_zh": caveats,
           "method_note": "odds scaled by exp(b × ΔSD) from the published MR estimate; assumes the exposure GWAS was in SD units and that the population causal effect applies to this member",
           "computed_at": now_iso()}
    path = ws / "work" / "insights" / "projections.json"
    allp = load_json(path)["projections"] if path.exists() else []
    allp = [x for x in allp if x["ref"] != rec["ref"]] + [rec]
    write_json(path, {"projections": allp})
    pubdata._cache(ws, "proj", rec["ref"], "local computation", [rec])
    return rec
