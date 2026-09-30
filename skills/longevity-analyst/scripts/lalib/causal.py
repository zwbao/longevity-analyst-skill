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
from typing import Any, Dict

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


def project(st: Dict[str, Any], ws: Path, mr_ref: str, analyte: str, target: float, baseline_readout: str) -> Dict[str, Any]:
    refs = pubdata.known_refs(ws)
    if mr_ref not in refs or not mr_ref.startswith("mr:"):
        raise LAError(f"{mr_ref!r} is not an MR record retrieved in this workspace (`la.py insights mr` first)", EXIT_INPUT)
    mr = refs[mr_ref]
    m = st["member"]
    band = _band(int(m.get("age") or 0))
    stratum = data("ref_population.json")["labs_nhanes"]["analytes"][analyte]["strata"].get(f"{m.get('sex')}:{band}") if band else None
    if not stratum or not stratum.get("sd"):
        raise LAError(f"no population SD for {analyte} in {m.get('sex')} {band}", EXIT_INPUT)
    cur = _member_value(st, analyte)
    ro = {r["id"]: r for r in load_json(ws / "work" / "readouts.json")["readouts"]}
    op = ws / "work" / "organs" / "organ_readouts.json"
    if op.exists():
        ro.update({r["id"]: r for r in load_json(op)["readouts"]})
    if baseline_readout not in ro:
        raise LAError(f"baseline {baseline_readout!r} is not a readout (e.g. china-par-ascvd-risk.risk_10y_pct or organ.<o>.risk.<n>)", EXIT_INPUT)
    b = ro[baseline_readout]
    p0 = float(b["value"]) / (100.0 if b.get("unit") in ("%",) else 1.0)
    if not 0 < p0 < 1:
        raise LAError(f"baseline risk {p0} is not a probability", EXIT_INPUT)
    d_sd = (target - cur["value"]) / stratum["sd"]
    def after(beta: float) -> float:
        orr = math.exp(beta * d_sd)
        return p0 * orr / (1 - p0 + p0 * orr)
    lo_b, hi_b = pubdata.se_ci(mr["b"], mr["se"])
    p1, pa, pb = after(mr["b"]), after(lo_b), after(hi_b)
    rec = {"ref": f"proj:{mr_ref}:{analyte}:{target}", "mr_ref": mr_ref, "exposure": mr["exposure"], "outcome": mr["outcome"],
           "analyte": analyte, "member_value": cur["value"], "target": target, "unit": cur["unit"], "population_sd": stratum["sd"],
           "delta_sd": round(d_sd, 3), "baseline_readout": baseline_readout, "baseline_risk": round(p0, 4),
           "risk_after": round(p1, 4), "risk_after_ci": [round(min(pa, pb), 4), round(max(pa, pb), 4)],
           "absolute_change": round(p1 - p0, 4), "method_note": "odds scaled by exp(b × ΔSD) from the published MR estimate; assumes the exposure GWAS was in SD units and that the population causal effect applies to this member",
           "computed_at": now_iso()}
    path = ws / "work" / "insights" / "projections.json"
    allp = load_json(path)["projections"] if path.exists() else []
    allp = [x for x in allp if x["ref"] != rec["ref"]] + [rec]
    write_json(path, {"projections": allp})
    pubdata._cache(ws, "proj", rec["ref"], "local computation", [rec])
    return rec
