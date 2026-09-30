"""Where the member stands: lab percentiles in a same-sex, same-age-band reference population (NHANES), GMHI against
healthy and non-healthy metagenomes (overall and East Asian healthy), and the several "ages" side by side.

Deterministic lookups in data/ref_population.json. Which reference fits a Chinese member, and what a gap between
ages means, is written by the agent, not here.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import data, load_json, now_iso, skillkit, write_json


def _pct_of(v: float, pcts: List[float], vals: List[float]) -> Dict[str, Any]:
    ties = [i for i, x in enumerate(vals) if x == v]
    if len(ties) > 1:                  # many people sit at the same value (usually a detection limit): report the span
        lo, hi = pcts[ties[0]], pcts[ties[-1]]
        return {"pct": round((lo + hi) / 2, 1), "pct_low": lo, "pct_high": hi, "bound": "tie"}
    if v <= vals[0]:
        return {"pct": pcts[0], "bound": "below"}
    if v >= vals[-1]:
        return {"pct": pcts[-1], "bound": "above"}
    for i in range(1, len(vals)):
        if v <= vals[i]:
            lo, hi = vals[i - 1], vals[i]
            f = (v - lo) / (hi - lo) if hi > lo else 0.0
            return {"pct": round(pcts[i - 1] + f * (pcts[i] - pcts[i - 1]), 1), "bound": None}
    return {"pct": pcts[-1], "bound": "above"}


def _band(age: int) -> Optional[str]:
    for lo, hi in ((20, 29), (30, 39), (40, 49), (50, 59), (60, 69), (70, 79), (80, 150)):
        if lo <= age <= hi:
            return f"{lo}-{hi}"
    return None


def position(st: Dict[str, Any], ws: Path) -> Dict[str, Any]:
    from .labnames import clean_value, name_keys, usable_rows
    kit = skillkit()
    ref = data("ref_population.json")
    pcts = ref["pcts"]
    m = st["member"]
    sex, age = m.get("sex"), m.get("age")
    band = _band(int(age)) if age else None
    labs_out, readouts = [], []
    rows = usable_rows(st)
    for key, a in ref["labs_nhanes"]["analytes"].items():
        if not a["names"] or not sex or not band:
            continue
        stratum = a["strata"].get(f"{sex}:{band}")
        if not stratum:
            continue
        names = {kit.fold_name(n) for n in a["names"]}
        for r in rows:
            if not ((name_keys(r["marker"]) | (name_keys(r["maps_to"]) if r.get("maps_to") else set())) & names):
                continue
            unit_ok = kit.normalize_unit(r.get("unit", "")) == kit.normalize_unit(a["unit"]) or not a["unit"]
            if not unit_ok:
                continue                              # a unit the table was not built in is never guessed
            try:
                v = kit.parse_number(clean_value(r["value"]))
            except ValueError:
                continue
            p = _pct_of(v, pcts, stratum["pct"])
            labs_out.append({"analyte": key, "label_zh": a["label_zh"], "marker": r["marker"], "value": v, "unit": a["unit"],
                             "stratum": f"{sex} {band}", "n": stratum["n"], **p,
                             "median": stratum["pct"][pcts.index(50)], "p90": stratum["pct"][pcts.index(90)]})
            readouts.append({"id": f"ref.{key}.pct_nhanes", "label_zh": f"{a['label_zh']}在美国同性别同龄段人群中的百分位",
                             "value": p["pct"], "unit": "%", "kind": "population_position", "method": "native.population_position",
                             "reference": ref["labs_nhanes"]["source"], "systems": []})
            break
    gm = None
    ro = {r["id"]: r for r in load_json(ws / "work" / "readouts.json")["readouts"]}
    if "native.gut.gmhi" in ro:
        v = float(ro["native.gut.gmhi"]["value"])
        g = ref["gmhi"]
        gm = {"value": v, "healthy": _pct_of(v, g["pcts"], g["healthy"]["pct"]),
              "nonhealthy": _pct_of(v, g["pcts"], g["nonhealthy"]["pct"]), "n_healthy": g["healthy"]["n"],
              "n_nonhealthy": g["nonhealthy"]["n"], "source": g["source"], "doi": g["doi"]}
        if "healthy_east_asia" in g:
            gm["healthy_east_asia"] = _pct_of(v, g["pcts"], g["healthy_east_asia"]["pct"])
            gm["n_healthy_east_asia"] = g["healthy_east_asia"]["n"]
            readouts.append({"id": "ref.gmhi.pct_healthy_east_asia", "label_zh": "GMHI 在东亚健康人群中的百分位",
                             "value": gm["healthy_east_asia"]["pct"], "unit": "%", "kind": "population_position",
                             "method": "native.population_position", "reference": g["source"], "systems": []})
        readouts.append({"id": "ref.gmhi.pct_healthy", "label_zh": "GMHI 在全部健康人群中的百分位", "value": gm["healthy"]["pct"],
                         "unit": "%", "kind": "population_position", "method": "native.population_position", "reference": g["source"], "systems": []})
    # the several ages, side by side (values as computed; a clock flagged low-coverage keeps its flag)
    ages = [{"id": "member.age", "label_zh": "实际年龄", "value": age, "kind": "fact"}] if age else []
    op = ws / "work" / "organs" / "organ_readouts.json"
    organ = load_json(op)["readouts"] if op.exists() else []
    for r in list(ro.values()) + organ:
        if r.get("unit") == "a" and isinstance(r.get("value"), (int, float)) and "advance" not in r["id"] and r.get("kind") != "implausible":
            ages.append({"id": r["id"], "label_zh": r.get("label_zh"), "value": r["value"], "kind": r.get("kind"),
                         "minus_age": round(r["value"] - age, 1) if age else None, "low": r.get("low"), "high": r.get("high")})
    clocks = [a["value"] for a in ages if str(a["id"]).startswith("epiage.") and a.get("kind") in ("computed", "computed_low_coverage")]
    summary = {"epigenetic_clock_range": [min(clocks), max(clocks)] if clocks else None,
               "epigenetic_clock_median": sorted(clocks)[len(clocks) // 2] if clocks else None}
    out = {"generated_at": now_iso(), "labs": labs_out, "gmhi": gm, "ages": ages, "ages_summary": summary,
           "reference_note": ref["note"]}
    p = ws / "work" / "insights" / "positions.json"
    write_json(p, out)
    return {"path": str(p), "labs": len(labs_out), "gmhi": bool(gm), "ages": len(ages), "readouts": readouts}
