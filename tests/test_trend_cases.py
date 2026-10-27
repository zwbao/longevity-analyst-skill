"""The shared trend cases (fixtures/trend_cases.json, trend-cases/1) through `twin series`. The coach's noise.py
trend and the LongPi site read the same file; each keeps its own implementation."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))

from lalib import series  # noqa: E402

CASES = json.loads((Path(__file__).parent / "fixtures" / "trend_cases.json").read_text(encoding="utf-8"))
LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
pytestmark = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="LONGEVITY_SKILLS_HOME not set")


def _snapshots(case: dict):
    """One la-twin/2 snapshot per distinct date, observations in fixture order."""
    by = {}
    for o in case["series"]:
        by.setdefault(o["t"], []).append({"t": o["t"], "t_source": "row", "kind": "lab", "marker": case["marker"],
                                          "value": str(o["value"]), "unit": o.get("unit", case["unit"]),
                                          "provenance_uncertain": bool(o.get("provenance_uncertain"))})
    snaps = [{"schema": "la-twin/2", "member": {"id": "synthetic-m"}, "snapshot": {"t": t}, "observations": obs,
              "readouts": [], "exposures": []} for t, obs in sorted(by.items())]
    return snaps, [f"s{i}.json" for i in range(len(snaps))]


def _marker(case):
    res = series.run(*_snapshots(case))
    hit = [m for m in res["markers"] if m["marker_key"] == case["marker"]]
    assert len(hit) == 1, res["markers"]
    return hit[0]


@pytest.mark.parametrize("case", CASES["cases"], ids=lambda c: c["id"])
def test_trend_case(case):
    m = _marker(case)
    exp = case["expect"]["analyst"]
    assert [p["verdict"] for p in m["pairs"]] == exp["pairs"]
    assert m["baseline_vs_latest"]["verdict"] == exp["baseline_vs_latest"]
    assert m["trend"]["verdict"] == exp["trend"]
    assert m["summary_verdict"] == exp["summary_verdict"]
    if case["expect"].get("trend_not_judged_reason"):
        assert m["trend"]["why"] == case["expect"]["trend_not_judged_reason"]
    codes = {c["code"] for c in m["caveats"]}
    assert set(case["expect"].get("caveats", [])) <= codes
    n = case["numbers"]
    assert len(m["points"]) == n["points"]
    assert abs(series.band(series.match(series.noise(), case["marker"]), CASES["z"])[1] - n["single_rcv_up_pct"]) <= 0.1
    assert abs(m["baseline_vs_latest"]["change_pct"] - n["first_vs_latest_pct"]) <= 0.1
    if "slope_per_year" in n:
        assert abs(m["trend"]["slope_per_year"] - n["slope_per_year"]) <= 0.01
        assert abs(m["trend"]["z"] - n["z_slope"]) <= 0.1
        assert abs(m["trend"]["annual_change_pct"] - n["annual_change_pct"]) <= 0.1


def test_fixture_constants_match_the_library():
    N = series.noise()
    assert N["z"] == CASES["z"] and N["z_trend"] == CASES["z_trend"]
