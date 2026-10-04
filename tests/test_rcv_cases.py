"""The shared RCV cases (fixtures/rcv_cases.json) through twin compare. The coach's noise.py and LongPi's
reference.ts read the same file; each keeps its own implementation."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))

from lalib import twin  # noqa: E402

CASES = json.loads((Path(__file__).parent / "fixtures" / "rcv_cases.json").read_text(encoding="utf-8"))
LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="LONGEVITY_SKILLS_HOME not set")


def _compare(tmp: Path, case: dict) -> dict:
    base = {"schema": "la-twin/1", "member": {"id": "m"}, "longevity_skills_home": LIB, "interventions": [], "retest_plan": [], "readouts": []}
    for name, values, t in (("a", case["before"], case["prev_date"]), ("b", case["after"], case["cur_date"])):
        obs = [{"marker": case["marker"], "value": str(values[0]), "unit": case["unit"], "t": t}]
        (tmp / f"{name}.json").write_text(json.dumps({**base, "observations": obs}, ensure_ascii=False), encoding="utf-8")
    rows = twin.compare(tmp / "a.json", tmp / "b.json")["rows"]
    assert len(rows) == 1
    return rows[0]


@needs_lib
@pytest.mark.parametrize("case", [c for c in CASES["cases"] if c["expect"]["twin"] != "not_applicable"], ids=lambda c: c["id"])
def test_twin_compare_on_shared_rcv_cases(tmp_path, case):
    row = _compare(tmp_path, case)
    assert row["verdict"] == case["expect"]["twin"], row
    if "change_pct" in case and row["verdict"] != "not_judged":
        assert row["change_pct"] == case["change_pct"]
        assert [row["rcv_down_pct"], row["rcv_up_pct"]] == case["band_pct"]
