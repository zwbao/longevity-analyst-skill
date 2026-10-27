"""WS-H: Chinese reference table (labs_cn) and the selection rule, offline."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import reference  # noqa: E402
from lalib import common as C  # noqa: E402
from test_insights import _ws  # noqa: E402


def test_labs_cn_table_is_sourced():
    cn = C.data("ref_population.json")["labs_cn"]
    assert cn["doi"] == "10.1186/s12944-019-1000-0" and cn["license"] == "CC BY 4.0" and cn["statistic"] == "mean_sd"
    assert set(cn["analytes"]) == {"tc", "ldl", "hdl"}
    assert cn["analytes"]["ldl"]["strata"]["male:40-44"] == {"mean": 2.90, "sd": 0.77}       # Table 2, men 40-44 LDL
    assert cn["analytes"]["tc"]["strata"]["female:70-150"] == {"mean": 5.75, "sd": 1.07}     # Table 2, women >= 70 total
    assert all(len(a["strata"]) == 14 for a in cn["analytes"].values())


@pytest.mark.parametrize("age,band", [(40, "40-44"), (49, "45-49"), (58, "55-59"), (70, "70-150"), (88, "70-150"), (39, None)])
def test_cn_bands(age, band):
    assert reference._band_cn(age) == band


def test_normal_percentile():
    assert reference._pct_normal(2.90, 2.90, 0.77) == {"pct": 50.0, "bound": None}
    assert reference._pct_normal(2.90 + reference.Z90 * 0.77, 2.90, 0.77)["pct"] == 90.0
    assert reference._pct_normal(0.1, 2.90, 0.77)["bound"] == "below"


def test_cn_resident_uses_labs_cn_and_falls_back_to_nhanes(tmp_path):
    ws, st = _ws(tmp_path)                                  # male 58, LDL-C 3.8 mmol/L
    st["labs"].append(dict(st["labs"][0], marker="肌酐", value="80", unit="μmol/L", confirm=dict(st["labs"][0]["confirm"])))
    from lalib import labnames
    for k, r in labnames.keyed(st).items():
        r["confirm"]["row_key"] = k
    st["member"]["answers"]["residence"] = "CN"
    pos = C.load_json(Path(reference.position(st, ws)["path"]))
    by = {x["analyte"]: x for x in pos["labs"]}
    assert by["ldl"]["reference"] == "cn" and by["ldl"]["stratum"] == "male 55-59" and by["ldl"]["median"] == 3.12
    if "creatinine" in by:
        assert by["creatinine"]["reference"] == "nhanes" and "回退" in by["creatinine"]["reference_label_zh"]
    st["member"]["answers"]["residence"] = "US"
    pos = C.load_json(Path(reference.position(st, ws)["path"]))
    assert {x["analyte"]: x for x in pos["labs"]}["ldl"]["reference"] == "nhanes"


def test_phenoage_note(tmp_path):
    ws, st = _ws(tmp_path)
    C.write_json(ws / "work" / "readouts.json", {"readouts": [{"id": "accelerated-biological-aging-risk.phenoage", "label_zh": "表型年龄",
                                                             "value": 55.0, "unit": "a", "kind": "computed"}]})
    pos = C.load_json(Path(reference.position(st, ws)["path"]))
    assert pos["phenoage_note_zh"] == "这个模型是用美国人群数据拟合的。"
