"""la-exposures/1 (WS-C.1): validation against the shared schema, storage in the twin, interventions filled from it.
All exposure records here are SYNTHETIC."""
from __future__ import annotations

import copy
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import common as C  # noqa: E402
from lalib import report, twin  # noqa: E402

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")
SCHEMA = json.loads((ROOT / "skills" / "longevity-analyst" / "schema" / "la-exposures-1.schema.json").read_text())

GOOD = {"schema": "la-exposures/1", "member": "m", "source": "coach_member_file", "generated_at": "2026-09-20T10:00:00+08:00",
        "items": [{"plan_item": "I3", "origin_snapshot": "2025-09-02", "label_zh": "工作日晚饭用限盐勺", "category": "diet",
                   "targets": ["收缩压"], "status": "adopted", "adopted_at": "2026-08-15", "status_since": "2026-08-15",
                   "weekly_checkins": [{"week_start": "2026-08-17", "count": 4}, {"week_start": "2026-08-24", "count": 5}],
                   "cumulative": 21},
                  {"plan_item": "I4", "label_zh": "每周运动", "category": "exercise", "targets": ["ALT"], "status": "graduated",
                   "cumulative": None}]}


def _bad(**patch):
    d = copy.deepcopy(GOOD)
    for path, v in patch.items():
        tgt = d
        keys = path.split("__")
        for k in keys[:-1]:
            tgt = tgt[int(k)] if k.isdigit() else tgt[k]
        tgt[keys[-1]] = v
    return d


BAD = [_bad(extra="x"), _bad(items__0__note="free text"), _bad(items__0__status="doing"), _bad(items__0__label_zh="长" * 41),
       _bad(items__0__weekly_checkins=[{"week_start": "2026-08-17", "count": -1}]), _bad(generated_at="2026-09-20"),
       _bad(source="wechat"), _bad(items__0__adopted_at="2026/08/15"), _bad(items__0__cumulative=True)]


@pytest.mark.parametrize("doc", BAD)
def test_schema_rejects(doc):
    assert twin._schema_errors(doc, SCHEMA)


def test_schema_accepts_and_agrees_with_jsonschema():
    assert not twin._schema_errors(GOOD, SCHEMA)
    jsonschema = pytest.importorskip("jsonschema")
    v = jsonschema.Draft202012Validator(SCHEMA)
    assert not list(v.iter_errors(GOOD))
    for d in BAD:
        assert list(v.iter_errors(d)), d


def _ws(tmp_path: Path):
    ws = tmp_path / "ws"
    (ws / "work" / "intervene").mkdir(parents=True)
    (ws / "work" / "intervene" / "plan.json").write_text(json.dumps({"items": [
        {"id": "I3", "category": "diet", "action_zh": "限盐", "targets": ["收缩压"], "executor": "member", "retest": {"what": "收缩压", "after_weeks": 6}},
        {"id": "I4", "category": "exercise", "action_zh": "运动", "targets": ["ALT"], "executor": "member", "retest": {"what": "ALT", "after_weeks": 12}},
        {"id": "I5", "category": "sleep", "action_zh": "早睡", "targets": ["收缩压"], "executor": "member"}]}, ensure_ascii=False))
    st = C.new_state({"id": "m", "age": 50, "sex": "male", "sample_date": "2026-09-01", "answers": {}}, tmp_path)
    st["labs"] = [{"marker": "收缩压", "value": "131", "unit": "mmHg", "source_file": "F001",
                   "confirm": None, "maps_to": "收缩压"}]
    return ws, st


@needs_lib
def test_build_without_exposures_is_unchanged(tmp_path):
    ws, st = _ws(tmp_path)
    twin.build(st, ws)
    tw = json.loads((ws / "work" / "twin" / "twin.json").read_text())
    assert tw["exposures"] == [] and all(i["status"] == "proposed" and i["since"] is None for i in tw["interventions"])
    assert "exposures" not in st and not (ws / "work" / "twin" / "exposures.json").exists()


@needs_lib
def test_build_with_exposures(tmp_path):
    ws, st = _ws(tmp_path)
    f = tmp_path / "exp.json"
    f.write_text(json.dumps(GOOD, ensure_ascii=False))
    res = twin.build(st, ws, f)
    tw = json.loads((ws / "work" / "twin" / "twin.json").read_text())
    assert tw["exposures"] == GOOD["items"] and res["exposures"] == 2 and "warnings" not in res
    st_i = {i["id"]: (i["status"], i["since"]) for i in tw["interventions"]}
    assert st_i == {"I3": ("adopted", "2026-08-15"), "I4": ("graduated", None), "I5": ("proposed", None)}
    bound = ws / "work" / "twin" / "exposures.json"
    assert json.loads(bound.read_text()) == GOOD and st["exposures"]["sha256"] == C.sha256_file(bound, limit=None)
    assert bound in report.bound_files(ws, st)                         # joins the review trace hash list
    assert {r["intervention"]: r.get("judgeable") for r in tw["retest_plan"]} == {"I3": True, "I4": False}
    res = twin.build(st, ws)                                           # a rebuild keeps them
    assert res["exposures"] == 2 and "kept the exposures" in res["warnings"][0]
    bound.write_text(bound.read_text() + " ")                          # an edited copy is not kept
    assert twin.build(st, ws)["exposures"] == 0


@needs_lib
def test_build_rejects_bad_exposures_and_warns_on_member(tmp_path):
    ws, st = _ws(tmp_path)
    f = tmp_path / "exp.json"
    f.write_text(json.dumps(_bad(items__0__note="x"), ensure_ascii=False))
    with pytest.raises(C.LAError) as e:
        twin.build(st, ws, f)
    assert e.value.code == C.EXIT_INPUT and "unknown field 'note'" in str(e.value)
    f.write_text(json.dumps(_bad(member="pi-laowang"), ensure_ascii=False))
    assert "check it is the same person" in twin.build(st, ws, f)["warnings"][0]
