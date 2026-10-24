"""Personal response prediction (PLAN WS-L.2, appendix B.7) against the shared prediction-cases/1 fixture."""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from clinical_helpers import ROOT, SALT_ITEM, make_ws, series, write_plan
from lalib import common as C, labels, predict

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
pytestmark = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")
CASES = json.loads((ROOT / "tests" / "fixtures" / "prediction_cases.json").read_text(encoding="utf-8"))
SCHEMA_DIR = ROOT / "skills" / "longevity-analyst" / "schema"


def _validator():
    import jsonschema
    from referencing import Registry, Resource
    lab = json.loads((SCHEMA_DIR / "la-label-1.schema.json").read_text())
    reg = Registry().with_resource("la-prediction/la-label-1.schema.json", Resource.from_contents(lab))   # $ref relative to $id la-prediction/1
    return jsonschema.Draft202012Validator(json.loads((SCHEMA_DIR / "la-prediction-1.schema.json").read_text()), registry=reg)


@pytest.mark.parametrize("case", CASES["cases"], ids=[c["id"] for c in CASES["cases"]])
def test_prediction_cases(case):
    effects, bv = predict._library()
    row = effects[case["effects_id"]]
    assert row["effect"]["value"] == case["effect"]["value"] and row["effect"]["ci"] == case["effect"]["ci"]
    m = next(x for x in bv["markers"] if x["key"] == case["marker"])
    assert (m["better"] == "lower") == case["lower_is_better"]
    k = case["repeats_each_side"]
    r = predict.model(row["effect"]["value"], tuple(row["effect"]["ci"]), case["baseline"], m["cvi_pct"], m.get("cva_pct"), k, k,
                      lower_better=case["lower_is_better"], a=case["adherence"], tau=case["tau_pct"], z=CASES["z"],
                      personal=[(p["change_pct"], p["repeats_each_side"], p["repeats_each_side"]) for p in case["personal"]])
    e = case["expect"]
    for key in ("mu0_pct", "s0_pct", "band_pct", "posterior_mu_pct", "posterior_s_pct"):
        assert r[key] == pytest.approx(e[key], abs=0.1), key
    for key in ("p_detect_favorable", "p_true_favorable"):
        assert r[key] == pytest.approx(e[key], abs=0.01), key
    assert r["interval_80_pct"] == pytest.approx(e["interval_80_pct"], abs=0.1)


def test_plan_item_gets_statistical_prediction_and_personal_update(tmp_path):
    ws, st = make_ws(tmp_path, labs=[("收缩压", 138, "mmHg")])
    write_plan(ws, [SALT_ITEM, dict(SALT_ITEM, id="I2", action_zh="多吃全谷物", evidence=[], targets=["收缩压"])])
    res = predict.run(st, ws, repeats=14)
    assert res["statistical"] == 1 and [x["plan_item"] for x in res["ai_requests"]] == ["I2"]
    p = json.loads((ws / "work" / "predict" / "predictions.json").read_text())["items"][0]
    _validator().validate(p)
    assert p["plan_item"] == "I1" and p["p_detect_favorable"] == pytest.approx(0.49, abs=0.01)
    assert p["p_true_favorable"] == pytest.approx(0.84, abs=0.01)
    assert p["label"]["basis_kind"] == "statistical_model" and p["label"]["label_zh"] == "统计模型估计"
    assert any("个人差异默认与平均效应同量级" in a for a in p["label"]["assumptions"])
    assert any("依从比例按 1 计" in a for a in p["label"]["assumptions"])
    assert {s["type"] for s in p["label"]["sources"]} >= {"effects", "biovar"}
    assert labels.check(st, ws) == []
    # the member ran the experiment once: 138 -> 131 over 14 + 14 home readings while salt reduction ran
    series(ws, [{"marker": "收缩压", "unit": "mmHg", "points": [{"t": "2026-07-01", "value": 138, "k": 14}, {"t": "2026-09-01", "value": 131, "k": 14}],
                 "pairs": [{"from": "2026-07-01", "to": "2026-09-01", "change_pct": 100 * (131 - 138) / 138,
                            "concurrent_exposures": [{"plan_item": "I1", "label_zh": "减盐"}]}]}])
    assert any("inputs that have changed" in b for b in labels.check(st, ws))      # stale until re-run
    predict.run(st, ws, repeats=14)
    p = json.loads((ws / "work" / "predict" / "predictions.json").read_text())["items"][0]
    assert p["p_detect_favorable"] == pytest.approx(0.764, abs=0.01) and p["p_true_favorable"] > 0.99
    assert p["interval_80_pct"] == pytest.approx([-7.4, -1.9], abs=0.1)
    assert any(s["type"] == "personal_experiment" for s in p["label"]["sources"])
    assert labels.check(st, ws) == []


def test_adherence_scaling_and_site_without_checkins():
    a, note = predict.adherence({"status": "adopted", "adopted_at": "2026-07-20", "cumulative": 21}, 7, "2026-08-31")
    assert a == pytest.approx(0.5) and "线性缩放" in note
    a, note = predict.adherence({"status": "adopted", "adopted_at": "2026-07-20", "cumulative": None}, 7, "2026-08-31")
    assert a == 1.0 and "不记录打卡次数" in note
    assert predict.adherence({"status": "adopted", "adopted_at": "2026-07-20", "cumulative": 500}, 7, "2026-08-31")[0] == 1.0
    assert predict.adherence(None, 7, "2026-08-31") == (1.0, "尚未开始执行，依从比例按 1 计")


def test_adherence_from_exposures_scales_prior(tmp_path):
    ws, st = make_ws(tmp_path, labs=[("收缩压", 138, "mmHg")])
    write_plan(ws, [SALT_ITEM])
    today = C.now_iso()[:10]
    import datetime as dt
    start = (dt.date.fromisoformat(today) - dt.timedelta(weeks=6)).isoformat()
    exp = {"schema": "la-exposures/1", "member": "synthetic", "source": "coach_member_file", "generated_at": today + "T09:00:00+08:00",
           "items": [{"plan_item": "I1", "label_zh": "减盐", "category": "diet", "targets": ["收缩压"], "status": "adopted",
                      "adopted_at": start, "cumulative": 21}]}
    (tmp_path / "exp.json").write_text(json.dumps(exp, ensure_ascii=False))
    predict.run(st, ws, repeats=14, exposures_file=tmp_path / "exp.json")
    p = json.loads((ws / "work" / "predict" / "predictions.json").read_text())["items"][0]
    assert p["design"]["adherence"] == pytest.approx(0.5, abs=0.01)
    assert p["p_detect_favorable"] == pytest.approx(0.236, abs=0.01)


def test_ai_fallback_is_interval_only(tmp_path, monkeypatch):
    ws, st = make_ws(tmp_path, labs=[("收缩压", 138, "mmHg")])
    write_plan(ws, [dict(SALT_ITEM, id="I2", action_zh="多吃全谷物", evidence=[])])
    predict.run(st, ws)
    C.write_json(ws / "work" / "evidence" / "pubmed_x.json", {"items": [{"pmid": "23558162"}]})
    import lalib.evidence as E
    monkeypatch.setattr(E, "verify_pmids", lambda ids: {i: {} for i in ids})
    good = {"plan_item": "I2", "intervention_zh": "多吃全谷物", "marker": "sbp", "marker_zh": "收缩压", "p_true_favorable": [0.4, 0.8],
            "confidence": "low", "pmids": ["23558162"], "basis_zh": "以全谷物干预试验为锚，证据间接。"}
    for bad in (dict(good, p_true_favorable=0.6), dict(good, confidence="moderate"), dict(good, pmids=["999"]),
                dict(good, basis_zh="约三成人有效"), dict(good, plan_item="I9")):
        (tmp_path / "ai.json").write_text(json.dumps({"items": [bad]}, ensure_ascii=False))
        with pytest.raises(C.LAError):
            predict.run(st, ws, ai_file=tmp_path / "ai.json")
    (tmp_path / "ai.json").write_text(json.dumps({"items": [good]}, ensure_ascii=False))
    res = predict.run(st, ws, ai_file=tmp_path / "ai.json")
    assert res["ai"] == 1
    p = json.loads((ws / "work" / "predict" / "predictions.json").read_text())["items"][0]
    _validator().validate(p)
    assert p["label"]["basis_kind"] == "ai_estimate" and p["label"]["confidence"] == {"level": "low"}
    assert "expected_change_pct" not in p and labels.check(st, ws) == []
    assert "AI 判断" in labels.render_item("p", p) and "40%–80%" in labels.render_item("p", p)
