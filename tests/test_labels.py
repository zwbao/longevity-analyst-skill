"""la-label/1 and the {{p:}} {{dx:}} {{rx:}} placeholders (PLAN WS-L.1, WS-L.5). Synthetic data only."""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from clinical_helpers import ROOT, SALT_ITEM, make_ws, write_plan
from lalib import common as C, labels as L
from lalib.labels import Bound
from lalib.report import substitute, trace_text

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")
LABEL_SCHEMA = json.loads((ROOT / "skills" / "longevity-analyst" / "schema" / "la-label-1.schema.json").read_text())


def _validate(lab):
    import jsonschema
    jsonschema.Draft202012Validator(LABEL_SCHEMA).validate(lab)


def test_label_zh_is_rendered_from_basis_kind():
    assert L.make_label("measurement", [{"type": "lab", "ref": "空腹血糖"}])["label_zh"] == "测量"
    g = L.make_label("guideline_criteria", [{"type": "guideline", "ref": "x", "quote": "q"}], guideline_zh="《ADA 2026》")
    assert g["label_zh"] == "按《ADA 2026》的标准"
    assert L.make_label("guideline_recommendation", [{"type": "rule", "ref": "R"}], guideline_zh="中国血脂管理指南（2023年）")["label_zh"] \
        == "《中国血脂管理指南（2023年）》推荐"
    assert L.make_label("cpic", [{"type": "cpic", "ref": "SLCO1B1"}])["label_zh"] == "CPIC 指南"
    ai = L.make_label("ai_estimate", [{"type": "pmid", "ref": "1"}], confidence={"level": "very_low", "range": [0.1, 0.3]})
    st = L.make_label("statistical_model", [{"type": "effects", "ref": "e"}], assumptions=["a"], confidence={"interval_80": [-7.3, -1.9], "unit": "%"})
    pub = L.make_label("published_model", [{"type": "readout", "ref": "china-par"}])
    for lab in (g, ai, st, pub):
        _validate(lab)
        assert L.label_problems(lab) == []
    assert L.render_label(ai) == "〔AI 判断，置信度很低〕" and L.render_label(st) == "〔统计模型估计〕"
    assert L.render_sources(st) == "来源：效应表 e"


@pytest.mark.parametrize("bad,why", [
    (None, "missing label"),
    ({"schema": "la-label/1", "basis_kind": "measurement", "label_zh": "测量", "sources": [], "generated_at": "t"}, "missing sources"),
    ({"schema": "la-label/1", "basis_kind": "ai_estimate", "label_zh": "AI 判断", "sources": [{"type": "pmid", "ref": "1"}],
      "confidence": {"level": "moderate"}, "generated_at": "t"}, "low or very_low"),
    ({"schema": "la-label/1", "basis_kind": "ai_estimate", "label_zh": "按《指南》的标准", "sources": [{"type": "pmid", "ref": "1"}],
      "confidence": {"level": "low"}, "generated_at": "t"}, "rendered tag"),
    ({"schema": "la-label/1", "basis_kind": "statistical_model", "label_zh": "统计模型估计", "sources": [{"type": "effects", "ref": "e"}],
      "confidence": {"interval_80": [1, 2]}, "generated_at": "t"}, "assumptions"),
])
def test_label_rejections(bad, why):
    assert any(why in p for p in L.label_problems(bad))
    with pytest.raises(C.LAError):
        L.make_label("ai_estimate", [{"type": "pmid", "ref": "1"}], confidence={"level": "high"})


def _registered_prediction(tmp_path):
    ws, st = make_ws(tmp_path, labs=[("收缩压", 138, "mmHg")])
    write_plan(ws, [SALT_ITEM])
    from lalib import predict
    predict.run(st, ws, repeats=14)
    return ws, st


@needs_lib
def test_placeholder_renders_tag_sentence_and_sources(tmp_path):
    ws, st = _registered_prediction(tmp_path)
    ro = Bound({})
    ro.labelled = L.registry(ws)
    out = substitute("减盐对你可能有用：{{p:P1}}", ro)
    assert out.startswith("减盐对你可能有用：〔统计模型估计〕减盐对收缩压：") and "约 49%" in out and "来源：效应表 salt-reduction-sbp" in out
    assert "[缺失条目 p:P9]" in substitute("{{p:P9}}", ro)
    assert trace_text("减盐对你可能有用：{{p:P1}}。") == []                 # the item's numbers are its own, checked in labels.check


def test_unregistered_probability_or_dose_in_prose_is_refused():
    assert trace_text("减盐有七成把握能测出来。")
    assert trace_text("能测出下降的概率约 76%。")
    assert trace_text("建议阿托伐他汀 20 mg 每晚一次。")
    assert trace_text("他汀每天吃两片。")
    assert trace_text("{{rx:RX1|20 mg}}")                                  # a placeholder cannot carry its own words


@needs_lib
def test_check_refuses_missing_label_source_and_ai_point(tmp_path):
    ws, st = _registered_prediction(tmp_path)
    p = ws / "work" / "predict" / "predictions.json"
    d = json.loads(p.read_text())
    good = json.loads(json.dumps(d))
    d["items"][0].pop("label")
    C.write_json(p, d)
    bad = L.check(st, ws)
    assert any("does not match" in b for b in bad)                         # edited by hand
    L.record(st, ws, "p", st["predict"]["inputs"])                         # even if re-bound, the content is refused
    assert any("missing label" in b for b in L.check(st, ws))
    d = json.loads(json.dumps(good))
    d["items"][0]["label"]["sources"] = []
    C.write_json(p, d)
    L.record(st, ws, "p", st["predict"]["inputs"])
    assert any("missing sources" in b for b in L.check(st, ws))
    ai = {"schema": "la-prediction/1", "id": "P2", "plan_item": "I1", "intervention_zh": "减盐", "marker": "收缩压", "p_detect_favorable": 0.5,
          "label": L.make_label("ai_estimate", [{"type": "pmid", "ref": "1"}], confidence={"level": "low"})}
    d = json.loads(json.dumps(good))
    d["items"].append(ai)
    C.write_json(p, d)
    L.record(st, ws, "p", st["predict"]["inputs"])
    assert any("never a point" in b for b in L.check(st, ws))
    ai["label"]["confidence"]["level"] = "moderate"
    C.write_json(p, d)
    L.record(st, ws, "p", st["predict"]["inputs"])
    assert any("low or very_low" in b for b in L.check(st, ws))


@needs_lib
def test_plan_medication_items(tmp_path, monkeypatch):
    import shutil
    from lalib import clinical, evidence
    ws, st = make_ws(tmp_path, labs=[("低密度脂蛋白胆固醇(LDL-C)", 5.0, "mmol/L")])
    (ws / "work" / "pgx").mkdir(parents=True)
    shutil.copy(ROOT / "tests" / "fixtures" / "pgx_slco1b1_synthetic.json", ws / "work" / "pgx" / "pgx.json")
    clinical.dx(st, ws, "cn")
    clinical.rx(st, ws)
    plan_p = ws / "work" / "intervene" / "draft.json"

    def check(items):
        C.write_json(plan_p, {"items": items})
        return evidence.check_plan(st, ws, plan_p)
    med = {"id": "I2", "category": "medication", "action_zh": "请医生评估：{{rx:RX1}}", "rationale_zh": "低密度脂蛋白胆固醇偏高。",
           "targets": ["低密度脂蛋白胆固醇(LDL-C)"], "evidence": [], "evidence_grade": "guideline", "executor": "physician"}
    assert check([med])["items"] == 1
    for bad, why in ((dict(med, action_zh="请医生评估他汀"), "{{rx:<id>}}"), (dict(med, action_zh="{{rx:RX99}}"), "not a registered"),
                     (dict(med, action_zh="{{rx:RX1}}，阿托伐他汀 20 mg"), "dose"),
                     (dict(SALT_ITEM, action_zh="减盐，{{rx:RX1}}", targets=["低密度脂蛋白胆固醇(LDL-C)"]), "belongs in a `medication` item"),
                     (dict(SALT_ITEM, action_zh="减盐 {{p:P1}}", targets=["低密度脂蛋白胆固醇(LDL-C)"]), "attached to plan items automatically")):
        with pytest.raises(C.LAError) as e:
            check([bad])
        assert why in str(e.value), why
    menu = {c["id"]: c for c in C.data("intervention_menu.json")["categories"]}
    assert menu["medication"]["executors"] == ["physician", "member"]
