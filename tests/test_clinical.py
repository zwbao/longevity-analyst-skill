"""Diagnostic conclusions and medication suggestions (PLAN WS-L.3, WS-L.4). Synthetic members only."""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from clinical_helpers import ROOT, make_ws, series
from lalib import clinical, common as C, labels

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
pytestmark = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")
GLU = [{"marker": "空腹血糖", "unit": "mmol/L", "points": [{"t": t, "value": v, "k": 1} for t, v in
        (("2023-09-01", 4.6), ("2024-09-03", 5.0), ("2025-09-02", 5.3), ("2026-09-01", 5.6))]}]


def _dx(tmp_path, gset, labs=(("空腹血糖", 5.6, "mmol/L", "2026-09-01"),), readouts=(), age=49, ser=GLU):
    ws, st = make_ws(tmp_path, labs=labs, readouts=readouts, age=age)
    if ser:
        series(ws, ser)
    res = clinical.dx(st, ws, gset)
    return ws, st, res, json.loads((ws / "work" / "clinical" / "diagnoses.json").read_text())


def test_case_c_fasting_glucose_5_6_ada_yes_cn_no(tmp_path):
    ws, st, res, d = _dx(tmp_path / "ada", "ada")
    assert res["guideline_set"] == "ada" and "DX-ifg-ada" in res["met"]
    it = next(x for x in d["items"] if x["rule_id"] == "DX-ifg-ada")
    assert it["conclusion_zh"] == "按《ADA Standards of Care in Diabetes—2026》的标准，符合空腹血糖受损（糖尿病前期）的诊断标准。"
    assert it["label"]["basis_kind"] == "guideline_criteria" and it["label"]["label_zh"] == "按《ADA Standards of Care in Diabetes—2026》的标准"
    g = next(s for s in it["label"]["sources"] if s["type"] == "guideline")
    assert g["quote"] == "FPG 100 mg/dL (5.6 mmol/L) to 125 mg/dL (6.9 mmol/L) (IFG)" and g["url"] == "https://doi.org/10.2337/dc26-S002"
    assert labels.label_problems(it["label"]) == [] and labels.check(st, ws) == []
    ws, st, res, d = _dx(tmp_path / "cn", "cn")
    assert res["guideline_set"] == "cn" and not res["met"] and "DX-ifg-cn" in res["not_met"]
    assert d["guidelines_used_zh"] and "Guidelines on primary healthcare for type 2 diabetes in China, 2025" in d["guidelines_used_zh"]
    from lalib import report_clinical
    md = dict(report_clinical.sections(ws, {}))["诊断性结论"]
    assert "所用指南集：中国指南" in md and "已核对、未达到标准：空腹血糖受损" in md and "拨打 120" in md


def test_unit_conversion_and_distinct_days(tmp_path):
    _, _, res, _ = _dx(tmp_path / "a", "ada", labs=[("空腹血糖", 101, "mg/dL")], ser=None)
    assert "DX-ifg-ada" in res["met"]
    _, _, res, _ = _dx(tmp_path / "b", "cn", labs=[("空腹血糖", 110, "mg/dL")], ser=None)       # 6.1 mmol/L
    assert "DX-ifg-cn" in res["met"]
    one = [{"marker": "空腹血糖", "unit": "mmol/L", "points": [{"t": "2026-09-01", "value": 7.3}]}]
    _, _, res, _ = _dx(tmp_path / "c", "cn", labs=[("空腹血糖", 7.3, "mmol/L", "2026-09-01")], ser=one)
    assert "DX-dm-fpg-cn" in res["not_met"]                   # one day only: retest on another day
    two = [{"marker": "空腹血糖", "unit": "mmol/L", "points": [{"t": "2026-08-01", "value": 7.4}, {"t": "2026-09-01", "value": 7.3}]}]
    _, _, res, _ = _dx(tmp_path / "d", "cn", labs=[("空腹血糖", 7.3, "mmol/L", "2026-09-01")], ser=two)
    assert "DX-dm-fpg-cn" in res["met"]
    _, _, res, _ = _dx(tmp_path / "e", "ada", labs=[("空腹血糖", 7.3, "mmol/L"), ("糖化血红蛋白", 6.7, "%")], ser=None)
    assert "DX-dm-ada" in res["met"]                           # two different tests at the same time


def test_kdigo_and_fib4_as_conclusions(tmp_path):
    ro = [{"id": "native.organ.egfr", "label_zh": "eGFR", "value": 52.0, "unit": "mL/min/1.73m²", "kind": "computed"},
          {"id": "native.organ.fib4", "label_zh": "FIB-4", "value": 1.1, "unit": "", "kind": "computed"}]
    _, _, res, d = _dx(tmp_path / "a", "cn", labs=[], readouts=ro, ser=None)
    assert {"DX-gfr-g3a", "DX-fib4-low"} <= set(res["met"]) and "DX-ckd-gfr" in res["no_data"]
    g3 = next(x for x in d["items"] if x["rule_id"] == "DX-gfr-g3a")
    assert "KDIGO 2024" in g3["conclusion_zh"] and "单次测量" in g3["conclusion_zh"]
    ro[1]["value"] = 1.8
    _, _, res, _ = _dx(tmp_path / "b", "cn", labs=[], readouts=ro, ser=None, age=70)   # over 65 the lower cut-off is 2.0
    assert "DX-fib4-low" in res["met"]
    _, _, res, _ = _dx(tmp_path / "c", "cn", labs=[], readouts=ro, ser=None, age=50)
    assert "DX-fib4-raised" in res["met"]
    cr = [{"marker": "肌酐", "unit": "umol/L", "points": [{"t": "2026-03-01", "value": 150}, {"t": "2026-09-01", "value": 155}]}]
    _, _, res, d = _dx(tmp_path / "d", "cn", labs=[("肌酐", 155, "umol/L", "2026-09-01")], ser=cr, age=70)
    assert "DX-ckd-gfr" in res["met"]
    cr[0]["points"][0]["t"] = "2026-07-01"                      # under three months apart
    _, _, res, _ = _dx(tmp_path / "e", "cn", labs=[("肌酐", 155, "umol/L", "2026-09-01")], ser=cr, age=70)
    assert "DX-ckd-gfr" in res["not_met"]


@pytest.fixture
def pubmed_ok(monkeypatch):
    import lalib.evidence as E
    monkeypatch.setattr(E, "verify_pmids", lambda ids: {i: {} for i in ids})


def test_ai_diagnosis_separate_and_interval_only(tmp_path, pubmed_ok):
    ws, st = make_ws(tmp_path, labs=[("空腹血糖", 5.6, "mmol/L")])
    C.write_json(ws / "work" / "evidence" / "pubmed_x.json", {"items": [{"pmid": "31234567"}]})
    good = {"condition_zh": "胰岛素抵抗", "probability_range": [0.2, 0.5], "confidence": "low", "basis": ["空腹血糖"],
            "pmids": ["31234567"], "basis_zh": "空腹血糖在正常高限，结合年龄。"}
    for bad in (dict(good, probability_range=0.3), dict(good, confidence="moderate"), dict(good, basis_zh="空腹血糖 5.6"),
                dict(good, basis=["别人的化验"]), dict(good, pmids=["1"])):
        (tmp_path / "ai.json").write_text(json.dumps({"items": [bad]}, ensure_ascii=False))
        with pytest.raises(C.LAError):
            clinical.dx(st, ws, "cn", ai_file=tmp_path / "ai.json")
    (tmp_path / "ai.json").write_text(json.dumps({"items": [good]}, ensure_ascii=False))
    clinical.dx(st, ws, "cn", ai_file=tmp_path / "ai.json")
    d = json.loads((ws / "work" / "clinical" / "diagnoses.json").read_text())
    ai = [x for x in d["items"] if x["kind"] == "ai_estimate"][0]
    assert ai["conclusion_zh"] == "AI 判断：可能存在胰岛素抵抗（概率区间 20%–50%，依据空腹血糖在正常高限，结合年龄。）。"
    assert ai["label"]["label_zh"] == "AI 判断" and ai["label"]["confidence"]["level"] == "low"
    assert labels.check(st, ws) == []
    from lalib import report_clinical
    md = dict(report_clinical.sections(ws, {}))["诊断性结论"]
    assert "AI 判断（单独列出" in md


def _rx_ws(tmp_path, ldl=5.0, meds=None):
    ws, st = make_ws(tmp_path, labs=[("低密度脂蛋白胆固醇(LDL-C)", ldl, "mmol/L")],
                     answers={"current_medications": meds} if meds else None)
    (ws / "work" / "pgx").mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "tests" / "fixtures" / "pgx_slco1b1_synthetic.json", ws / "work" / "pgx" / "pgx.json")
    clinical.dx(st, ws, "cn")
    return ws, st


def test_rx_guideline_rule_and_cpic_verbatim(tmp_path):
    ws, st = _rx_ws(tmp_path, meds="辛伐他汀")
    res = clinical.rx(st, ws)
    assert res["guideline"] == ["RX-lipid-ldl49-cn"] and set(res["cpic"]) == {"simvastatin", "atorvastatin", "rosuvastatin"}
    d = json.loads((ws / "work" / "clinical" / "medications.json").read_text())
    g = next(x for x in d["items"] if x["kind"] == "guideline_recommendation")
    rule = next(r for r in C.data("medication_rules.json")["rules"] if r["id"] == "RX-lipid-ldl49-cn")
    assert g["action"] == "consider_start" and g["dose_text_zh"] == rule["dose_text_zh"] and "Atorvastatin 10–20 mg" in g["dose_text_zh"]
    assert g["label"]["label_zh"] == "《中国血脂管理指南（2023年）》推荐" and g["currently_taking"] is True
    pgx = json.loads((ROOT / "tests" / "fixtures" / "pgx_slco1b1_synthetic.json").read_text())
    for x in (x for x in d["items"] if x["kind"] == "cpic"):
        src = next(s for s in pgx["results"][0]["drugs"] if s["drug"] == x["drug"])
        assert x["text_zh"] == src["recommendation_text"] and x["label"]["label_zh"] == "CPIC 指南"
        assert x["label"]["sources"][0]["quote"] == src["recommendation_text"]
    act = {x["drug"]: x["action"] for x in d["items"] if x["kind"] == "cpic"}
    assert act == {"simvastatin": "consider_switch", "atorvastatin": "adjust", "rosuvastatin": "monitor"}
    assert labels.check(st, ws) == []
    d["items"][0]["dose_text_zh"] = "Atorvastatin 80 mg"               # a dose that is not the rule's own words
    C.write_json(ws / "work" / "clinical" / "medications.json", d)
    assert labels.check(st, ws)                                       # hash and dose both refused
    clinical.rx(st, ws)
    assert labels.check(st, ws) == []
    st["member"]["answers"]["current_medications"] = "无"                 # inputs changed after rx
    assert any("inputs that have changed" in b for b in labels.check(st, ws))


def test_rx_dm_rule_needs_dx_hit_and_ldl_target(tmp_path):
    ws, st = make_ws(tmp_path, labs=[("低密度脂蛋白胆固醇", 2.0, "mmol/L"), ("糖化血红蛋白", 6.8, "%")], age=52)
    clinical.dx(st, ws, "cn")
    assert clinical.rx(st, ws)["guideline"] == ["RX-lipid-dm40-cn"]
    ws, st = make_ws(tmp_path / "y", labs=[("低密度脂蛋白胆固醇", 2.0, "mmol/L"), ("糖化血红蛋白", 6.8, "%")], age=38)
    clinical.dx(st, ws, "cn")
    assert clinical.rx(st, ws)["guideline"] == []


def test_ai_medication_always_labelled_ai(tmp_path, pubmed_ok):
    ws, st = _rx_ws(tmp_path, ldl=3.0)
    C.write_json(ws / "work" / "evidence" / "pubmed_x.json", {"items": [{"pmid": "31234567"}]})
    good = {"drug_or_class_zh": "鱼油补充剂", "action": "monitor", "text_zh": "与他汀同用时留意肌肉酸痛。", "dose_text_zh": "每日 1 g",
            "confidence": "very_low", "basis": ["低密度脂蛋白胆固醇(LDL-C)"], "pmids": ["31234567"], "basis_zh": "文献提示相互作用很少。"}
    for bad in (dict(good, action="take"), dict(good, text_zh="每天吃 2 粒"), dict(good, confidence="high")):
        (tmp_path / "ai.json").write_text(json.dumps({"items": [bad]}, ensure_ascii=False))
        with pytest.raises(C.LAError):
            clinical.rx(st, ws, ai_file=tmp_path / "ai.json")
    (tmp_path / "ai.json").write_text(json.dumps({"items": [good]}, ensure_ascii=False))
    res = clinical.rx(st, ws, ai_file=tmp_path / "ai.json")
    assert res["ai"] == 1 and res["guideline"] == []
    ai = [x for x in json.loads((ws / "work" / "clinical" / "medications.json").read_text())["items"] if x["kind"] == "ai_estimate"][0]
    assert ai["label"]["basis_kind"] == "ai_estimate" and ai["label"]["confidence"] == {"level": "very_low"}
    assert labels.check(st, ws) == []
    txt = labels.render_item("rx", ai)
    assert txt.startswith("〔AI 判断，置信度很低〕") and "每日 1 g" in txt
    from lalib import report_clinical
    md = dict(report_clinical.sections(ws, {}))["用药建议"]
    assert "剂量也是 AI 判断" in md and "CPIC 指南" in md
