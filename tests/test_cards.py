"""deliver/cards.json (pi-card/1): every number resolves to the export it was made from. Synthetic data only."""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from clinical_helpers import ROOT
from lalib import cards, common as C, labels as L

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")


def _export_v2():
    pred = {"schema": "la-prediction/1", "id": "P1", "plan_item": "I3", "intervention_zh": "减盐", "marker": "收缩压",
            "design": {"repeats_each_side": 14, "weeks": 6, "adherence": 1.0}, "p_detect_favorable": 0.764, "p_true_favorable": 0.999,
            "expected_change_pct": -4.6, "interval_80_pct": [-7.4, -1.9],
            "label": L.make_label("statistical_model", [{"type": "effects", "ref": "salt-reduction-sbp"}], assumptions=["a"],
                                  confidence={"interval_80": [-7.4, -1.9], "unit": "%"})}
    dx = {"id": "DX1", "kind": "guideline_criteria", "rule_id": "DX-ifg-ada", "label_zh": "空腹血糖受损（糖尿病前期）",
          "conclusion_zh": "按《ADA Standards of Care in Diabetes—2026》的标准，符合空腹血糖受损（糖尿病前期）的诊断标准。",
          "caveats_zh": ["这是 2026-09-01 一次测量得出的结论，需要复查确认。"],
          "evidence": [{"type": "lab", "ref": "空腹血糖", "value": 5.6, "unit": "mmol/L", "t": "2026-09-01"}], "executor": ["physician"],
          "label": L.make_label("guideline_criteria", [{"type": "guideline", "ref": "ADA", "quote": "q"}], guideline_zh="ADA Standards of Care in Diabetes—2026")}
    aidx = {"id": "DX2", "kind": "ai_estimate", "label_zh": "胰岛素抵抗", "probability_range": [0.2, 0.5],
            "conclusion_zh": "AI 判断：可能存在胰岛素抵抗（概率区间 20%–50%，依据空腹血糖。）。",
            "label": L.make_label("ai_estimate", [{"type": "pmid", "ref": "1"}], confidence={"level": "low", "range": [0.2, 0.5]})}
    cp = {"id": "RX1", "kind": "cpic", "gene": "SLCO1B1", "phenotype": "Decreased Function", "drug": "atorvastatin", "drug_or_class_zh": "阿托伐他汀",
          "action": "adjust", "label_zh": "SLCO1B1 与 阿托伐他汀", "text_zh": "Prescribe ≤40mg as a starting dose", "executor": ["physician", "member"],
          "label": L.make_label("cpic", [{"type": "cpic", "ref": "SLCO1B1, ABCG2, CYP2C9, and Statins", "quote": "Prescribe ≤40mg as a starting dose"}])}
    med = dict(cp, id="RX2", kind="guideline_recommendation", drug_or_class_zh="他汀类降脂药", action="consider_start", label_zh="他汀：考虑开始",
               label=L.make_label("guideline_recommendation", [{"type": "rule", "ref": "RX-lipid-ldl49-cn"}], guideline_zh="中国血脂管理指南（2023年）"))
    return {"schema": "la-export/2", "generated_at": "2026-09-30T10:00:00+08:00",
            "trends": [{"marker": "空腹血糖", "unit": "mmol/L", "trend": {"verdict": "rising_beyond_noise", "annual_change_pct": 6.4},
                        "baseline_vs_latest": {"change_pct": 21.7, "verdict": "increase_beyond_noise"}, "bv_source": "10.1515/cclm-2021-0370"}],
            "concordance": [{"analyte": "urate", "label_zh": "尿酸", "genetic_pct": 86, "class": "concordant_adverse", "wording_zh": "尿酸偏高的基因倾向已经表达。"}],
            "surveillance": [{"id": "S1", "marker": "脂蛋白(a)", "interval_weeks": None, "once": True, "done": True, "reason": "guideline",
                              "source": {"title": "某指南", "doi": "10.1/x", "quote": "q"}},
                             {"id": "S2", "marker": "同型半胱氨酸", "interval_weeks": 26, "once": False, "done": False, "reason": "genetic",
                              "source": {"title": "某指南", "url": "https://example.org/g", "quote": "q"}}],
            "predictions": [pred], "diagnoses": [dx, aidx], "medications": [cp, med],
            "retests": [{"item": "I3", "what": "收缩压", "after_weeks": 6, "due_at": "2026-11-10T08:00:00+08:00", "reason": "intervention", "judgeable": True}]}


def test_cards_cover_every_kind_and_numbers_resolve(tmp_path):
    ex = tmp_path / "deliver" / "la-export.json"
    C.write_json(ex, _export_v2())
    cs = cards.build(tmp_path, ex)
    kinds = {c["kind"] for c in cs}
    assert kinds == {"trend", "change", "prediction", "diagnosis", "medication", "pgx", "concordance", "retest_due", "surveillance"}
    assert cards.verify(cs, ex, tmp_path) == []
    sha = C.sha256_file(ex, limit=None)
    for c in cs:
        assert c["schema"] == "pi-card/1" and c["source"] == {"export": "la-export/2", "sha256": sha}
        assert c["deep_link"] == f"longpi://cards/{c['id']}" and c["actions"] and c["title_zh"] and c["summary_zh"]
        assert c["kind"] == "retest_due" or L.label_problems(c["label"]) == []
    dxc = next(c for c in cs if c["id"].endswith("dx-DX1"))
    assert dxc["summary_zh"].endswith("这是 2026-09-01 一次测量得出的结论，需要复查确认。")          # caveats as in the report
    tr = next(c for c in cs if c["kind"] == "trend")
    assert tr["numbers"] == [{"ref": "trends[0].trend.annual_change_pct", "value": 6.4, "unit": "%/年", "label": "测量"}]
    pr = next(c for c in cs if c["kind"] == "prediction")
    assert {n["ref"]: n["value"] for n in pr["numbers"]}["predictions[0].p_detect_favorable"] == 0.764
    assert pr["label"]["label_zh"] == "统计模型估计"
    ai = next(c for c in cs if c["id"].endswith("DX2"))
    assert ai["label"]["label_zh"] == "AI 判断" and ai["numbers"][0]["value"] == [0.2, 0.5]
    assert next(c for c in cs if c["kind"] == "pgx")["label"]["label_zh"] == "CPIC 指南"
    # a card number that does not match the export is caught
    cs[0]["numbers"][0]["value"] = 9.9
    assert cards.verify(cs, ex, tmp_path)


def test_cards_read_handover_files_when_export_v1(tmp_path):
    ex = tmp_path / "deliver" / "la-export.json"
    v2 = _export_v2()
    C.write_json(tmp_path / "work" / "predict" / "predictions.json", {"schema": "la-predictions/1", "items": v2["predictions"]})
    C.write_json(ex, {"schema": "la-export/1", "generated_at": "2026-09-30T10:00:00+08:00", "retests": v2["retests"]})
    cs = cards.build(tmp_path, ex)
    pr = next(c for c in cs if c["kind"] == "prediction")
    assert pr["numbers"][0]["ref"].startswith("file:work/predict/predictions.json#items[0].")
    assert cards.verify(cs, ex, tmp_path) == []


def ok(la, *args):
    assert la.main(list(args)) == 0, args


@needs_lib
def test_chain_case_a_labelled_items_reach_report_export_and_cards(tmp_path, capsys, monkeypatch):
    """CLI chain on examples/case-A (public data) with dx, rx, a medication plan item, an auto prediction and placeholders."""
    import la
    ws = tmp_path / "ws"
    ok(la, "init", str(ROOT / "examples" / "case-A" / "raw"), str(ws), "--member-id", "caseA", "--age", "67", "--sex", "f",
       "--sample-date", "2026-09-10", "--longevity-skills", LIB)
    ok(la, "intake", str(ws))
    st = json.loads((ws / "state.json").read_text())
    beta = next(f["id"] for f in st["files"] if f["kind"] == "methylation_beta")
    ok(la, "assign", str(ws), "--file", beta, "--what", "tissue", "--value", "whole_blood", "--reason", "GEO: whole blood")
    ok(la, "assign", str(ws), "--what", "identity", "--value", "consistent", "--reason", "fixture")
    capsys.readouterr()
    ok(la, "labs", "candidates", str(ws))
    pend = json.loads(capsys.readouterr().out)["pending"]
    (ws / "a.json").write_text(json.dumps({"answers": [{"row_key": p["row_key"], "answer": "yes", "why": "fixture"} for p in pend]}))
    ok(la, "labs", "confirm", str(ws), "--answers", str(ws / "a.json"), "--reason", "fixture")
    ok(la, "preflight", str(ws), "--no-network")
    ok(la, "member", str(ws), "smoker=no", "treated=no", "diabetes=no", "north=no", "genetic_disclosure=yes",
       "current_medications=阿托伐他汀", "--source", "fixture")
    ok(la, "methods", "plan", str(ws))
    ok(la, "methods", "run", str(ws))
    ok(la, "integrate", "bundle", str(ws))
    st = json.loads((ws / "state.json").read_text())
    ana = ws / "work" / "integrate" / "analyses"
    ana.mkdir(parents=True, exist_ok=True)
    for b in st["integrate"]["bundles"]:
        rid = json.loads(Path(b["bundle"]).read_text())["readouts"][0]["id"]
        (ana / f"{b['system']}.md").write_text(f"## 结论\n见读数表。\n\n## 证据\n- {{{{r:{rid}|label}}}} = {{{{r:{rid}}}}}\n\n"
                                              "## 不确定性\n单次测量。\n\n## 建议复测\n按计划复测。\n", encoding="utf-8")
        ok(la, "integrate", "register", str(ws), "--system", b["system"])
    ok(la, "organ", "bundle", str(ws))
    for organ in C.data("organs.json")["organs"]:
        ok(la, "organ", "skip", str(ws), "--organ", organ, "--reason", "fixture")
    ok(la, "insights", "skip-genomics", str(ws), "--reason", "offline test run")
    ok(la, "insights", "position", str(ws))
    bdir = ws / "work" / "insights" / "board"
    bdir.mkdir(parents=True, exist_ok=True)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": f"问题{'甲乙丙'[i - 1]}", "hypothesis_zh": "待检验。", "basis": ["epiage.dnam_hannum"], "why_zh": "本人数据。"}
                        for i in (1, 2, 3)]}
    (bdir / "q.json").write_text(json.dumps(qs, ensure_ascii=False), encoding="utf-8")
    ok(la, "board", "questions", str(ws), "--file", str(bdir / "q.json"))
    for i in (1, 2, 3):
        ok(la, "board", "skip", str(ws), "--id", f"Q{i}", "--reason", "offline test")
    (ws / "work" / "pgx").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "pgx" / "pgx.json").write_text((ROOT / "tests" / "fixtures" / "pgx_slco1b1_synthetic.json").read_text())
    capsys.readouterr()
    ok(la, "clinical", "dx", str(ws), "--guideline-set", "ada")
    dx = json.loads(capsys.readouterr().out)
    assert "DX-pre-a1c-ada" in dx["met"]                         # HbA1c 5.7 %
    ok(la, "clinical", "rx", str(ws))
    rx = json.loads(capsys.readouterr().out)
    assert set(rx["cpic"]) == {"simvastatin", "atorvastatin", "rosuvastatin"}
    meds = json.loads((ws / "work" / "clinical" / "medications.json").read_text())
    ator = next(x for x in meds["items"] if x.get("drug") == "atorvastatin")
    assert ator["currently_taking"] and ator["action"] == "adjust"
    plan = {"items": [dict(json.loads(json.dumps(__import__("clinical_helpers").SALT_ITEM)), action_zh="减盐，用限盐勺"),
                      {"id": "I2", "category": "medication", "action_zh": f"请医生评估：{{{{rx:{ator['id']}}}}}", "rationale_zh": "用药基因提示。",
                       "targets": ["低密度脂蛋白胆固醇(LDL-C)"], "evidence": [], "evidence_grade": "guideline", "executor": "physician"}]}
    (ws / "work" / "intervene").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "intervene" / "plan.json").write_text(json.dumps(plan, ensure_ascii=False), encoding="utf-8")
    ok(la, "intervene", "register", str(ws), "--plan", str(ws / "work" / "intervene" / "plan.json"))
    reg = json.loads(capsys.readouterr().out)
    assert reg["predictions"] == {"statistical": 1, "ai": 0}           # attached automatically
    ok(la, "twin", "build", str(ws))
    (ws / "work" / "report").mkdir(parents=True, exist_ok=True)
    summ = ws / "work" / "report" / "summary.md"
    summ.write_text("减盐：{{p:P1}}\n\n血糖：{{dx:DX1}}\n\n用药基因：{{rx:%s}}\n" % ator["id"], encoding="utf-8")
    capsys.readouterr()
    ok(la, "review", "trace", str(ws))
    trace = json.loads(capsys.readouterr().out)
    summ.write_text(summ.read_text() + "\n减盐有效的概率约 80%。\n", encoding="utf-8")       # an unregistered probability
    assert la.main(["review", "trace", str(ws)]) == C.EXIT_BLOCKED
    summ.write_text(summ.read_text().replace("\n减盐有效的概率约 80%。\n", "\n"), encoding="utf-8")
    capsys.readouterr()
    ok(la, "review", "trace", str(ws))
    trace = json.loads(capsys.readouterr().out)
    rev = ws / "work" / "review" / "reviewer-1.json"
    rev.write_text(json.dumps({"trace_id": trace["trace_id"], "verdict": "pass", "findings": []}))
    ok(la, "review", "record", str(ws), "--verdict", "pass", "--findings", str(rev))
    ok(la, "report", str(ws))
    md = (ws / "deliver" / "report.md").read_text()
    assert "〔统计模型估计〕减盐对收缩压" in md and "〔按《ADA Standards of Care in Diabetes—2026》的标准〕" in md and "〔CPIC 指南〕" in md
    assert "Prescribe ≤40mg as a starting dose" in md
    cs = json.loads((ws / "deliver" / "cards.json").read_text())
    assert cards.verify(cs, ws / "deliver" / "la-export.json", ws) == []
    assert {"prediction", "diagnosis", "pgx", "retest_due"} <= {c["kind"] for c in cs}
    st = json.loads((ws / "state.json").read_text())
    assert "cards.json" in st["report"]["deliver_sha256"]
    ok(la, "export", str(ws))
    from lalib import report_clinical
    secs = dict(report_clinical.sections(ws, {}))
    assert set(secs) == {"个人反应预测", "诊断性结论", "用药建议"}
    # a labelled file changed after the report: the delivery is no longer current
    ok(la, "clinical", "dx", str(ws), "--guideline-set", "cn")
    assert la.main(["export", str(ws)]) == C.EXIT_BLOCKED
