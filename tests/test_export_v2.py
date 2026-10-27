"""la-export/2 and the privacy rules for what leaves the workspace (WS-B.2, WS-K, PLAN A.8).

The chain test runs case-A (public data) through the CLI offline, adds SYNTHETIC handover files (an earlier visit's
twin, exposures, concordance, fingerprint, surveillance plan) and checks the export, the report section and the
delivered twin. Every value in the synthetic files is made up for the test.
"""
from __future__ import annotations

import copy
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
import la  # noqa: E402
from lalib import common as C  # noqa: E402
from lalib import export as X  # noqa: E402

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")

CALLS = {"rs4149056": "TT", "rs429358": "TT", "rs7412": "CC", "rs1801133": "AG", "rs671": "GG"}   # synthetic


def ok(*args):
    rc = la.main(list(args))
    assert rc == 0, args


# ------------------------------------------------------------ privacy, unit level
@pytest.mark.parametrize("doc", [
    {"x": [{"rsid": "rs4149056", "genotype": "TT"}]},
    {"x": {"rs4149056": "CT"}},
    {"note": "你的 rs4149056 (C;T) 提示……"},
    {"note": "rs671:GA"},
    {"member": {"name": "王某"}},
    {"contact": "请拨打 138 1234 5678"},
    {"text": "姓名：王某某"},
    {"genotype_fingerprint": {"sites": 3}},
    {"twin": {"calls": {"rs1": "AA"}}},
])
def test_privacy_rejects(doc):
    assert X.privacy_problems(doc)


@pytest.mark.parametrize("doc", [
    {"report": {"sha256": "a13812345678b" * 4}, "generation": "e1381234567890abcdef"},
    {"label_zh": "APOE 基因型", "value": "ε3/ε4"},
    {"note": "SLCO1B1 功能降低，按 CPIC 指南", "ref": "PMID 34554658"},
    {"rsid_list": ["rs4149056"]},
])
def test_privacy_allows(doc):
    assert not X.privacy_problems(doc), X.privacy_problems(doc)


def test_as_v1_drops_v2_fields():
    d = {"schema": "la-export/2", "trends": [1], "concordance": [], "concordance_note_zh": None, "surveillance": [],
         "predictions": [], "diagnoses": [], "medications": [], "identity_alerts": [], "concordance_transitions": [],
         "retests": [{"item": "I1", "what": "x", "after_weeks": 4, "due": "d", "from": "f", "reason": "intervention",
                      "judgeable": True, "due_at": "z"},
                     {"item": "S1", "what": "y", "after_weeks": 52, "due": "d", "from": "f", "reason": "genetic",
                      "judgeable": False, "due_at": "z"}]}
    v1 = X.as_v1(d)
    assert v1["schema"] == "la-export/1" and set(v1) == {"schema", "retests"}
    assert v1["retests"] == [{"item": "I1", "what": "x", "after_weeks": 4, "due": "d", "from": "f"}]


# ------------------------------------------------------------ the chain
def _chain_to_intervene(ws: Path, capsys, monkeypatch):
    raw = ROOT / "examples" / "case-A" / "raw"
    ok("init", str(raw), str(ws), "--member-id", "caseA", "--age", "67", "--sex", "f", "--sample-date", "2026-09-10",
       "--longevity-skills", LIB)
    ok("intake", str(ws))
    st = json.loads((ws / "state.json").read_text())
    beta = next(f["id"] for f in st["files"] if f["kind"] == "methylation_beta")
    ok("assign", str(ws), "--file", beta, "--what", "tissue", "--value", "whole_blood", "--reason", "GEO: whole blood")
    ok("assign", str(ws), "--what", "identity", "--value", "consistent", "--reason", "test fixture treated as one person")
    capsys.readouterr()
    ok("labs", "candidates", str(ws))
    pend = json.loads(capsys.readouterr().out)["pending"]
    assert all(p["collected_at"] == "2026-09-10" for p in pend)            # no date column: the workspace date
    ans = ws / "lab_answers.json"
    ans.write_text(json.dumps({"answers": [{"row_key": p["row_key"], "answer": "yes", "why": "fixture"} for p in pend]}))
    ok("labs", "confirm", str(ws), "--answers", str(ans), "--reason", "agent reviewed each row")
    ok("preflight", str(ws), "--no-network")
    ok("member", str(ws), "smoker=no", "treated=no", "diabetes=no", "north=no", "genetic_disclosure=yes", "--source", "fixture")
    ok("methods", "plan", str(ws))
    ok("methods", "run", str(ws))
    ok("integrate", "bundle", str(ws))
    st = json.loads((ws / "state.json").read_text())
    ana = ws / "work" / "integrate" / "analyses"
    ana.mkdir(parents=True, exist_ok=True)
    for b in st["integrate"]["bundles"]:
        rid = json.loads(Path(b["bundle"]).read_text())["readouts"][0]["id"]
        (ana / f"{b['system']}.md").write_text(
            f"## 结论\n主要读数见下。\n\n## 证据\n- {{{{r:{rid}|label}}}} = {{{{r:{rid}}}}}\n\n"
            "## 不确定性\n单次测量。\n\n## 建议复测\n{{n:12 周后}}复测。\n", encoding="utf-8")
        ok("integrate", "register", str(ws), "--system", b["system"])
    ok("organ", "bundle", str(ws))
    for organ in C.data("organs.json")["organs"]:
        ok("organ", "skip", str(ws), "--organ", organ, "--reason", "fixture: no relevant signal")
    ok("insights", "skip-genomics", str(ws), "--reason", "offline test run")
    ok("insights", "position", str(ws))
    bdir = ws / "work" / "insights" / "board"
    bdir.mkdir(parents=True, exist_ok=True)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": t, "hypothesis_zh": "待检验的假设。", "basis": [b], "why_zh": "来自本人数据。"}
                        for i, (t, b) in enumerate([("表观年龄与血液年龄为何不一致", "epiage.dnam_hannum"),
                                                    ("血脂偏高是否有遗传因素", "高密度脂蛋白胆固醇(HDL-C)"),
                                                    ("肠道菌群是否偏离健康人群", "native.gut.gmhi")], 1)]}
    (bdir / "q.json").write_text(json.dumps(qs, ensure_ascii=False), encoding="utf-8")
    ok("board", "questions", str(ws), "--file", str(bdir / "q.json"))
    for i in (1, 2, 3):
        ok("board", "skip", str(ws), "--id", f"Q{i}", "--reason", "offline test")
    plan = {"items": [
        {"id": "I1", "category": "diet", "action_zh": "在营养师指导下尝试适度热量限制", "targets": ["epiage.dnam_stoch"],
         "evidence": [{"type": "effects", "ref": "calerie-cr-dunedinpace"}], "evidence_grade": "human_rct",
         "executor": "nutritionist", "retest": {"what": "甲基化时钟", "after_weeks": 52}},
        {"id": "I2", "category": "diet", "action_zh": "工作日晚饭少放盐", "targets": ["收缩压"],
         "evidence": [{"type": "effects", "ref": "salt-reduction-sbp"}], "evidence_grade": "human_rct",
         "executor": "nutritionist", "retest": {"what": "收缩压", "after_weeks": 2}}]}
    (ws / "work" / "intervene").mkdir(parents=True, exist_ok=True)
    pp = ws / "work" / "intervene" / "plan.json"
    pp.write_text(json.dumps(plan, ensure_ascii=False), encoding="utf-8")
    capsys.readouterr()
    ok("intervene", "register", str(ws), "--plan", str(pp))
    return json.loads(capsys.readouterr().out)


def _earlier_twin(path: Path, t: str, glucose: str, sbp: str, calls: dict, conc_class: str):
    """A SYNTHETIC earlier visit of the same member."""
    doc = {"schema": "la-twin/2", "member": {"id": "caseA", "age": 66, "sex": "female"}, "longevity_skills_home": LIB,
           "snapshot": {"t": t}, "readouts": [], "exposures": [],
           "observations": [{"t": t, "t_source": "row", "kind": "lab", "marker": "空腹血糖(GLU)", "value": glucose,
                             "unit": "mmol/L", "ref_range": "3.9-6.1", "provenance_uncertain": False},
                            {"t": t, "t_source": "row", "kind": "lab", "marker": "收缩压", "value": sbp, "unit": "mmHg",
                             "provenance_uncertain": False}],
           "concordance": [{"t": t, "analyte": "glucose", "genetic_band": "elevated", "genetic_pct": 85, "lab_state": "in_range",
                            "class": conc_class, "gwas_retrieved_on": "2026-09-30", "trait_map_version": "la-traits/1@synthetic"}],
           "genotype_fingerprint": {"sites": len(calls), "sha256": "synthetic", "calls": calls}}
    path.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
    return path


@needs_lib
def test_chain_export_v2_series_exposures_and_privacy(tmp_path, capsys, monkeypatch):
    ws = tmp_path / "ws"
    reg = _chain_to_intervene(ws, capsys, monkeypatch)
    rt = {r["item"]: r for r in reg["retests"]}
    assert rt["I1"]["judgeable"] is False and rt["I2"]["judgeable"] is True          # clock without noise table; SBP has one
    assert any("I2" in w and "sooner" in w for w in reg["warnings"])                   # 2 weeks < 28 days
    assert any("I1" in w and "no within-person noise model" in w for w in reg["warnings"])
    # SYNTHETIC handover files from the genomics lane
    ins = ws / "work" / "insights"
    (ins / "fingerprint.json").write_text(json.dumps({"schema": "la-fingerprint/1", "sites": len(CALLS), "sha256": "synthetic", "calls": CALLS}))
    (ins / "concordance.json").write_text(json.dumps({"schema": "la-concordance/1", "rules": {"band_high_pct": 80, "band_low_pct": 20},
        "retrieved_on": "2026-09-30", "trait_map_version": "la-traits/1@synthetic",
        "items": [{"analyte": "glucose", "adverse_direction": "high", "genetic_pct": 85, "genetic_band": "elevated", "lab_flag": "high",
                   "lab_state": "adverse", "class": "concordant_adverse", "near_band_edge": False, "clinvar_override": False,
                   "wording_zh": "这一项有较明显的先天倾向，而且已经表现出来", "basis": ["work/insights/genotype_phenotype.json#glucose"]}]},
        ensure_ascii=False), encoding="utf-8")
    (ws / "work" / "surveillance").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "surveillance" / "plan.json").write_text(json.dumps({"schema": "la-surveillance-plan/1", "items": [
        {"id": "SV1", "marker": "空腹血糖", "interval_weeks": 52, "once": False, "done": False, "reason": "genetic",
         "rule_ids": ["S-synthetic"], "executor": "physician", "note_zh": "合成测试规则", "source": {"title": "synthetic"}},
        {"id": "SV2", "marker": "脂蛋白(a)", "interval_weeks": None, "once": True, "done": True, "reason": "guideline",
         "rule_ids": ["S-lpa-once"], "executor": "physician", "note_zh": "合成测试规则", "source": {"title": "synthetic"}}]},
        ensure_ascii=False), encoding="utf-8")
    exp = {"schema": "la-exposures/1", "member": "caseA", "source": "longpi_site", "generated_at": "2026-09-09T10:00:00+08:00",
           "items": [{"plan_item": "I2", "origin_snapshot": "2025-09-10", "label_zh": "工作日晚饭少放盐", "category": "diet",
                      "targets": ["收缩压"], "status": "adopted", "adopted_at": "2026-03-02", "status_since": "2026-03-02",
                      "weekly_checkins": [{"week_start": "2026-03-02", "count": 4}, {"week_start": "2026-03-09", "count": 5}],
                      "cumulative": 9}]}
    ef = tmp_path / "exposures.json"
    ef.write_text(json.dumps(exp, ensure_ascii=False), encoding="utf-8")
    bad = dict(exp, items=[dict(exp["items"][0], note="free text")])
    bf = tmp_path / "bad.json"
    bf.write_text(json.dumps(bad, ensure_ascii=False), encoding="utf-8")
    assert la.main(["twin", "build", str(ws), "--exposures", str(bf)]) == C.EXIT_INPUT       # unknown field
    ok("twin", "build", str(ws), "--exposures", str(ef))
    tw = json.loads((ws / "work" / "twin" / "twin.json").read_text())
    assert tw["schema"] == "la-twin/2" and tw["genotype_fingerprint"]["calls"] == CALLS
    assert {i["id"]: i["status"] for i in tw["interventions"]} == {"I1": "proposed", "I2": "adopted"}
    assert next(i for i in tw["interventions"] if i["id"] == "I2")["since"] == "2026-03-02"
    assert all(o["t"] == "2026-09-10" and o["t_source"] == "workspace" for o in tw["observations"])
    assert {r["intervention"]: r["judgeable"] for r in tw["retest_plan"]} == {"I1": False, "I2": True}
    # three SYNTHETIC earlier visits; the 2024 one is another person's sample (three fingerprint sites differ)
    swapped = dict(CALLS, rs4149056="CC", rs1801133="GG", rs671="AA")
    snaps = [_earlier_twin(tmp_path / "t2023.json", "2023-09-10", "4.4", "150", CALLS, "latent"),
             _earlier_twin(tmp_path / "t2024.json", "2024-09-10", "6.9", "150", swapped, "latent"),
             _earlier_twin(tmp_path / "t2024b.json", "2024-09-20", "4.7", "152", CALLS, "latent"),
             _earlier_twin(tmp_path / "t2025.json", "2025-09-10", "5.1", "150", CALLS, "latent")]
    capsys.readouterr()
    ok("twin", "series", "--snapshots", *map(str, snaps), "--ws", str(ws))
    ser = json.loads(capsys.readouterr().out)
    assert ser["written"].endswith("series.json")
    assert [a["snapshot"] for a in ser["identity_alerts"]] == [str(snaps[1])]
    assert [s["excluded"] for s in ser["snapshots"]] == [False, True, False, False, False]
    g = next(m for m in ser["markers"] if m["marker_key"] == "glucose")
    assert [p["value"] for p in g["points"]] == [4.4, 4.7, 5.1, 5.4] and g["trend"]["verdict"] == "rising_beyond_noise"
    sbp = next(m for m in ser["markers"] if m["marker_key"] == "sbp")
    assert sbp["summary_verdict"] == "change_beyond_noise"
    assert sbp["concurrent_exposures"] == [{"plan_item": "I2", "label_zh": "工作日晚饭少放盐", "status": "adopted",
                                            "adopted_at": "2026-03-02", "cumulative_in_window": 9, "cumulative_total": 9}]
    assert [(t["analyte"], t["kind"]) for t in ser["transitions"]] == [("glucose", "tendency_expressed")]
    assert la.main(["review", "trace", str(ws)]) == C.EXIT_BLOCKED       # predictions predate the exposures and the series
    ok("predict", "response", str(ws))
    # trace, review, report, export
    (ws / "work" / "report").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "report" / "summary.md").write_text("{{r:accelerated-biological-aging-risk.phenoage|label}}为 "
                                                      "{{r:accelerated-biological-aging-risk.phenoage}}。", encoding="utf-8")
    capsys.readouterr()
    ok("review", "trace", str(ws))
    trace = json.loads(capsys.readouterr().out)
    assert any(k.endswith("series.json") for k in trace["hashes"]) and any(k.endswith("exposures.json") for k in trace["hashes"])
    rev = ws / "work" / "review" / "reviewer-1.json"
    rev.write_text(json.dumps({"trace_id": trace["trace_id"], "verdict": "pass", "findings": []}))
    ok("review", "record", str(ws), "--verdict", "pass", "--findings", str(rev))
    ok("report", str(ws))
    md = (ws / "deliver" / "report.md").read_text()
    assert "## 复测对比与趋势" in md and "出现了超出你个人正常波动的缓慢上升，下一次复测可以确认。" in md
    assert "这段时间你同期在做：工作日晚饭少放盐（累计 9 次）。" in md
    assert "这一项复测时只能并排看两次结果，不能判断是不是真的变化。" in md
    assert "两次检测的胚系基因型不一致" in md
    assert md.index("## 复测对比与趋势") < md.index("## 边界")
    for word in ("因为", "由于", "证明", "恶化", "变差"):
        assert word not in md.split("## 复测对比与趋势")[1].split("## 边界")[0]
    dtw = json.loads((ws / "deliver" / "twin.json").read_text())
    assert "calls" not in dtw["genotype_fingerprint"] and dtw["genotype_fingerprint"]["sites"] == len(CALLS)
    ex = json.loads((ws / "deliver" / "la-export.json").read_text())
    assert ex["schema"] == "la-export/2" and not X.privacy_problems(ex)
    assert "rs4149056" not in json.dumps(ex)
    tr = {t["marker_key"]: t for t in ex["trends"]}
    assert tr["glucose"]["summary_verdict"] == "trend_beyond_noise" and tr["glucose"]["observed_at"] == "2026-09-10"
    assert tr["sbp"]["concurrent_zh"] == "这段时间你同期在做：工作日晚饭少放盐（累计 9 次）。"
    assert ex["concordance"][0]["class"] == "concordant_adverse" and "basis" not in ex["concordance"][0]
    assert ex["concordance_note_zh"] == "基因分值只能解释这一项差异的一部分。"
    assert ex["concordance_transitions"][0]["kind"] == "tendency_expressed" and ex["identity_alerts"]
    assert [s["id"] for s in ex["surveillance"]] == ["SV1", "SV2"]
    assert [(p["plan_item"], p["label"]["basis_kind"]) for p in ex["predictions"]] == [("I2", "statistical_model")]
    assert ex["diagnoses"] == [] and ex["medications"] == []
    rts = {r["item"]: r for r in ex["retests"]}
    assert rts["I1"]["judgeable"] is False and rts["I1"]["reason"] == "intervention"
    assert rts["I1"]["due_at"] == "2027-09-09T00:00:00+08:00"                     # sample date + 52 weeks, China time
    assert rts["I2"]["id"] == "I2" and rts["I2"]["marker"] == "收缩压" and rts["SV1"]["marker"] == "空腹血糖"
    assert ex["concordance"][0]["transition"] == {"from": "latent", "to": "concordant_adverse", "observed_at": "2026-09-10",
                                                  "kind": "tendency_expressed"}
    assert {x["source_item"]: x["targets"] for x in ex["plan"]["items"]} == {"I1": ["epiage.dnam_stoch"], "I2": ["收缩压"]}
    assert rts["SV1"]["reason"] == "genetic" and rts["SV1"]["judgeable"] is True and "SV2" not in rts
    capsys.readouterr()
    ok("export", str(ws))
    assert json.loads(capsys.readouterr().out)["trends"] == len(ex["trends"])
    ok("export", str(ws), "--schema", "la-export/1")
    v1 = json.loads((ws / "work" / "export" / "la-export-1.json").read_text())
    assert v1["schema"] == "la-export/1" and "trends" not in v1 and set(v1["retests"][0]) == {"item", "what", "after_weeks", "due", "from"}
    # an edit after delivery to a bound handover file makes the delivery stale
    sp = ws / "work" / "surveillance" / "plan.json"
    sp.write_text(sp.read_text() + " ")
    assert la.main(["export", str(ws)]) == C.EXIT_BLOCKED
    # stage clock
    st = json.loads((ws / "state.json").read_text())
    assert st["timing"]["intake"]["started_at"] and st["timing"]["report"]["ended_at"]
    capsys.readouterr()
    ok("observe", str(ws))
    assert "report" in json.loads(capsys.readouterr().out)["timing"]


@needs_lib
def test_export_refuses_a_genotype_in_a_handover_file(tmp_path):
    ws = tmp_path
    (ws / "work" / "clinical").mkdir(parents=True)
    (ws / "work" / "clinical" / "medications.json").write_text(json.dumps({"schema": "la-medications/1", "items": [
        {"id": "RX1", "kind": "cpic", "text_zh": "rs4149056 TT：按常规剂量", "label": {"basis_kind": "cpic"}}]}, ensure_ascii=False))
    doc = X._v2(ws, {})
    assert X.privacy_problems(doc)
