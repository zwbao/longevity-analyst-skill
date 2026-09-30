"""Whole chain through the CLI on real public data (examples/case-A), no network.

Catches wiring bugs unit tests miss (stage bookkeeping vs. content binding, derived files, placeholders).
The plan uses an effects-table id instead of a PubMed id so the test runs offline.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
import la  # noqa: E402
from lalib import common as C  # noqa: E402

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
pytestmark = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")


def ok(*args):
    rc = la.main(list(args))
    assert rc == 0, args


def test_full_chain_case_a(tmp_path, capsys, monkeypatch):
    ws = tmp_path / "ws"
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
    assert pend and any(p["candidates"] for p in pend)                          # every numeric row is asked
    ans = ws / "lab_answers.json"
    ans.write_text(json.dumps({"answers": [{"row_key": p["row_key"], "answer": "yes", "why": "fixture: serum, fasting"} for p in pend]}))
    assert la.main(["preflight", str(ws), "--no-network"]) != 0                    # intake not done until rows are confirmed
    ok("labs", "confirm", str(ws), "--answers", str(ans), "--reason", "agent reviewed each row")
    ok("preflight", str(ws), "--no-network")
    ok("member", str(ws), "smoker=no", "treated=no", "diabetes=no", "north=no", "genetic_disclosure=yes", "--source", "fixture")
    ok("methods", "plan", str(ws))
    ok("methods", "run", str(ws))
    readouts = json.loads((ws / "work" / "readouts.json").read_text())["readouts"]
    ids = {r["id"] for r in readouts}
    assert "accelerated-biological-aging-risk.phenoage" in ids and "epiage.dnam_hannum" in ids and "native.apoe.genotype" in ids
    assert not any("grimage" in i or "dunedin" in i for i in ids)            # commercial mode
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
    (ws / "work" / "evidence").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "evidence" / "pubmed_fixture.json").write_text(json.dumps({"items": [{"pmid": "34554658"}]}))
    monkeypatch.setattr(la.evidence, "verify_pmids", lambda ids: set(ids))
    sha = {b["organ"]: b["bundle_sha256"] for b in json.loads((ws / "state.json").read_text())["organs"]["bundles"]}
    est = {"organ": "kidney", "method": "llm_estimate", "confidence": "low", "bundle_sha256": sha["kidney"], "age_estimate": None,
           "age_estimate_null_reason": "无肾脏相关读数",
           "disease_risks": [{"disease": "慢性肾脏病 3 期及以上", "horizon_years": 10, "low": 0.04, "point": 0.08, "high": 0.16,
                              "basis": ["member.age"], "evidence": [{"type": "pubmed", "ref": "34554658"}],
                              "rationale_zh": "以队列基线发病率为锚，按同龄女性估计。"}]}
    (ws / "work" / "organs" / "estimates").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "organs" / "estimates" / "kidney.json").write_text(json.dumps(est, ensure_ascii=False), encoding="utf-8")
    assert la.main(["intervene", "register", str(ws), "--plan", "x"]) != 0          # organs stage not done yet
    bad = dict(est, disease_risks=[dict(est["disease_risks"][0], rationale_zh="按{{n:67 岁}}女性估计，{{r:0.08}}")])
    (ws / "work" / "organs" / "estimates" / "kidney.json").write_text(json.dumps(bad, ensure_ascii=False), encoding="utf-8")
    assert la.main(["organ", "register", str(ws), "--organ", "kidney"]) == C.EXIT_INPUT   # own age as literal, value as id
    (ws / "work" / "organs" / "estimates" / "kidney.json").write_text(json.dumps(est, ensure_ascii=False), encoding="utf-8")
    ok("organ", "register", str(ws), "--organ", "kidney")
    for organ in C.data("organs.json")["organs"]:
        if organ != "kidney":
            ok("organ", "skip", str(ws), "--organ", organ, "--reason", "fixture: no relevant signal")
    # insights: genomics needs live GWAS/ClinVar (tested in test_insights_live.py); position is offline
    assert la.main(["intervene", "register", str(ws), "--plan", "x"]) != 0          # insights stage not done yet
    ok("insights", "skip-genomics", str(ws), "--reason", "offline test run")
    ok("insights", "position", str(ws))
    pos = json.loads((ws / "work" / "insights" / "positions.json").read_text())
    assert any(x["analyte"] == "creatinine" for x in pos["labs"]) and len(pos["ages"]) >= 3
    bdir = ws / "work" / "insights" / "board"
    bdir.mkdir(parents=True, exist_ok=True)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": t, "hypothesis_zh": "待检验的假设。", "basis": [b], "why_zh": "来自本人数据。"}
                        for i, (t, b) in enumerate([("表观年龄与血液年龄为何不一致", "epiage.dnam_hannum"),
                                                    ("血脂偏高是否有遗传因素", "高密度脂蛋白胆固醇(HDL-C)"),
                                                    ("肠道菌群是否偏离健康人群", "native.gut.gmhi")], 1)]}
    (bdir / "q.json").write_text(json.dumps(qs, ensure_ascii=False), encoding="utf-8")
    bad = dict(qs, questions=qs["questions"][:2])
    (bdir / "q2.json").write_text(json.dumps(bad, ensure_ascii=False), encoding="utf-8")
    assert la.main(["board", "questions", str(ws), "--file", str(bdir / "q2.json")]) == C.EXIT_INPUT   # 3-10 questions
    ok("board", "questions", str(ws), "--file", str(bdir / "q.json"))
    for i in (1, 2):
        (bdir / f"Q{i}.json").write_text(json.dumps({"id": f"Q{i}", "verdict": "insufficient", "confidence": "low",
            "member_evidence": [qs["questions"][i - 1]["basis"][0]], "public_evidence": [],
            "summary_zh": "现有数据不足以判断。", "next_step_zh": "复测后再评估。", "limitations_zh": "只有一次测量。"}, ensure_ascii=False), encoding="utf-8")
        ok("board", "finding", str(ws), "--id", f"Q{i}")
    (bdir / "Q3.json").write_text(json.dumps({"id": "Q3", "verdict": "supported", "confidence": "high", "member_evidence": ["native.gut.gmhi"],
        "public_evidence": ["gwas:made-up"], "summary_zh": "你的指数低于健康人群 3 个百分位。", "next_step_zh": "x", "limitations_zh": "x"}, ensure_ascii=False), encoding="utf-8")
    assert la.main(["board", "finding", str(ws), "--id", "Q3"]) == C.EXIT_INPUT     # confidence, unretrieved ref, bare digit
    ok("board", "skip", str(ws), "--id", "Q3", "--reason", "offline test")
    assert json.loads((ws / "state.json").read_text())["stages"]["insights"] == "done"
    plan = {"items": [{"id": "I1", "category": "diet", "action_zh": "在营养师指导下尝试适度热量限制",
                       "targets": ["epiage.dnam_hannum"], "evidence": [{"type": "effects", "ref": "calerie-cr-dunedinpace"}],
                       "evidence_grade": "human_rct", "executor": "nutritionist", "retest": {"what": "甲基化时钟", "after_weeks": 52}}]}
    (ws / "work" / "intervene").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "intervene" / "plan.json").write_text(json.dumps(plan, ensure_ascii=False), encoding="utf-8")
    ok("intervene", "register", str(ws), "--plan", str(ws / "work" / "intervene" / "plan.json"))
    ok("twin", "build", str(ws))
    (ws / "work" / "report").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "report" / "summary.md").write_text(
        "{{r:accelerated-biological-aging-risk.phenoage|label}}为 {{r:accelerated-biological-aging-risk.phenoage}}。"
        "{{r:organ.kidney.risk.1|label}}：{{r:organ.kidney.risk.1}}。", encoding="utf-8")
    capsys.readouterr()
    ok("review", "trace", str(ws))
    trace = json.loads(capsys.readouterr().out)
    rev = ws / "work" / "review" / "reviewer-1.json"
    rev.write_text(json.dumps({"trace_id": trace["trace_id"], "verdict": "pass", "findings": []}))
    ok("review", "record", str(ws), "--verdict", "pass", "--findings", str(rev))
    ok("report", str(ws))
    ok("validate", str(ws))
    md = (ws / "deliver" / "report.md").read_text()
    assert "表型年龄" in md and "12 周后" in md and "{{" not in md
    assert "## 你在同龄人群中的位置" in md and "## 问题看板" in md and "证据不足" in md
    assert "## 器官体检表" in md and "AI 估计" in md and "8%（AI 估计，区间 4%–16%，10 年）" in md
    tw = json.loads((ws / "deliver" / "twin.json").read_text())
    assert [e["id"] for e in tw["organ_estimates"]] == ["organ.kidney.risk.1"]
    # tampering after delivery is detected
    (ws / "deliver" / "report.md").write_text(md + "x")
    assert la.main(["validate", str(ws)]) == C.EXIT_INPUT
    # a later upstream change moves the stale deliverable aside
    ok("member", str(ws), "smoker=yes", "--source", "fixture correction")
    assert not (ws / "deliver").exists() and not (ws / ".stale").exists()     # removed, not kept
