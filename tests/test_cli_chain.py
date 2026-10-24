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
    ex = json.loads((ws / "deliver" / "la-export.json").read_text())
    assert ex["schema"] == "la-export/2" and ex["member"]["age"]
    assert all("{{" not in (x.get("title") or "") + (x.get("detail") or "") for x in ex["plan"]["items"])
    assert all(x["category"] in ("diet", "exercise", "sleep", "supplement", "behavior", "other") for x in ex["plan"]["items"])
    assert ex["readouts"] and ex["report"]["html"].endswith("report.html")
    assert all("overrides" in o for o in ex["organs"]) and all("](" not in x["title"] + x["detail"] for x in ex["plan"]["items"])
    assert all(len(x["title"]) <= 60 for x in ex["plan"]["items"])
    confirmed = {r["marker"] for r in __import__("lalib.labnames", fromlist=["usable_rows"]).usable_rows(json.loads((ws / "state.json").read_text()))}
    assert all(m in confirmed for x in ex["plan"]["items"] for m in x["markers"])      # retests tracked by the member's lab names
    capsys.readouterr()
    ok("export", str(ws))
    assert f'"plan_items": {len(ex["plan"]["items"])}' in capsys.readouterr().out
    plan_p = ws / "work" / "intervene" / "plan.json"
    saved = plan_p.read_text()
    plan_p.write_text(saved.replace('"I1"', '"I1" ', 1))                   # an input changed after the report
    assert la.main(["export", str(ws)]) == C.EXIT_BLOCKED
    plan_p.write_text(saved)
    assert [e["id"] for e in tw["organ_estimates"]] == ["organ.kidney.risk.1"]
    # tampering after delivery is detected
    (ws / "deliver" / "report.md").write_text(md + "x")
    assert la.main(["validate", str(ws)]) == C.EXIT_INPUT
    # a later upstream change moves the stale deliverable aside
    ok("member", str(ws), "smoker=yes", "--source", "fixture correction")
    assert not (ws / "deliver").exists() and not (ws / ".stale").exists()     # removed, not kept


# ================================================================== case-C: four visits, one member (PLAN WS-J.1, appendix C)
# Labs, home blood pressure, plans and exposures are SYNTHETIC (examples/case-C); the genotypes are real 1000 Genomes
# CHB/CHS samples, small regions only (examples/README.md). Public-database answers (GWAS Catalog, myvariant.info,
# Ensembl, ClinVar) were recorded once from the live services into tests/fixtures/case_c/pubdata (see its README) and
# are served through pubdata's on-machine cache; nothing in this test touches the network.
CASE_C = ROOT / "examples" / "case-C"
CASE_C_SWAP = ROOT / "examples" / "case-C-swap"
PUBDATA = ROOT / "tests" / "fixtures" / "case_c" / "pubdata"
LPA_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC9639807/"
LPA_TITLE = ("Lipoprotein(a) in atherosclerotic cardiovascular disease and aortic stenosis: a European Atherosclerosis "
             "Society consensus statement - PMC")                   # the page title as fetched on 2026-10-24
VISITS = {"2023-09": ("2023-09-04", 49), "2024-09": ("2024-09-09", 50), "2025-09": ("2025-09-15", 51), "2026-09": ("2026-09-15", 52)}


class _Page:
    status = 200

    def __init__(self, body: bytes):
        self.body = body

    def read(self, n=-1):
        return self.body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _offline(monkeypatch, record: Path = None):
    """Public data from the recorded cache only; the one guideline page `evidence fetch` reads, served with its title."""
    from lalib import pubdata
    import urllib.request
    monkeypatch.setattr(pubdata, "SHARED", record or PUBDATA)
    monkeypatch.setattr(pubdata, "SHARED_TTL_DAYS", 10 ** 6)
    if record is None:
        def no_net(url, *a, **k):
            raise AssertionError(f"case-C test went to the network: {url}")
        monkeypatch.setattr(pubdata, "_http", no_net)
    real = urllib.request.urlopen

    def page(req, *a, **k):
        url = getattr(req, "full_url", req)
        if url == LPA_URL:
            return _Page(f"<html><head><title>{LPA_TITLE}</title></head></html>".encode())
        if record is not None:
            return real(req, *a, **k)
        raise AssertionError(f"case-C test went to the network: {url}")
    monkeypatch.setattr(urllib.request, "urlopen", page)


def run_json(capsys, *args, rc=0):
    capsys.readouterr()
    got = la.main(list(args))
    cap = capsys.readouterr()
    if got != rc:
        sys.__stderr__.write(cap.out[-6000:] + cap.err[-3000:])
    assert got == rc, (args, got)
    out = cap.out
    return json.loads(out) if out.strip() else None


def _visit(tmp: Path, ym: str, capsys, *, raw: Path = None, prev: Path = None, history=(), plan: Path = None,
           exposures: Path = None, dx=("ada",), summary=None, until: str = "report") -> Path:
    """One visit of the case-C member through the CLI, no pipelines. Returns the workspace."""
    date, age = VISITS[ym]
    ws = tmp / ym
    raw = raw or CASE_C / ym / "raw"
    run_json(capsys, "init", str(raw), str(ws), "--member-id", "caseC", "--age", str(age), "--sex", "m", "--sample-date", date,
             "--longevity-skills", LIB, *(["--reuse-germline-from", str(prev)] if prev else []))
    ok("intake", str(ws))
    ok("assign", str(ws), "--what", "identity", "--value", "consistent", "--reason", "case-C: one member's files (fixture)")
    pend = run_json(capsys, "labs", "candidates", str(ws))["pending"]
    ans = ws / "lab_answers.json"
    ans.write_text(json.dumps({"answers": [{"row_key": p["row_key"], "answer": "yes", "why": "fixture: serum, fasting; home BP mornings"}
                                           for p in pend]}))
    ok("labs", "confirm", str(ws), "--answers", str(ans), "--reason", "agent reviewed each row")
    ok("preflight", str(ws), "--no-network")
    ok("member", str(ws), "smoker=no", "treated=no", "diabetes=no", "north=no", "genetic_disclosure=yes", "residence=CN",
       "--source", "fixture: member questionnaire")
    ok("methods", "plan", str(ws))
    ok("methods", "run", str(ws))
    ok("integrate", "bundle", str(ws))
    st = json.loads((ws / "state.json").read_text())
    ana = ws / "work" / "integrate" / "analyses"
    ana.mkdir(parents=True, exist_ok=True)
    for b in st["integrate"]["bundles"]:
        rid = json.loads(Path(b["bundle"]).read_text())["readouts"][0]["id"]
        (ana / f"{b['system']}.md").write_text(
            f"## 结论\n这一组本期的读数见下；化验的多次变化见「复测对比与趋势」。\n\n## 证据\n- {{{{r:{rid}|label}}}} = {{{{r:{rid}}}}}\n\n"
            "## 不确定性\n这些读数只来自本期数据；化验有没有真实变化，以复测对比为准。\n\n## 建议复测\n按方案复测。\n", encoding="utf-8")
        ok("integrate", "register", str(ws), "--system", b["system"])
    ok("organ", "bundle", str(ws))
    for organ in C.data("organs.json")["organs"]:
        ok("organ", "skip", str(ws), "--organ", organ, "--reason", "本期以复查化验为主，没有做器官 AI 估计")
    run_json(capsys, "insights", "genomics", str(ws))
    run_json(capsys, "insights", "pgx", str(ws))
    ok("insights", "position", str(ws))
    bdir = ws / "work" / "insights" / "board"
    bdir.mkdir(parents=True, exist_ok=True)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": t, "hypothesis_zh": "待检验的假设。", "basis": [b], "why_zh": "来自本人数据。"}
                        for i, (t, b) in enumerate([("空腹血糖是否在慢慢升高", "空腹血糖(GLU)"), ("尿酸偏高与先天倾向的关系", "尿酸(UA)"),
                                                    ("转氨酶的变化与运动有没有同期关系", "谷丙转氨酶(ALT)")], 1)]}
    (bdir / "q.json").write_text(json.dumps(qs, ensure_ascii=False), encoding="utf-8")
    ok("board", "questions", str(ws), "--file", str(bdir / "q.json"))
    for i in (1, 2, 3):
        ok("board", "skip", str(ws), "--id", f"Q{i}", "--reason", "本期没有安排研究员，留到下次")
    ok("evidence", "fetch", str(ws), "--url", LPA_URL)
    ok("intervene", "register", str(ws), "--plan", str(plan or CASE_C / ym / "plan.json"))
    ok("surveillance", "plan", str(ws), *(["--history", *map(str, history)] if history else []))
    ok("twin", "build", str(ws), *(["--exposures", str(exposures)] if exposures else []))
    if until == "twin":
        return ws
    if history:
        ok("twin", "series", "--snapshots", *map(str, history), "--ws", str(ws))
    ok("predict", "response", str(ws), "--repeats", "14")
    runs = {gs: run_json(capsys, "clinical", "dx", str(ws), "--guideline-set", gs) for gs in dx}
    (tmp / f"{ym}-dx-runs.json").write_text(json.dumps(runs, ensure_ascii=False), encoding="utf-8")
    ok("clinical", "rx", str(ws))
    (ws / "work" / "report").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "report" / "summary.md").write_text(summary(ws) if summary else "本次报告的读数和条目见下文。", encoding="utf-8")
    trace = run_json(capsys, "review", "trace", str(ws))
    rev = ws / "work" / "review" / "reviewer-1.json"
    rev.write_text(json.dumps({"trace_id": trace["trace_id"], "verdict": "pass", "findings": []}))
    ok("review", "record", str(ws), "--verdict", "pass", "--findings", str(rev))
    ok("report", str(ws))
    ok("export", str(ws))
    assert run_json(capsys, "validate", str(ws))["ok"] is True
    return ws


def _twin(ws: Path) -> Path:
    return ws / "work" / "twin" / "twin.json"


def _load(ws: Path, *parts) -> dict:
    return json.loads(ws.joinpath("work", *parts).read_text(encoding="utf-8"))


def _summary_v4(ws: Path) -> str:
    """What an agent would write for visit 4: every member number through a placeholder (a SYNTHETIC member)."""
    pr = next(p for p in _load(ws, "predict", "predictions.json")["items"] if p["plan_item"] == "I3")
    dxi = _load(ws, "clinical", "diagnoses.json")["items"]
    rx = next(x for x in _load(ws, "clinical", "medications.json")["items"] if x["kind"] == "cpic" and x["gene"] == "SLCO1B1")
    lines = ["这几年的体检连起来看：空腹血糖每次都在参考范围内，但几次结果连起来出现了超出你个人正常波动的缓慢上升；低密度脂蛋白胆固醇一直在波动范围内。",
             "",
             "家庭血压在减盐期间下降，超出了你的正常波动；这是同期发生的变化，不能说明是减盐造成的。下一次同样设计的减盐实验："
             f"{{{{p:{pr['id']}}}}}。"]
    if dxi:
        lines += ["", f"诊断性结论：{{{{dx:{dxi[0]['id']}}}}}。请带着这一条去看医生。"]
    lines += ["", f"用药基因（他汀类降脂药）：{{{{rx:{rx['id']}}}}}。这一条只在医生考虑用这类药时才用得上。",
              "", "尿酸这一项有先天倾向，从上一年那次检测起已经表现出来，这次仍高于参考范围，已列入报告末尾的问医生的单子。"]
    return "\n".join(lines)


def _case_c_chain(tmp: Path, capsys):
    vcf1 = next((CASE_C / "2023-09" / "raw").glob("*.vcf.gz"))
    v1 = _visit(tmp, "2023-09", capsys)
    v2 = _visit(tmp, "2024-09", capsys, prev=v1, history=[_twin(v1)], exposures=CASE_C / "2024-09" / "exposures.json")
    v3 = _visit(tmp, "2025-09", capsys, prev=v2, history=[_twin(v1), _twin(v2)], exposures=CASE_C / "2025-09" / "exposures.json")
    v4 = _visit(tmp, "2026-09", capsys, prev=v3, history=[_twin(v1), _twin(v2), _twin(v3)],
                exposures=CASE_C / "2026-09" / "exposures.json", dx=("cn", "ada"), summary=_summary_v4)
    return vcf1, v1, v2, v3, v4


# ================================================================== case-B: PDF labs, MS proteome, unknown reads (release gate 4)
CASE_B_PDF_ROWS = """marker,value,unit,ref_range,page
白蛋白 ALB,39.2,g/L,40-55,1
肌酐 CREA,1.12,mg/dL,0.7-1.3,1
葡萄糖 GLU,6.3,mmol/L,3.9-6.1,1
C反应蛋白 CRP,4.6,mg/L,0-10,1
淋巴细胞% LYM%,24.0,%,20-50,1
平均红细胞体积 MCV,95.4,fL,82-100,1
红细胞分布宽度 RDW-CV,14.2,%,11.5-14.5,1
碱性磷酸酶 ALP,96,U/L,45-125,1
白细胞 WBC,6.9,10^9/L,3.5-9.5,1
总胆固醇 TC,4.9,mmol/L,<5.2,1
高密度脂蛋白 HDL-C,1.05,mmol/L,>1.0,1
甘油三酯 TG,2.1,mmol/L,<1.7,1
糖化血红蛋白 HbA1c,6.1,%,4-6,1
收缩压,136,mmHg,,1
舒张压,84,mmHg,,1
腰围,94,cm,,1
"""      # page 1 of examples/case-B/raw/2026年度体检报告_检验部分.pdf, transcribed as printed (synthetic report)


def test_full_chain_case_b(tmp_path, capsys):
    ws = tmp_path / "ws"
    ok("init", str(ROOT / "examples" / "case-B" / "raw"), str(ws), "--member-id", "caseB", "--age", "73", "--sex", "m",
       "--sample-date", "2026-09-10", "--longevity-skills", LIB)
    ok("intake", str(ws))
    st = json.loads((ws / "state.json").read_text())
    fid = {f["name"]: f["id"] for f in st["files"]}
    csv = tmp_path / "pdf_rows.csv"
    csv.write_text(CASE_B_PDF_ROWS, encoding="utf-8")
    ok("labs", "add", str(ws), "--csv", str(csv), "--file", fid["2026年度体检报告_检验部分.pdf"], "--reason", "transcribed page 1 as printed")
    prot = fid["plasma_proteome_DIA_matrix.tsv"]
    ok("assign", str(ws), "--file", prot, "--what", "platform", "--value", "ms_dia", "--reason", "file header: DIA-MS plasma")
    ok("assign", str(ws), "--file", prot, "--what", "value_scale", "--value", "log2_intensity", "--reason", "file header: log2 LFQ")
    ok("assign", str(ws), "--file", fid["unknown_run_R1.fastq.gz"], "--what", "modality", "--value", "unknown",
       "--reason", "no delivery note; the reads are a wiring fixture")
    ok("assign", str(ws), "--file", fid["甲基化_beta值.csv"], "--what", "tissue", "--value", "whole_blood", "--reason", "GEO GSM989833: whole blood")
    ok("assign", str(ws), "--what", "identity", "--value", "consistent", "--reason", "test fixture treated as one person")
    pend = run_json(capsys, "labs", "candidates", str(ws))["pending"]
    ans = ws / "lab_answers.json"
    ans.write_text(json.dumps({"answers": [{"row_key": p["row_key"], "answer": "yes", "why": "fixture: serum, fasting"} for p in pend]}))
    ok("labs", "confirm", str(ws), "--answers", str(ans), "--reason", "agent reviewed each row")
    ok("preflight", str(ws), "--no-network")
    ok("member", str(ws), "smoker=no", "treated=no", "diabetes=no", "north=no", "genetic_disclosure=yes", "--source", "fixture")
    ok("methods", "plan", str(ws))
    items = {i["method"]: i for i in json.loads((ws / "state.json").read_text())["methods"]["items"]}
    assert not any(i["status"] == "ready" for m, i in items.items() if "prote" in m and not m.startswith("native."))   # platform gate
    ok("methods", "run", str(ws))
    ids = {r["id"] for r in json.loads((ws / "work" / "readouts.json").read_text())["readouts"]}
    assert "native.apoe.genotype" in ids and not any("grimage" in i or "dunedin" in i for i in ids)
    ok("integrate", "bundle", str(ws))
    st = json.loads((ws / "state.json").read_text())
    ana = ws / "work" / "integrate" / "analyses"
    ana.mkdir(parents=True, exist_ok=True)
    for b in st["integrate"]["bundles"]:
        rid = json.loads(Path(b["bundle"]).read_text())["readouts"][0]["id"]
        (ana / f"{b['system']}.md").write_text(f"## 结论\n见下。\n\n## 证据\n- {{{{r:{rid}|label}}}} = {{{{r:{rid}}}}}\n\n"
                                               "## 不确定性\n单次测量。\n\n## 建议复测\n按方案复测。\n", encoding="utf-8")
        ok("integrate", "register", str(ws), "--system", b["system"])
    ok("organ", "bundle", str(ws))
    for organ in C.data("organs.json")["organs"]:
        ok("organ", "skip", str(ws), "--organ", organ, "--reason", "fixture: no organ estimates in this test")
    ok("insights", "skip-genomics", str(ws), "--reason", "offline test run")
    ok("insights", "pgx", str(ws))
    ok("insights", "position", str(ws))
    bdir = ws / "work" / "insights" / "board"
    bdir.mkdir(parents=True, exist_ok=True)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": t, "hypothesis_zh": "待检验的假设。", "basis": [b], "why_zh": "来自本人数据。"}
                        for i, (t, b) in enumerate([("血糖为何略高", "葡萄糖 GLU"), ("甘油三酯为何偏高", "甘油三酯 TG"),
                                                    ("白蛋白为何略低", "白蛋白 ALB")], 1)]}
    (bdir / "q.json").write_text(json.dumps(qs, ensure_ascii=False), encoding="utf-8")
    ok("board", "questions", str(ws), "--file", str(bdir / "q.json"))
    for i in (1, 2, 3):
        ok("board", "skip", str(ws), "--id", f"Q{i}", "--reason", "offline test")
    plan = {"items": [{"id": "I1", "category": "diet", "action_zh": "晚饭减少精米白面，增加蔬菜和豆类", "targets": ["甘油三酯 TG"],
                       "evidence": [{"type": "effects", "ref": "oat-glucose"}], "evidence_grade": "human_rct", "executor": "nutritionist",
                       "retest": {"what": "甘油三酯", "after_weeks": 12}}]}
    pp = ws / "work" / "intervene" / "plan.json"
    pp.parent.mkdir(parents=True, exist_ok=True)
    pp.write_text(json.dumps(plan, ensure_ascii=False), encoding="utf-8")
    ok("intervene", "register", str(ws), "--plan", str(pp))
    ok("twin", "build", str(ws))
    ok("clinical", "dx", str(ws))
    ok("clinical", "rx", str(ws))
    (ws / "work" / "report").mkdir(parents=True, exist_ok=True)
    (ws / "work" / "report" / "summary.md").write_text("{{r:native.apoe.genotype|label}}为 {{r:native.apoe.genotype}}。", encoding="utf-8")
    trace = run_json(capsys, "review", "trace", str(ws))
    rev = ws / "work" / "review" / "reviewer-1.json"
    rev.write_text(json.dumps({"trace_id": trace["trace_id"], "verdict": "pass", "findings": []}))
    ok("review", "record", str(ws), "--verdict", "pass", "--findings", str(rev))
    ok("report", str(ws))
    ok("export", str(ws))
    assert run_json(capsys, "validate", str(ws))["ok"] is True
    md = (ws / "deliver" / "report.md").read_text()
    assert "{{" not in md and "## 边界" in md


@pytest.fixture
def case_c_net(monkeypatch):
    rec = os.environ.get("LA_CASE_C_RECORD")          # set only when re-recording the public-data fixture (network)
    _offline(monkeypatch, Path(rec) if rec else None)


def test_full_chain_case_c(tmp_path, capsys, case_c_net):
    vcf1, v1, v2, v3, v4 = _case_c_chain(tmp_path, capsys)
    for ws in (v1, v2, v3, v4):
        assert run_json(capsys, "validate", str(ws))["ok"] is True
    # germline: one VCF at the first visit, reused afterwards (APOE readouts, genes vs labs, fingerprint, pgx)
    for ws in (v2, v3, v4):
        assert json.loads((ws / "state.json").read_text())["reuse_germline"]["fingerprint_check"] == "no_new_genotype_data"
        ro = {r["id"]: r for r in _load(ws, "readouts.json")["readouts"]}
        assert ro["native.apoe.genotype"]["reused_from"] and ro["native.apoe.genotype"]["value"] == \
            {r["id"]: r for r in _load(v1, "readouts.json")["readouts"]}["native.apoe.genotype"]["value"]
        assert any(r["id"] == "native.apoe.genotype" for r in json.loads(_twin(ws).read_text())["readouts"])
    ser = _load(v4, "twin", "series.json")
    m = {s["marker_key"]: s for s in ser["markers"] if s["marker_key"]}
    # slow drift: every adjacent pair within noise, first vs latest beyond, trend rising
    g = m["glucose"]
    assert [p["value"] for p in g["points"]] == [4.6, 5.0, 5.3, 5.6]
    assert [p["verdict"] for p in g["pairs"]] == ["within_noise"] * 3
    assert g["baseline_vs_latest"]["verdict"] == "increase_beyond_noise" and g["trend"]["verdict"] == "rising_beyond_noise"
    # stable: nothing at any level
    ldl = m["ldl"]
    assert all(p["verdict"] == "within_noise" for p in ldl["pairs"]) and ldl["baseline_vs_latest"]["verdict"] == "within_noise"
    assert ldl["trend"]["verdict"] == "no_trend_detected" and ldl["summary_verdict"] == "within_noise"
    # home blood pressure: 14 mornings each side, salt reduction concurrent
    sbp = m["sbp"]
    before, after = sbp["points"][-2], sbp["points"][-1]
    assert (before["k"], after["k"], before["value"], after["value"]) == (14, 14, 138, 131)
    last = sbp["pairs"][-1]
    assert last["verdict"] == "decrease_beyond_noise" and last["band_pct"] == [-3.1, 3.1] and last["change_pct"] == -5.1
    salt = next(e for e in sbp["concurrent_exposures"] if e["plan_item"] == "I3")
    assert salt["label_zh"] == "做饭用限盐勺" and salt["status"] == "adopted" and salt["cumulative_in_window"] > 0
    # concordance transitions: urate tendency expressed (doctor list); ALT acquired -> typical with exercise concurrent
    tr = {(t["analyte"], t["from"], t["to"]): t for t in ser["transitions"]}
    ua = tr[("urate", "latent", "concordant_adverse")]
    assert ua["kind"] == "tendency_expressed" and ua["ask_doctor"] is True and ua["to_t"] == "2025-09-15"
    alt = tr[("alt", "acquired_adverse", "concordant_typical")]
    assert alt["kind"] == "back_in_range" and alt["real_change"] is True
    assert [e["plan_item"] for e in alt["concurrent_exposures"]] == ["I2"]
    assert m["alt"]["summary_verdict"] == "change_beyond_noise"
    assert not ser["identity_alerts"] and not any(s["excluded"] for s in ser["snapshots"])
    # Lp(a): once at visit 1, then done and out of the retests
    for ws in (v1, v2, v3, v4):
        sv = {x["id"]: x for x in _load(ws, "surveillance", "plan.json")["items"]}
        assert sv["S-lpa-once"]["once"] and sv["S-lpa-once"]["done"]
        ex = json.loads((ws / "deliver" / "la-export.json").read_text())
        assert not any("脂蛋白(a)" in str(r.get("what")) for r in ex["retests"])
    # homocysteine as a retest target is judgeable now (biovar first batch landed)
    rt = {r["intervention"]: r for r in json.loads(_twin(v4).read_text())["retest_plan"]}
    assert rt["I4"]["what"] == "同型半胱氨酸" and rt["I4"]["judgeable"] is True
    # pharmacogenomics: SLCO1B1 recommendations labelled CPIC, CPIC's words verbatim
    pg = _load(v1, "pgx", "pgx.json")
    slco = next(r for r in pg["results"] if r["gene"] == "SLCO1B1")
    rx = [x for x in _load(v4, "clinical", "medications.json")["items"] if x["kind"] == "cpic"]
    assert rx and all(x["label"]["basis_kind"] == "cpic" and x["text_zh"] == x["label"]["sources"][0]["quote"] for x in rx)
    texts = {d["drug"]: d["recommendation_text"] for d in slco["drugs"]}
    assert all(x["text_zh"] == texts[x["drug"]] for x in rx if x["gene"] == "SLCO1B1")
    # personal response prediction for salt reduction: before the experiment and after one
    p3 = next(p for p in _load(v3, "predict", "predictions.json")["items"] if p["plan_item"] == "I3")
    p4 = next(p for p in _load(v4, "predict", "predictions.json")["items"] if p["plan_item"] == "I3")
    assert p3["p_detect_favorable"] == pytest.approx(0.49, abs=0.02) and p3["label"]["basis_kind"] == "statistical_model"
    assert p4["p_detect_favorable"] == pytest.approx(0.76, abs=0.02) and p4["p_true_favorable"] > 0.99
    assert any(s["type"] == "personal_experiment" for s in p4["label"]["sources"])
    # diagnostic conclusion for fasting glucose 5.6: ADA impaired fasting glucose yes, China no; the guideline is named
    cn = json.loads((tmp_path / "2026-09-dx-runs.json").read_text())["cn"]
    assert not cn["met"] and "DX-ifg-cn" in cn["not_met"]
    dx = _load(v4, "clinical", "diagnoses.json")
    assert dx["guideline_set"] == "ada" and [x["rule_id"] for x in dx["items"]] == ["DX-ifg-ada"]
    md = (v4 / "deliver" / "report.md").read_text()
    assert "{{" not in md and "## 复测对比与趋势" in md and "做饭用限盐勺" in md and "约 超过" not in md
    doctor = md.split("## 问医生的单子")[1].split("\n## ")[0]
    assert "尿酸" in doctor and "空腹血糖受损" in doctor and "只按研究结果、不算你自己的前后对比时" in md
    assert "不报告" in md and "LinAge" not in md.split("## 本期没有给出结果的项目")[0]       # mostly-missing LinAge withheld
    ex = json.loads((v4 / "deliver" / "la-export.json").read_text())
    assert ex["diagnoses"][0]["rule"] == "DX-ifg-ada" and ex["diagnoses"][0]["hit"] is True and ex["diagnoses"][0]["source_zh"]
    assert all("drug" in x and "drug_zh" in x and isinstance(x["routine"], bool) for x in ex["medications"])
    assert "rs4149056" not in json.dumps(ex) and json.loads((v4 / "deliver" / "cards.json").read_text())
    assert "中国" in md.split("## 你在同龄人群中的位置")[1][:400]          # LDL against the Chinese lipid reference


def test_case_c_swap_identity_alert(tmp_path, capsys, case_c_net):
    """case-C-swap: the 2025-09 folder brings another person's VCF (a real 1000 Genomes sample)."""
    v1 = _visit(tmp_path, "2023-09", capsys, until="twin")
    v2 = _visit(tmp_path, "2024-09", capsys, prev=v1, history=[_twin(v1)], exposures=CASE_C / "2024-09" / "exposures.json", until="twin")
    v3 = _visit(tmp_path, "2025-09", capsys, raw=CASE_C_SWAP / "2025-09" / "raw", prev=v2, history=[_twin(v1), _twin(v2)],
                exposures=CASE_C / "2025-09" / "exposures.json", until="twin")
    rec = json.loads((v3 / "state.json").read_text())["reuse_germline"]
    assert rec["reused"] is False and "different VCF" in rec["reason"]                 # genomics ran again on the new file
    v4 = _visit(tmp_path, "2026-09", capsys, prev=v2, history=[_twin(v1), _twin(v2), _twin(v3)],
                exposures=CASE_C / "2026-09" / "exposures.json", until="twin")
    ser = run_json(capsys, "twin", "series", "--snapshots", *(str(_twin(w)) for w in (v1, v2, v3, v4)))
    alerts = [a for a in ser["identity_alerts"] if a.get("basis") == "fingerprint"]
    assert [a["snapshot"] for a in alerts] == [str(_twin(v3))] and alerts[0]["mismatches"] >= 2
    assert [s["excluded"] for s in ser["snapshots"]] == [False, False, True, False]
    g = next(m for m in ser["markers"] if m["marker_key"] == "glucose")
    assert [p["value"] for p in g["points"]] == [4.6, 5.0, 5.6]                         # the swapped visit is not a point
    assert any(c["code"] == "snapshot_excluded_identity" for c in g["caveats"])
