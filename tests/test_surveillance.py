"""WS-E surveillance rules (la-surveillance/1) and plan (la-surveillance-plan/1), offline."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import genomics, surveillance  # noqa: E402
from lalib import common as C  # noqa: E402

LPA_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC9639807/"


def _ws(tmp_path, monkeypatch, labs=("lpa",), retests=(), fetched=(LPA_URL,), variants=True):
    ws = tmp_path / "ws"
    for u in fetched:                                   # what `la.py evidence fetch` records
        C.write_json(ws / "work" / "evidence" / f"url_{hashlib.sha1(u.encode()).hexdigest()[:10]}.json",
                     {"url": u, "status": 200, "title": "t", "retrieved_at": "2026-10-24"})
    items = [{"id": f"I{i + 1}", "category": "diet", "action_zh": "x", "targets": [], "evidence": [], "executor": "member",
              "retest": {"what": w, "after_weeks": wk}} for i, (w, wk) in enumerate(retests)]
    C.write_json(ws / "work" / "intervene" / "plan.json", {"items": items})
    monkeypatch.setattr(genomics, "member_labs", lambda st: {k: {"flag": "in_range"} for k in labs})
    st = {"member": {"age": 49, "sex": "male", "answers": {}}, "labs": [], "files": [],
          "processed": {"variants": [{"path": "x.vcf"}] if variants else []}}
    return ws, st


def _by(doc):
    return {x["id"]: x for x in doc["items"]}


def test_rules_file_is_valid_and_quoted():
    rs = surveillance.rules()
    assert surveillance.check_rules(rs) == []
    lpa = next(r for r in rs if r["id"] == "S-lpa-once")
    assert lpa["action"]["interval"] == "once" and lpa["source"]["url"] == LPA_URL and lpa["source"]["doi"] == "10.1093/eurheartj/ehac361"
    assert lpa["source"]["quote"].startswith("Lp(a) should be measured at least once in adults")
    assert not any("apoe" in str(r).lower() for r in rs)          # no guideline sets a screening interval by APOE genotype


def test_once_measured_lpa_drops_out_of_retests(tmp_path, monkeypatch):
    ws, st = _ws(tmp_path, monkeypatch, retests=[("空腹血糖", 12)])
    doc = surveillance.plan(st, ws)
    lpa = _by(doc)["S-lpa-once"]
    assert lpa["once"] and lpa["done"] and lpa["interval_weeks"] is None and lpa["reason"] == "guideline"
    assert lpa["source"]["quote"] and lpa["executor"] == "physician"
    assert "S-lpa-once" not in {x["id"] for x in surveillance.retests(doc)}
    assert {x["id"] for x in surveillance.retests(doc)} == {"R-I1"}
    assert _by(doc)["S-germline"]["done"] and _by(doc)["S-germline"]["once"]
    assert C.load_json(ws / "work" / "surveillance" / "plan.json")["schema"] == "la-surveillance-plan/1"


def test_lpa_measured_in_an_earlier_snapshot_counts(tmp_path, monkeypatch):
    ws, st = _ws(tmp_path, monkeypatch, labs=())
    assert not _by(surveillance.plan(st, ws))["S-lpa-once"]["done"]       # never measured: measure once
    tw = tmp_path / "twin_2023.json"
    C.write_json(tw, {"observations": [{"kind": "lab", "marker": "脂蛋白(a)", "value": "35", "unit": "mg/dL"}]})
    doc = surveillance.plan(st, ws, [tw])
    assert _by(doc)["S-lpa-once"]["done"] and "S-lpa-once" not in {x["id"] for x in surveillance.retests(doc)}


def test_no_rules_keeps_only_intervention_retests(tmp_path, monkeypatch):
    monkeypatch.setattr(surveillance, "rules", lambda: [])
    ws, st = _ws(tmp_path, monkeypatch, retests=[("空腹血糖", 12), ("低密度脂蛋白胆固醇", 24)], variants=False)
    doc = surveillance.plan(st, ws)
    assert [(x["id"], x["interval_weeks"], x["reason"]) for x in doc["items"]] == [("R-I1", 12, "intervention"), ("R-I2", 24, "intervention")]


def test_rule_needs_its_source_fetched(tmp_path, monkeypatch):
    ws, st = _ws(tmp_path, monkeypatch, fetched=())
    with pytest.raises(C.LAError, match="evidence fetch"):
        surveillance.plan(st, ws)
    assert surveillance.unfetched(ws, surveillance.rules()) == ["S-lpa-once"]


def test_merge_keeps_shorter_interval_and_both_reasons(tmp_path, monkeypatch):
    rule = {"id": "S-x", "when": {}, "action": {"analyte": "ldl", "markers": ["低密度脂蛋白胆固醇"], "interval": 52},
            "executor": "physician", "reason": "guideline", "source": {"title": "t", "quote": "q", "url": LPA_URL}}
    monkeypatch.setattr(surveillance, "rules", lambda: [rule])
    ws, st = _ws(tmp_path, monkeypatch, retests=[("复查 LDL-C", 12)])
    x = _by(surveillance.plan(st, ws))["S-x"]
    assert x["interval_weeks"] == 12 and x["reasons"] == ["guideline", "intervention"] and x["interventions"] == ["I1"]
    rule["action"]["interval"] = 8
    x = _by(surveillance.plan(st, ws))["S-x"]
    assert x["interval_weeks"] == 8 and x["reason"] == "guideline" and set(x["reasons"]) == {"guideline", "intervention"}
    assert "R-I1" not in _by(surveillance.plan(st, ws))                    # deduplicated into the rule item


def test_conditions(tmp_path, monkeypatch):
    rule = {"id": "S-c", "when": {"concordance_class": {"analyte": "urate", "class": "concordant_adverse"}},
            "action": {"analyte": "urate", "markers": ["尿酸"], "interval": 26}, "executor": "physician", "reason": "genetic",
            "source": {"title": "t", "quote": "q", "url": LPA_URL}}
    monkeypatch.setattr(surveillance, "rules", lambda: [rule])
    ws, st = _ws(tmp_path, monkeypatch)
    assert "S-c" not in _by(surveillance.plan(st, ws))
    C.write_json(ws / "work" / "insights" / "concordance.json", {"items": [{"analyte": "urate", "class": "concordant_adverse"}]})
    assert _by(surveillance.plan(st, ws))["S-c"]["interval_weeks"] == 26
    bad = dict(rule, when={"apoe": "e4"})
    assert surveillance.check_rules([bad]) and surveillance.check_rules([dict(rule, source={"title": "t"})])


def test_report_section_quotes_the_guideline(tmp_path, monkeypatch):
    from lalib import report_genomics
    ws, st = _ws(tmp_path, monkeypatch, retests=[("空腹血糖", 12)])
    surveillance.plan(st, ws)
    st["surveillance"] = {"sha256": C.sha256_file(ws / "work" / "surveillance" / "plan.json", limit=None)}
    secs = dict(report_genomics.sections(ws, {"st": st}))
    md = secs["监测计划"]
    assert "测一次就够，不需要复测" in md and "12 周后" in md
    assert "Lp(a) should be measured at least once in adults" in md and "10.1093/eurheartj/ehac361" in md and LPA_URL in md
