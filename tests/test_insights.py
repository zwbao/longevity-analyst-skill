"""Insights layer (offline): effect parsing, reference ranges, VCF genotype reads, percentiles, board validation,
MR projection arithmetic, wearable summaries."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import board, causal, genomics, pubdata, reference, wearable  # noqa: E402
from lalib import common as C  # noqa: E402


def test_beta_direction_parsing():
    assert pubdata._beta("1.7 ln nmol/L increase") == 1.7
    assert pubdata._beta("0.3 unit decrease") == -0.3
    assert pubdata._beta("0.3 unit") is None and pubdata._beta(None) is None


@pytest.mark.parametrize("text,lo,hi", [("3.9-6.1", 3.9, 6.1), ("<5.2", None, 5.2), (">1.0", 1.0, None), ("≤30", None, 30.0),
                                        ("5～15", 5.0, 15.0), ("阴性", None, None)])
def test_ref_bounds(text, lo, hi):
    assert genomics._ref_bounds(text) == (lo, hi)


def _vcf(tmp_path, lines):
    p = tmp_path / "m.vcf"
    hdr = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    p.write_text("\n".join(hdr + lines) + "\n")
    return genomics.Vcf(p, {"min_gq": 20, "min_dp": 10})


def test_vcf_dosage_calls_blocks_and_gaps(tmp_path):
    v = _vcf(tmp_path, ["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30",
                        "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t1/1:40:30",
                        "chr1\t300\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t0/1:5:30",
                        "chr1\t400\t.\tC\t<NON_REF>\t.\tPASS\tEND=500\tGT:GQ:DP\t0/0:40:25",
                        "chr1\t600\t.\tC\t<NON_REF>\t.\tPASS\tEND=700\tGT:GQ:DP\t0/0:10:25",
                        "chr1\t800\t.\tC\tT\t50\tLowQual\t.\tGT:GQ:DP\t0/1:40:30"])
    assert v.alt_dosage("1", 100, "A", "G") == {"dosage": 1, "status": "called"}
    assert v.alt_dosage("chr1", 200, "C", "T")["dosage"] == 2
    assert v.alt_dosage("1", 300, "C", "T")["dosage"] is None                      # GQ 5
    assert v.alt_dosage("1", 450, "C", "T") == {"dosage": 0, "status": "ref_block"}
    assert v.alt_dosage("1", 650, "C", "T")["dosage"] is None                      # low-quality block
    assert v.alt_dosage("1", 900, "C", "T") == {"dosage": None, "status": "not_called"}   # never assumed reference
    assert v.alt_dosage("1", 800, "C", "T")["status"].startswith("filtered")


def test_percentile_interpolation():
    pcts, vals = [1, 50, 99], [1.0, 5.0, 9.0]
    assert reference._pct_of(5.0, pcts, vals)["pct"] == 50
    assert reference._pct_of(3.0, pcts, vals)["pct"] == 25.5
    assert reference._pct_of(0.5, pcts, vals)["bound"] == "below" and reference._pct_of(10, pcts, vals)["bound"] == "above"


def _ws(tmp_path):
    ws = tmp_path / "ws"
    (ws / "work" / "evidence").mkdir(parents=True)
    (ws / "work" / "insights" / "board").mkdir(parents=True)
    C.write_json(ws / "work" / "readouts.json", {"readouts": [{"id": "native.gut.gmhi", "value": -1.8, "unit": ""},
                                                             {"id": "china-par-ascvd-risk.risk_10y_pct", "value": 6.0, "unit": "%"}]})
    st = {"member": {"age": 58, "sex": "male", "answers": {"smoker": "no"}},
          "labs": [{"marker": "低密度脂蛋白胆固醇", "value": "3.8", "unit": "mmol/L", "confirm": {"answer": "yes", "row_key": "", "why": "t"}}],
          "files": [], "stages": {}}
    from lalib import labnames
    for k, r in labnames.keyed(st).items():
        r["confirm"]["row_key"] = k
    return ws, st


def test_board_questions_and_findings(tmp_path):
    ws, st = _ws(tmp_path)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": "问题", "hypothesis_zh": "假设", "basis": ["native.gut.gmhi"], "why_zh": "因为"} for i in (1, 2, 3)]}
    C.write_json(ws / "q.json", qs)
    board.register_questions(st, ws, ws / "q.json")
    bad = {"questions": qs["questions"][:2]}
    C.write_json(ws / "b.json", bad)
    with pytest.raises(C.LAError, match="3-10"):
        board.register_questions(st, ws, ws / "b.json")
    bad2 = {"questions": [dict(q, basis=["made.up"]) for q in qs["questions"]]}
    C.write_json(ws / "b2.json", bad2)
    with pytest.raises(C.LAError, match="not this member"):
        board.register_questions(st, ws, ws / "b2.json")
    pubdata._cache(ws, "mr", "k", "u", [{"ref": "mr:x->y:IVW", "b": 0.4, "se": 0.04, "exposure": "LDL", "outcome": "CHD"}])
    ok = {"id": "Q1", "verdict": "supported", "confidence": "moderate", "member_evidence": ["低密度脂蛋白胆固醇", "member.answers.smoker"],
          "public_evidence": ["mr:x->y:IVW"], "summary_zh": "证据一致。", "next_step_zh": "复查血脂。", "limitations_zh": "单次测量。"}
    C.write_json(ws / "work" / "insights" / "board" / "Q1.json", ok)
    board.register_finding(st, ws, "Q1")
    for bad_f, msg in ((dict(ok, confidence="high"), "confidence"), (dict(ok, public_evidence=["gwas:none"]), "not retrieved"),
                       (dict(ok, public_evidence=[]), "at least one"), (dict(ok, summary_zh="高出 3 倍"), "number rules"),
                       (dict(ok, verdict="maybe"), "verdict"), (dict(ok, member_evidence=["x.y"]), "not this member")):
        C.write_json(ws / "work" / "insights" / "board" / "Q2.json", dict(bad_f, id="Q2"))
        with pytest.raises(C.LAError, match=msg):
            board.register_finding(st, ws, "Q2")
    assert not board.complete(st)
    board.skip(st, "Q2", "no data")
    C.write_json(ws / "work" / "insights" / "board" / "Q3.json", dict(ok, id="Q3", verdict="insufficient", public_evidence=[]))
    board.register_finding(st, ws, "Q3")
    assert board.complete(st)


def test_mr_projection_arithmetic(tmp_path):
    ws, st = _ws(tmp_path)
    pubdata._cache(ws, "mr", "k", "u", [{"ref": "mr:ldl->chd:IVW", "b": math.log(1.5), "se": 0.05, "exposure": "LDL", "outcome": "CHD"}])
    rec = causal.project(st, ws, "mr:ldl->chd:IVW", "ldl", 2.8, "china-par-ascvd-risk.risk_10y_pct")
    sd = C.data("ref_population.json")["labs_nhanes"]["analytes"]["ldl"]["strata"]["male:50-59"]["sd"]
    d = (2.8 - 3.8) / sd
    orr = 1.5 ** d
    exp = 0.06 * orr / (1 - 0.06 + 0.06 * orr)
    assert abs(rec["risk_after"] - round(exp, 4)) < 1e-4 and rec["risk_after"] < 0.06
    with pytest.raises(C.LAError, match="not an MR record"):
        causal.project(st, ws, "mr:made-up", "ldl", 2.8, "china-par-ascvd-risk.risk_10y_pct")


def test_wearable_summary(tmp_path):
    ws, _ = _ws(tmp_path)
    f = tmp_path / "w.csv"
    rows = ["日期,步数,静息心率"] + [f"2026-09-{d:02d},{5000 + d * 10},{70}" for d in range(1, 29)]
    f.write_text("\n".join(rows), encoding="utf-8")
    st = {"files": [{"id": "F009", "path": str(f), "name": "w.csv"}]}
    res = wearable.summarize(st, ws, "F009", {"date": "日期", "steps": "步数", "rhr": "静息心率"})
    ids = {r["id"]: r["value"] for r in res["readouts"]}
    assert ids["wear.rhr.mean30"] == 70 and ids["wear.steps.trend30"] > 0
    with pytest.raises(C.LAError, match="not in"):
        wearable.summarize(st, ws, "F009", {"date": "日期", "steps": "步数X"})
