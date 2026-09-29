"""Organ layer: published calculators, AI-estimate registration rules, report rendering."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import common as C, organs, report  # noqa: E402


def test_egfr_ckd_epi_2021_reference_values():
    # female 60 y, Scr 0.8 mg/dL -> 84 (NKF calculator); male 50 y, Scr 1.0 -> 92
    assert round(organs.egfr_ckd_epi_2021(0.8 * 88.4, 60, "female")) == 84
    assert round(organs.egfr_ckd_epi_2021(1.0 * 88.4, 50, "male")) == 92
    assert organs.kdigo_g(84) == "G2" and organs.kdigo_g(95) == "G1" and organs.kdigo_g(50) == "G3a" and organs.kdigo_g(10) == "G5"


def test_fib4_and_tiers():
    v = organs.fib4(50, 30, 25, 250)
    assert abs(v - 1.2) < 0.01
    assert organs.fib4_tier(1.0, 50).startswith("低于") and organs.fib4_tier(2.0, 50).startswith("界值之间")
    assert organs.fib4_tier(3.0, 50).startswith("高于") and organs.fib4_tier(1.8, 70).startswith("低于")   # age >=65 cut-off 2.0
    assert organs.fib4_tier(1.8, 50).startswith("界值之间") and organs.fib4_tier(1.2951, 50).startswith("低于")
    with pytest.raises(ValueError):
        organs.egfr_ckd_epi_2021(70, 60, "f")
    assert organs.fib4_tier(1.0, 30) is None                                            # not validated under 35


def _ws(tmp_path, labs=None):
    ws = tmp_path / "ws"
    (ws / "work" / "organs" / "estimates").mkdir(parents=True)
    (ws / "work" / "evidence").mkdir(parents=True)
    C.write_json(ws / "work" / "readouts.json", {"readouts": [{"id": "native.organ.egfr", "value": 84.0, "unit": "mL/min/1.73m²",
                                                               "label_zh": "eGFR", "kind": "computed", "method": "native.organ_indices"}]})
    C.write_json(ws / "work" / "evidence" / "pubmed_x.json", {"items": [{"pmid": "111"}]})
    st = {"member": {"age": 60, "sex": "female"}, "labs": labs or [{"marker": "肌酐", "value": 70.7, "unit": "umol/L"}],
          "files": [], "organs": {}}
    _conf(st["labs"])
    organs.bundle(st, ws)
    return ws, st


def _sha(st, organ="kidney"):
    return next(b["bundle_sha256"] for b in st["organs"]["bundles"] if b["organ"] == organ)


def _est(**kw):
    r = {"disease": "慢性肾脏病 3 期及以上", "horizon_years": 10, "low": 0.04, "point": 0.08, "high": 0.16,
         "basis": ["native.organ.egfr", "肌酐"], "evidence": [{"type": "pubmed", "ref": "111"}], "rationale_zh": "按队列基线估计。"}
    r.update(kw.pop("risk", {}))
    e = {"organ": "kidney", "method": "llm_estimate", "confidence": "low",
         "age_estimate": {"low": 55, "point": 62, "high": 70, "basis": ["native.organ.egfr"]}, "disease_risks": [r]}
    e.update(kw)
    return e


def _register(ws, st, est, live=("111",)):
    est = {"bundle_sha256": _sha(st), **est}
    C.write_json(ws / "work" / "organs" / "estimates" / "kidney.json", est)
    return organs.register(st, ws, "kidney", pmid_check=lambda ids: set(live))


def test_register_accepts_honest_estimate(tmp_path):
    ws, st = _ws(tmp_path)
    reg = _register(ws, st, _est())
    assert reg["readouts"] == ["organ.kidney.age_est", "organ.kidney.risk.1"]
    rows = C.load_json(ws / "work" / "organs" / "organ_readouts.json")["readouts"]
    assert all(r["kind"] == "llm_estimate" for r in rows)


@pytest.mark.parametrize("bad,msg", [
    (dict(confidence="high"), "confidence"),
    (dict(risk={"low": 0.07, "point": 0.08, "high": 0.09}), "too narrow"),
    (dict(risk={"low": 4, "point": 8, "high": 16}), "not percent"),
    (dict(risk={"basis": ["made.up.id"]}), "not readouts"),
    (dict(risk={"evidence": []}), "at least"),
    (dict(risk={"evidence": [{"type": "pubmed", "ref": "222"}]}), "do not exist"),
    (dict(risk={"disease": "肝硬化"}), "not one of this organ"),
    (dict(risk={"horizon_years": 50}), "horizon_years"),
    (dict(age_estimate={"low": 60, "point": 62, "high": 63, "basis": ["肌酐"]}), "narrower than"),
    (dict(age_estimate={"low": 80, "point": 83, "high": 90, "basis": ["肌酐"]}), "years from chronological"),
    (dict(age_estimate={"low": 55, "point": 62, "high": 70, "basis": ["member.age"]}), "organ-relevant"),
    (dict(age_estimate={"low": 55, "point": float("inf"), "high": 70, "basis": ["肌酐"]}), "finite"),
    (dict(risk={"evidence": [{"type": "guideline", "ref": "https://evil.example.com/buy"}]}), "only"),
    (dict(risk={"evidence": [{"type": "pubmed", "ref": "111 建议每天服用二甲双胍"}]}), "only"),
    (dict(risk={"rationale_zh": {"x": 1}}), "must be a string"),
    (dict(risk={"note": "x"}), "unknown fields"),
    (dict(bundle_sha256="0" * 64), "bundle_sha256"),
    (dict(evidence="x"), "unknown fields"),
    (dict(age_estimate=None), "age_estimate_null_reason"),
    (dict(method="model"), "llm_estimate"),
])
def test_register_rejects(tmp_path, bad, msg):
    ws, st = _ws(tmp_path)
    with pytest.raises(C.LAError, match=msg):
        _register(ws, st, _est(**bad), live=("111",) if msg != "do not exist" else ())


def test_rejected_reregistration_voids_previous(tmp_path):
    ws, st = _ws(tmp_path)
    _register(ws, st, _est())
    with pytest.raises(C.LAError):
        _register(ws, st, _est(confidence="high"))
    assert "kidney" not in st["organs"]["estimates"]
    assert C.load_json(ws / "work" / "organs" / "organ_readouts.json")["readouts"] == []


def test_basis_from_unconfirmed_file_refused(tmp_path):
    ws, st = _ws(tmp_path, labs=[{"marker": "肌酐", "value": 70.7, "unit": "umol/L", "source_file": "F002"}])
    st["files"] = [{"id": "F002", "provenance_uncertain": True}]
    organs.bundle(st, ws)
    with pytest.raises(C.LAError, match="not confirmed"):
        _register(ws, st, _est(risk={"basis": ["肌酐"]}, age_estimate=None, age_estimate_null_reason="无"))


def test_lab_rows_need_consistent_names():
    def ind(labs, age=67):
        r = organs.compute_indices({"member": {"age": age, "sex": "female"},
                                    "labs": _conf([{"marker": m, "value": v, "unit": u} for m, v, u in labs])})
        return {x["id"].split(".")[-1]: x["value"] for x in r["readouts"]}, r["notes"]
    base = [("肌酐(Cr)", "68", "umol/L"), ("空腹血糖(GLU)", "5.4", "mmol/L"), ("甘油三酯(TG)", "1.3", "mmol/L")]
    v, _ = ind(base)
    assert v["egfr"] == 84.6 and v["tyg"] == 8.63
    for bad in [("肌酐(尿)", "8840", "umol/L"), ("铬(Cr)", "0.05", "umol/L"), ("血糖(餐后2小时)", "11.8", "mmol/L")]:
        rows = _conf([{"marker": bad[0], "value": bad[1], "unit": bad[2]}], yes=False) + _conf([{"marker": m, "value": x, "unit": u} for m, x, u in base])
        r = organs.compute_indices({"member": {"age": 67, "sex": "female"}, "labs": rows})
        assert {x["id"].split(".")[-1]: x["value"] for x in r["readouts"]} == v
    v3, notes = ind([("肌酐(Cr)", "150", "umol/L")] + base)
    assert "egfr" not in v3 and any("多行" in n for n in notes)
    for val in ("NaN", "inf", "0.77"):
        assert "egfr" not in ind([("肌酐(Cr)", val, "umol/L")] + base[1:])[0]
    assert "egfr" not in ind(base, age=12)[0]


def test_pmid_must_be_retrieved_this_run(tmp_path):
    ws, st = _ws(tmp_path)
    with pytest.raises(C.LAError, match="not retrieved"):
        _register(ws, st, _est(risk={"evidence": [{"type": "pubmed", "ref": "999"}]}), live=("999",))


def test_skip_after_register_drops_rows(tmp_path):
    ws, st = _ws(tmp_path)
    _register(ws, st, _est())
    organs.skip(st, ws, "kidney", "member declined")
    assert C.load_json(ws / "work" / "organs" / "organ_readouts.json")["readouts"] == []
    assert "kidney" not in st["organs"]["estimates"]
    _register(ws, st, _est())
    assert "kidney" not in st["organs"]["skipped"]


def test_ai_estimate_never_renders_as_measurement(tmp_path):
    ws, st = _ws(tmp_path)
    _register(ws, st, _est())
    ro = report._readouts(ws)
    assert report.substitute("{{r:organ.kidney.age_est}}", ro) == "62 岁（AI 估计，区间 55–70 岁）"
    assert "AI 估计" in report.substitute("{{r:organ.kidney.risk.1}}", ro)
    md = "\n".join(report.organ_section(ws, ro))
    assert "## 器官体检表" in md and "不是测量" in md and "62 岁（55–70，置信度低）" in md and "8%（4%–16%，10 年）" in md


def test_organ_rationale_digits_are_refused_at_register(tmp_path):
    ws, st = _ws(tmp_path)
    with pytest.raises(C.LAError, match="number rules"):
        _register(ws, st, _est(risk={"rationale_zh": "风险约为同龄人的 2 倍"}))


def _conf(labs, yes=True, st=None):
    """Mark rows answered; keys are computed over the list the rows will live in (st['labs'] or labs itself)."""
    from lalib import labnames
    whole = st["labs"] if st else labs
    keys = {id(r): k for k, r in labnames.keyed({"labs": whole}).items()}
    for r in labs:
        r["confirm"] = {"row_key": keys[id(r)], "answer": "yes" if yes else "no", "why": "test"}
    return labs


def test_unconfirmed_candidate_rows_never_feed_indices():
    from lalib import labnames
    labs = [{"marker": m, "value": v, "unit": u} for m, v, u in [("肌酐(酶法)", "68", "umol/L"), ("葡萄糖(GLU-2h)", "11.8", "mmol/L"),
                                                               ("甘油三酯(TG)", "1.3", "mmol/L")]]
    st = {"member": {"age": 67, "sex": "female"}, "labs": labs}
    assert {p["marker"] for p in labnames.pending(st)} == {"肌酐(酶法)", "葡萄糖(GLU-2h)", "甘油三酯(TG)"}
    assert organs.compute_indices(st)["readouts"] == []
    _conf(labs)                                           # agent: all three are serum / fasting values
    labs[1]["confirm"]["answer"] = "no"                   # ... except the 2-hour glucose
    ids = {r["id"].split(".")[-1]: r["value"] for r in organs.compute_indices(st)["readouts"]}
    assert ids.get("egfr") == 84.6 and "tyg" not in ids
    labs[0]["value"] = "150"                              # editing a row voids its confirmation
    assert "egfr" not in {r["id"].split(".")[-1] for r in organs.compute_indices(st)["readouts"]}
    assert labnames.pending(st)[0]["marker"] == "肌酐(酶法)"


def test_bundle_is_idempotent_and_keeps_registrations(tmp_path):
    ws, st = _ws(tmp_path)
    _register(ws, st, _est())
    before = [b["bundle_sha256"] for b in st["organs"]["bundles"]]
    organs.bundle(st, ws)
    assert [b["bundle_sha256"] for b in st["organs"]["bundles"]] == before
    assert "kidney" in st["organs"]["estimates"]
    st["labs"][0]["value"] = 250                                     # data changed -> bundle changes -> estimate voided
    organs.bundle(st, ws)
    assert "kidney" not in st["organs"]["estimates"]


def test_calibrated_override_blocks_contradicting_estimate(tmp_path):
    ws, st = _ws(tmp_path)
    ro = C.load_json(ws / "work" / "readouts.json")
    ro["readouts"][0]["value"] = 17.7                                 # eGFR G4
    C.write_json(ws / "work" / "readouts.json", ro)
    organs.bundle(st, ws)
    with pytest.raises(C.LAError, match="not estimated"):
        _register(ws, st, _est(risk={"low": 0.01, "point": 0.03, "high": 0.06}))
    md = "\n".join(report.organ_section(ws, report._readouts(ws)))
    assert "达到 3 期的 GFR 阈值" in md
    est = _est(disease_risks=[])                                       # kidney's only disease is overridden: age alone is fine
    _register(ws, st, est)


def test_overrides_for_diabetes_and_fibrosis(tmp_path):
    ws, st = _ws(tmp_path, labs=[{"marker": "空腹血糖", "value": "8.4", "unit": "mmol/L"},
                                  {"marker": "糖化血红蛋白", "value": "6.0", "unit": "%"}])
    assert [o["disease"] for o in organs.active_overrides(ws, "metabolic", st)] == ["2 型糖尿病"]
    st["labs"][0]["value"] = "5.2"
    _conf(st["labs"])
    assert organs.active_overrides(ws, "metabolic", st) == []
    st["labs"][1]["value"] = "7.4"
    _conf(st["labs"])
    assert [o["lab"] for o in organs.active_overrides(ws, "metabolic", st)] == ["hba1c"]
    ro = C.load_json(ws / "work" / "readouts.json")
    ro["readouts"].append({"id": "native.organ.fib4", "value": 9.42, "unit": "", "label_zh": "FIB-4", "kind": "computed"})
    C.write_json(ws / "work" / "readouts.json", ro)
    assert organs.active_overrides(ws, "liver", st)[0]["disease"] == "进展期肝纤维化或肝硬化"


def test_verify8_no_unconfirmed_row_reaches_anything(tmp_path):
    from lalib import labnames, methods
    labs = [{"marker": m, "value": v, "unit": u} for m, v, u in
            [("葡萄糖 餐后2h", "11.8", "mmol/L"), ("餐后2小时血糖 GLU", "11.8", "mmol/L"), ("肌酐　尿", "18", "mg/dL"),
             ("白蛋白(ALB)", "44.1", "g/L"), ("白蛋白(ALB)", "44.1", "g/L")]]
    st = {"member": {"age": 67, "sex": "female"}, "labs": labs, "files": []}
    pend = labnames.pending(st)
    assert len(pend) == 5 and len({p["row_key"] for p in pend}) == 5                # spaced names and duplicates all asked
    assert methods._lab_rows(st) == []
    keys = labnames.keyed(st)
    for k, r in keys.items():
        r["confirm"] = {"row_key": k, "answer": "yes" if r["marker"].startswith("白蛋白") else "no", "why": "t"}
    assert [r["marker"] for r in methods._lab_rows(st)] == ["白蛋白(ALB)", "白蛋白(ALB)"]
    assert labnames.pending(st) == []
    assert "肌酐　尿" in organs._unusable_basis(st, [])                          # a refused row cannot be a basis


def test_verify9_spaced_abbrev_names_are_candidates_and_overrides_robust(tmp_path):
    from lalib import labnames
    assert labnames.analytes({"marker": "肌酐 Cr"}) and labnames.analytes({"marker": "空腹血糖 GLU"})
    assert labnames.candidates_zh(labnames.analytes({"marker": "葡萄糖(GLU)"})) == ["空腹血糖（静脉血浆/血清；不是餐后、OGTT、随机或尿糖）"]
    ws, st = _ws(tmp_path, labs=[{"marker": "空腹血糖 GLU", "value": "7.2", "unit": "mmol/L"},
                                  {"marker": "空腹血糖 GLU", "value": "8.1", "unit": "mmol/L"}])
    assert organs.active_overrides(ws, "metabolic", st)                   # two differing rows, both past the threshold
    ws2, st2 = _ws(tmp_path / "b", labs=[{"marker": "空腹血糖", "value": "126", "unit": "mg/dL"}])
    assert organs.active_overrides(ws2, "metabolic", st2)                 # 126 mg/dL is 7.0 mmol/L
    ws3, st3 = _ws(tmp_path / "c", labs=[{"marker": "糖化血红蛋白", "value": "57", "unit": "mmol/mol"}])
    assert organs.active_overrides(ws3, "metabolic", st3)                 # 57 mmol/mol = 7.4 %
    ws4, st4 = _ws(tmp_path / "d", labs=[{"marker": "HbA1c", "value": "7.4", "unit": ""}])
    assert organs.active_overrides(ws4, "metabolic", st4)
    ro = C.load_json(ws / "work" / "readouts.json")
    ro["readouts"].append({"id": "native.organ.fib4", "value": 2.67, "value_raw": 2.674, "label_zh": "FIB-4", "kind": "computed"})
    C.write_json(ws / "work" / "readouts.json", ro)
    assert organs.active_overrides(ws, "liver", st)                       # the unrounded value decides


def test_verify10_flags_and_blank_units_still_trigger_overrides(tmp_path):
    from lalib.labnames import clean_value
    assert [clean_value(v) for v in ("8.4↑", "8.4 H", "250*", "6.8（偏高）")] == ["8.4", "8.4", "250", "6.8"]
    ws, st = _ws(tmp_path, labs=[{"marker": "空腹血糖", "value": "8.4↑", "unit": ""}])
    assert organs.active_overrides(ws, "metabolic", st)
    ws2, st2 = _ws(tmp_path / "b", labs=[{"marker": "肌酐", "value": "250↑", "unit": ""},
                                           {"marker": "肌酐", "value": "240", "unit": "umol/L"}])
    C.write_json(ws2 / "work" / "readouts.json", {"readouts": []})       # duplicates: no eGFR readout was computed
    assert organs.active_overrides(ws2, "kidney", st2)                   # still overridden from the rows themselves


def test_verify11_unreadable_threshold_row_blocks_registration(tmp_path):
    from lalib.labnames import clean_value
    assert [clean_value(v) for v in ("8.4H", "8.4 HH", "8.4(↑)", "8.4 偏高", "8.4▲", "7.2%↑")] == ["8.4"] * 5 + ["7.2"]
    ws, st = _ws(tmp_path, labs=[{"marker": "肌酐", "value": "70.7", "unit": "umol/L"},
                                  {"marker": "空腹血糖", "value": "8.4", "unit": "mM"}])
    assert organs.unreadable_check_rows(st, "metabolic") == ["空腹血糖"]


def test_verify12_blank_unit_us_value_is_unreadable(tmp_path):
    ws, st = _ws(tmp_path, labs=[{"marker": "肌酐", "value": "70.7", "unit": "umol/L"},
                                  {"marker": "空腹血糖", "value": "151", "unit": ""}])
    assert organs.unreadable_check_rows(st, "metabolic") == ["空腹血糖"]


def test_twin_keeps_provenance_flag_and_compare_skips_it(tmp_path):
    from lalib import twin
    prev = {"observations": [{"t": "2025-09-10", "marker": "肌酐", "value": "68", "unit": "umol/L", "provenance_uncertain": True}],
            "readouts": [], "identity": {}}
    cur = {"observations": [{"t": "2026-09-10", "marker": "肌酐", "value": "90", "unit": "umol/L"}], "readouts": [], "identity": {}}
    C.write_json(tmp_path / "p.json", prev)
    C.write_json(tmp_path / "c.json", cur)
    res = twin.compare(tmp_path / "p.json", tmp_path / "c.json")
    rows = res.get("rows") or res.get("changes") or []
    assert any(r.get("verdict") == "not_judged" and "not confirmed" in r.get("why", "") for r in rows)
