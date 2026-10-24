"""twin series beyond the shared trend cases: identity, readout noise, concordance transitions, la-twin/1 input,
per-row sampling dates (WS-A.1), retest judgeability (WS-B.2), germline reuse and the stage clock (WS-I).
All member data here is SYNTHETIC."""
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
from lalib import intake, labnames, series, twin  # noqa: E402

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
pytestmark = pytest.mark.skipif(not (LIB and (Path(LIB) / "catalog.json").exists()), reason="needs longevity-skills")


def snap(t, obs=(), readouts=(), schema="la-twin/2", member="m", **kw):
    return {"schema": schema, "member": {"id": member}, "snapshot": {"t": t}, "longevity_skills_home": LIB,
            "observations": [{"t": t, "marker": m, "value": v, "unit": u, **({"t": t} if schema == "la-twin/2" else {})}
                             for m, v, u in obs],
            "readouts": [{"id": i, "label_zh": l, "value": v, "unit": "岁"} for i, l, v in readouts], **kw}


def run(*snaps):
    return series.run(list(snaps), [f"s{i}" for i in range(len(snaps))])


def test_refuses_different_members():
    with pytest.raises(C.LAError) as e:
        run(snap("2025-01-01", member="a"), snap("2026-01-01", member="b"))
    assert e.value.code == C.EXIT_INPUT


def test_twin_1_input_runs_with_a_warning_and_snapshot_dates():
    a = snap("2025-01-10", [("空腹血糖", "4.6", "mmol/L")], schema="la-twin/1")
    a["observations"][0]["t"] = "1999-01-01"                       # /1: the snapshot date is what counts
    b = snap("2026-01-10", [("空腹血糖", "5.6", "mmol/L")])
    res = run(a, b)
    assert res["warnings"] and "la-twin/1" in res["warnings"][0]
    g = res["markers"][0]
    assert [p["t"] for p in g["points"]] == ["2025-01-10", "2026-01-10"] and g["summary_verdict"] == "change_beyond_noise"


def test_series_orders_snapshots_and_merges_names_of_one_analyte():
    res = run(snap("2026-01-10", [("葡萄糖(GLU)", "5.6", "mmol/L")]), snap("2024-01-10", [("空腹血糖", "4.6", "mmol/L")]),
              snap("2025-01-10", [("血糖", "5.0", "mmol/L")]))
    g = [m for m in res["markers"] if m["marker_key"] == "glucose"]
    assert len(g) == 1 and [p["value"] for p in g[0]["points"]] == [4.6, 5.0, 5.6] and g[0]["marker"] == "葡萄糖(GLU)"
    assert res["observed_at"] == "2026-01-10"


def test_no_model_marker_is_side_by_side_only():
    res = run(snap("2025-01-10", [("握力", "38", "kg")]), snap("2026-01-10", [("握力", "30", "kg")]))
    m = res["markers"][0]
    assert m["summary_verdict"] == "not_judged" and m["marker_key"] is None and len(m["points"]) == 2


def test_fingerprint_identity_alert_and_threshold():
    calls = {f"rs{i}": "AG" for i in range(100, 160)}
    other = dict(calls, rs100="AA", rs101="GG")
    a = snap("2024-01-10", [("空腹血糖", "4.6", "mmol/L")], genotype_fingerprint={"calls": calls})
    b = snap("2025-01-10", [("空腹血糖", "9.9", "mmol/L")], genotype_fingerprint={"calls": other})
    c = snap("2026-01-10", [("空腹血糖", "4.7", "mmol/L")], genotype_fingerprint={"calls": calls})
    res = run(a, b, c)
    assert [x["snapshot"] for x in res["identity_alerts"]] == ["s1"] and res["identity_alerts"][0]["mismatches"] == 2
    assert res["identity_alerts"][0]["wording_zh"] == series.IDENTITY_ZH
    g = res["markers"][0]
    assert [p["value"] for p in g["points"]] == [4.6, 4.7]                    # the swapped visit leaves every comparison
    res = series.run([a, b, c], ["s0", "s1", "s2"], max_mismatch=3)            # configurable: 2 of 60 (3.3 %) still > 1 %
    assert res["identity_alerts"]
    assert series.fingerprint_mismatch({f"rs{i}": "AG" for i in range(300)}, dict({f"rs{i}": "AG" for i in range(300)}, rs1="AA"))["alert"] is False
    assert series.fingerprint_mismatch({"rs1": "AG"}, {"rs1": "GA"})["mismatches"] == 0     # unphased order


def test_apoe_change_is_an_identity_alert():
    a = snap("2025-01-10", readouts=[("native.apoe.genotype", "APOE", "ε3/ε3")])
    b = snap("2026-01-10", readouts=[("native.apoe.genotype", "APOE", "ε4/ε4")])
    res = run(a, b)
    assert res["identity_alerts"] and res["identity_alerts"][0]["basis"] == "apoe"


def test_readouts_not_judged_without_noise_table_and_judged_with_one(tmp_path, monkeypatch):
    a = snap("2025-01-10", readouts=[("epiage.dnam_horvath", "Horvath 时钟", 50.0)])
    b = snap("2026-01-10", readouts=[("epiage.dnam_horvath", "Horvath 时钟", 58.0)])
    monkeypatch.setattr(series, "readout_noise_path", lambda: tmp_path / "missing.json")
    r = run(a, b)["readouts"][0]
    assert r["baseline_vs_latest"]["verdict"] == "not_judged"
    rn = {"schema": "longevity-readout-noise/1", "z": 1.96, "readouts": [
        {"key": "horvath", "match": ["epiage.dnam_horvath", "Horvath"], "unit": "岁", "noise_kind": "technical_replicate",
         "sd_abs": 2.0, "source": {"title": "SYNTHETIC test value, not from a paper"}}]}
    (tmp_path / "rn.json").write_text(json.dumps(rn))
    monkeypatch.setattr(series, "readout_noise_path", lambda: tmp_path / "rn.json")
    r = run(a, b)["readouts"][0]
    b2 = r["baseline_vs_latest"]
    assert b2["verdict"] == "increase_beyond_noise" and b2["band_abs"] == [-5.54, 5.54]   # 1.96 × √2 × 2.0
    assert b2["wording_zh"] == "超出测量误差" and "不等于生物学上的真实变化" in b2["note_zh"]
    small = snap("2026-01-10", readouts=[("epiage.dnam_horvath", "Horvath 时钟", 53.0)])
    assert run(a, small)["readouts"][0]["baseline_vs_latest"]["verdict"] == "within_noise"
    res = twin.retest_check(series.noise(), "甲基化时钟", ["epiage.dnam_horvath"], 52)
    assert res["judgeable"] and res["noise"]["table"] == "readout_noise"


def _conc(cls, band="elevated", retrieved="2026-01-01", tmv="la-traits/1@x", analyte="ldl"):
    return [{"analyte": analyte, "genetic_band": band, "class": cls, "gwas_retrieved_on": retrieved, "trait_map_version": tmv}]


def test_transitions():
    ldl = lambda v: [("低密度脂蛋白胆固醇", v, "mmol/L")]
    # latent -> concordant_adverse
    res = run(snap("2025-01-10", ldl("3.0"), concordance=_conc("latent")), snap("2026-01-10", ldl("4.2"), concordance=_conc("concordant_adverse")))
    t = res["transitions"][0]
    assert t["kind"] == "tendency_expressed" and t["ask_doctor"] and t["marker"] == "低密度脂蛋白胆固醇"
    # concordant_adverse -> latent within noise / beyond noise
    res = run(snap("2025-01-10", ldl("3.5"), concordance=_conc("concordant_adverse")), snap("2026-01-10", ldl("3.3"), concordance=_conc("latent")))
    assert res["transitions"][0]["real_change"] is False and "正常波动之内" in res["transitions"][0]["wording_zh"]
    res = run(snap("2025-01-10", ldl("4.6"), concordance=_conc("concordant_adverse")), snap("2026-01-10", ldl("3.2"), concordance=_conc("latent")))
    assert res["transitions"][0]["real_change"] is True
    # acquired_adverse -> concordant_typical with a concurrent exposure (no check-in counts: the site)
    ex = [{"plan_item": "I3", "label_zh": "每周运动", "category": "exercise", "targets": ["LDL-C"], "status": "adopted",
           "adopted_at": "2025-03-01", "cumulative": None}]
    res = run(snap("2025-01-10", ldl("4.6"), concordance=_conc("acquired_adverse", band="typical")),
              snap("2026-01-10", ldl("3.2"), concordance=_conc("concordant_typical", band="typical"), exposures=ex))
    t = res["transitions"][0]
    assert t["kind"] == "back_in_range" and t["concurrent_exposures"][0]["cumulative_in_window"] is None
    # band changed because the public data were re-retrieved: not a member change
    res = run(snap("2025-01-10", ldl("3.0"), concordance=_conc("latent")),
              snap("2026-01-10", ldl("3.0"), concordance=_conc("concordant_typical", band="typical", retrieved="2026-06-01")))
    assert res["transitions"][0]["kind"] == "public_data_updated" and not res["identity_alerts"]
    # band changed with the same retrieval: the genotype changed, so it is a sample problem
    res = run(snap("2025-01-10", ldl("3.0"), concordance=_conc("latent")),
              snap("2026-01-10", ldl("3.0"), concordance=_conc("concordant_typical", band="typical")))
    assert not res["transitions"] and res["identity_alerts"][0]["basis"] == "concordance"


def test_concurrent_exposure_counts_only_inside_the_window():
    ex = [{"plan_item": "I3", "label_zh": "工作日晚饭用限盐勺", "category": "diet", "targets": ["收缩压"], "status": "adopted",
           "adopted_at": "2026-07-20", "weekly_checkins": [{"week_start": "2026-06-01", "count": 3},
                                                          {"week_start": "2026-07-20", "count": 4}, {"week_start": "2026-07-27", "count": 5}],
           "cumulative": 12},
          {"plan_item": "I4", "label_zh": "散步", "category": "exercise", "targets": ["收缩压"], "status": "dropped"}]
    res = run(snap("2026-07-01", [("收缩压", "150", "mmHg")]), snap("2026-08-31", [("收缩压", "128", "mmHg")], exposures=ex))
    m = res["markers"][0]
    assert m["summary_verdict"] == "change_beyond_noise"
    assert [(c["plan_item"], c["cumulative_in_window"]) for c in m["concurrent_exposures"]] == [("I3", 9)]


def test_compare_is_the_two_point_case_and_converts_units(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps(snap("2025-01-10", [("空腹血糖", "83", "mg/dL")]), ensure_ascii=False))
    b.write_text(json.dumps(snap("2026-01-10", [("空腹血糖", "5.6", "mmol/L")]), ensure_ascii=False))
    row = twin.compare(a, b)["rows"][0]
    assert row["verdict"] == "increase_beyond_noise" and row["rcv_up_pct"] == 14.6
    pr = series.run([json.loads(a.read_text()), json.loads(b.read_text())], ["a", "b"])["markers"][0]["pairs"][0]
    assert (pr["verdict"], pr["change_pct"]) == (row["verdict"], row["change_pct"])


def test_compare_with_history_rows_uses_each_snapshots_latest_date(tmp_path):
    a = snap("2025-01-10", [("空腹血糖", "5.0", "mmol/L")])
    b = snap("2026-01-10", [("空腹血糖", "5.1", "mmol/L")])
    b["observations"].insert(0, {"t": "2024-01-10", "marker": "空腹血糖", "value": "4.0", "unit": "mmol/L"})
    (tmp_path / "a.json").write_text(json.dumps(a, ensure_ascii=False))
    (tmp_path / "b.json").write_text(json.dumps(b, ensure_ascii=False))
    row = twin.compare(tmp_path / "a.json", tmp_path / "b.json")["rows"][0]
    assert row["cur"] == "5.1" and row["verdict"] == "within_noise"


# ------------------------------------------------------------ WS-A.1 sampling dates
@pytest.mark.parametrize("raw,iso", [("2024-09-03", "2024-09-03"), ("2024/9/3", "2024-09-03"), ("2024.09.03 08:10", "2024-09-03"),
                                     ("20240903", "2024-09-03"), ("2024年9月3日", "2024-09-03"), ("9/3/2024", None), ("昨天", None),
                                     ("2024-13-01", None)])
def test_parse_date(raw, iso):
    assert labnames.parse_date(raw) == iso


def _state(raw: Path, sample_date=None):
    return C.new_state({"id": "m", "age": 50, "sex": "male", "sample_date": sample_date, "answers": {}}, raw)


def test_row_dates_file_dates_and_pending(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "mirobody_labs_all.csv").write_text("marker,value,unit,ref_range,date\n空腹血糖,4.6,mmol/L,3.9-6.1,2023-09-01\n"
                                               "空腹血糖,5.0,mmol/L,3.9-6.1,2024-09-03\n空腹血糖,5.6,mmol/L,3.9-6.1,2026-09-01\n", encoding="utf-8")
    (raw / "mixed.csv").write_text("项目,结果,单位,采样日期\n白蛋白,44,g/L,2026-09-01\n肌酐,68,umol/L,\n", encoding="utf-8")
    st = _state(raw)
    intake.run_intake(st, raw)
    rows = {(r["marker"], r["value"]): r for r in st["labs"]}
    assert rows[("空腹血糖", "4.6")]["collected_at"] == "2023-09-01" and rows[("空腹血糖", "4.6")]["collected_at_source"] == "row"
    fid = next(f["id"] for f in st["files"] if f["name"] == "mixed.csv")
    assert any(j["kind"] == "collected_at" and j["target"] == fid for j in C.pending_judgments(st))     # never defaulted
    with pytest.raises(C.LAError):
        intake.set_file_date(st, fid, "2026-09-01", "", "r")                     # the quote is required
    res = intake.set_file_date(st, fid, "2026-09-01", "报告单抬头：采样日期 2026-09-01", "read the header")
    assert res == {"file": fid, "date": "2026-09-01", "rows_with_own_date": 1, "rows_dated_by_file": 1}
    assert not [j for j in C.pending_judgments(st) if j["kind"] == "collected_at"]
    for r in st["labs"]:
        r["confirm"] = {"row_key": labnames.row_key(r, st["labs"]), "answer": "yes", "why": "t"}
    assert labnames.visit_date(st) == "2026-09-01"                              # no sample date: the latest row date
    assert {r["value"] for r in labnames.usable_rows(st) if r["marker"] == "空腹血糖"} == {"5.6"}       # methods: this visit
    assert {r["value"] for r in labnames.usable_rows(st, history=True) if r["marker"] == "空腹血糖"} == {"4.6", "5.0", "5.6"}
    snapd = twin.provisional(st, tmp_path)
    assert sorted(o["t"] for o in snapd["observations"] if o["marker"] == "空腹血糖") == ["2023-09-01", "2024-09-03", "2026-09-01"]
    g = next(m for m in series.run([snapd], ["p"])["markers"] if m["marker_key"] == "glucose")
    assert g["trend"]["verdict"] == "rising_beyond_noise"                       # one workspace, three visits


def test_undated_file_in_a_dated_workspace_without_sample_date_is_asked(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "a.csv").write_text("marker,value,unit,date\n空腹血糖,4.6,mmol/L,2023-09-01\n", encoding="utf-8")
    (raw / "b.csv").write_text("marker,value,unit\n空腹血糖,5.6,mmol/L\n", encoding="utf-8")
    st = _state(raw)
    intake.run_intake(st, raw)
    assert [j["target"] for j in C.pending_judgments(st) if j["kind"] == "collected_at"] == ["F002"]
    st2 = _state(raw, "2026-09-01")
    intake.run_intake(st2, raw)
    assert not [j for j in C.pending_judgments(st2) if j["kind"] == "collected_at"]
    assert next(r for r in st2["labs"] if r["value"] == "5.6")["collected_at_source"] == "workspace"


# ------------------------------------------------------------ WS-B.2 retest judgeability
@pytest.mark.parametrize("what,targets,weeks,judgeable,warn", [
    ("空腹血糖", [], 12, True, None),
    ("同型半胱氨酸", [], 12, False, "no within-person noise model"),
    ("低密度脂蛋白胆固醇", [], 2, True, "sooner"),
    ("{{r:x|label}}", ["低密度脂蛋白胆固醇(LDL-C)"], 12, True, None),
    ("甲基化时钟", ["epiage.dnam_hannum"], 52, False, "no within-person noise model"),
])
def test_retest_check(what, targets, weeks, judgeable, warn, monkeypatch, tmp_path):
    monkeypatch.setattr(series, "readout_noise_path", lambda: tmp_path / "none.json")
    res = twin.retest_check(series.noise(), what, targets, weeks, {"x": "x"})
    assert res["judgeable"] is judgeable
    assert (warn is None and not res["warnings"]) or (warn and any(warn in w for w in res["warnings"]))


# ------------------------------------------------------------ WS-I germline reuse and stage clock
def _prev_ws(tmp_path, raw):
    prev = tmp_path / "prev"
    (prev / "work" / "insights").mkdir(parents=True)
    (prev / "work" / "pgx").mkdir(parents=True)
    vcf_sha = C.sha256_file(raw / "x.vcf") if (raw / "x.vcf").exists() else "deadbeef"
    st = {"schema": C.SCHEMA, "member": {"id": "m"},
          "files": [{"id": "F001", "name": "x.vcf", "kind": "variants_vcf", "sha256": vcf_sha}]}
    (prev / "state.json").write_text(json.dumps(st))
    (prev / "work" / "insights" / "fingerprint.json").write_text(json.dumps({"schema": "la-fingerprint/1", "sites": 1, "calls": {"rs1": "AA"}}))
    (prev / "work" / "insights" / "genotype_phenotype.json").write_text("{}")
    (prev / "work" / "insights" / "concordance.json").write_text(json.dumps({"trait_map_version": intake.trait_map_version(), "items": []}))
    (prev / "work" / "pgx" / "pgx.json").write_text(json.dumps({"schema": "la-pgx/1", "status": "ran"}))
    (prev / "work" / "readouts.json").write_text(json.dumps({"readouts": [{"id": "native.apoe.genotype", "label_zh": "APOE", "value": "ε3/ε3"}]}))
    return prev


def test_reuse_germline(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "x.vcf").write_text("##fileformat=VCFv4.2\n")
    prev = _prev_ws(tmp_path, raw)
    st = _state(raw)
    ws = tmp_path / "new"
    rec = intake.reuse_germline(st, ws, prev, raw)                               # same VCF bytes
    assert rec["reused"] and rec["fingerprint_check"] == "identical_vcf_sha256"
    assert set(rec["files"]) == {"fingerprint.json", "pgx.json", "genotype_phenotype.json", "concordance.json", "apoe.json"}
    assert (ws / "work" / "reuse" / "germline.json").exists()
    raw2 = tmp_path / "raw2"
    raw2.mkdir()
    (raw2 / "labs.csv").write_text("marker,value,unit\n空腹血糖,5.6,mmol/L\n")
    (prev / "work" / "insights" / "concordance.json").write_text(json.dumps({"trait_map_version": "la-traits/1@old", "items": []}))
    rec = intake.reuse_germline(_state(raw2), tmp_path / "new2", prev, raw2)     # labs only; trait map changed since
    assert rec["reused"] and rec["fingerprint_check"] == "no_new_genotype_data" and "genotype_phenotype.json" in rec["skipped"]
    assert "genotype_phenotype.json" not in rec["files"] and "fingerprint.json" in rec["files"]
    raw3 = tmp_path / "raw3"
    raw3.mkdir()
    (raw3 / "y.vcf").write_text("##fileformat=VCFv4.2\n##other\n")
    rec = intake.reuse_germline(_state(raw3), tmp_path / "new3", prev, raw3)     # a different VCF: genomics runs again
    assert not rec["reused"] and "different VCF" in rec["reason"]
    other = _state(raw)
    other["member"]["id"] = "someone-else"
    with pytest.raises(C.LAError):
        intake.reuse_germline(other, tmp_path / "new4", prev, raw)


def test_init_reuse_flag_and_fast_route(tmp_path, capsys):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "labs.csv").write_text("marker,value,unit\n空腹血糖,5.6,mmol/L\n")
    prev = _prev_ws(tmp_path, raw)
    ws = tmp_path / "ws"
    capsys.readouterr()
    assert la.main(["init", str(raw), str(ws), "--member-id", "m", "--sample-date", "2026-09-01", "--longevity-skills", LIB,
                    "--reuse-germline-from", str(prev)]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["reuse_germline"]["reused"] and "fingerprint.json" in out["reuse_germline"]["files"]
    assert la.main(["intake", str(ws)]) == 0
    capsys.readouterr()
    assert la.main(["observe", str(ws)]) == 0
    obs = json.loads(capsys.readouterr().out)
    assert "fast_route" in obs and "series verdicts beyond noise" in obs["fast_route"]
    assert "intake" in obs["timing"]
    st = json.loads((ws / "state.json").read_text())
    assert st["timing"]["intake"]["started_at"]
