"""Harness tests. Every gate carries a negative test; the adversarial-review repros are kept as regressions."""
from __future__ import annotations

import csv
import gzip
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "skills" / "longevity-analyst" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import la  # noqa: E402
from lalib import common as C  # noqa: E402
from lalib import evidence, intake, methods, native, pipelines, preflight, report, twin  # noqa: E402

FIX = Path(__file__).parent / "fixtures"
LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
HAS_LIB = bool(LIB) and (Path(LIB) / "catalog.json").exists()
needs_lib = pytest.mark.skipif(not HAS_LIB, reason="LONGEVITY_SKILLS_HOME not set")


def _state(tmp: Path, **member):
    m = {"id": "t1", "age": 60, "sex": "female", "mode": "commercial", "answers": {}}
    m.update(member)
    st = C.new_state(m, tmp)
    st["workspace"] = str(tmp)
    return st


def _vcf(path: Path, rows, sample="S1", contig_len="248956422", extra_header=""):
    """rows: (chrom, pos, ref, alt, gt[, info[, fmt_extra dict[, filter]]])"""
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID=chr1,length={contig_len}>\n{extra_header}")
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for row in rows:
            chrom, pos, ref, alt, gt = row[:5]
            info = row[5] if len(row) > 5 else "."
            extra = row[6] if len(row) > 6 else {}
            flt = row[7] if len(row) > 7 else "PASS"
            fmt = ":".join(["GT", *extra.keys()])
            val = ":".join([gt, *[str(v) for v in extra.values()]])
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{flt}\t{info}\t{fmt}\t{val}\n")


def _ingest(tmp, files, **member):
    raw = tmp / "raw"
    raw.mkdir(exist_ok=True)
    for name, content in files.items():
        p = raw / name
        if isinstance(content, bytes):
            p.write_bytes(content)
        else:
            p.write_text(content, encoding="utf-8")
    st = _state(tmp, **member)
    intake.run_intake(st, raw)
    return st, {f["name"]: f for f in st["files"]}


BETAS = "cpg,beta\n" + "\n".join(f"cg{i:08d},0.{i % 9 + 1}" for i in range(60))


# ================================================================ intake
def test_intake_classifies_and_routes(tmp_path):
    st, by = _ingest(tmp_path, {
        "betas.csv": BETAS,
        "mvals.csv": "cpg,m\n" + "\n".join(f"cg{i:08d},{i - 25}.5" for i in range(50)),
        "gut.txt": "#mpa_v20\n#clade_name\trelative_abundance\ns__Akkermansia_muciniphila\t30.2\ns__Bacteroides_vulgatus\t69.8\n",
        "prot.tsv": "Protein\tS1\nP02768\t25.1\nP02741\t20.3\nP01023\t19\n",
        "labs.csv": "项目,结果,单位\n白蛋白,44,g/L\n肌酐,70,umol/L\n",
        "run_R1.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n"),
        "run_R2.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n"),
        "report.pdf": b"%PDF-1.4\n"})
    _vcf(tmp_path / "raw" / "v.vcf.gz", [("chr19", 44908684, "T", "C", "0/1")])
    intake.run_intake(st, tmp_path / "raw")
    by = {f["name"]: f for f in st["files"]}
    assert by["betas.csv"]["status"] == "needs_judgment"
    assert by["mvals.csv"]["status"] == "rejected"
    assert by["gut.txt"]["status"] == "ready" and st["processed"]["gut_profile"][0]["provenance"]["mpa_version"] == "mpa_v20"
    assert by["prot.tsv"]["status"] == "needs_judgment"
    assert len(st["labs"]) == 2
    assert by["v.vcf.gz"]["assembly"] == "GRCh38" and by["v.vcf.gz"]["status"] == "ready"
    kinds = {j["kind"] for j in C.pending_judgments(st)}
    assert kinds == {"tissue", "platform", "value_scale", "modality", "transcribe", "identity"}
    intake.assign(st, by["run_R1.fastq.gz"]["id"], "modality", "wgs", "lab note says WGS")
    assert by["run_R2.fastq.gz"]["status"] == "needs_pipeline"


def test_assign_requires_reason_and_valid_option(tmp_path):
    st, by = _ingest(tmp_path, {"prot.tsv": "Protein\tS1\nP02768\t25.1\nP02741\t20.3\n"})
    fid = by["prot.tsv"]["id"]
    with pytest.raises(C.LAError):
        intake.assign(st, fid, "platform", "olink", "  ")
    with pytest.raises(C.LAError):
        intake.assign(st, fid, "platform", "luminex", "a reason")


def test_resolving_a_judgment_that_does_not_exist_fails(tmp_path):
    with pytest.raises(C.LAError):
        C.resolve_judgment(_state(tmp_path), "identity", "member", "consistent", "looked")


def test_multisample_beta_asks_column_and_ignores_pvalue_columns(tmp_path):
    rows = "\n".join(f"cg{i:08d},0.0{i % 5},0.{i % 9 + 1},0.8" for i in range(60))
    st, by = _ingest(tmp_path, {"b.csv": "ID_REF,Detection Pval,MEMBER,OTHER\n" + rows})
    f = by["b.csv"]
    assert f["sample_columns"] == ["MEMBER", "OTHER"] and f["pvalue_columns"] == ["Detection Pval"]
    assert any(j["kind"] == "sample_column" for j in C.pending_judgments(st))
    intake.assign(st, f["id"], "sample_column", "MEMBER", "lab sample sheet lists MEMBER as this person")
    intake.assign(st, f["id"], "tissue", "whole_blood", "blood draw")
    derived = Path(st["processed"]["methylation"][0]["path"])
    first = gzip.open(derived, "rt").read().splitlines()[1]
    assert first == "cg00000000,0.1"                     # MEMBER column, not the p-value or OTHER


def test_merged_metaphlan_asks_column_and_uses_it(tmp_path):
    st, by = _ingest(tmp_path, {"m.tsv": "#mpa_v20\nclade_name\tA\tB\ns__X\t60\t10\ns__Y\t40\t90\n"})
    intake.assign(st, by["m.tsv"]["id"], "sample_column", "B", "member is sample B")
    sp = native.read_metaphlan(Path(st["processed"]["gut_profile"][0]["path"]))
    assert sp == {"s__X": 10.0, "s__Y": 90.0}


def test_gbk_lab_table_with_title_row_and_synonym_header(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "化验.csv").write_bytes("某医院检验报告\n检测项目,检测结果,单位\n白蛋白,44,g/L\n".encode("gb18030"))
    st = _state(tmp_path)
    intake.run_intake(st, raw)
    assert st["labs"] and st["labs"][0]["marker"] == "白蛋白"


def test_hidden_symlink_and_corrupt_files_do_not_crash(tmp_path):
    raw = tmp_path / "raw"
    (raw / ".hidden").mkdir(parents=True)
    (raw / ".hidden" / "x.csv").write_text(BETAS)
    (raw / "real.csv").write_text(BETAS)
    os.symlink(raw / "real.csv", raw / "link.csv")
    (raw / "bad.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00garbage")
    st = _state(tmp_path)
    intake.run_intake(st, raw)
    names = {f["name"]: f for f in st["files"]}
    assert "x.csv" not in names and len([n for n in names if n in ("real.csv", "link.csv")]) == 1
    assert names["bad.fastq.gz"]["status"] == "unreadable"


def test_vcf_without_assembly_asks(tmp_path):
    st, by = _ingest(tmp_path, {})
    _vcf(tmp_path / "raw" / "v.vcf.gz", [("19", 45411941, "T", "C", "0/1")], contig_len="59128983")
    intake.run_intake(st, tmp_path / "raw")
    assert any(j["kind"] == "assembly" for j in C.pending_judgments(st))


def test_unknown_modality_parks_file_and_pair(tmp_path):
    st, by = _ingest(tmp_path, {"x_R1.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n"),
                                "x_R2.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n")})
    intake.assign(st, by["x_R1.fastq.gz"]["id"], "modality", "unknown", "user and lab note silent")
    assert {f["status"] for f in st["files"]} == {"deferred_ask_lab"}
    assert preflight.plan_needed(st) == []


def test_identity_exclusion_removes_file_and_derivatives(tmp_path):
    st, by = _ingest(tmp_path, {"labs.csv": "项目,结果,单位\n白蛋白,44,g/L\n"})
    _vcf(tmp_path / "raw" / "v.vcf.gz", [("chr19", 44908684, "T", "C", "0/1")], sample="HG00096")
    intake.run_intake(st, tmp_path / "raw")
    vid = next(f["id"] for f in st["files"] if f["kind"] == "variants_vcf")
    st["processed"]["variants"].append({"file_id": None, "source_file_ids": [vid], "path": "/x", "provenance": {}})
    intake.exclude_files(st, [vid], "VCF sample is 1000 Genomes HG00096")
    assert not st["processed"].get("variants")
    with pytest.raises(C.LAError):
        intake.add_labs(st, tmp_path / "raw" / "labs.csv", vid, "try to re-import")


def test_lab_fix_and_drop_keep_history(tmp_path):
    st = _state(tmp_path)
    st["labs"] = [{"marker": "葡萄糖", "value": "6.8↑", "unit": "mmol/L"}, {"marker": "X", "value": "1", "unit": ""}]
    intake.fix_lab(st, "葡萄糖", "fix", "6.8", None, "arrow is the out-of-range flag, not part of the value")
    assert st["labs"][0]["value"] == "6.8" and st["labs"][0]["history"][0]["value"] == "6.8↑"
    intake.fix_lab(st, "X", "drop", None, None, "not a lab")
    assert [l["marker"] for l in st["labs"]] == ["葡萄糖"]


# ================================================================ native
def test_gmhi_reproduces_published_output(tmp_path):
    rows = list(csv.reader(open(FIX / "gmhi_example_input.csv", encoding="utf-8-sig")))
    gold = {r[0]: float(r[1]) for r in list(csv.reader(open(FIX / "GMHI_output.csv")))[1:]}
    worst = 0.0
    for j, s in enumerate(rows[0][1:], 1):
        p = tmp_path / f"{s}.tsv"
        p.write_text("".join(f"{r[0]}\t{r[j]}\n" for r in rows[1:]))
        worst = max(worst, abs(native.gmhi(p)["readouts"][0]["value"] - gold[s]))
    assert len(gold) == 25 and worst < 1e-3


@pytest.mark.parametrize("gt429,gt7412,expect,e4", [
    ("0/0", "0/0", "ε3/ε3", 0), ("0/1", "0/0", "ε3/ε4", 1), ("1/1", "0/0", "ε4/ε4", 2),
    ("0/0", "0/1", "ε2/ε3", 0), ("0/0", "1/1", "ε2/ε2", 0), ("0/1", "0/1", "ε2/ε4", None),
])
def test_apoe_table(tmp_path, gt429, gt7412, expect, e4):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "T", "C", gt429), ("chr19", 44908822, "C", "T", gt7412)])
    ro = {x["id"]: x["value"] for x in native.apoe(v, "GRCh38")["readouts"]}
    assert ro["native.apoe.genotype"].startswith(expect) and ro["native.apoe.e4_count"] == e4


def test_apoe_split_multiallelic_records_keep_the_carrier_call(tmp_path):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "T", "C", "0/1"), ("chr19", 44908684, "T", "G", "0/0"),
             ("chr19", 44908822, "C", "T", "0/0")])
    ro = {x["id"]: x["value"] for x in native.apoe(v, "GRCh38")["readouts"]}
    assert ro["native.apoe.genotype"] == "ε3/ε4"


def test_apoe_gvcf_confident_reference_block_is_a_call(tmp_path):
    v = tmp_path / "a.g.vcf.gz"
    _vcf(v, [("chr19", 44908600, "A", "<NON_REF>", "0/0", "END=44908900", {"GQ": 50, "MIN_DP": 30})])
    ro = {x["id"]: x["value"] for x in native.apoe(v, "GRCh38")["readouts"]}
    assert ro["native.apoe.genotype"] == "ε3/ε3"


@pytest.mark.parametrize("extra,flt", [({"GQ": 0, "MIN_DP": 0}, "PASS"), ({}, "PASS"), ({"GQ": 1, "MIN_DP": 3}, "RefCall")])
def test_apoe_low_quality_or_uncovered_reference_is_not_a_call(tmp_path, extra, flt):
    v = tmp_path / "a.g.vcf.gz"
    _vcf(v, [("chr19", 44908600, "A", "<NON_REF>", "0/0", "END=44908900", extra, flt)])
    assert native.apoe(v, "GRCh38")["status"] == "not_determined"


def test_apoe_low_quality_variant_call_is_not_a_call(tmp_path):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "T", "C", "1/1", ".", {"GQ": 5, "DP": 2}, "LowQual"), ("chr19", 44908822, "C", "T", "0/0")])
    assert native.apoe(v, "GRCh38")["status"] == "not_determined"


def test_detection_pvalue_matrix_without_keyword_is_refused(tmp_path):
    rows = "\n".join(f"cg{i:08d},{1e-10 * (i + 1):.3g}" for i in range(80))
    st, by = _ingest(tmp_path, {"p.csv": "cpg,S1\n" + rows})
    intake.assign(st, by["p.csv"]["id"], "tissue", "whole_blood", "stated")
    assert not st["processed"].get("methylation") and by["p.csv"]["status"] == "rejected"
    assert "detection p-values" in by["p.csv"]["reason"]


def test_duplicate_probe_ids_are_refused(tmp_path):
    st, by = _ingest(tmp_path, {"d.csv": BETAS + "\ncg00000001,0.9"})
    intake.assign(st, by["d.csv"]["id"], "tissue", "whole_blood", "stated")
    assert not st["processed"].get("methylation") and "appears twice" in by["d.csv"]["reason"]


def test_apoe_missing_site_in_variant_only_vcf_is_not_called(tmp_path):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "T", "C", "0/1")])
    assert native.apoe(v, "GRCh38")["status"] == "not_determined"


def test_apoe_wrong_assembly_refused(tmp_path):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "G", "A", "0/1"), ("chr19", 44908822, "C", "T", "0/0")])
    with pytest.raises(C.LAError):
        native.apoe(v, "GRCh38")


# ================================================================ gates
def _prot_state(tmp, platform, scale):
    base = {"log2_intensity": 20, "npx": 3, "cohort_z": 0, "raw_intensity": 1e6}[scale]
    st, by = _ingest(tmp, {"prot.tsv": "# header\nProtein\tS1\n" + "\n".join(f"P0{2700 + i}\t{base + (i % 7) * 0.1}" for i in range(40))})
    fid = by["prot.tsv"]["id"]
    intake.assign(st, fid, "platform", platform, "stated by lab")
    intake.assign(st, fid, "value_scale", scale, "stated by lab")
    return st


@needs_lib
def test_ms_proteomics_blocked_from_affinity_clocks(tmp_path):
    by = {i["method"]: i for i in methods.plan(_prot_state(tmp_path, "ms_dia", "log2_intensity"), tmp_path)["items"]}
    assert by["plasma-proteomics-brain-immune"]["status"] == "blocked_platform"
    assert by["proteomic-aging-clock"]["status"] == "blocked_platform"


@needs_lib
def test_olink_raw_npx_still_blocked_where_cohort_z_is_required(tmp_path):
    by = {i["method"]: i for i in methods.plan(_prot_state(tmp_path, "olink", "npx"), tmp_path)["items"]}
    assert by["plasma-proteomics-brain-immune"]["status"] == "blocked_platform"      # needs cohort z-scores
    assert by["proteomic-aging-clock"]["status"] != "blocked_platform"               # declared for NPX
    assert by["plasma-proteomic-cellular-aging"]["status"] == "blocked_platform"     # SomaScan-only


REAL_BETAS = ROOT / "examples" / "case-A" / "raw" / "methylation_450k_betas.csv.gz"


def _beta_state(tmp, mode="commercial", tissue="whole_blood", transform=None):
    raw = tmp / "raw"
    raw.mkdir(exist_ok=True)
    lines = gzip.open(REAL_BETAS, "rt").read().splitlines()
    if transform:
        lines = [lines[0]] + [f"{l.split(',')[0]},{transform(float(l.split(',')[1]))}" for l in lines[1:] if l.split(",")[1]]
    (raw / "b.csv").write_text("\n".join(lines))
    st = _state(tmp, mode=mode)
    intake.run_intake(st, raw)
    intake.assign(st, st["files"][0]["id"], "tissue", tissue, "stated")
    return st


@needs_lib
@pytest.mark.parametrize("transform", [lambda v: round(1 - v, 4), lambda v: 0.5])
def test_non_blood_beta_values_are_refused_by_reference_profile(tmp_path, transform):
    by = {i["method"]: i for i in methods.plan(_beta_state(tmp_path, transform=transform), tmp_path)["items"]}
    assert by["(all methylation methods)"]["status"] == "blocked_platform" and by["epiage"]["status"] == "blocked_platform"


@needs_lib
@pytest.mark.parametrize("mode,has_grim,has_dunedin", [("commercial", False, False), ("research", True, True)])
def test_license_mode_controls_clocks(tmp_path, mode, has_grim, has_dunedin):
    ep = next(i for i in methods.plan(_beta_state(tmp_path, mode), tmp_path)["items"] if i["method"] == "epiage")
    assert ("grimagev2" in ep["run"]["clocks"]) is has_grim
    assert ("dunedinpace" in ep["run"]["clocks"]) is has_dunedin


@needs_lib
def test_pyaging_clock_menu_is_never_auto_run(tmp_path):
    st = _beta_state(tmp_path)
    st["member"]["answers"]["matrix"] = "/x"
    by = {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}
    assert by["pyaging"]["status"] == "blocked_license"


@needs_lib
def test_blood_only_clocks_refuse_saliva(tmp_path):
    ep = next(i for i in methods.plan(_beta_state(tmp_path, tissue="saliva"), tmp_path)["items"] if i["method"] == "epiage")
    assert ep["status"] == "blocked_platform"


@needs_lib
def test_gmhi_blocked_on_non_metaphlan2_profiles(tmp_path):
    st, by = _ingest(tmp_path, {"g.txt": "#mpa_vJan21_CHOCOPhlAnSGB_202103\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
                                         "k__Bacteria|s__Akkermansia_muciniphila\t2|74217\t100\t\n"})
    by = {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}
    assert by["native.gmhi"]["status"] == "blocked_platform" and by["native.gut_diversity"]["status"] == "ready"


@needs_lib
def test_labs_map_to_protein_or_cpg_key_refused():
    with pytest.raises(C.LAError):
        la._check_map_target("cg00921350")
    la._check_map_target("albumin_gL")


PHENO_ROWS = [("白蛋白 ALB", "39.2", "g/L"), ("肌酐 CREA", "1.12", "mg/dL"), ("葡萄糖 GLU", "6.3", "mmol/L"),
              ("C反应蛋白 CRP", "4.6", "mg/L"), ("淋巴细胞% LYM%", "24.0", "%"), ("平均红细胞体积 MCV", "95.4", "fL"),
              ("红细胞分布宽度 RDW-CV", "14.2", "%"), ("碱性磷酸酶 ALP", "96", "U/L"), ("白细胞 WBC", "6.9", "10^9/L")]


@needs_lib
def test_pdf_style_names_reach_phenoage(tmp_path):
    st = _state(tmp_path, age=73, sex="male")
    st["labs"] = [{"marker": m, "value": v, "unit": u} for m, v, u in PHENO_ROWS]
    assert methods.plan(st, tmp_path)["items"] and \
        {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}["accelerated-biological-aging-risk"]["status"] == "needs_answers"
    _confirm_all(st)                                                   # nothing reaches a method before the agent's yes
    by = {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}
    assert by["accelerated-biological-aging-risk"]["status"] == "ready"


@needs_lib
def test_unmatched_lab_names_surface_instead_of_vanishing(tmp_path):
    st = _state(tmp_path, age=73, sex="male")
    st["labs"] = [{"marker": f"项目{i}", "value": "1", "unit": ""} for i in range(9)]
    _confirm_all(st)
    by = {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}
    assert by["accelerated-biological-aging-risk"]["status"] == "input_problem"


@needs_lib
def test_optional_answers_are_passed(tmp_path):
    st = _state(tmp_path, age=73, sex="male")
    st["labs"] = [{"marker": m, "value": v, "unit": u} for m, v, u in
                  [("收缩压", "136", "mmHg"), ("总胆固醇", "4.9", "mmol/L"), ("高密度脂蛋白胆固醇", "1.05", "mmol/L"), ("腰围", "94", "cm")]]
    _confirm_all(st)
    st["member"]["answers"] = {"treated": "no", "smoker": "no", "diabetes": "no", "north": "no", "urban": "yes"}
    cp = next(i for i in methods.plan(st, tmp_path)["items"] if i["method"] == "china-par-ascvd-risk")
    assert cp["status"] == "ready" and cp["run"]["args"].get("--urban") == "yes"


@needs_lib
def test_methods_run_survives_a_crashing_method(tmp_path):
    st = _state(tmp_path)
    st["methods"] = {"mode": "commercial", "items": [
        {"method": "native.gut_diversity", "status": "ready", "run": {"kind": "native", "fn": "gut_diversity", "path": "/nope"}},
        {"method": "native.vcf_summary", "status": "ready", "run": {"kind": "native", "fn": "vcf_summary", "path": "/nope"}}]}
    saves = []
    res = methods.run(st, tmp_path, checkpoint=lambda: saves.append(1))
    assert [r["status"] for r in res["ran"]] == ["failed", "failed"] and len(saves) == 2
    assert (tmp_path / "work" / "readouts.json").exists()


# ================================================================ pipelines
def _planned(st, tier, status="planned"):
    st["pipelines"] = [{"run_id": "x", "pipeline": "sarek", "revision": "3.10.0", "tier": tier, "optional": False,
                        "file_ids": [], "cmd": ["nextflow", "run", "nf-core/sarek", "-r", "3.10.0"], "dir": "/tmp/none",
                        "blocked": None, "status": status}]


@pytest.mark.parametrize("tier,consent,accept", [("green", "", False), ("green", "ok", False), ("red", "yes run it", True),
                                                 ("yellow", "yes please", False)])
def test_launch_gates(tmp_path, tier, consent, accept):
    st = _state(tmp_path)
    _planned(st, tier)
    with pytest.raises(C.LAError) as e:
        pipelines.launch(st, "x", consent, force_yellow=accept)
    assert e.value.code == C.EXIT_BLOCKED


def test_skipped_run_needs_reopen_and_stub_needs_dev(tmp_path, monkeypatch):
    st = _state(tmp_path)
    _planned(st, "green", status="skipped")
    with pytest.raises(C.LAError, match="reopen"):
        pipelines.launch(st, "x", "用户说可以跑了")
    _planned(st, "red")
    monkeypatch.delenv("LA_DEV", raising=False)
    with pytest.raises(C.LAError, match="LA_DEV"):
        pipelines.launch(st, "x", "ok, developer test", stub=True)


def test_pipeline_skip_needs_reason(tmp_path):
    st = _state(tmp_path)
    _planned(st, "red")
    with pytest.raises(C.LAError):
        pipelines.skip(st, "x", " ")
    assert pipelines.skip(st, "x", "用户说：先不跑")["status"] == "skipped"


def test_preflight_detects_java_stub(monkeypatch):
    monkeypatch.setattr(preflight, "_run", lambda cmd, **k: {"rc": 1, "out": "Unable to locate a Java Runtime."})
    r = preflight.java_status()
    assert r["ok"] is False and "fix" in r


def test_preflight_detects_docker_without_daemon(monkeypatch):
    monkeypatch.setattr(preflight.shutil, "which", lambda x: "/usr/local/bin/docker")
    monkeypatch.setattr(preflight, "_run", lambda cmd, **k: {"rc": 1, "out": "Cannot connect to the Docker daemon"})
    assert preflight.docker_status()["ok"] is False


def test_deferred_fastq_gets_what_if_estimates(tmp_path):
    st, by = _ingest(tmp_path, {"u_R1.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n")})
    intake.assign(st, by["u_R1.fastq.gz"]["id"], "modality", "unknown", "lab silent")
    pf = preflight.run_preflight(st, tmp_path / "w", check_network=False)
    assert {o["pipeline"] for o in pf["deferred_files"][0]["what_if"]} >= {"sarek", "taxprofiler"}


# ================================================================ trace / review / render
@pytest.mark.parametrize("text,ok", [
    # verify5 (R5-1/R5-2/R5-11): numerals containing an ordinary word, measure words, split runs, disguised links
    ("你的表型年龄是七十一岁", False), ("收缩压一百四十一", False), ("超重十一公斤", False), ("风险高出一万倍", False),
    ("十分之三", False), ("高一倍", False), ("三个月后复查", False), ("每周两次", False), ("七\u200b十\u200b岁", False),
    ("陆 拾 柒", False), ("eGFR 为 八 四", False), ("见 evil.example.com/x", False), ("http\u200b://evil.example", False),
    ("hxxps://evil", False), ("每日 叁 粒", False),
    ("千万别熬夜", True), ("一次性说清", True), ("每天一次散步", True), ("三三两两", True), ("万一", True),
    ("这一期间应保持规律作息", True), ("逐一克服困难", True), ("这是唯一支持该结论的研究", True), ("五十肩", True), ("吃饭七八分饱", True),
    ("你的表型年龄是 seventy-one 岁", False), ("比实际年龄小半年", False), ("第九十分位", False), ("活到一百岁", False),
    ("一三七", False), ("相差一年", False), ("褪黑素{{n:3}}{{r:x|毫克}}", False), ("维生素D {{n:2000单位}}", False),
    ("二甲双胍{{n:1 g}}", False), ("每天{{n:2克}}鱼油", False), ("褪黑素{{n:500 毫 克}}", False),
    ("hs-CRP 高于{{n:3 mg/L}}提示炎症", True), ("尿白蛋白肌酐比{{n:30 mg/g}}", True), ("百岁老人", True), ("三七粉", True),
    ("体重减少{{n:5}}千克", True), ("比实际年龄少了一年", False), ("年轻了十来岁", False), ("六旬出头", False),
    ("每天{{n:3克}}三七粉", False), ("富含维生素C的水果每天{{n:200克}}左右", True), ("素食者每天豆类{{n:50克}}", True),
    ("淋巴细胞百分比偏低", True), ("阿司匹林{{n:0.1g}}", False), ("司美格鲁肽每周{{n:一支}}", False),
    ("坚果（约{{n:10 克}}），富含维生素E", True), ("发酵食品{{n:100–200 克}}，其中含益生菌", True),
    ("维生素C {{n:1 g}}", False), ("肌酸 {{n:5 g}}", False), ("阿卡波糖 {{n:0.05 g}}", False),
    ("每天补充{{n:1 克}}维生素C", False), ("每天冲服{{n:1 袋}}益生菌", False), ("每天吃{{n:200克}}水果，富含维生素C", True),
    ("补充{{n:300毫升}}水分", True),
    ("第二步是复查", True), ("三步走", False),
    ("建议每天食盐不超过{{n:5克}}", True), ("每天喝{{n:一杯}}牛奶", True), ("每天喝一杯牛奶", False), ("维生素D {{n:5000IU}}", False),
    ("你的表型年龄是 {{r:x.phenoage}}。", True),
    ("{{n:IL-6}}、{{n:GLP-1}}、{{n:GRCh38}} 与 {{n:每周 3~5 次}} 都写成字面量", True),
    ("参见 {{pmid:29897866}}", True),
    ("这一点十分重要，一起努力", True),
    ("你的表型年龄是 58 岁。", False),
    ("IL-6 是炎症因子", False),                      # non-member digits must be marked {{n:…}}
    ("你的表型年龄是五十八岁", False),
    ("大九岁", False),
    ("年龄&#x36;&#x37;岁", False),
    ("见 https://example.com", False),
    ("{{r:x.r|风险 58 分}}", False),
    ("PMID:67 岁", False),
])
def test_trace_strict_placeholders(text, ok):
    assert (not report.trace_text(text, [])) == ok, report.trace_text(text, [])


def test_substitute_placeholders():
    ro = {"x.age": {"label_zh": "表型年龄", "value": 57.684, "unit": "a"}, "x.r": {"label_zh": "风险", "value": 5.1, "unit": "%"}}
    assert report.substitute("{{r:x.age|label}} {{r:x.age}}", ro) == "表型年龄 57.68 岁"
    assert report.substitute("{{r:x.r|心血管风险}} = {{r:x.r}}", ro) == "心血管风险 = 5.1 %"
    assert report.substitute("{{n:IL-6}} {{pmid:1}}", ro) == "IL-6 [PMID 1](https://pubmed.ncbi.nlm.nih.gov/1/)"


def _rev(st, verdict="pass", findings=None):
    return json.dumps({"trace_id": st["review"]["trace"]["trace_id"], "verdict": verdict, "findings": findings or []}, ensure_ascii=False)


def _review_ready(tmp):
    st = _state(tmp)
    C.write_json(tmp / "work" / "readouts.json", {"readouts": [{"id": "x.a", "label_zh": "A", "value": 1.0, "unit": ""}]})
    st["methods"] = {"items": [], "readouts_sha256": C.sha256_file(tmp / "work" / "readouts.json", limit=None)}
    (tmp / "work" / "report").mkdir(parents=True, exist_ok=True)
    (tmp / "work" / "report" / "summary.md").write_text("结论见 {{r:x.a}}。", encoding="utf-8")
    assert report.trace(st, tmp)["ok"]
    return st


def test_edit_after_trace_blocks_review(tmp_path):
    st = _review_ready(tmp_path)
    (tmp_path / "work" / "report" / "summary.md").write_text("你是 88 岁。", encoding="utf-8")
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    with pytest.raises(C.LAError, match="changed"):
        report.record_review(st, tmp_path, "pass", f)


def test_edit_after_review_blocks_render(tmp_path):
    st = _review_ready(tmp_path)
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    report.record_review(st, tmp_path, "pass", f)
    (tmp_path / "work" / "report" / "summary.md").write_text("你是 88 岁。", encoding="utf-8")
    with pytest.raises(C.LAError, match="changed"):
        report.render(st, tmp_path)


def test_review_verdict_must_match_file_and_p0_must_block(tmp_path):
    st = _review_ready(tmp_path)
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st, "block"))
    with pytest.raises(C.LAError, match="does not match"):
        report.record_review(st, tmp_path, "pass", f)
    for fd in ({"severity": "P0"}, {"severity": "p0 "}, {"severity": "P０"}, {"level": "P0"}, {"severity": "", "level": "P0"},
               {"Severity": "P0-阻断"}, {"severity": ["P0"]}, {"severity": "critical"}, {"severity": "P1"}):
        f.write_text(_rev(st, "pass", [fd]))
        with pytest.raises(C.LAError, match="P0/P1|severity"):
            report.record_review(st, tmp_path, "pass", f)
    f.write_text(_rev(st, "pass", [{"severity": "P2", "note": "wording"}]))
    report.record_review(st, tmp_path, "pass", f)                    # a P2 alone may pass
    f.write_text("not json")
    with pytest.raises(C.LAError):
        report.record_review(st, tmp_path, "pass", f)
    outside = tmp_path / "rev.json"
    outside.write_text(_rev(st))
    f.write_text(json.dumps({"trace_id": "old", "verdict": "pass", "findings": []}))
    with pytest.raises(C.LAError, match="trace_id"):
        report.record_review(st, tmp_path, "pass", f)
    with pytest.raises(C.LAError, match="work/review"):
        report.record_review(st, tmp_path, "pass", outside)


def test_same_verdict_text_is_fine_for_a_new_trace_but_old_trace_id_is_not(tmp_path):
    st = _review_ready(tmp_path)
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    report.record_review(st, tmp_path, "pass", f)
    (tmp_path / "work" / "report" / "summary.md").write_text("改写后的结论见 {{r:x.a}}。", encoding="utf-8")
    old = f.read_text()
    assert report.trace(st, tmp_path)["ok"]
    f.write_text(old)                                   # the earlier review of different content
    with pytest.raises(C.LAError, match="trace_id"):
        report.record_review(st, tmp_path, "pass", f)
    f.write_text(_rev(st))                              # a fresh review of this trace, same verdict text
    report.record_review(st, tmp_path, "pass", f)


def test_render_escapes_html_in_file_names(tmp_path):
    st = _review_ready(tmp_path)
    st["files"] = [{"id": "F001", "name": "<img src=x onerror=alert(1)>.csv", "kind": "unknown", "status": "unsupported", "sha256": "x"}]
    assert report.trace_text("![](evil) [点这里](javascript:alert('x'))", [])   # links never pass the trace
    assert report.trace(st, tmp_path)["ok"]
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    report.record_review(st, tmp_path, "pass", f)
    report.render(st, tmp_path)
    html_out = (tmp_path / "deliver" / "report.html").read_text()
    assert "<img" not in html_out and "javascript:" not in html_out


# ================================================================ plan checks
@needs_lib
def test_plan_check_rejects_doses_bad_executor_and_forged_pmid(tmp_path, monkeypatch):
    st = _state(tmp_path)
    st["labs"] = [{"marker": "甘油三酯", "value": "2.1", "unit": "mmol/L"}]
    C.write_json(tmp_path / "work" / "readouts.json", {"readouts": [{"id": "native.gut.gmhi", "value": -1}]})
    C.write_json(tmp_path / "work" / "evidence" / "pubmed_fake.json", {"items": [{"pmid": "99999999"}]})
    monkeypatch.setattr(evidence, "verify_pmids", lambda ids: {})          # PubMed says: no such article
    plan = {"items": [
        {"id": "I1", "category": "supplement", "action_zh": "每日两粒鱼油", "targets": ["甘油三酯"],
         "evidence": [{"type": "pubmed", "ref": "99999999"}], "executor": "member", "retest": {"what": "TG", "after_weeks": 12}},
        {"id": "I2", "category": "diet", "action_zh": "饮食", "rationale_zh": "每天 1000mg", "targets": ["甘油三酯"],
         "evidence": [{"type": "guideline", "ref": "https://example.com/g"}], "executor": "nutritionist", "retest": {"what": "TG", "after_weeks": 12}},
    ]}
    p = tmp_path / "plan.json"
    p.write_text(json.dumps(plan, ensure_ascii=False))
    with pytest.raises(C.LAError) as e:
        evidence.check_plan(st, tmp_path, p)
    msg = str(e.value)
    assert msg.count("dose") >= 2 and "executor" in msg and "not retrieved and verified" in msg


@needs_lib
def test_plan_check_accepts_valid_and_test_items_without_retest(tmp_path):
    st = _state(tmp_path)
    st["files"] = [{"id": "F005", "name": "x"}]
    C.write_json(tmp_path / "work" / "readouts.json", {"readouts": [{"id": "epiage.dnam_hannum", "value": 60}]})
    plan = {"items": [
        {"id": "I1", "category": "diet", "action_zh": "在营养师指导下尝试适度热量限制", "targets": ["epiage.dnam_hannum"],
         "evidence": [{"type": "effects", "ref": "calerie-cr-dunedinpace"}], "executor": "nutritionist",
         "retest": {"what": "甲基化时钟", "after_weeks": 52}},
        {"id": "I2", "category": "test", "action_zh": "请机构说明 F005 是什么检测", "targets": ["file:F005"], "evidence": [],
         "executor": "physician"}]}
    p = tmp_path / "plan.json"
    p.write_text(json.dumps(plan, ensure_ascii=False))
    assert evidence.check_plan(st, tmp_path, p)["items"] == 2


# ================================================================ twin
def test_rcv_symmetric_and_lognormal():
    r = twin.rcv({"cvi_pct": 2.5, "cva_pct": None, "log_normal": False}, 1.96)
    assert abs(r["up_pct"] - 7.75) < 0.01 and abs(r["down_pct"] + 7.75) < 0.01
    ln = twin.rcv({"cvi_pct": 29.4, "cva_pct": None, "log_normal": True}, 1.96)
    assert ln["up_pct"] > -ln["down_pct"]


def _twins(tmp, prev_obs, cur_obs, prev_ro=None, cur_ro=None, t0="2026-01-01", t1="2026-07-01"):
    base = {"schema": "la-twin/1", "member": {"id": "m"}, "longevity_skills_home": LIB, "interventions": [], "retest_plan": []}
    a, b = tmp / "a.json", tmp / "b.json"
    a.write_text(json.dumps({**base, "observations": [{"t": t0, **o} for o in prev_obs], "readouts": prev_ro or []}, ensure_ascii=False))
    b.write_text(json.dumps({**base, "observations": [{"t": t1, **o} for o in cur_obs], "readouts": cur_ro or []}, ensure_ascii=False))
    return twin.compare(a, b)


@needs_lib
def test_twin_compare_judges_labs_not_readouts_and_alerts_on_genotype_change(tmp_path):
    res = _twins(tmp_path,
                 [{"marker": "白蛋白", "value": "40", "unit": "g/L"}, {"marker": "C反应蛋白", "value": "<0.5", "unit": "mg/L"}],
                 [{"marker": "白蛋白", "value": "45", "unit": "g/L"}, {"marker": "C反应蛋白", "value": "3", "unit": "mg/L"}],
                 [{"id": "native.apoe.genotype", "label_zh": "APOE", "value": "ε3/ε3"}, {"id": "epiage.dnam_hannum", "label_zh": "H", "value": 60}],
                 [{"id": "native.apoe.genotype", "label_zh": "APOE", "value": "ε4/ε4"}, {"id": "epiage.dnam_hannum", "label_zh": "H", "value": 58}])
    rows = {r.get("marker") or r.get("readout"): r for r in res["rows"]}
    assert rows["白蛋白"]["verdict"] == "increase_beyond_noise"
    assert rows["C反应蛋白"]["verdict"] == "not_judged"                 # censored baseline
    assert rows["epiage.dnam_hannum"]["verdict"] == "not_judged"
    assert rows["native.apoe.genotype"]["verdict"] == "identity_alert" and res["alerts"]


def test_twin_compare_refuses_other_member(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps({"member": {"id": "m1"}, "observations": [], "readouts": []}))
    b.write_text(json.dumps({"member": {"id": "m2"}, "observations": [], "readouts": []}))
    with pytest.raises(C.LAError):
        twin.compare(a, b)


def _confirm_all(st):
    from lalib import labnames
    for k, r in labnames.keyed(st).items():
        r["confirm"] = {"row_key": k, "answer": "yes", "why": "test fixture"}


# ================================================================ CLI: order, invalidation, init, lock
@needs_lib
def test_cli_stage_order_invalidation_and_force(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "labs.csv").write_text("项目,结果,单位\n白蛋白,44,g/L\n", encoding="utf-8")
    ws = tmp_path / "ws"
    assert la.main(["init", str(raw), str(ws), "--member-id", "a", "--age", "60", "--sex", "f", "--longevity-skills", LIB]) == 0
    assert la.main(["init", str(raw), str(ws), "--member-id", "b", "--age", "60", "--sex", "f", "--longevity-skills", LIB]) == C.EXIT_USAGE
    assert la.main(["intake", str(ws)]) == 0
    assert la.main(["methods", "plan", str(ws)]) == C.EXIT_BLOCKED            # identity + preflight not done
    assert la.main(["assign", str(ws), "--what", "identity", "--value", "consistent", "--reason", "single lab sheet, no ids"]) == 0
    assert la.main(["preflight", str(ws), "--no-network"]) == C.EXIT_BLOCKED     # lab rows not confirmed yet
    from lalib import labnames
    stj = json.loads((ws / "state.json").read_text())
    (ws / "a.json").write_text(json.dumps({"answers": [{"row_key": x["row_key"], "answer": "yes", "why": "serum"} for x in labnames.pending(stj)]}))
    assert la.main(["labs", "confirm", str(ws), "--answers", str(ws / "a.json"), "--reason", "checked"]) == 0
    assert la.main(["preflight", str(ws), "--no-network"]) == 0
    assert la.main(["methods", "plan", str(ws)]) == 0
    assert la.main(["report", str(ws)]) == C.EXIT_BLOCKED
    assert la.main(["member", str(ws), "smoker=no", "--source", "user"]) == 0
    st = json.loads((ws / "state.json").read_text())
    assert st["stages"]["methods"] == "todo" and st["methods"] is None     # an answer voids methods and everything after
    assert la.main(["member", str(ws), "id=someone_else", "--source", "user"]) == C.EXIT_USAGE
    assert la.main(["member", str(ws), "clocks=grimagev2", "--source", "user"]) == C.EXIT_USAGE
    assert la.main(["member", str(ws), "age=300", "--source", "user"]) == C.EXIT_USAGE
    assert la.main(["init", str(raw), str(ws), "--member-id", "b", "--age", "61", "--sex", "m", "--force", "--longevity-skills", LIB]) == 0
    assert any(p.name.startswith("ws.old-") for p in tmp_path.iterdir())
    assert json.loads((ws / "state.json").read_text())["member"]["id"] == "b" and not (ws / "work").exists()


def test_workspace_lock_is_exclusive(tmp_path):
    a = C.Workspace(tmp_path)
    a.lock()
    with pytest.raises(C.LAError):
        C.Workspace(tmp_path).lock(wait=0.3)
    a.unlock()
    C.Workspace(tmp_path).lock(wait=0.3)


def test_epiage_coverage_rows_are_parsed(tmp_path):
    p = tmp_path / "report.md"
    p.write_text("| Hannum 甲基化年龄（`hannum`） | 74.12 岁 | +3.12 岁 | 97%（缺 2/71） | 2 | 0 |\n"
                 "| Horvath（`horvath`） | 66.99 岁 | -4.01 岁 | 80%（缺 70/353） | 70 | 3 |\n", encoding="utf-8")
    cov = methods._epiage_coverage(p)
    assert cov["dnam_hannum"] == {"coverage_pct": 97.0, "missing_cpgs": 2, "model_cpgs": 71}
    assert cov["dnam_horvath"]["coverage_pct"] == 80.0


def test_vcf_summary_skips_gvcf_reference_blocks(tmp_path):
    v = tmp_path / "g.vcf.gz"
    _vcf(v, [("chr1", 100, "A", "<NON_REF>", "0/0", "END=200"), ("chr1", 300, "A", "G,<NON_REF>", "0/1")])
    ro = {r["id"]: r["value"] for r in native.vcf_summary(v)["readouts"]}
    assert ro["native.vcf.records"] == 1


# ================================================================ re-attack regressions (round 3)
def test_excluded_inputs_void_their_pipeline_runs(tmp_path):
    st, by = _ingest(tmp_path, {"a_R1.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n"),
                                "a_R2.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n")})
    f1 = by["a_R1.fastq.gz"]["id"]
    intake.assign(st, f1, "modality", "wgs", "lab note")
    st["pipelines"] = [{"run_id": "sarek-x", "pipeline": "sarek", "revision": "3.10.0", "tier": "green", "optional": False,
                        "file_ids": sorted(f["id"] for f in st["files"]), "source_file_ids": sorted(f["id"] for f in st["files"]),
                        "cmd": ["nextflow"], "dir": str(tmp_path), "blocked": None, "status": "planned"}]
    intake.exclude_files(st, [f1], "someone else's reads")
    assert st["pipelines"][0]["status"] == "void_excluded"
    with pytest.raises(C.LAError):
        pipelines.launch(st, "sarek-x", "用户同意运行")


def test_verify_refuses_never_launched_runs(tmp_path):
    st = _state(tmp_path)
    _planned(st, "green")
    st["pipelines"][0]["dir"] = str(tmp_path)
    with pytest.raises(C.LAError, match="never launched"):
        pipelines.verify(st, "x")


@needs_lib
def test_init_force_refuses_while_a_pipeline_process_is_alive(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "labs.csv").write_text("项目,结果,单位\n白蛋白,44,g/L\n", encoding="utf-8")
    ws = tmp_path / "ws"
    assert la.main(["init", str(raw), str(ws), "--member-id", "a", "--age", "60", "--sex", "f", "--longevity-skills", LIB]) == 0
    st = json.loads((ws / "state.json").read_text())
    st["pipelines"] = [{"run_id": "r", "pid": os.getpid(), "status": "running"}]
    (ws / "state.json").write_text(json.dumps(st))
    assert la.main(["init", str(raw), str(ws), "--member-id", "b", "--force", "--longevity-skills", LIB]) == C.EXIT_BLOCKED


def test_state_text_change_after_review_blocks_render(tmp_path):
    st = _review_ready(tmp_path)
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    report.record_review(st, tmp_path, "pass", f)
    st["pipelines"] = [{"run_id": "p", "pipeline": "sarek", "revision": "3", "status": "skipped", "tier": "red",
                        "skip_reason": "用户 LDL 6.2，已用阿托伐他汀"}]
    with pytest.raises(C.LAError, match="changed"):
        report.render(st, tmp_path)


def test_free_text_reasons_are_never_printed(tmp_path):
    st = _review_ready(tmp_path)
    st["pipelines"] = [{"run_id": "p", "pipeline": "sarek", "revision": "3", "status": "skipped", "tier": "red",
                        "skip_reason": "用户 LDL 6.2，吃阿托伐他汀 20mg"}]
    assert report.trace(st, tmp_path)["ok"]
    f = tmp_path / "work" / "review" / "rev.json"
    f.write_text(_rev(st))
    report.record_review(st, tmp_path, "pass", f)
    report.render(st, tmp_path)
    assert "阿托伐他汀" not in (tmp_path / "deliver" / "report.md").read_text()


def test_readouts_edited_before_trace_are_refused(tmp_path):
    st = _review_ready(tmp_path)
    C.write_json(tmp_path / "work" / "readouts.json", {"readouts": [{"id": "x.a", "label_zh": "A", "value": 99.0, "unit": ""}]})
    with pytest.raises(C.LAError, match="does not match"):
        report.trace(st, tmp_path)


def test_invalidation_moves_a_stale_deliverable_aside(tmp_path):
    st = _state(tmp_path)
    (tmp_path / "deliver").mkdir()
    (tmp_path / "deliver" / "report.html").write_text("old")
    st["report"] = {"md": "x"}
    st["stages"]["report"] = "done"
    C.invalidate_after(st, "intake", "new files")
    C.Workspace(tmp_path).save(st)
    assert not (tmp_path / "deliver").exists()


def test_ch_probes_are_kept_in_derived_beta(tmp_path):
    st, by = _ingest(tmp_path, {"b.csv": BETAS + "\nch.1.1234567R,0.4\nch.2.7654321F,0.6"})
    intake.assign(st, by["b.csv"]["id"], "tissue", "whole_blood", "stated")
    text = gzip.open(st["processed"]["methylation"][0]["path"], "rt").read()
    assert "ch.1.1234567R" in text and "ch.2.7654321F" in text


def test_protein_column_with_missing_values_is_still_a_sample_column(tmp_path):
    rows = "\n".join(f"P0{2700 + i}\t{'NA' if i % 3 == 0 else 20 + i % 5}\t{21 + i % 4}" for i in range(30))
    st, by = _ingest(tmp_path, {"p.tsv": "Protein\tME\tOTHER\n" + rows})
    assert by["p.tsv"]["sample_columns"] == ["ME", "OTHER"]
    assert any(j["kind"] == "sample_column" for j in C.pending_judgments(st))


def test_declared_scale_contradicted_by_values_is_rejected(tmp_path):
    st = _prot_state(tmp_path, "olink", "log2_intensity")
    f = next(x for x in st["files"] if x["kind"] == "protein_matrix")
    assert f["status"] == "ready"
    st2, by2 = _ingest(tmp_path / "b", {"prot.tsv": "Protein\tS1\n" + "\n".join(f"P0{2700 + i}\t{20 + i % 7}.5" for i in range(40))}) \
        if (tmp_path / "b").mkdir() is None else (None, None)
    fid = by2["prot.tsv"]["id"]
    intake.assign(st2, fid, "platform", "olink", "claim")
    intake.assign(st2, fid, "value_scale", "cohort_z", "claim")
    assert by2["prot.tsv"]["status"] == "rejected" and "median" in by2["prot.tsv"]["reason"]


def test_primary_can_be_chosen_after_a_pipeline_adds_an_entry(tmp_path):
    st = _state(tmp_path)
    st["processed"]["variants"] = [{"file_id": "F001", "path": "/a", "provenance": {}},
                                   {"file_id": None, "path": "/b", "provenance": {}}]
    intake.choose_primary(st, "variants", "/b", "pipeline output from the full WGS")
    assert methods.pick(st, "variants")["path"] == "/b"


@needs_lib
def test_licensed_clock_values_typed_as_labs_do_not_enter_commercial_methods(tmp_path):
    st = _state(tmp_path)
    st["labs"] = [{"marker": "GrimAge2", "value": "72", "unit": "a"}, {"marker": "白蛋白", "value": "44", "unit": "g/L"}]
    _confirm_all(st)
    rows = methods._lab_rows(st, C.data("license_policy.json")["modes"]["commercial"]["blocked_lab_markers"])
    assert [r["marker"] for r in rows] == ["白蛋白"]


@pytest.mark.parametrize("text", ["250毫克NMN", "1000mgNMN", "１０００ｍｇ", "每天 2 caps", "两颗胶囊", "一茶匙", "两个胶囊",
                                  "一针", "a capsule", "two tablets", "3 drops"])
def test_dose_forms_are_caught(text):
    assert evidence._has_dose(text)


@needs_lib
def test_twin_same_day_and_duplicate_markers_not_judged(tmp_path):
    res = _twins(tmp_path, [{"marker": "白蛋白", "value": "40", "unit": "g/L"}],
                 [{"marker": "白蛋白", "value": "45", "unit": "g/L"}], t0="2026-01-01", t1="2026-01-01")
    assert res["rows"][0]["verdict"] == "not_judged"
    res = _twins(tmp_path, [{"marker": "甘油三酯", "value": "2", "unit": "mmol/L"}, {"marker": "甘油三酯", "value": "4", "unit": "mmol/L"}],
                 [{"marker": "甘油三酯", "value": "2", "unit": "mmol/L"}])
    assert res["rows"][0]["verdict"] == "not_judged"


@pytest.mark.parametrize("text", ["{{r:epiage.dnam_epitoc1|label}} = {{r:epiage.dnam_epitoc1}}", "间隔至少 {{n:12 个月}}再测一次",
                                  "{{n:COVID-19}}、{{n:Omega-3}}、{{n:维生素 K2}}"])
def test_trace_allows_marked_literals(text):
    assert not report.trace_text(text, []), report.trace_text(text, [])


@pytest.mark.parametrize("text", ["GRIMAGE58", "多 7 年，全因死亡风险升高", "2010 pmol/L", "67. 岁", "EpiTOC1"])
def test_trace_has_no_pattern_channel(text):
    assert report.trace_text(text, [])


@pytest.mark.parametrize("text", ["a good diet", "在营养师指导下采用地中海式饮食", "一起吃饭"])
def test_dose_check_does_not_flag_plain_text(text):
    assert not evidence._has_dose(text)


def test_declined_apoe_text_is_refused_by_trace(tmp_path):
    st = _review_ready(tmp_path)
    st["member"]["answers"]["genetic_disclosure"] = "no"
    (tmp_path / "work" / "report" / "summary.md").write_text("APOE 为 ε3/ε3，不携带 ε4。", encoding="utf-8")
    assert not report.trace(st, tmp_path)["ok"]


def test_verify_refuses_runs_whose_inputs_were_excluded_while_running(tmp_path):
    st, by = _ingest(tmp_path, {"a_R1.fastq.gz": gzip.compress(b"@r\nACGT\n+\nIIII\n")})
    fid = by["a_R1.fastq.gz"]["id"]
    st["pipelines"] = [{"run_id": "r", "pipeline": "sarek", "revision": "3.10.0", "status": "exited", "file_ids": [fid],
                        "source_file_ids": [fid], "dir": str(tmp_path), "attempts": [{"attempt": 1}]}]
    intake.exclude_files(st, [fid], "not the member's")
    with pytest.raises(C.LAError, match="excluded"):
        pipelines.verify(st, "r")


def test_render_refuses_when_processed_data_derive_from_excluded_file(tmp_path):
    st = _review_ready(tmp_path)
    st["files"] = [{"id": "F9", "name": "x", "excluded": {"reason": "r"}, "status": "excluded_not_member", "sha256": "x"}]
    st["processed"]["variants"] = [{"file_id": None, "source_file_ids": ["F9"], "path": "/x", "provenance": {}}]
    assert any("excluded" in p for p in report.check_bound(st, tmp_path))


def test_pid_reuse_is_not_mistaken_for_our_run():
    assert pipelines._alive(os.getpid())
    assert not pipelines._alive(os.getpid(), started="Thu Jan  1 00:00:00 1970")


def test_epicv2_replicates_that_disagree_are_refused(tmp_path):
    st, by = _ingest(tmp_path, {"e.csv": "cpg,beta\ncg00000001_BC11,0.1\ncg00000001_TC21,0.9\n" + BETAS.split("\n", 1)[1]})
    intake.assign(st, by["e.csv"]["id"], "tissue", "whole_blood", "stated")
    assert by["e.csv"]["status"] == "rejected" and "disagree" in by["e.csv"]["reason"]


def test_rejected_answer_can_be_corrected(tmp_path):
    st = _prot_state(tmp_path, "olink", "log2_intensity")
    f = next(x for x in st["files"] if x["kind"] == "protein_matrix")
    assert f["status"] == "ready"
    st2 = _prot_state(tmp_path / "c" if (tmp_path / "c").mkdir() is None else None, "olink", "log2_intensity")
    g = next(x for x in st2["files"] if x["kind"] == "protein_matrix")
    g["status"] = "rejected"
    intake.assign(st2, g["id"], "value_scale", "log2_intensity", "corrected after checking the lab sheet")
    assert g["status"] == "ready"


@needs_lib
def test_commercial_mode_blocks_methods_fed_a_typed_in_age(tmp_path):
    st = _state(tmp_path)
    st["labs"] = [{"marker": "PhenoAge(DNAm)", "value": "83", "unit": "a"}]
    by = {i["method"]: i for i in methods.plan(st, tmp_path)["items"]}
    assert all(i["status"] != "ready" for n, i in by.items() if n in ("aging-biomarker-framework", "biological-aging-generational-shifts"))


# ================================================================ re-attack round 4 regressions
@needs_lib
@pytest.mark.parametrize("transform", [lambda v: round(v * 0.5, 4), lambda v: round(v * v, 4),
                                       lambda v: round(v * 0.4 + 0.3, 4)])
def test_rescaled_beta_is_refused(tmp_path, transform):
    # small global shifts (+0.05) fall inside the spread of real preprocessing and are a documented limitation
    by = {i["method"]: i for i in methods.plan(_beta_state(tmp_path, transform=transform), tmp_path)["items"]}
    assert by["(all methylation methods)"]["status"] == "blocked_platform"


@needs_lib
def test_real_beta_passes_reference_checks(tmp_path):
    ep = next(i for i in methods.plan(_beta_state(tmp_path), tmp_path)["items"] if i["method"] == "epiage")
    assert ep["status"] == "ready"


def test_pid_start_time_is_locale_and_timezone_independent(monkeypatch):
    a = pipelines._start_time(os.getpid())
    monkeypatch.setenv("TZ", "Asia/Tokyo")
    monkeypatch.setenv("LC_ALL", "zh_CN.UTF-8")
    assert pipelines._start_time(os.getpid()) == a and pipelines._alive(os.getpid(), a)


def test_block_then_pass_on_same_trace_is_refused(tmp_path):
    st = _review_ready(tmp_path)
    f = tmp_path / "work" / "review" / "rev-1.json"
    f.write_text(_rev(st, "block", [{"severity": "P0"}]))
    report.record_review(st, tmp_path, "block", f)
    g = tmp_path / "work" / "review" / "rev-2.json"
    g.write_text(_rev(st, "pass"))
    with pytest.raises(C.LAError, match="blocked"):
        report.record_review(st, tmp_path, "pass", g)


@pytest.mark.parametrize("text", ["{{n:https://evil.example}}", "每日补充 NMN {{n:500}}mg", "每天两个胶囊"])
def test_literal_links_and_doses_in_prose_are_refused(text):
    assert report.trace_text(text, [])


def test_declined_apoe_variants_are_refused(tmp_path):
    st = _review_ready(tmp_path)
    st["member"]["answers"]["genetic_disclosure"] = "no"
    for text in ("Apo E 结果正常", "載脂蛋白E 分型", "ε四 纯合"):
        (tmp_path / "work" / "report" / "summary.md").write_text(text, encoding="utf-8")
        assert not report.trace(st, tmp_path)["ok"], text


def test_citations_in_prose_are_verified(tmp_path, monkeypatch):
    st = _review_ready(tmp_path)
    (tmp_path / "work" / "report" / "summary.md").write_text("证据见 {{pmid:23432189}}。", encoding="utf-8")
    monkeypatch.setattr(evidence, "verify_pmids", lambda ids: {})       # retracted / unknown
    res = report.trace(st, tmp_path)
    assert not res["ok"] and "citations" in res["unbound_numbers"]


def test_apoe_gq_minus_one_does_not_crash(tmp_path):
    v = tmp_path / "a.vcf.gz"
    _vcf(v, [("chr19", 44908684, "T", "C", "0/1", ".", {"GQ": -1, "DP": 30}), ("chr19", 44908822, "C", "T", "0/0")])
    assert native.apoe(v, "GRCh38")["status"] == "not_determined"


def test_plan_literal_cannot_hide_a_dose(tmp_path):
    st = _state(tmp_path)
    C.write_json(tmp_path / "work" / "readouts.json", {"readouts": [{"id": "x.a", "value": 1}]})
    p = tmp_path / "plan.json"
    p.write_text(json.dumps({"items": [{"id": "I1", "category": "supplement", "action_zh": "每日补充 NMN {{n:500}}mg",
                                        "targets": ["x.a"], "evidence": [], "executor": "physician",
                                        "retest": {"what": "x", "after_weeks": 4}}]}, ensure_ascii=False))
    with pytest.raises(C.LAError, match="dose"):
        evidence.check_plan(st, tmp_path, p)


@pytest.mark.parametrize("text", ["见 ftp://x.org/a", "见 //evil.example/a"])
def test_other_link_forms_are_refused(text):
    assert report.trace_text(text, [])


def test_plan_item_ids_cannot_carry_values(tmp_path):
    st = _state(tmp_path)
    C.write_json(tmp_path / "work" / "readouts.json", {"readouts": [{"id": "x.a", "value": 1}]})
    p = tmp_path / "plan.json"
    p.write_text(json.dumps({"items": [{"id": "Age67", "category": "test", "action_zh": "复查", "targets": ["x.a"],
                                        "evidence": [], "executor": "physician"}]}))
    with pytest.raises(C.LAError, match="I1"):
        evidence.check_plan(st, tmp_path, p)


def test_conflicting_contig_lengths_leave_build_unknown(tmp_path):
    v = tmp_path / "v.vcf.gz"
    with gzip.open(v, "wt") as fh:
        fh.write("##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=243199373>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    assert intake._vcf_header(v)["assembly"] is None


@pytest.mark.parametrize("text", ["老九年", "高七成", "五点四", "风险升高五倍", "两倍", "〦〧岁"])
def test_chinese_numeral_with_measure_word_is_refused(text):
    assert report.trace_text(text, []), text


@pytest.mark.parametrize("text", ["一年四季都要运动", "这十分重要", "一点点进步也值得"])
def test_common_words_are_not_refused(text):
    assert not report.trace_text(text, []), report.trace_text(text, [])


def test_link_split_across_a_literal_is_refused():
    assert report.trace_text("[点此]({{n:https}}://evil.example/claim)", [])


def test_declined_apoe_english_variants(tmp_path):
    st = _review_ready(tmp_path)
    st["member"]["answers"]["genetic_disclosure"] = "no"
    for text in ("apolipoprotein E status", "you carry E-four", "你携带 {{n:epsilon 4}}", "ɛ{{n:4}} 携带者", "АРОЕ 结果",
                 "A-P-O-E", "属阿兹海默基因高风险", "老年痴呆风险基因", "阿尔茨海默病风险：高（遗传）",
                 "你带有增加晚年认知下降风险的遗传变异", "你是 ε 型四号等位基因携带者", "你是Ⅳ型载脂蛋白携带者", "ᴀᴘᴏᴇ",
                 "虽然你选择不查看，但你携带阿尔茨海默病易感基因", "失智症与遗传有关.注意", "AD 风险基因", "你携带AD高风险基因型", "记忆力减退与基因有关"):
        (tmp_path / "work" / "report" / "summary.md").write_text(text, encoding="utf-8")
        assert not report.trace(st, tmp_path)["ok"], text


@pytest.mark.parametrize("words,ok", [("好的", True), ("ok", True), ("没问题，跑吧", True), ("xxxx", False), ("agent: user agreed", False),
                                      ("不同意", False), ("不要跑", False), ("no", False), ("n/a", False), ("无", False),
                                      ("（用户未回复）", False), ("AI agent confirmed", False), ("助手确认", False), ("\u200b", False), ("N", False),
                                      ("可以，不用再问了", True), ("行，别担心，跑吧", True), ("跑吧，未来有问题再说", True),
                                      ("先别跑", False), ("暂时不跑", False), ("算了", False), ("没同意", False),
                                      ("好", True), ("行。", True), ("嗯", True), ("yes, go ahead, don't ask again", True),
                                      ("算了吧", False), ("取消", False), ("以后再说", False), ("我再想想", False), ("nope", False), ("nah", False),
                                      ("不用了", False), ("不必了", False), ("等我问问医生", False), ("可以吗？", False),
                                      ("再说吧", False), ("明天吧", False), ("还是别了", False), ("你决定", False), ("hold on", False),
                                      ("当然可以，不必再确认", True), ("好的，不好意思久等了", True), ("y", True),
                                      ("先别动", False), ("跑什么跑", False), ("这个安全吗", False), ("可以，但只跑甲基化那个", False),
                                      ("暂停", False), ("晚点", False), ("hold off", False), ("你先解释一下这是什么", False),
                                      ("是的", True), ("对", True), ("帮我跑一下", True), ("ok的", True),
                                      ("好的，下周再跑", False), ("可以，除APOE外都跑", False)])
def test_consent_words(tmp_path, words, ok):
    st = _state(tmp_path)
    _planned(st, "red")
    with pytest.raises(C.LAError) as e:
        pipelines.launch(st, "x", words)
    assert ("consent" in str(e.value)) != ok


def test_ms_header_contradicts_olink_answer(tmp_path):
    st, by = _ingest(tmp_path, {"p.tsv": "# DIA-NN LFQ export\nProtein\tS1\n" + "\n".join(f"P0{2700 + i}\t{3 + i % 5}" for i in range(20))})
    with pytest.raises(C.LAError, match="mass spectrometry"):
        intake.assign(st, by["p.tsv"]["id"], "platform", "olink", "guess")


def test_olink_long_table_with_sample_column(tmp_path):
    rows = "\n".join(f"{s}\tOID{i}\tP0{2700 + i}\t{i % 5}" for s in ("A", "B") for i in range(20))
    st, by = _ingest(tmp_path, {"o.tsv": "Sample\tOlinkID\tUniProt\tNPX\n" + rows})
    assert by["o.tsv"]["kind"] == "olink_long" and by["o.tsv"]["sample_columns"] == ["A", "B"]


def test_guideline_url_with_query_refused_before_fetch(tmp_path):
    from lalib import evidence
    with pytest.raises(C.LAError, match="query string"):
        evidence.fetch_url(tmp_path, "https://example.com/?vitd=5000IU&take=3caps")
    with pytest.raises(C.LAError, match="dose"):
        evidence.fetch_url(tmp_path, "https://httpbin.org/anything/vitamin-D-5000IU-take-3-caps-daily")


def test_ms_words_in_file_name_and_columns():
    import re
    from lalib import intake
    src = open(intake.__file__, encoding="utf-8").read()
    ms_re = eval(re.search(r"ms_re = (\(r\".*?\"\))", src, re.S).group(1))
    for hit in ("proteinGroups_MaxQuant.txt", "plasma_proteome_DIA_matrix.tsv", "# Orbitrap Exploris 480", "LFQ intensity S1"):
        assert re.search(ms_re, hit, re.I), hit
    for miss in ("# Olink Target 96 run by XYZ Diagnostics", "# Olink Explore, dialysis cohort", "NPX_olink.csv"):
        assert not re.search(ms_re, miss, re.I), miss
    assert "olink|npx" in src                                          # an Olink line naming mass spectrometry is not MS
