"""la.py age-evidence: calculators for short questions on pasted measurements (synthetic inputs only)."""
import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "skills" / "longevity-analyst"
sys.path.insert(0, str(ROOT / "scripts"))

from lalib import age_evidence as AE  # noqa: E402

PANEL = ("Demographic information: {sex}. Body measurements: BMI 23.0 kg/m2. Blood panel: Albumin (g/dL) {alb}; "
         "Alkaline phosphatase (U/L) {alp}; Blood urea nitrogen (mg/dL) {bun}; Creatinine (mg/dL) {cr}; "
         "Serum glucose (mmol/L) {glu}; Total cholesterol (mg/dL) 190; Red cell distribution width (%) {rdw}; "
         "Lymphocyte percent (%) 32; White blood cell count (1000 cells/uL) 6.1; Glycohemoglobin (%) {a1c}.")
YOUNG = dict(alb=4.8, alp=60, bun=11, cr=0.8, glu=4.8, rdw=12.4, a1c=5.0)
OLD = dict(alb=3.9, alp=95, bun=24, cr=1.3, glu=6.4, rdw=14.8, a1c=6.3)


def _ref(name):
    with gzip.open(ROOT / "data" / "reference" / name, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_labs_single_person_gives_an_age_with_interval():
    out = AE.run("labs", PANEL.format(sex="Male", **OLD))
    ref = out["people"][0]["lab_age_reference"]
    assert ref["status"] == "ok"
    lo, hi = ref["outputs"]["interval_90"]
    assert 18 <= lo <= ref["outputs"]["age_estimate"] <= hi <= 80
    assert ref["outputs"]["n_values_used"] >= 5


def test_labs_two_people_older_profile_scores_older():
    text = "Participant A:\n" + PANEL.format(sex="Female", **OLD) + "\nParticipant B:\n" + PANEL.format(sex="Female", **YOUNG)
    out = AE.run("labs", text)
    assert len(out["people"]) == 2
    assert out["p_A_older_by_lab_age"] > 0.5


def test_labs_with_age_given_uses_survival_not_lab_age():
    text = "Demographics: 71-year-old male. " + PANEL.format(sex="Male", **OLD)
    p = AE.run("labs", text)["people"][0]
    assert p["lab_age_reference"]["status"] == "skipped"
    s = p["survival_cox"]["survival_probability"]
    assert 0 < s["15y"] <= s["10y"] <= s["5y"] <= 1
    assert "survival_band_probability" not in p["survival_cox"]


def test_proteome_pair_vote_follows_age_direction():
    rows = [r for r in _ref("olink_age_direction_ckb.csv.gz") if float(r["p"]) < 1e-6 and float(r["sd_per_year"]) > 0][:25]
    text = "; ".join(f"{r['gene']} (A: 1.5, B: 0.2)" for r in rows)
    vote = AE.run("proteome", text)["age_direction_vote"]
    assert vote["status"] == "ok" and vote["outputs"]["looks_older"] == "A"


def test_methylation_pair_vote_follows_age_direction():
    rows = [r for r in _ref("dnam_age_direction_gs.csv.gz") if float(r["z"]) > 10][:25]
    text = "\n".join(f"{r['cpg']} (chr1:1, +, OpenSea). A: 0.40, B: 0.70, B much higher." for r in rows)
    vote = AE.run("methylation", text)["age_direction_vote"]
    assert vote["status"] == "ok" and vote["outputs"]["looks_older"] == "B"


def test_expression_and_bad_modality_say_so():
    single = AE.run("expression", "GENE1 (3.2)")
    assert single.get("status") == "no_calibrated_estimator"
    assert AE.run("bogus", "x")["status"] == "bad_input"


def test_expression_top_list_vote_runs():
    rises = ["CDKN1A", "C1QA", "C1QB", "IL6", "SERPINE1", "GDF15", "TIMP1", "COL1A1"]
    falls = ["PRM1", "PRM2", "TNP1", "ND1", "COX1"]
    # A looks older: more senescence/inflammation genes early; fewer germline markers
    a_genes = rises + falls[2:] + [f"FILL{i}" for i in range(40)]
    b_genes = falls + rises[4:] + [f"FILL{i}" for i in range(40)]
    text = (
        "Which of the following two Testis samples comes from an older individual?\n"
        "Options: A. Sample A  B. Sample B\n"
        "Sample a: Demographic information: male\n"
        "Top-400 expressed genes: " + "; ".join(a_genes) + "\n"
        "Sample b: Demographic information: male\n"
        "Top-400 expressed genes: " + "; ".join(b_genes) + "\n"
        "GSEA results:\n"
        "Sample-A shows high activity in Complement Cascade (NES: 0.40, genes: C3, C1QA), "
        "while Sample-B shows low activity (NES: -0.10, genes: C3, C1QA). "
        "Sample-A shows low activity in Mitochondrial Translation (NES: -0.20, genes: MRPL14), "
        "while Sample-B shows high activity (NES: 0.35, genes: MRPL14).\n"
    )
    out = AE.run("expression", text)
    vote = out["age_direction_vote"]
    assert vote["status"] == "ok"
    assert vote["outputs"]["looks_older"] == "A"
    assert vote["outputs"]["gsea_pathways_used"] >= 1
    assert vote["outputs"]["fusion"] == "fused"


def test_cli_reads_text_and_emits_json():
    p = subprocess.run([sys.executable, str(ROOT / "scripts" / "la.py"), "age-evidence", "--modality", "labs",
                        "--text", PANEL.format(sex="Male", **YOUNG)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    out = json.loads(p.stdout)
    assert out["input"]["source"] == "--text" and out["people"][0]["lab_age_reference"]["status"] == "ok"
