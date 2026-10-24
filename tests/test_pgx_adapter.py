"""WS-F: VCF -> 23andMe adapter for cpic-pharmacogenomics and nutrigenetic-variant-panel.

Fixture (real data, one region only): tests/fixtures/pgx/NA18526.chr12_SLCO1B1.vcf.gz is sample NA18526 (CHB) from the
1000 Genomes 30x phased panel (20220422, 1kGP_high_coverage_Illumina.chr12.filtered.SNV_INDEL_SV_phased_panel.vcf.gz,
ftp.1000genomes.ebi.ac.uk), chr12:21177615-21179615 (GRCh38), read with pysam over tabix on 2026-10-24. At rs4149056
(chr12:21178615 T>C, SLCO1B1 c.521T>C) the panel calls 1|0. The equivalent 23andMe file is written independently from
Ensembl REST /variation/human/rs4149056?genotypes=1 (1000GENOMES:phase_3:NA18526 = C|T).
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import genomics, methods, pgx  # noqa: E402
from lalib import common as C  # noqa: E402

LIB = os.environ.get("LONGEVITY_SKILLS_HOME", "")
needs_lib = pytest.mark.skipif(not (LIB and (Path(LIB) / "skills" / "cpic-pharmacogenomics").exists()), reason="needs longevity-skills")
FIX = ROOT / "tests" / "fixtures" / "pgx"
VCF = FIX / "NA18526.chr12_SLCO1B1.vcf.gz"
RULES = {"min_gq": 20, "min_dp": 10, "min_alt_reads": 3}


def _st(tmp_path, vcf=VCF, consent="yes", assembly="GRCh38", meds=None):
    ans = {"genetic_disclosure": consent}
    if meds:
        ans["medications"] = meds
    st = {"member": {"id": "t", "age": 50, "sex": "female", "answers": ans}, "files": [], "labs": [], "stages": {},
          "processed": {"variants": [{"path": str(vcf), "provenance": {"assembly": assembly}}] if vcf else []}}
    ws = tmp_path / "ws"
    ws.mkdir()
    return st, ws


def test_fixture_is_the_recorded_sample():
    with gzip.open(VCF, "rt") as fh:
        head = [l for l in fh if l.startswith("#CHROM")][0].split()
    assert head[-1] == "NA18526"


def test_whitelist_is_separate_from_auto_tiers():
    assert methods.GENOTYPE_ADAPTER_SKILLS == {"cpic-pharmacogenomics", "nutrigenetic-variant-panel"}
    assert "B" not in methods.AUTO_TIERS


@needs_lib
def test_sites_come_from_the_skills_and_sit_in_their_genes():
    s = {x["rsid"]: x for x in pgx.sites()}
    assert (s["rs4149056"]["chrom"], s["rs4149056"]["pos"], s["rs4149056"]["ref"], s["rs4149056"]["alts"]) == ("12", 21178615, "T", ["C"])
    assert s["rs1142345"]["alts"] == ["C", "G"]                                  # two alts in the CPIC HGVS
    assert s["rs1801133"]["skill"] == "nutrigenetic-variant-panel" and s["rs1801133"]["chrom"] == "1"
    assert "rs9923231" in s                                                      # warfarin VKORC1 site


def _vcf(tmp_path, lines):
    p = tmp_path / "m.vcf"
    p.write_text("\n".join(["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"] + lines) + "\n")
    return p


SITES = [{"rsid": "rsA", "chrom": "1", "pos": 100, "ref": "T", "alts": ["C"]},
         {"rsid": "rsB", "chrom": "1", "pos": 200, "ref": "C", "alts": ["T"]},
         {"rsid": "rsC", "chrom": "1", "pos": 300, "ref": "T", "alts": ["C", "G"]},
         {"rsid": "rsD", "chrom": "1", "pos": 450, "ref": "A", "alts": ["G"]},
         {"rsid": "rsE", "chrom": "1", "pos": 900, "ref": "A", "alts": ["G"]},
         {"rsid": "rsF", "chrom": "1", "pos": 600, "ref": "G", "alts": ["A"]}]


def test_read_sites_follows_vcf_quality_rules(tmp_path):
    p = _vcf(tmp_path, ["1\t100\t.\tT\tC\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30",
                        "1\t200\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t1/1:5:30",                  # GQ 5: no call
                        "1\t300\t.\tT\tC,G\t50\tPASS\t.\tGT:GQ:DP\t1/2:40:30",
                        "1\t400\t.\tA\t<NON_REF>\t.\tPASS\tEND=500\tGT:GQ:DP\t0/0:40:25",
                        "1\t600\t.\tG\tA\t50\tLowQual\t.\tGT:GQ:DP\t0/1:40:30"])
    c = pgx.read_sites(genomics.Vcf(p, RULES), SITES)
    assert c["rsA"] == {"genotype": "CT", "status": "called"}
    assert c["rsB"]["genotype"] is None and c["rsB"]["status"] == "low_gq"
    assert c["rsC"]["genotype"] == "CG"
    assert c["rsD"] == {"genotype": "AA", "status": "ref_block"}
    assert c["rsE"] == {"genotype": None, "status": "not_called"}                 # absent is never reference by default
    assert c["rsF"]["genotype"] is None and c["rsF"]["status"].startswith("filtered")
    c2 = pgx.read_sites(genomics.Vcf(p, RULES, absent_is_ref=True), SITES)
    assert c2["rsE"] == {"genotype": "AA", "status": "assumed_ref"}


def test_23andme_writer_marks_no_calls(tmp_path):
    out = tmp_path / "g.txt"
    pgx.write_23andme(out, SITES[:2], {"rsA": {"genotype": "CT"}, "rsB": {"genotype": None}})
    rows = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
    assert rows == [["rsA", "1", "100", "CT"], ["rsB", "1", "200", "--"]]


@pytest.mark.parametrize("kw,why", [({"consent": "no"}, "consent"), ({"vcf": None}, "no_vcf"), ({"assembly": "GRCh37"}, "grch37")])
def test_not_run_with_reason(tmp_path, kw, why):
    st, ws = _st(tmp_path, **kw)
    doc = pgx.run(st, ws)
    assert doc["status"] == "not_run" and doc["reason"] == why and doc["reason_zh"]
    assert json.loads((ws / "work" / "pgx" / "pgx.json").read_text())["status"] == "not_run"
    assert not (ws / "work" / "pgx" / "genotype.txt").exists()


def _skill(genotype: Path, out: Path):
    script = Path(LIB) / "skills" / "cpic-pharmacogenomics" / "scripts" / "personal_report.py"
    subprocess.run([sys.executable, str(script), "--genotype", str(genotype), "--out", str(out)], check=True, capture_output=True)
    return json.loads((out / "result.json").read_text())


@needs_lib
def test_adapter_equals_skill_on_equivalent_23andme_file(tmp_path):
    st, ws = _st(tmp_path, meds="辛伐他汀")
    doc = pgx.run(st, ws)
    assert doc["status"] == "ran" and doc["assembly"] == "GRCh38"
    g = {l.split("\t")[0]: l.split("\t")[3] for l in (ws / "work" / "pgx" / "genotype.txt").read_text().splitlines() if not l.startswith("#")}
    assert g["rs4149056"] == "CT"
    assert g["rs4244285"] == "--"                                              # outside the fixture region: not called, not reference
    direct = _skill(FIX / "NA18526_rs4149056_ensembl.23andme.txt", tmp_path / "direct")
    adapted = json.loads((ws / "work" / "pgx" / "cpic-pharmacogenomics" / "result.json").read_text())
    assert adapted["outputs"] == direct["outputs"]                             # every CPIC output identical
    slco = next(r for r in doc["results"] if r["gene"] == "SLCO1B1")
    assert slco["phenotype"] == direct["outputs"]["slco1b1_phenotype"]["value"] == "Decreased Function"
    sim = next(d for d in slco["drugs"] if d["drug"] == "simvastatin")
    assert sim["member_takes"] and sim["guideline_url"].startswith("https://") and sim["skill_text_zh"]
    if doc["cpic_text"] == "gap":                                                # the skill cannot print CPIC text yet
        assert sim["recommendation_text"] is None and doc["cpic_text_gap_zh"]
    else:
        assert sim["recommendation_text"]
    assert "rs4149056" not in json.dumps(doc)                                    # pgx.json carries no rsid+genotype pairs
    ng = {x["key"]: x for x in doc["nutrigenetic"]}
    assert set(ng) == {"aldh2", "mthfr", "lct"} and all(x["copies"] is None and x["kind"] == "descriptive" for x in ng.values())


def test_cpic_text_rows_are_read_from_the_result(tmp_path):
    rows = [{"drug": "simvastatin", "recommendation": "verbatim", "classification": "Strong", "id": 1}]
    assert pgx._cpic_text_rows({"outputs": {"cpic_recommendations": rows}}) == rows
    assert pgx._cpic_text_rows({"cpic_text": rows}) == rows
    assert pgx._cpic_text_rows({"outputs": {}}) == []
