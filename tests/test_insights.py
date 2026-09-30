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
    v = genomics.Vcf(p, RULES)
    v.prefetch({("1", x) for x in range(0, 2000)}, [("1", 0, 5000)])
    return v


RULES = {"min_gq": 20, "min_dp": 10, "min_alt_reads": 3}


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
    pubdata._cache(ws, "mr", "k", "u", [{"ref": "mr:ldl->chd:IVW", "b": math.log(1.5), "se": 0.05, "exposure": "LDL cholesterol", "outcome": "CHD"}])
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


# ---------------------------------------------------------------- v0.6 adversarial fixes
def test_split_multiallelic_rows_are_all_read(tmp_path):
    v = _vcf(tmp_path, ["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30",
                        "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30"])
    assert v.alt_dosage("1", 100, "A", "T")["dosage"] == 1          # the second row, not "0 from the first row"
    assert v.alt_dosage("1", 100, "A", "C") == {"dosage": 0, "status": "called_other_alt"}


def test_pick_hit_multiallelic():
    h = lambda r, a, f=0.2: {"pos": 1, "ref_allele": r, "alt_allele": a, "af_eas": f, "chrom": "1", "ref": "x"}
    hit, how = genomics.pick_hit([h("A", "G"), h("A", "T", 0.05)], "T")
    assert how == "effect_is_alt" and hit["p_effect"] == 0.05
    hit, how = genomics.pick_hit([h("A", "G"), h("A", "T", 0.05)], "A")          # REF effect: 1 - every alt
    assert how == "effect_is_ref" and abs(hit["p_effect"] - 0.75) < 1e-9 and hit["alt_allele"] is None
    assert genomics.pick_hit([h("A", "G")], "C") == (None, "effect_allele_not_in_record")
    assert genomics.pick_hit([h("A", "G", None)], "G") == (None, "no_eas_frequency")
    assert genomics.pick_hit([h("A", "G"), dict(h("C", "T"), pos=9)], "G") == (None, "rsid_maps_to_several_positions")


def test_ref_dosage_joined_and_split(tmp_path):
    v = _vcf(tmp_path, ["chr1\t100\t.\tA\tG,T\t50\tPASS\t.\tGT:GQ:DP\t1/2:40:30",
                        "chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30",
                        "chr1\t200\t.\tA\tT\t50\tPASS\t.\tGT:GQ:DP\t0/0:40:30",
                        "chr1\t300\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:5:30"])
    assert v.ref_dosage("1", 100, "A")["dosage"] == 0
    assert v.ref_dosage("1", 200, "A")["dosage"] == 1
    assert v.ref_dosage("1", 300, "A")["dosage"] is None
    assert v.ref_dosage("1", 999, "A") == {"dosage": None, "status": "not_called"}


def test_carried_only_alleles_present_and_quality(tmp_path):
    v = _vcf(tmp_path, ["chr1\t100\t.\tA\tG,T\t50\tPASS\t.\tGT:AD:GQ:DP\t0/2:10,0,12:40:30",
                        "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:28,2:40:30",
                        "chr1\t300\t.\tC\tT\t50\tLowQual\t.\tGT:AD:GQ:DP\t1/1:0,30:40:30",
                        "chr1\t400\t.\tCA\tC\t50\tPASS\t.\tGT:AD:GQ:DP\t1/1:0,30:40:30"])
    c = {(x["pos"], x["alt"]): x for x in v.carried("1", 1, 1000)}
    assert (100, "G") not in c and c[(100, "T")]["zygosity"] == "het"            # 0/2 carries T only
    assert c[(200, "T")]["quality_ok"] is False and c[(200, "T")]["quality"] == "few_alt_reads"
    assert c[(300, "T")]["quality_ok"] is False
    assert c[(400, "C")]["zygosity"] == "hom" and c[(400, "C")]["quality_ok"]


def test_no_tabix_binary_falls_back_to_streaming(tmp_path, monkeypatch):
    import gzip
    p = tmp_path / "m.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t1/1:40:30\n")
    (tmp_path / "m.vcf.gz.tbi").write_bytes(b"x")
    monkeypatch.setattr(genomics.shutil, "which", lambda _: None)
    v = genomics.Vcf(p, RULES)
    assert v._tabix is False
    v.prefetch({("1", 100)}, [])
    assert v.alt_dosage("1", 100, "A", "G")["dosage"] == 2


def test_absent_as_ref_is_explicit(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
    v = genomics.Vcf(p, RULES, absent_is_ref=True)
    v.prefetch(set(), [])
    assert v.alt_dosage("1", 5, "A", "G") == {"dosage": 0, "status": "assumed_ref"}


def test_clinvar_tier():
    ok = {"classification": "Pathogenic/Likely pathogenic", "review_status": "criteria provided, multiple submitters, no conflicts"}
    assert genomics.clinvar_tier(ok) == "plp"
    assert genomics.clinvar_tier(dict(ok, review_status="no assertion criteria provided")) == "plp_not_unanimous"
    assert genomics.clinvar_tier(dict(ok, classification="Conflicting classifications of pathogenicity",
                                      review_status="criteria provided, conflicting classifications")) == "plp_not_unanimous"
    assert genomics.clinvar_tier(dict(ok, classification="Pathogenic; risk factor")) == "plp_not_unanimous"
    assert genomics.clinvar_tier(dict(ok, review_status="reviewed by expert panel")) == "plp"
    assert genomics.clinvar_tier(dict(ok, classification="Benign")) is None


def _explain_setup(tmp_path, monkeypatch, vcf_lines, gene_ok=True, clinvar_ok=True, clinvar=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    ws, st = _ws(tmp_path)
    st["member"]["answers"]["genetic_disclosure"] = "yes"
    p = tmp_path / "g.vcf"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + "\n".join(vcf_lines) + "\n")
    from lalib import methods
    monkeypatch.setattr(methods, "pick", lambda st, k: {"path": str(p)})
    monkeypatch.setattr(genomics, "member_labs", lambda st: {"ldl": {"analyte": "ldl", "label_zh": "LDL", "marker": "LDL", "value": 3.8,
                                                                     "unit": "mmol/L", "ref_range": "<3.4", "flag": "high"}})
    loci = [{"rsid": f"rs{i}", "p_value": 1e-10, "effect_allele": "G", "direction": 1.0, "mapped_genes": [], "beta_text": "",
             "ref": f"gwas:x:rs{i}-G", "pubmed_id": "1"} for i in range(10)]
    monkeypatch.setattr(pubdata, "gwas_trait", lambda ws, efo, pmax: {"records": loci})
    monkeypatch.setattr(pubdata, "myvariant_rsids", lambda ws, rs, assembly: {
        f"rs{i}": [{"pos": 1_000_000 * (i + 1), "chrom": "1", "ref_allele": "A", "alt_allele": "G", "af_eas": 0.3, "ref": f"myvariant:{i}"}] for i in range(10)})

    def region(ws, g, asm):
        if not gene_ok:
            return None
        return {"chrom": "19", "start": 100, "end": 900} if g == "LDLR" else {"chrom": "2", "start": 100, "end": 900}
    monkeypatch.setattr(pubdata, "gene_region", region)

    def cv(ws, gene):
        if not clinvar_ok:
            raise C.LAError("eutils down", 3)
        if gene != "LDLR" or not clinvar:
            return []
        return [dict({"ref": "clinvar:VCV1", "gene": gene, "title": "LDLR c.1G>A", "spdi": "NC_000019.10:199:G:A", "deleted": "G",
                      "inserted": "A", "loc": {"GRCh38": {"chr": "19", "start": 200, "stop": 200}}, "conditions": []}, **clinvar),
                dict({"ref": "clinvar:VCV2", "gene": gene, "title": "LDLR c.9dup", "spdi": "NC_000019.10:300:G:GG", "deleted": "G",
                      "inserted": "GG", "loc": {"GRCh38": {"chr": "19", "start": 301, "stop": 301}}, "conditions": []}, **clinvar)]
    monkeypatch.setattr(pubdata, "clinvar_gene", cv)
    monkeypatch.setattr(pubdata, "canonical_spdi", lambda ws, asm, c, p, r, a: {"spdi": "NC_000019.10:300:G:GG" if (p, r, a) == (299, "C", "CG") else None, "warning": None})
    return ws, st


def test_explain_failure_is_not_scanned_not_zero(tmp_path, monkeypatch):
    line = ["chr19\t200\t.\tG\tA\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:15,15:40:30"]
    ws, st = _explain_setup(tmp_path, monkeypatch, line, clinvar_ok=False)
    res = genomics.explain(st, ws, absent_as_ref="交付说明：WGS 30x 联合分型，列出全部非参考位点")
    assert "LDLR" in res["analytes"]["ldl"]["not_scanned"]            # ClinVar down: not scanned, never "0 found"
    assert not any(r["id"] == "gen.ldl.clinvar_plp" for r in res["readouts"])
    ws, st = _explain_setup(tmp_path / "b", monkeypatch, line, gene_ok=False)
    res = genomics.explain(st, ws)
    assert "LDLR" in res["analytes"]["ldl"]["not_scanned"] and not any(r["id"] == "gen.ldl.clinvar_plp" for r in res["readouts"])


def test_explain_variant_only_vcf_needs_judgment(tmp_path, monkeypatch):
    ws, st = _explain_setup(tmp_path, monkeypatch, [])
    res = genomics.explain(st, ws)
    assert res["analytes"]["ldl"]["pct_eas"] is None and "absent-as-ref" in res["analytes"]["ldl"]["not_computed"]
    res2 = genomics.explain(st, ws, absent_as_ref="交付说明：WGS 30x 联合分型，列出全部非参考位点")
    assert res2["analytes"]["ldl"]["pct_eas"] is not None
    with pytest.raises(C.LAError, match="reason"):
        genomics.explain(st, ws, absent_as_ref="yes")


def test_explain_recessive_het_is_carrier(tmp_path, monkeypatch):
    tm = C.data("trait_map.json")
    ldlr_mode = tm["gene_modes"]["LDLR"]
    assert ldlr_mode == "AD"
    ws, st = _explain_setup(tmp_path, monkeypatch, ["chr19\t200\t.\tG\tA\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:15,15:40:30"],
                            clinvar={"classification": "Pathogenic", "review_status": "criteria provided, single submitter"})
    res = genomics.explain(st, ws, absent_as_ref="交付说明：WGS 30x 联合分型，列出全部非参考位点")
    gp = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]
    assert gp["pathogenic_or_likely"][0]["carrier_only"] is False and gp["pathogenic_or_likely"][0]["zygosity"] == "het"
    assert next(r for r in res["readouts"] if r["id"] == "gen.ldl.clinvar_plp")["value"] == 1
    real = genomics.data

    def as_recessive(name):
        d = real(name)
        if name == "trait_map.json":
            d["gene_modes"]["LDLR"] = "AR"                              # the same het in a recessive gene is a carrier
        return d
    monkeypatch.setattr(genomics, "data", as_recessive)
    AR = "交付说明：WGS 30x 联合分型，列出全部非参考位点"
    genomics.explain(st, ws, absent_as_ref=AR)
    gp = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]
    assert gp["pathogenic_or_likely"][0]["carrier_only"] is True
    assert next(r for r in genomics.explain(st, ws, absent_as_ref=AR)["readouts"] if r["id"] == "gen.ldl.clinvar_plp")["value"] == 0


def test_board_edit_after_registration_is_caught(tmp_path):
    ws, st = _ws(tmp_path)
    qs = {"questions": [{"id": f"Q{i}", "title_zh": "问题", "hypothesis_zh": "假设", "basis": ["native.gut.gmhi"], "why_zh": "因为"} for i in (1, 2, 3)]}
    C.write_json(ws / "q.json", qs)
    board.register_questions(st, ws, ws / "q.json")
    C.write_json(ws / "work" / "insights" / "board" / "Q1.json", {"id": "Q1", "verdict": "insufficient", "confidence": "low",
                 "member_evidence": ["native.gut.gmhi"], "public_evidence": [], "summary_zh": "证据不足。", "next_step_zh": "复查。", "limitations_zh": "单次。"})
    board.register_finding(st, ws, "Q1")
    assert board.integrity(ws, st) == []
    f = C.load_json(ws / "work" / "insights" / "board" / "Q1.json")
    C.write_json(ws / "work" / "insights" / "board" / "Q1.json", dict(f, verdict="supported"))
    assert any("Q1" in x for x in board.integrity(ws, st))
    q = C.load_json(ws / "work" / "insights" / "board" / "questions.json")
    q["questions"][0]["title_zh"] = "改过"
    C.write_json(ws / "work" / "insights" / "board" / "questions.json", q)
    assert any("questions.json" in x for x in board.integrity(ws, st))


def test_board_basis_only_usable_rows(tmp_path):
    ws, st = _ws(tmp_path)
    st["labs"].append({"marker": "被拒绝的行", "value": "1", "unit": "", "confirm": {"answer": "no", "row_key": "x", "why": "t"}})
    assert "被拒绝的行" not in board.member_ids(st, ws) and "低密度脂蛋白胆固醇" in board.member_ids(st, ws)


def test_projection_baseline_and_exposure_guards(tmp_path):
    ws, st = _ws(tmp_path)
    C.write_json(ws / "work" / "readouts.json", {"readouts": [{"id": "china-par-ascvd-risk.risk_10y_pct", "value": 6.0, "unit": "%"},
                                                             {"id": "x.score", "value": 0.4, "unit": ""}]})
    pubdata._cache(ws, "mr", "k", "u", [{"ref": "mr:a->b:IVW", "b": 0.4, "se": 0.04, "exposure": "LDL cholesterol", "outcome": "CHD"},
                                        {"ref": "mr:c->b:IVW", "b": 0.4, "se": 0.04, "exposure": "Apolipoprotein B", "outcome": "CHD"}])
    with pytest.raises(C.LAError, match="absolute risk"):
        causal.project(st, ws, "mr:a->b:IVW", "ldl", 2.6, "x.score")
    with pytest.raises(C.LAError, match="exposure-match"):
        causal.project(st, ws, "mr:c->b:IVW", "ldl", 2.6, "china-par-ascvd-risk.risk_10y_pct")
    rec = causal.project(st, ws, "mr:c->b:IVW", "ldl", 2.6, "china-par-ascvd-risk.risk_10y_pct", exposure_match="ApoB 与 LDL-C 高度共线，作为同一脂蛋白负荷的替代")
    assert rec["ref"].endswith(":china-par-ascvd-risk.risk_10y_pct") and len(rec["caveats_zh"]) == 2
    with pytest.raises(C.LAError, match="no population reference"):
        causal.project(st, ws, "mr:a->b:IVW", "lpa", 30, "china-par-ascvd-risk.risk_10y_pct")


def test_percentile_ties_report_span():
    r = reference._pct_of(0.2, [1, 5, 10, 25, 50], [0.2, 0.2, 0.2, 0.5, 1.0])
    assert r["bound"] == "tie" and r["pct_low"] == 1 and r["pct_high"] == 10


def test_dup_matched_by_spdi_and_variant_only_vcf_not_scanned(tmp_path, monkeypatch):
    lines = ["chr19\t299\t.\tC\tCG\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:15,15:40:30"]
    ws, st = _explain_setup(tmp_path, monkeypatch, lines, clinvar={"classification": "Pathogenic", "review_status": "criteria provided, single submitter"})
    res = genomics.explain(st, ws)                                  # variant-only VCF, no judgment recorded
    ms = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]
    assert [x["clinvar_ref"] for x in ms["pathogenic_or_likely"]] == ["clinvar:VCV2"]   # the dup is found (P0-1)
    assert "LDLR" in ms["not_scanned"]                               # but the gene is not called fully scanned (P0-2)
    assert not any(r["id"] == "gen.ldl.clinvar_plp" for r in res["readouts"])


def test_gvcf_covered_region_counts_as_scanned(tmp_path, monkeypatch):
    lines = ["chr19\t1\t.\tN\t<NON_REF>\t.\tPASS\tEND=5000\tGT:GQ:DP\t0/0:40:30"] + \
            [f"chr2\t1\t.\tN\t<NON_REF>\t.\tPASS\tEND=5000\tGT:GQ:DP\t0/0:40:30"]
    ws, st = _explain_setup(tmp_path, monkeypatch, lines, clinvar={"classification": "Pathogenic", "review_status": "criteria provided, single submitter"})
    res = genomics.explain(st, ws)
    assert res["analytes"]["ldl"]["not_scanned"] == []
    assert next(r for r in res["readouts"] if r["id"] == "gen.ldl.clinvar_plp")["value"] == 0


def test_absent_as_ref_never_overrides_a_record(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 "chr1\t100\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30\n"
                 "chr1\t200\t.\tA\tT\t50\tLowQual\t.\tGT:GQ:DP\t0/1:40:30\n"
                 "chr1\t290\t.\tATTTT\tA\t50\tPASS\t.\tGT:GQ:DP\t1/1:40:30\n")
    v = genomics.Vcf(p, RULES, absent_is_ref=True)
    v.prefetch({("1", 100), ("1", 200), ("1", 292), ("1", 500)}, [])
    assert v.alt_dosage("1", 100, "A", "G")["status"] == "ref_mismatch"          # REF differs: not "assumed reference"
    assert v.alt_dosage("1", 200, "A", "G")["dosage"] is None                    # another allele's call failed here
    assert v.ref_dosage("1", 292, "T")["status"] == "within_carried_deletion"
    assert v.alt_dosage("1", 500, "A", "G") == {"dosage": 0, "status": "assumed_ref"}


def test_hemizygous_and_half_calls(tmp_path):
    v = _vcf(tmp_path, ["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t1:40:30",
                        "chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t./1:40:30"])
    c = {x["pos"]: x for x in v.carried("1", 1, 1000)}
    assert c[100]["zygosity"] == "hemi" and c[100]["quality_ok"]
    assert c[200]["quality_ok"] is False and c[200]["quality"] == "half_call"


def test_gvcf_blocks_kept_only_where_wanted(tmp_path):
    lines = [f"chr1\t{i * 100 + 1}\t.\tN\t<NON_REF>\t.\tPASS\tEND={i * 100 + 100}\tGT:GQ:DP\t0/0:40:30" for i in range(2000)]
    p = tmp_path / "g.vcf"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + "\n".join(lines) + "\n")
    v = genomics.Vcf(p, RULES)
    v.prefetch({("1", 150)}, [("1", 1000, 1300)])
    assert len(v._blocks["1"]) == 5                                   # 1 block for the position + 4 for the region
    assert v.alt_dosage("1", 150, "N", "G")["status"] == "ref_block"
    assert v.region_coverage("1", 1001, 1300) == 1.0 and v.region_coverage("1", 1001, 1600) < 0.9


def test_palindromic_near_half_is_skipped():
    h = {"pos": 1, "ref_allele": "A", "alt_allele": "T", "af_eas": 0.45, "chrom": "1", "ref": "x"}
    assert genomics.pick_hit([h], "T") == (None, "palindromic_ambiguous")
    assert genomics.pick_hit([dict(h, af_eas=0.1)], "T")[1] == "effect_is_alt"


def test_forged_cache_file_is_ignored(tmp_path):
    ws, _ = _ws(tmp_path)
    C.write_json(ws / "work" / "evidence" / "pub_mr_forged.json", {"records": [{"ref": "mr:fake->x:IVW", "b": 5, "se": 0.1}]})
    assert "mr:fake->x:IVW" not in pubdata.known_refs(ws)
    pubdata._cache(ws, "mr", "real", "u", [{"ref": "mr:real->x:IVW", "b": 0.4, "se": 0.1}])
    assert "mr:real->x:IVW" in pubdata.known_refs(ws)
    p = next(x for x in (ws / "work" / "evidence").glob("pub_mr_*.json") if "forged" not in x.name)
    d = C.load_json(p)
    d["records"][0]["b"] = 9
    C.write_json(p, d)
    assert "mr:real->x:IVW" not in pubdata.known_refs(ws)            # edited after the harness wrote it


def test_projection_target_and_organ_baseline(tmp_path):
    ws, st = _ws(tmp_path)
    C.write_json(ws / "work" / "readouts.json", {"readouts": [{"id": "china-par-ascvd-risk.risk_10y_pct", "value": 6.0, "unit": "%"},
                                                             {"id": "lab.hba1c_pct", "value": 5.8, "unit": "%"}]})
    (ws / "work" / "organs").mkdir(parents=True)
    C.write_json(ws / "work" / "organs" / "organ_readouts.json", {"readouts": [{"id": "organ.heart.risk.1", "value": 0.08, "unit": "概率"}]})
    pubdata._cache(ws, "mr", "k2", "u", [{"ref": "mr:a->b:IVW", "b": 0.4, "se": 0.04, "exposure": "LDL cholesterol", "outcome": "CHD"}])
    for bad in (float("nan"), float("inf"), -50.0, 99.0):
        with pytest.raises(C.LAError, match="target"):
            causal.project(st, ws, "mr:a->b:IVW", "ldl", bad, "china-par-ascvd-risk.risk_10y_pct")
    with pytest.raises(C.LAError, match="absolute risk"):
        causal.project(st, ws, "mr:a->b:IVW", "ldl", 2.6, "lab.hba1c_pct")
    rec = causal.project(st, ws, "mr:a->b:IVW", "ldl", 2.6, "organ.heart.risk.1")
    assert rec["baseline_risk"] == 0.08 and rec["risk_after"] < 0.08


def test_percentile_edge_is_not_below():
    assert reference._pct_of(1.0, [1, 50, 99], [1.0, 5.0, 9.0]) == {"pct": 1, "bound": None}


def test_board_rejects_non_string_items(tmp_path):
    ws, st = _ws(tmp_path)
    C.write_json(ws / "q.json", {"questions": [{"id": f"Q{i}", "title_zh": "问题", "hypothesis_zh": "假设", "basis": [{"x": 1}], "why_zh": "因为"} for i in (1, 2, 3)]})
    with pytest.raises(C.LAError, match="id strings"):
        board.register_questions(st, ws, ws / "q.json")


def test_xlr_male_hemizygous_is_not_carrier_female_het_is(tmp_path, monkeypatch):
    AR = "交付说明：WGS 30x 联合分型，列出全部非参考位点"
    real = genomics.data

    def xlr(name):
        d = real(name)
        if name == "trait_map.json":
            d["gene_modes"]["LDLR"] = "XLR"
        return d
    for sex, gt, carrier in (("male", "1", False), ("female", "0/1", True)):
        ws, st = _explain_setup(tmp_path / sex, monkeypatch, [f"chr19\t200\t.\tG\tA\t50\tPASS\t.\tGT:AD:GQ:DP\t{gt}:15,15:40:30"],
                                clinvar={"classification": "Pathogenic", "review_status": "criteria provided, single submitter"})
        st["member"]["sex"] = sex
        monkeypatch.setattr(genomics, "data", xlr)
        genomics.explain(st, ws, absent_as_ref=AR)
        x = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]["pathogenic_or_likely"][0]
        assert x["carrier_only"] is carrier, sex


def test_not_unanimous_pathogenic_is_listed_not_dropped(tmp_path, monkeypatch):
    ws, st = _explain_setup(tmp_path, monkeypatch, ["chr19\t200\t.\tG\tA\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:15,15:40:30"],
                            clinvar={"classification": "Conflicting classifications of pathogenicity", "review_status": "criteria provided, conflicting classifications"})
    genomics.explain(st, ws, absent_as_ref="交付说明：WGS 30x 联合分型，列出全部非参考位点")
    ms = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]
    assert not ms["pathogenic_or_likely"] and ms["plp_not_unanimous"][0]["clinvar_ref"] == "clinvar:VCV1"


def test_insight_files_locked_and_reasons_traced(tmp_path, monkeypatch):
    from lalib import report
    ws, st = _explain_setup(tmp_path, monkeypatch, [], clinvar=None)
    genomics.explain(st, ws, absent_as_ref="交付说明写明 99.9% 位点均已检出")
    gp = ws / "work" / "insights" / "genotype_phenotype.json"
    assert report.insight_file_problems(ws, st) == ["genotype_phenotype.json was not written by `la.py insights`; run that step again"]
    st.setdefault("insights", {})["files"] = {"genotype_phenotype.json": C.sha256_file(gp, limit=None)}
    assert report.insight_file_problems(ws, st) == []
    d = C.load_json(gp)
    d["analytes"]["ldl"]["monogenic_scan"]["pathogenic_or_likely"] = [{"x": 1}]
    C.write_json(gp, d)
    assert report.insight_file_problems(ws, st)
    reasons = dict(report._authored(ws, st))
    assert "99.9%" in reasons["absent-as-ref reason"] and report.trace_text(reasons["absent-as-ref reason"])


def test_clinvar_oversized_batches_are_split(tmp_path, monkeypatch):
    ws, _ = _ws(tmp_path)
    monkeypatch.setattr(pubdata, "_shared_get", lambda k: None)
    monkeypatch.setattr(pubdata, "_shared_put", lambda k, d: None)
    import urllib.parse as up

    def fake(url, data=None, need="result"):
        if need == "esearchresult":
            return {"esearchresult": {"idlist": ["1", "2", "3", "4"]}}
        ids = up.parse_qs(data.decode())["id"][0].split(",")
        if "3" in ids:                                   # record 3 alone is too big to transform
            raise C.LAError("NCBI E-utilities returned no result: Input XML size is 26306057 bytes, and cannot be transformed", 3)
        return {"result": {"uids": ids, **{i: {"accession": f"VCV{i}", "title": "t", "germline_classification": {"description": "Pathogenic"},
                                              "variation_set": [{"canonical_spdi": f"NC_000019.10:{i}:G:A", "variation_loc": []}]} for i in ids}}}
    monkeypatch.setattr(pubdata, "_ncbi", fake)
    recs = pubdata.clinvar_gene(ws, "LDLR")
    assert sorted(r["ref"] for r in recs) == ["clinvar:VCV1", "clinvar:VCV2", "clinvar:VCV4"]
