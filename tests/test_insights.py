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


@pytest.mark.parametrize("pos,ref,alt,want", [
    (100, "A", "G", "chr1:g.100A>G"), (100, "AT", "A", "chr1:g.101del"), (100, "ATTC", "A", "chr1:g.101_103del"),
    (100, "A", "AGG", "chr1:g.100_101insGG"), (100, "AT", "GC", "chr1:g.100_101delinsGC"), (100, "CAT", "CG", "chr1:g.101_102delinsG"),
    (100, "A", "<DEL>", None)])
def test_hgvs_ids_for_indels(pos, ref, alt, want):
    assert genomics.hgvs_id("1", pos, ref, alt) == want


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


def test_plp_requires_unanimous_reviewed():
    ok = {"clinvar_significance": ["Pathogenic", "Likely pathogenic"], "review_status": ["criteria provided, multiple submitters, no conflicts"]}
    assert genomics.plp_unanimous(ok)
    assert not genomics.plp_unanimous(dict(ok, clinvar_significance=["Pathogenic", "Uncertain significance"]))
    assert not genomics.plp_unanimous(dict(ok, review_status=["no assertion criteria provided"]))
    assert not genomics.plp_unanimous(dict(ok, review_status=["criteria provided, conflicting classifications"]))
    assert genomics.plp_unanimous(dict(ok, review_status=["reviewed by expert panel"]))


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

    def hg(ws, ids, assembly):
        if not clinvar_ok:
            raise C.LAError("myvariant down", 3)
        return {h: dict(clinvar or {}, hgvs=h, ref=f"clinvar:{h}", conditions=[], af_eas=None) for h in ids}
    monkeypatch.setattr(pubdata, "myvariant_hgvs", hg)
    return ws, st


def test_explain_failure_is_not_scanned_not_zero(tmp_path, monkeypatch):
    line = ["chr19\t200\t.\tG\tA\t50\tPASS\t.\tGT:AD:GQ:DP\t0/1:15,15:40:30"]
    ws, st = _explain_setup(tmp_path, monkeypatch, line, clinvar_ok=False)
    res = genomics.explain(st, ws)
    assert res["analytes"]["ldl"]["not_scanned"] == ["LDLR"]           # ClinVar down: not scanned, never "0 found"
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
                            clinvar={"clinvar_significance": ["Pathogenic"], "review_status": ["criteria provided, single submitter"]})
    res = genomics.explain(st, ws)
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
    genomics.explain(st, ws)
    gp = C.load_json(ws / "work" / "insights" / "genotype_phenotype.json")["analytes"]["ldl"]["monogenic_scan"]
    assert gp["pathogenic_or_likely"][0]["carrier_only"] is True
    assert next(r for r in genomics.explain(st, ws)["readouts"] if r["id"] == "gen.ldl.clinvar_plp")["value"] == 0


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
    with pytest.raises(C.LAError, match="percent"):
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
