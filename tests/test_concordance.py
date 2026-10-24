"""WS-D.1 genetic-lab concordance (la-concordance/1) and WS-D.3 fingerprint, offline."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import genomics  # noqa: E402
from lalib import common as C  # noqa: E402

LEAD = "基因分值只能解释这一项差异的一部分。"


def _res(flag, pct, plp=(), not_computed=None):
    r = {"member_lab": {"flag": flag}, "monogenic_scan": {"pathogenic_or_likely": list(plp)}}
    if pct is None:
        r["score_not_computed"] = not_computed or "only 3 lead loci vary in East Asians (need 8)"
    else:
        r["pct_eas"] = pct
    return r


def _one(key, res, pgx=None):
    doc = genomics.concordance({"generated_at": "2026-10-24T10:00:00+08:00", "analytes": {key: res}}, pgx)
    return doc, doc["items"][0]


# ldl: adverse high; hdl: adverse low
CASES = [
    ("ldl", "high", 85, "elevated", "adverse", "concordant_adverse"),
    ("ldl", "high", 50, "typical", "adverse", "acquired_adverse"),
    ("ldl", "high", 10, "reduced", "adverse", "acquired_adverse"),
    ("ldl", "in_range", 90, "elevated", "in_range", "latent"),
    ("ldl", "in_range", 40, "typical", "in_range", "concordant_typical"),
    ("hdl", "low", 12, "elevated", "adverse", "concordant_adverse"),
    ("hdl", "low", 60, "typical", "adverse", "acquired_adverse"),
    ("hdl", "low", 95, "reduced", "adverse", "acquired_adverse"),
    ("hdl", "in_range", 5, "elevated", "in_range", "latent"),
    ("hdl", "in_range", 50, "typical", "in_range", "concordant_typical"),
]


@pytest.mark.parametrize("key,flag,pct,band,state,cls", CASES)
def test_four_classes_both_directions(key, flag, pct, band, state, cls):
    doc, it = _one(key, _res(flag, pct))
    assert (it["genetic_band"], it["lab_state"], it["class"]) == (band, state, cls)
    assert it["wording_zh"] == LEAD + genomics.CLASS_ZH[cls] + "。"
    assert it["clinvar_override"] is False and it["basis"] == [f"work/insights/genotype_phenotype.json#{key}"]
    assert doc["schema"] == "la-concordance/1" and doc["rules"]["band_high_pct"] == 80 and doc["rules"]["band_low_pct"] == 20
    assert doc["trait_map_version"].startswith("la-traits/1@") and len(doc["trait_map_version"].split("@")[1]) == 12
    assert doc["retrieved_on"] == "2026-10-24"


def test_fixed_wordings_are_the_plan_table():
    assert genomics.CLASS_ZH == {
        "concordant_adverse": "这一项有较明显的先天倾向，而且已经表现出来",
        "acquired_adverse": "先天倾向不突出，这一项主要受后天因素影响，生活方式调整的空间较大",
        "latent": "有先天倾向，但目前还没有表现出来",
        "concordant_typical": "先天倾向不突出，目前也在范围内",
    }
    assert genomics.LEAD_ZH == LEAD and genomics.OVERRIDE_ZH == "需要临床确认和遗传咨询"


def test_wrong_side_or_no_range_is_not_classified():
    assert _one("ldl", _res("low", 90))[1]["class"] == "not_classified"          # low LDL is not the adverse side
    assert _one("hdl", _res("no_range", 10))[1]["class"] == "not_classified"
    assert _one("hdl", _res("no_range", 10))[1]["wording_zh"] is None


def test_not_computed():
    _, it = _one("urate", _res("high", None))
    assert it["genetic_band"] == "not_computed" and it["class"] == "not_computed" and it["wording_zh"] is None
    assert "loci" in it["not_computed_reason"] and it["near_band_edge"] is False


def test_clinvar_plp_overrides_class():
    hit = {"gene": "LDLR", "carrier_only": False}
    _, it = _one("ldl", _res("high", 50, plp=[hit]))
    assert it["clinvar_override"] is True and it["class"] == "clinvar_override"
    assert it["wording_zh"] == "需要临床确认和遗传咨询" and it["class_without_override"] == "acquired_adverse"
    _, it2 = _one("ldl", _res("high", 50, plp=[dict(hit, carrier_only=True)]))     # a carrier finding does not explain the value
    assert it2["clinvar_override"] is False and it2["class"] == "acquired_adverse"


@pytest.mark.parametrize("pct,near", [(75, True), (80, True), (85, True), (86, False), (50, False), (14, False), (15, True), (25, True)])
def test_near_band_edge(pct, near):
    assert _one("ldl", _res("high", pct))[1]["near_band_edge"] is near


@pytest.mark.parametrize("flag,pct,band", [("high", 85, "elevated"), ("high", 10, "reduced"), ("low", 10, "elevated"),
                                           ("low", 90, "reduced"), ("in_range", 90, "elevated"), ("in_range", 10, "elevated"),
                                           ("in_range", 50, "typical")])
def test_both_direction_rule(flag, pct, band):
    rules = {"band_high_pct": 80, "band_low_pct": 20}
    assert genomics.genetic_band(pct, "both", flag, rules) == band
    assert genomics.lab_state(flag, "both") == ("in_range" if flag == "in_range" else "adverse")


def test_every_lab_item_has_an_adverse_direction():
    tm = C.data("trait_map.json")
    assert all(a.get("adverse_direction") in ("high", "low", "both") for a in tm["analytes"].values())
    assert tm["analytes"]["hdl"]["adverse_direction"] == "low" and tm["analytes"]["vitd"]["adverse_direction"] == "low"
    assert tm["score_rules"]["band_high_pct"] == 80 and tm["score_rules"]["band_low_pct"] == 20


def test_mthfr_677tt_feeds_homocysteine():
    pgx = {"nutrigenetic": [{"key": "mthfr", "copies": 2, "caveat_zh": "ACMG 提醒"}]}
    _, it = _one("hcy", _res("high", 50), pgx)
    assert it["genetic_band"] == "elevated" and it["class"] == "concordant_adverse" and it["band_basis"] == "mthfr_677tt"
    assert it["acmg_note_zh"] == "ACMG 提醒" and "work/pgx/pgx.json#nutrigenetic.mthfr" in it["basis"]
    _, het = _one("hcy", _res("high", 50), {"nutrigenetic": [{"key": "mthfr", "copies": 1, "caveat_zh": "x"}]})
    assert het["class"] == "acquired_adverse" and "band_basis" not in het


def test_write_concordance_and_readouts(tmp_path):
    ws = tmp_path / "ws"
    C.write_json(ws / "work" / "insights" / "genotype_phenotype.json",
                 {"generated_at": "2026-10-24T10:00:00", "analytes": {"ldl": _res("in_range", 90), "urate": _res("high", None)}})
    doc, rows = genomics.write_concordance(ws)
    assert json.loads((ws / "work" / "insights" / "concordance.json").read_text(encoding="utf-8")) == doc
    assert [(r["id"], r["value"], r["kind"]) for r in rows] == [("gen.ldl.concordance", "latent", "genetic_concordance")]


def test_fingerprint(tmp_path):
    fp = genomics.write_fingerprint(tmp_path, {"rs2": "CT", "rs1": "AA"})
    assert fp["schema"] == "la-fingerprint/1" and fp["sites"] == 2 and list(fp["calls"]) == ["rs1", "rs2"]
    assert fp["sha256"] == hashlib.sha256(b"rs1:AA\nrs2:CT").hexdigest()
    assert json.loads((tmp_path / "work" / "insights" / "fingerprint.json").read_text()) == fp


def test_explain_writes_fingerprint_of_real_calls_only(tmp_path, monkeypatch):
    from test_insights import _explain_setup
    ws, st = _explain_setup(tmp_path, monkeypatch, ["chr1\t1000000\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30",
                                                     "chr1\t2000000\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t1/1:5:30"])
    genomics.explain(st, ws, absent_as_ref="交付说明：WGS 30x 联合分型，列出全部非参考位点")
    fp = C.load_json(ws / "work" / "insights" / "fingerprint.json")
    assert fp["calls"] == {"rs0": "AG"}                     # low GQ and assumed-reference sites are not fingerprint calls
    doc, rows = genomics.write_concordance(ws)
    assert doc["items"][0]["analyte"] == "ldl" and doc["items"][0]["lab_state"] == "adverse"


def test_report_section_uses_fixed_wordings_and_refuses_edits(tmp_path):
    from lalib import report_genomics
    ws = tmp_path / "ws"
    C.write_json(ws / "work" / "insights" / "genotype_phenotype.json",
                 {"generated_at": "2026-10-24", "analytes": {"ldl": _res("in_range", 82), "hdl": _res("low", 50)}})
    genomics.write_concordance(ws)
    sha = C.sha256_file(ws / "work" / "insights" / "concordance.json", limit=None)
    st = {"member": {"answers": {"genetic_disclosure": "yes"}}, "insights": {"files": {"concordance.json": sha}}}
    (heading, md), = report_genomics.sections(ws, {"st": st})
    assert heading == "基因与化验的一致性"
    assert LEAD + genomics.CLASS_ZH["latent"] in md and "接近分档边界" in md and genomics.CLASS_ZH["acquired_adverse"] in md
    assert "诊断" in md and "注定" not in md and "一定会" not in md
    assert report_genomics.sections(ws, {"st": {"member": {"answers": {"genetic_disclosure": "no"}}, "insights": st["insights"]}}) == []
    (ws / "work" / "insights" / "concordance.json").write_text("{}")
    with pytest.raises(C.LAError, match="not written by la.py"):
        report_genomics.sections(ws, {"st": st})
    assert report_genomics.sections(ws, {"st": dict(st, insights={})}) == []          # unregistered: not rendered
