"""Pharmacogenomics and nutrigenetics from the member's VCF (WS-F).

The library skills `cpic-pharmacogenomics` and `nutrigenetic-variant-panel` read consumer raw-data files (23andMe
format, rsID matched, plus strand). This adapter reads the member's genotype at exactly the sites those skills list,
with the same quality rules as genomics.Vcf (an absent site is never reference unless a gVCF block covers it or the
recorded --absent-as-ref judgment says so), writes `work/pgx/genotype.txt` in that format, and runs the two skills
through the method runner. CPIC text is stored as the skill prints it, never rewritten. The genotype file stays in
the workspace: it never goes to deliver/ or a host.
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import EXIT_INPUT, LAError, data, load_json, longevity_skills_home, now_iso, write_json

SKILLS = ("cpic-pharmacogenomics", "nutrigenetic-variant-panel")
CPIC_TEXT_FLAG = "--include-cpic-text"
# Gene -> (chromosome, start, end) on GRCh38, from Ensembl REST POST /lookup/symbol/homo_sapiens, retrieved 2026-10-24
# (sources/genomics/ensembl-gene-lookup/lookup.json). Every site position is checked against its gene's span (+/- GENE_PAD).
GENE_CHROM = {
    "CYP2C19": ("10", 94762536, 94922264), "CYP2C9": ("10", 94938544, 94990149), "TPMT": ("6", 18128311, 18155348),
    "NUDT15": ("13", 48037514, 48054196), "DPYD": ("1", 97077743, 97995000), "CYP3A5": ("7", 99647442, 99680043),
    "SLCO1B1": ("12", 21127862, 21298065), "ABCG2": ("4", 88088127, 88231818), "VKORC1": ("16", 31090260, 31095980),
    "CYP4F2": ("19", 15877893, 15900528), "ALDH2": ("12", 111766887, 111817532), "MTHFR": ("1", 11785723, 11806455),
    "LCT": ("2", 135787850, 135837184), "MCM6": ("2", 135839626, 135876484),
}
HGVS = re.compile(r"g\.(\d+)([ACGT])>([ACGT])")
PANEL_POS = re.compile(r"^\s*(\w+):(\d+)\s+([ACGT])/([ACGT])")
NOT_RUN_ZH = {
    "consent": "会员没有同意在报告中查看遗传结果（genetic_disclosure 不是 yes），用药基因没有分析。",
    "no_vcf": "没有可用的基因组数据（VCF），用药基因没有分析。",
    "grch37": "这份基因组数据是 GRCh37 坐标，本版本不做坐标转换，用药基因没有分析。",
}


def _presets(home: Path):
    p = home / "skills" / "nutrigenetic-variant-panel" / "scripts" / "presets.py"
    spec = importlib.util.spec_from_file_location("la_nutri_presets", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)                     # plain data module: constants only
    return mod


def _in_gene(gene: str, chrom: str, pos: int) -> bool:
    """Inside the gene span plus the promoter/UTR margin genomics uses (CYP2C19*17 and VKORC1 -1639 are upstream)."""
    from .genomics import GENE_PAD
    g = GENE_CHROM.get(gene)
    return bool(g) and g[0] == chrom and g[1] - GENE_PAD <= pos <= g[2] + GENE_PAD


def sites(home: Optional[Path] = None) -> List[Dict[str, Any]]:
    """The sites the two skills read: CPIC genes[].sites and warfarin_sites (GRCh38 HGVS) plus the panel's VARIANTS."""
    home = home or longevity_skills_home()
    cp = load_json(home / "skills" / "cpic-pharmacogenomics" / "data" / "cpic_tables.json")
    out: List[Dict[str, Any]] = []
    pairs = [(g, s) for g, x in cp["genes"].items() for s in x["sites"]] + list(cp.get("warfarin_sites", {}).items())
    for gene, s in pairs:
        hits = HGVS.findall(s["grch38"])
        if not hits or gene not in GENE_CHROM:
            raise LAError(f"cpic site {s['rsid']} ({gene}): no GRCh38 position or chromosome", EXIT_INPUT)
        pos, ref = int(hits[0][0]), hits[0][1]
        chrom = GENE_CHROM[gene][0]
        if not _in_gene(gene, chrom, pos) or any(int(h[0]) != pos or h[1] != ref for h in hits):
            raise LAError(f"cpic site {s['rsid']} at {pos} is not inside {gene}", EXIT_INPUT)
        out.append({"rsid": s["rsid"], "gene": gene, "chrom": chrom, "pos": pos, "ref": ref,
                    "alts": sorted({h[2] for h in hits}), "skill": "cpic-pharmacogenomics"})
    for key, v in _presets(home).VARIANTS.items():
        m = PANEL_POS.match(v["grch38"])
        gene = v["gene"].split("（")[0]
        if not m or m.group(3) != v["reference"] or m.group(4) != v["effect"]:
            raise LAError(f"nutrigenetic site {v['rsid']}: position not readable", EXIT_INPUT)
        chrom, pos = m.group(1), int(m.group(2))
        if not (_in_gene(gene, chrom, pos) or (gene == "LCT" and GENE_CHROM["MCM6"][0] == chrom)):
            raise LAError(f"nutrigenetic site {v['rsid']} at {chrom}:{pos} is not on {gene}'s chromosome span", EXIT_INPUT)
        out.append({"rsid": v["rsid"], "gene": gene, "chrom": chrom, "pos": pos, "ref": v["reference"], "alts": [v["effect"]],
                    "skill": "nutrigenetic-variant-panel", "key": key})
    return out


def read_sites(vcf, site_list: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """rsid -> {genotype: sorted plus-strand bases or None, status}. Uses genomics.Vcf rules: FILTER/GQ/DP, gVCF blocks,
    deletions spanning the site, absent = not called unless the recorded judgment says reference."""
    from .genomics import _norm_chrom
    vcf.prefetch({(s["chrom"], s["pos"]) for s in site_list}, [])
    out: Dict[str, Dict[str, Any]] = {}
    for s in site_list:
        r = vcf.ref_dosage(s["chrom"], s["pos"], s["ref"])
        if r["dosage"] is None:
            out[s["rsid"]] = {"genotype": None, "status": r["status"]}
            continue
        if r["dosage"] == 2:
            out[s["rsid"]] = {"genotype": s["ref"] * 2, "status": r["status"]}
            continue
        rows, _, _ = vcf._site(_norm_chrom(s["chrom"]), s["pos"])
        alts = sorted({a.upper() for f in rows if f[3].upper() == s["ref"] for a in f[4].split(",")})
        bases, status = s["ref"] * r["dosage"], r["status"]
        for a in alts:
            d = vcf.alt_dosage(s["chrom"], s["pos"], s["ref"], a)
            if d["dosage"]:
                if len(a) != 1 or a not in "ACGT":
                    status = "indel_at_site"
                bases += a * d["dosage"]
        if len(bases) != 2 or status == "indel_at_site":
            out[s["rsid"]] = {"genotype": None, "status": status if status == "indel_at_site" else "unresolved"}
        else:
            out[s["rsid"]] = {"genotype": "".join(sorted(bases)), "status": status}
    return out


def write_23andme(path: Path, site_list: List[Dict[str, Any]], calls: Dict[str, Dict[str, Any]]) -> None:
    lines = ["# Generated by longevity-analyst from the member's VCF (GRCh38, plus strand) at the sites of "
             "cpic-pharmacogenomics and nutrigenetic-variant-panel. -- = no usable call.",
             "# rsid\tchromosome\tposition\tgenotype"]
    seen = set()
    for s in site_list:
        if s["rsid"] in seen:
            continue
        seen.add(s["rsid"])
        lines.append(f"{s['rsid']}\t{s['chrom']}\t{s['pos']}\t{calls[s['rsid']]['genotype'] or '--'}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def has_cpic_text_flag(home: Path) -> bool:
    """The tooling lane adds --include-cpic-text to the skill; detect it in the manifest or the entry script."""
    sd = home / "skills" / "cpic-pharmacogenomics"
    man = load_json(sd / "skill.json")
    if any(i.get("flag") == CPIC_TEXT_FLAG for i in man.get("inputs", [])) or CPIC_TEXT_FLAG in str(man.get("entry", {})):
        return True
    script = sd / man["entry"]["script"]
    return script.exists() and CPIC_TEXT_FLAG in script.read_text(encoding="utf-8")


MED_SPLIT = re.compile(r"[,，、;；\n]+")


def medications(st: Dict[str, Any]) -> List[str]:
    a = (st["member"].get("answers") or {})
    raw = a.get("medications") or a.get("current_medications") or ""
    return [m.strip() for m in MED_SPLIT.split(str(raw)) if m.strip() and m.strip() not in ("无", "没有", "none", "no")]


def _not_run(ws: Path, why: str, assembly: Optional[str] = None) -> Dict[str, Any]:
    doc = {"schema": "la-pgx/1", "status": "not_run", "reason": why, "reason_zh": NOT_RUN_ZH[why], "assembly": assembly,
           "generated_at": now_iso(), "results": [], "nutrigenetic": []}
    write_json(ws / "work" / "pgx" / "pgx.json", doc)
    return doc


def _drug_blocks(report_md: str) -> Dict[str, List[List[str]]]:
    """drug (English generic name) -> each block of lines the skill printed for it, verbatim (the current-medication
    section repeats a drug with 'name →' in front)."""
    out: Dict[str, List[List[str]]] = {}
    cur = None
    for line in report_md.splitlines():
        m = re.match(r"^- (?:.*→ )?\*\*[^*]+\*\*（([a-z][a-z \-]+)）", line)
        if m:
            cur = [line]
            out.setdefault(m.group(1), []).append(cur)
        elif cur is not None and line.startswith("  - "):
            cur.append(line)
        else:
            cur = None
    return out


def _gene_rows(report_md: str) -> Dict[str, tuple]:
    """gene -> (inferred diplotype, phenotype in the skill's Chinese wording) from its gene table."""
    out = {}
    for line in report_md.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 4 and re.match(r"^[A-Z0-9-]+（", cells[0]):
            out[cells[0].split("（")[0]] = (cells[2], cells[3])
    return out


def _val(v: Any) -> Any:
    return v.get("value") if isinstance(v, dict) else v


def _cpic_text_rows(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    for k in ("cpic_recommendations", "cpic_text"):
        v = result.get(k) if isinstance(result.get(k), list) else (result.get("outputs") or {}).get(k)
        if isinstance(v, list):
            return [r for r in v if isinstance(r, dict)]
    return []


def run(st: Dict[str, Any], ws: Path, absent_as_ref: Optional[str] = None) -> Dict[str, Any]:
    from . import genomics, methods
    if (st["member"].get("answers") or {}).get("genetic_disclosure") != "yes":
        return _not_run(ws, "consent")
    variants = methods.pick(st, "variants")
    if not variants and genomics.reused_germline(ws, "pgx.json"):
        doc = load_json(genomics.reused_germline(ws, "pgx.json"))      # germline does not change (WS-I.1)
        doc.update(reused_from=load_json(ws / "work" / "reuse" / "germline.json").get("from"), reused_at=now_iso(),
                   reuse_note_zh="用药基因沿用上一次的结果；当时登记的在用药物可能与现在不同，请把现在用的药告诉医生。")
        write_json(ws / "work" / "pgx" / "pgx.json", doc)
        return doc
    if not variants:
        return _not_run(ws, "no_vcf")
    asm = genomics._assembly(st)
    if asm != "hg38":
        return _not_run(ws, "grch37", "GRCh37")
    home = longevity_skills_home()
    sl = sites(home)
    absent = absent_as_ref or ((st.get("insights") or {}).get("absent_as_ref") or {}).get("reason")
    vcf = genomics.Vcf(Path(variants["path"]), data("trait_map.json")["score_rules"], absent_is_ref=bool(absent))
    calls = read_sites(vcf, sl)
    gfile = ws / "work" / "pgx" / "genotype.txt"
    write_23andme(gfile, sl, calls)
    meds = medications(st)
    mfile = ws / "work" / "pgx" / "medications.txt"
    mfile.write_text("\n".join(meds) + ("\n" if meds else ""), encoding="utf-8")
    text_flag = has_cpic_text_flag(home)
    ran: Dict[str, Dict[str, Any]] = {}
    for name in SKILLS:
        args: Dict[str, Any] = {"--genotype": str(gfile)}
        if name == "cpic-pharmacogenomics":
            if meds:
                args["--medications"] = str(mfile)
            if text_flag:
                args[CPIC_TEXT_FLAG] = True
        it = {"method": name, "run": {"kind": "skill", "rows": [], "args": args}}
        res = methods._run_skill(home, it, st["member"], ws / "work" / "pgx" / name)
        rj = ws / "work" / "pgx" / name / "result.json"
        md = ws / "work" / "pgx" / name / "report.md"
        ran[name] = {"exit": res["exit"], "result": load_json(rj) if rj.exists() else {},
                     "report_md": md.read_text(encoding="utf-8") if md.exists() else "", "readouts": res.get("readouts") or []}
    doc = {"schema": "la-pgx/1", "status": "ran", "reason_zh": None, "assembly": "GRCh38", "generated_at": now_iso(),
           "genotype_file": "work/pgx/genotype.txt", "absent_as_ref": absent,
           "sites": {"listed": len(calls), "called": sum(1 for c in calls.values() if c["genotype"]),
                     "not_called": sorted(rs for rs, c in calls.items() if not c["genotype"])},
           "medications": meds, "cpic_text": "included" if text_flag else "gap", "results": [], "nutrigenetic": []}
    if not text_flag:
        doc["cpic_text_gap_zh"] = ("方法库的 CPIC 技能还不能输出 CPIC 推荐原文；这里只保存技能打印的中文对照（类别和推荐强度）和指南链接，"
                                   "原文请看指南。")
    cp = load_json(home / "skills" / "cpic-pharmacogenomics" / "data" / "cpic_tables.json")
    c = ran["cpic-pharmacogenomics"]
    outs = c["result"].get("outputs", c["result"]) if isinstance(c["result"].get("outputs", c["result"]), dict) else {}
    blocks = _drug_blocks(c["report_md"])
    gene_rows = _gene_rows(c["report_md"])
    discuss = {d.strip() for d in str(_val(outs.get("drugs_to_discuss")) or "").split(",") if d.strip()}
    med_drugs = {d for d, bl in blocks.items() if any("→" in b[0] for b in bl)}
    text_rows = _cpic_text_rows(c["result"])
    for gene in [k[:-len("_phenotype")].upper() for k in outs if k.endswith("_phenotype")]:
        ph = _val(outs.get(gene.lower() + "_phenotype"))
        if not ph:
            continue
        drugs = []
        for drug, e in sorted(cp["drugs"].items()):
            if gene not in e.get("genes", []) or not (drug in discuss or drug in med_drugs):
                continue
            rows = [r for r in text_rows if r.get("drug") == drug]
            drugs.append({"drug": drug, "drug_zh": next((r.get("drug_zh") for r in rows if r.get("drug_zh")), None),
                          "category": next((r.get("category") for r in rows if r.get("category")), None),
                          "recommendation_text": "\n".join(str(r.get("recommendation") or r.get("recommendation_text") or "") for r in rows) or None,
                          "dose_text": next((r.get("dose_text") for r in rows if r.get("dose_text")), None),
                          "strength": next((r.get("classification") or r.get("strength") for r in rows if r.get("classification") or r.get("strength")), None),
                          "guideline": e.get("guideline"), "guideline_url": e.get("guideline_url"), "guideline_doi": e.get("guideline_doi"),
                          "member_takes": drug in med_drugs,
                          "skill_text_zh": next((b for b in blocks.get(drug, []) if "→" not in b[0]), (blocks.get(drug) or [[]])[0]),
                          "cpic_source": {"api": cp["source"].get("api"), "release": cp["source"].get("cpic_data_release"),
                                          "accessed": cp["source"].get("accessed"), "license": cp["source"].get("license"),
                                          "row_ids": [r.get("recommendation_id", r.get("id")) for r in rows
                                                      if r.get("recommendation_id", r.get("id")) is not None],
                                          "versions": sorted({r["version"] for r in rows if r.get("version") is not None})}})
        dip, ph_zh = gene_rows.get(gene, (None, None))
        doc["results"].append({"gene": gene, "diplotype_or_genotype": dip, "phenotype": ph, "phenotype_zh": ph_zh, "drugs": drugs})
    if c["exit"] != 0:
        doc["cpic_error"] = (c["report_md"] or "")[:500]
    n = ran["nutrigenetic-variant-panel"]
    nouts = n["result"].get("outputs", n["result"]) if isinstance(n["result"].get("outputs", n["result"]), dict) else {}
    pres = _presets(home)
    for key, v in pres.VARIANTS.items():
        copies = _val(nouts.get(f"{key}_{v['rsid']}_copies"))
        doc["nutrigenetic"].append({"key": key, "gene": v["gene"], "effect_name": v.get("effect_name"),
                                    "copies": copies, "kind": "descriptive",
                                    "reading_zh": v["reading"].get(copies) if copies is not None else None,
                                    "action_zh": v["action"].get(copies) if copies is not None else None,
                                    "caveat_zh": v.get("caveat") if copies else None, "doi": v.get("doi"), "cite": v.get("cite")})
    write_json(ws / "work" / "pgx" / "pgx.json", doc)
    return doc


def readouts(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    rows = [{"id": f"pgx.{r['gene'].lower()}.phenotype", "label_zh": f"{r['gene']} 表型（CPIC 用语）", "value": r["phenotype"], "unit": "",
             "kind": "pgx", "method": "cpic-pharmacogenomics", "source": "work/pgx/pgx.json", "systems": []} for r in doc.get("results", [])]
    rows += [{"id": f"nutrigenetic.{x['key']}.copies", "label_zh": f"{x['gene']} {x['effect_name']} 的份数", "value": x["copies"],
              "unit": "份", "kind": "descriptive", "method": "nutrigenetic-variant-panel", "source": "work/pgx/pgx.json", "systems": []}
             for x in doc.get("nutrigenetic", []) if x.get("copies") is not None]
    return rows
