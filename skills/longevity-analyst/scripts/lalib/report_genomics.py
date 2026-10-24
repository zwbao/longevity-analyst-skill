"""Report sections of the genomics lane: genes vs labs concordance (WS-D.1), the surveillance plan (WS-E) and
pharmacogenomics (WS-F). report.py calls `sections(ws_root, ctx)` after the body and before the boundary text.

Wordings are fixed (genomics.CLASS_ZH / LEAD_ZH / OVERRIDE_ZH); CPIC text is printed as stored, never rewritten, and
every CPIC line carries the CPIC label and the guideline link. Files must be the ones `la.py` wrote (hash in state).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import EXIT_BLOCKED, LAError, load_json, sha256_file

CPIC_TAG = "【CPIC】"
EXECUTOR_ZH = {"physician": "医生", "member": "你自己", "nutritionist": "营养师"}
REASON_ZH = {"genetic": "基因", "guideline": "指南", "intervention": "干预复测"}


def _cell(x: Any) -> str:
    return str("" if x is None else x).replace("|", "／").replace("\n", " ").strip()


def _load(ws: Path, st: Dict[str, Any], rel: str, sha: Optional[str]) -> Optional[Dict[str, Any]]:
    p = ws / rel
    if not p.exists() or sha is None:                   # not registered (or voided upstream): nothing is rendered
        return None
    if sha256_file(p, limit=None) != sha:
        raise LAError(f"{rel} was not written by la.py or changed afterwards; run that step again", EXIT_BLOCKED)
    return load_json(p)


def concordance_section(doc: Dict[str, Any]) -> List[str]:
    from .genomics import CLASS_ZH
    rows = [x for x in doc.get("items", []) if x["class"] in CLASS_ZH or x["class"] == "clinvar_override"]
    if not rows:
        return []
    r = doc.get("rules") or {}
    L = [f"> 遗传倾向按「偏向不利方向的程度」分档：东亚人群百分位达到第 {r.get('band_high_pct')} 或第 {r.get('band_low_pct')} 百分位"
         "（看哪一侧不利）算「倾向明显」。这个分档是约定，不是生物学界限。下面每一条都只是倾向，不是诊断。", "",
         "| 指标 | 化验 | 怎么看 |", "|---|---|---|"]
    lab_zh = {"adverse": "超出范围（不利的一侧）", "in_range": "在范围内"}
    for x in rows:
        text = x["wording_zh"]
        if x["class"] != "clinvar_override" and x.get("near_band_edge"):
            text += "（遗传分位接近分档边界，这个分类不太稳定。）"
        if x.get("band_basis") == "mthfr_677tt":
            text += "（这里的先天倾向来自 MTHFR 677TT。）"
        if x.get("acmg_note_zh"):
            text += x["acmg_note_zh"]
        L.append(f"| {_cell(x.get('label_zh') or x['analyte'])} | {lab_zh.get(x['lab_state'], '—')} | {_cell(text)} |")
    nc = [x for x in doc.get("items", []) if x["class"] == "not_computed"]
    if nc:
        L += ["", "没有算出遗传倾向、所以没有分类的指标：" + "、".join(_cell(x.get("label_zh") or x["analyte"]) for x in nc) + "。"]
    return L


def surveillance_section(doc: Dict[str, Any]) -> List[str]:
    items = doc.get("items") or []
    if not items:
        return []
    L = ["| 测什么 | 多久一次 | 为什么 | 谁来安排 |", "|---|---|---|---|"]
    for x in items:
        if x["once"]:
            mx = x.get("measured") or {}
            seen = (f"你在 {mx['t']} 测过（{_cell(mx.get('value'))} {_cell(mx.get('unit') or '')}"
                    + (f"，参考 {_cell(mx['ref_range'])}" if mx.get("ref_range") else "") + "）") if mx.get("t") else "你已经测过"
            when = f"成人至少测一次；{seen}，不放进复测计划" if x["done"] else "成人至少测一次"
        else:
            when = f"{int(x['interval_weeks'])} 周后" if x.get("interval_weeks") else "—"
        why = "、".join(REASON_ZH.get(r, r) for r in x.get("reasons") or [x["reason"]])
        if x.get("note_zh"):
            why += f"：{x['note_zh']}"
        L.append(f"| {_cell(x['marker'])} | {when} | {_cell(why)} | {EXECUTOR_ZH.get(x.get('executor'), _cell(x.get('executor')))} |")
    srcs = [x for x in items if x.get("source")]
    if srcs:
        L += ["", "指南依据（原文）："]
        for x in srcs:
            s = x["source"]
            loc = f"，{s['locator']}" if s.get("locator") else ""
            doi = f"，doi:{s['doi']}" if s.get("doi") else ""
            L.append(f"- {_cell(x['marker'])}：“{_cell(s['quote'])}”——{_cell(s['title'])}{doi}{loc}（{s['url']}）")
    return L


def pgx_section(doc: Dict[str, Any], hcy_measured: bool = False) -> List[str]:
    if doc.get("status") != "ran":
        return [_cell(doc.get("reason_zh") or "用药基因没有分析。")]
    L = ["> 下面是按 CPIC（临床药物基因组学实施联盟）指南对照的结果。用药、换药、调剂量都要由开药的医生或药师决定，"
         "请把这一节给他们看；不要自己停药或改药。", ""]
    if doc.get("reuse_note_zh"):
        L += [f"> {_cell(doc['reuse_note_zh'])}", ""]
    if doc.get("cpic_text") != "included":
        L += [f"> {_cell(doc.get('cpic_text_gap_zh'))}", ""]
    if not doc.get("results"):
        L.append("这次的基因数据里没有能用来对照 CPIC 的位点结果。")
    for r in doc.get("results") or []:
        dip = f"，推断型别 {_cell(r['diplotype_or_genotype'])}" if r.get("diplotype_or_genotype") not in (None, "—") else ""
        L.append(f"**{_cell(r['gene'])}**：{_cell(r.get('phenotype_zh') or r['phenotype'])}（CPIC 用语：{_cell(r['phenotype'])}{dip}）")
        if not r["drugs"]:
            L += ["", "- 没有需要先和医生核对的药。", ""]
            continue
        L.append("")
        for d in r["drugs"]:
            mine = "（你正在用）" if d.get("member_takes") else ""
            link = f"指南：{d['guideline_url']}" + (f"（doi:{d['guideline_doi']}）" if d.get("guideline_doi") else "")
            if d.get("recommendation_text"):
                strength = f"，推荐强度 {_cell(d['strength'])}" if d.get("strength") else ""
                ut = "、".join(d.get("untested_genes") or [])
                miss = f"{ut} 没有检测结果，这条推荐没有考虑它；用这个药之前，请医生确认{ut}。" if ut else ""
                name = f"{_cell(d['drug_zh'])}（{_cell(d['drug'])}）" if d.get("drug_zh") else _cell(d["drug"])
                if len(d.get("lookup") or {}) > 1:
                    name += "（这条推荐同时按 " + "、".join(f"{_cell(g)} {_cell(v)}" for g, v in d["lookup"].items()) + " 给出）"
                L.append(f"- {CPIC_TAG}{name}{mine}：{_cell(d['recommendation_text'])}（CPIC 原文{strength}）。{miss}{link}")
            else:
                head = (d.get("skill_text_zh") or [f"- {d['drug']}"])[0].lstrip("- ").strip()
                rest = [l.strip().lstrip("- ").strip() for l in (d.get("skill_text_zh") or [])[1:] if not l.strip().lstrip("- ").startswith("指南")]
                body = "；".join(_cell(t).rstrip("。") for t in [head + mine, *rest])
                L.append(f"- {CPIC_TAG}{body}。{link}")
        L.append("")
    ng = [x for x in doc.get("nutrigenetic") or [] if x.get("copies") is not None]
    if ng:
        L += ["**饮食相关的基因位点**（只描述，不是诊断，也不是营养处方）", ""]
        for x in ng:
            act = x.get("action_zh")
            if x.get("key") == "mthfr" and hcy_measured and act:     # the panel's generic advice, when the member has the value
                act = "你已经测过血同型半胱氨酸，结果见「基因与化验的一致性」和化验对比。"
            extra = "".join(f" {_cell(t)}" for t in (act, x.get("caveat_zh")) if t)
            L.append(f"- {_cell(x['gene'])}：{_cell(x.get('reading_zh'))}{extra}（依据：{_cell(x.get('cite'))}，doi:{_cell(x.get('doi'))}）")
        L.append("")
    return L


def sections(ws_root: Path, ctx: Optional[Dict[str, Any]] = None) -> List[Tuple[str, str]]:
    ws = Path(ws_root)
    sp = ws / "state.json"
    st = (ctx or {}).get("st") or (load_json(sp) if sp.exists() else {})
    consent = ((st.get("member") or {}).get("answers") or {}).get("genetic_disclosure") == "yes"
    files = (st.get("insights") or {}).get("files") or {}
    out: List[Tuple[str, str]] = []
    conc = _load(ws, st, "work/insights/concordance.json", files.get("concordance.json")) if consent else None
    if conc:
        L = concordance_section(conc)
        if L:
            out.append(("基因与化验的一致性", "\n".join(L)))
    sv = st.get("surveillance")
    plan = _load(ws, st, "work/surveillance/plan.json", (sv or {}).get("sha256"))
    pp = ws / "work" / "intervene" / "plan.json"
    if plan and (sha256_file(pp, limit=None) if pp.exists() else None) != (sv or {}).get("intervene_plan_sha256", "-"):
        raise LAError("the intervention plan changed after `la.py surveillance plan`; run it again", EXIT_BLOCKED)
    if plan:
        if not consent:                                 # a rule triggered by a genetic result stays out without consent
            plan = dict(plan, items=[x for x in plan["items"] if not (x["reason"] == "genetic" and x["rule_ids"])])
        L = surveillance_section(plan)
        if L:
            out.append(("监测计划", "\n".join(L)))
    pg = _load(ws, st, "work/pgx/pgx.json", files.get("pgx.json"))
    if pg and (consent or pg.get("status") != "ran"):
        cp = ws / "work" / "insights" / "concordance.json"
        hcy = any(i.get("analyte") == "hcy" and i.get("lab_flag") in ("in_range", "high", "low")
                  for i in ((load_json(cp).get("items") or []) if cp.exists() else []))
        out.append(("用药基因", "\n".join(pgx_section(pg, hcy_measured=hcy)).rstrip()))
    return out
