"""Report sections for labelled items: 个人反应预测, 诊断性结论, 用药建议.

Each section prints the category tag and the sources of every row; AI judgements are in their own table, like the
organ table's AI column. Called by report.py's section hook: sections(ws_root, ctx) -> [(heading, markdown)].
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

from . import labels as L
from .common import load_json

EMERGENCY = "出现急症（胸痛、呼吸困难、一侧肢体无力、意识不清等）时，先拨打 120，不要等报告或复诊。"
REVIEW = "下面的内容由医生复核后才执行；是否用药、怎么用，由你的医生决定。"
ACTION_ZH = {"consider_start": "考虑开始", "consider_stop": "考虑停用", "consider_switch": "考虑换药", "adjust": "调整剂量",
             "avoid": "避免使用", "monitor": "监测"}
SET_ZH = {"cn": "中国指南", "ada": "美国糖尿病协会（ADA）标准"}


def _cell(x: Any) -> str:
    from .report import _cell as c
    return c(x)


def _p(x: Any) -> str:
    return L._p(x)


def _load(ws: Path, kind: str) -> Dict[str, Any]:
    p = L.item_path(ws, kind)
    return load_json(p) if p.exists() else {}


def _src(lab: Dict[str, Any]) -> str:
    return _cell(L.render_sources(lab)[3:])


def predictions(ws: Path) -> List[str]:
    d = _load(ws, "p")
    items = d.get("items") or []
    if not items:
        return []
    stat = [x for x in items if x["label"]["basis_kind"] == "statistical_model"]
    ai = [x for x in items if x["label"]["basis_kind"] == "ai_estimate"]
    out = ["> 预测的是「这件事对你有没有用、下次复测能不能看出来」。「统计模型估计」用试验平均效应、人和人之间的差异、"
           "你本人做过的实验和你的测量波动算出；「AI 判断」是没有可用试验数据时由模型参照文献给出的区间，不是计算结果。", ""]
    if stat:
        out += ["| 方案项 | 做什么 | 指标 | 下次能测出有利变化 | 真实效应有利 | 预计变化（80% 区间） | 类别 | 来源 |",
                "|---|---|---|---|---|---|---|---|"]
        for x in stat:
            lo, hi = x["interval_80_pct"]
            out.append(f"| {_cell(x.get('plan_item'))} | {_cell(x['intervention_zh'])} | {_cell(x['marker'])} | "
                       f"约 {_p(x['p_detect_favorable'])}（前后各测 {x['design']['repeats_each_side']} 次） | 约 {_p(x['p_true_favorable'])} | "
                       f"{lo:.1f}% 到 {hi:.1f}% | {L.render_label(x['label'])} | {_src(x['label'])} |")
        out.append("")
        assum = []
        for x in stat:
            assum += [a for a in x["label"].get("assumptions", []) if a not in assum]
        out += ["**模型假设**", ""] + [f"- {_cell(a)}" for a in assum] + [""]
    if ai:
        out += ["**AI 判断（单独列出，不是计算结果）**", "", "| 方案项 | 做什么 | 指标 | 能测出有利变化（区间） | 真实效应有利（区间） | 置信度 | 依据 | 来源 |",
                "|---|---|---|---|---|---|---|---|"]
        for x in ai:
            out.append(f"| {_cell(x.get('plan_item'))} | {_cell(x['intervention_zh'])} | {_cell(x['marker'])} | "
                       f"{_p(x['p_detect_favorable']) if x.get('p_detect_favorable') else '—'} | {_p(x['p_true_favorable']) if x.get('p_true_favorable') else '—'} | "
                       f"{L.LEVEL_ZH[x['label']['confidence']['level']]} | {_cell('；'.join(x['label'].get('assumptions', [])))} | {_src(x['label'])} |")
        out.append("")
    return out


def diagnoses(ws: Path) -> List[str]:
    d = _load(ws, "dx")
    if not d:
        return []
    items = d.get("items") or []
    gl = [x for x in items if x["kind"] == "guideline_criteria"]
    ai = [x for x in items if x["kind"] == "ai_estimate"]
    used = "、".join(f"《{g}》" for g in d.get("guidelines_used_zh", []))
    out = [f"> {EMERGENCY}", "", f"所用指南集：{SET_ZH.get(d.get('guideline_set'), d.get('guideline_set'))}（{_cell(used)}）。{REVIEW}", ""]
    if gl:
        out += ["| 结论 | 你的数值 | 类别 | 指南原文 |", "|---|---|---|---|"]
        for x in gl:
            vals = "；".join(f"{_cell(e['ref'])} {e['value']:g} {_cell(e.get('unit', ''))}（{e.get('t', '')}）" if isinstance(e.get("value"), (int, float))
                            else f"{_cell(e['ref'])} {_cell(e.get('value'))}" for e in x.get("evidence", []))
            g = next(s for s in x["label"]["sources"] if s["type"] == "guideline")
            note = f"（{_cell(x['note_zh'])}）" if x.get("note_zh") else ""
            out.append(f"| {_cell(x['conclusion_zh'])}{note} | {vals or '—'} | {L.render_label(x['label'])} | {_cell(g.get('quote'))}（{_cell(g.get('ref'))}） |")
        out.append("")
    else:
        out += ["按所用指南的标准，本次数据没有符合任何一条诊断标准。", ""]
    nm = [e for e in d.get("evaluated", []) if e["result"] == "not_met"]
    if nm:
        out += ["已核对、未达到标准：" + "；".join(f"{_cell(e['label_zh'])}（《{_cell(e['guideline_zh'])}》）" for e in nm) + "。", ""]
    if ai:
        out += ["**AI 判断（单独列出，不是指南标准）**", "", "| AI 判断 | 概率区间 | 置信度 | 来源 |", "|---|---|---|---|"]
        for x in ai:
            out.append(f"| {_cell(x['conclusion_zh'])} | {_p(x['probability_range'])} | {L.LEVEL_ZH[x['label']['confidence']['level']]} | {_src(x['label'])} |")
        out.append("")
    return out


def medications(ws: Path) -> List[str]:
    d = _load(ws, "rx")
    if not d:
        return []
    items = d.get("items") or []
    main = [x for x in items if x["kind"] in ("guideline_recommendation", "cpic")]
    ai = [x for x in items if x["kind"] == "ai_estimate"]
    out = [f"> {EMERGENCY}", "", REVIEW, ""]
    if main:
        out += ["| 药物 | 建议 | 原文推荐 | 剂量（原文） | 类别 | 来源 |", "|---|---|---|---|---|---|"]
        for x in main:
            taking = "（你正在用）" if x.get("currently_taking") else ""
            who = f"{_cell(x.get('gene'))} {_cell(x.get('phenotype'))}：" if x["kind"] == "cpic" else ""
            out.append(f"| {_cell(x['drug_or_class_zh'])}{taking} | {ACTION_ZH[x['action']]} | {who}{_cell(x['text_zh'])} | "
                       f"{_cell(x.get('dose_text_zh') or '—')} | {L.render_label(x['label'])} | {_src(x['label'])} |")
        out.append("")
    elif not ai:
        out += ["本次没有按指南或用药基因得出的用药建议。" + (f"（用药基因：{_cell(d['cpic_not_used'])}）" if d.get("cpic_not_used") else ""), ""]
    if ai:
        out += ["**AI 判断（单独列出，不是指南或 CPIC 推荐；剂量也是 AI 判断）**", "", "| 药物 | 建议 | 说明 | 剂量（AI 判断） | 置信度 | 来源 |",
                "|---|---|---|---|---|---|"]
        for x in ai:
            taking = "（你正在用）" if x.get("currently_taking") else ""
            out.append(f"| {_cell(x['drug_or_class_zh'])}{taking} | {ACTION_ZH[x['action']]} | {_cell(x['text_zh'])} | {_cell(x.get('dose_text_zh') or '—')} | "
                       f"{L.LEVEL_ZH[x['label']['confidence']['level']]} | {_src(x['label'])} |")
        out.append("")
    return out


def sections(ws_root: Path, ctx: Dict[str, Any]) -> List[Tuple[str, str]]:
    out = []
    for head, fn in (("个人反应预测", predictions), ("诊断性结论", diagnoses), ("用药建议", medications)):
        body = fn(ws_root)
        if body:
            out.append((head, "\n".join(body).strip()))
    return out
