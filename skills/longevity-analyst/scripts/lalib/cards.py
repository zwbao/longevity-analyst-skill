"""deliver/cards.json: the health-page cards (pi-card/1, PLAN A.17) a host shows, built at `report` time.

Every card is made from the export written in the same render (deliver/la-export.json), so report, export and
cards say the same thing. Each number carries the path it was read from (`ref`, e.g. "predictions[0].p_detect_favorable")
and its label; `verify` resolves every ref again. While the export has no section for a kind (la-export/1), the card is
read from the registered handover file instead and its ref says so ("file:work/predict/predictions.json#items[0]...").
Kinds made here: trend, change, prediction, diagnosis, medication, pgx, concordance, retest_due, surveillance
(commitment, milestone and experiment cards come from the host).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import labels as L
from .common import load_json, sha256_file

SCHEMA = "pi-card/1"
CHAT = {"id": "a1", "label_zh": "去聊聊", "type": "open_chat"}
DOCTOR = {"id": "a2", "label_zh": "加入问医生的单子", "type": "add_to_doctor_list"}
TREND_ZH = {"rising_beyond_noise": "缓慢上升", "falling_beyond_noise": "缓慢下降"}
CHANGE_ZH = {"increase_beyond_noise": "升高，超出个人正常波动", "decrease_beyond_noise": "下降，超出个人正常波动"}
CLASS_ZH = {"concordant_adverse": "先天倾向已经表现出来", "latent": "有基因倾向，化验尚在范围内", "acquired_adverse": "化验异常，主要不是基因倾向",
            "concordant_typical": "基因倾向和化验都在常见范围", "protective_discordant": "化验异常，基因倾向反而有利"}
FILES = {"predictions": ("work/predict/predictions.json", "items"), "diagnoses": ("work/clinical/diagnoses.json", "items"),
         "medications": ("work/clinical/medications.json", "items"), "concordance": ("work/insights/concordance.json", "items"),
         "surveillance": ("work/surveillance/plan.json", "items")}


def _slug(s: Any) -> str:
    return re.sub(r"[^A-Za-z0-9_\-]+", "-", str(s)).strip("-")[:40] or "x"


def _section(export: Dict[str, Any], ws: Path, key: str) -> List[tuple]:
    """[(ref prefix, item)] from the export, or from the handover file while the export has no such section."""
    if isinstance(export.get(key), list):
        return [(f"{key}[{i}]", x) for i, x in enumerate(export[key])]
    if key in FILES:
        rel, k = FILES[key]
        p = ws / rel
        if p.exists():
            try:
                items = load_json(p).get(k) or []
            except (OSError, ValueError, AttributeError):
                return []
            return [(f"file:{rel}#{k}[{i}]", x) for i, x in enumerate(items) if isinstance(x, dict)]
    return []


def _num(ref: str, value: Any, unit: str, label: str) -> Dict[str, Any]:
    return {"ref": ref, "value": value, "unit": unit, "label": label}


def _measure(sources: List[Dict[str, Any]]) -> Dict[str, Any]:
    return L.make_label("measurement", sources or [{"type": "lab", "ref": "member"}])


def build(ws: Path, export_path: Path) -> List[Dict[str, Any]]:
    export = load_json(export_path)
    src = {"export": export.get("schema"), "sha256": sha256_file(export_path, limit=None)}
    at = export.get("generated_at")
    cards: List[Dict[str, Any]] = []

    def card(cid: str, kind: str, title: str, summary: str, numbers: List[Dict[str, Any]], label: Optional[Dict[str, Any]],
             actions: Optional[List[Dict[str, Any]]] = None, **extra: Any) -> None:
        cid = f"C-{kind}-{_slug(cid)}"
        cards.append({"schema": SCHEMA, "id": cid, "kind": kind, "title_zh": title, "summary_zh": summary, "numbers": numbers,
                      "label": label, "actions": actions or [CHAT], "deep_link": f"longpi://cards/{cid}", "updated_at": at,
                      "source": src, **extra})

    for ref, t in _section(export, ws, "trends"):
        name = t.get("marker") or t.get("label_zh") or "指标"
        lab = _measure([{"type": "lab", "ref": str(name)}] + ([{"type": "doi", "ref": t["bv_source"]}] if t.get("bv_source") else []))
        tr = t.get("trend") or {}
        if tr.get("verdict") in TREND_ZH:
            nums = [_num(f"{ref}.trend.annual_change_pct", tr["annual_change_pct"], "%/年", lab["label_zh"])] \
                if isinstance(tr.get("annual_change_pct"), (int, float)) else []
            card(f"trend-{name}", "trend", f"{name}：{TREND_ZH[tr['verdict']]}",
                 t.get("wording_zh") or f"几次检测下来，{name}出现超出你个人正常波动的{TREND_ZH[tr['verdict']][-2:]}趋势，下一次复查可以确认。",
                 nums, lab)
        bl = t.get("baseline_vs_latest") or {}
        if bl.get("verdict") in CHANGE_ZH:
            card(f"change-{name}", "change", f"{name}：{CHANGE_ZH[bl['verdict']][:2]}",
                 f"和第一次相比，{name}{CHANGE_ZH[bl['verdict']]}。" + ("这段时间你在做的事见详情。" if t.get("concurrent_exposures") else ""),
                 [_num(f"{ref}.baseline_vs_latest.change_pct", bl["change_pct"], "%", lab["label_zh"])]
                 if isinstance(bl.get("change_pct"), (int, float)) else [], lab)

    for ref, p in _section(export, ws, "predictions"):
        lab = p.get("label") or {}
        if L.label_problems(lab):
            continue
        tag = lab["label_zh"]
        nums = []
        own = any(x.get("type") == "personal_experiment" for x in lab.get("sources") or [])
        for k, u in (("p_detect_favorable", "概率"), ("p_true_favorable", "概率"), ("interval_80_pct", "%"), ("expected_change_pct", "%")):
            if p.get(k) is not None and not (own and k == "p_true_favorable"):     # see labels.item_text
                nums.append(_num(f"{ref}.{k}", p[k], u, tag))
        card(f"pred-{p.get('id')}", "prediction", f"{p.get('intervention_zh')}对{p.get('marker')}：有多大可能有用",
             L.item_text("p", p), nums, lab, plan_item=p.get("plan_item"))

    for ref, d in _section(export, ws, "diagnoses"):
        lab = d.get("label") or {}
        if L.label_problems(lab):
            continue
        if d.get("kind") == "ai_estimate":
            nums = [_num(f"{ref}.probability_range", d.get("probability_range"), "概率", lab["label_zh"])]
        else:
            nums = [_num(f"{ref}.evidence[{i}].value", e.get("value"), e.get("unit", ""), "测量")
                    for i, e in enumerate(d.get("evidence") or []) if isinstance(e.get("value"), (int, float))]
        card(f"dx-{d.get('id')}", "diagnosis", str(d.get("label_zh")),
             str(d.get("conclusion_zh")) + "".join(d.get("caveats_zh") or []), nums, lab, [CHAT, DOCTOR])   # caveats as in the report

    for ref, m in _section(export, ws, "medications"):
        lab = m.get("label") or {}
        if L.label_problems(lab):
            continue
        kind = "pgx" if m.get("kind") == "cpic" else "medication"
        title = str(m.get("label_zh") or m.get("drug_or_class_zh"))
        card(f"rx-{m.get('id')}", kind, title, L.item_text("rx", m), [], lab, [CHAT, DOCTOR], action=m.get("action"),
             currently_taking=bool(m.get("currently_taking")))

    for ref, c in _section(export, ws, "concordance"):
        if c.get("class") not in CLASS_ZH:
            continue
        name = c.get("label_zh") or c.get("analyte")
        lab = L.make_label("published_model", [{"type": "readout", "ref": f"genetic_score:{c.get('analyte')}"},
                                               {"type": "lab", "ref": str(c.get("analyte"))}])
        nums = [_num(f"{ref}.genetic_pct", c["genetic_pct"], "百分位", lab["label_zh"])] if isinstance(c.get("genetic_pct"), (int, float)) else []
        card(f"conc-{c.get('analyte')}", "concordance", f"{name}：{CLASS_ZH[c['class']]}",
             c.get("wording_zh") or f"{name}：{CLASS_ZH[c['class']]}。", nums, lab)

    for ref, r in _section(export, ws, "retests"):           # a plan date, not a judgement about the member: no label
        what = r.get("what") or r.get("marker") or "复测"
        nums = [_num(f"{ref}.after_weeks", r["after_weeks"], "周", "方案")] if isinstance(r.get("after_weeks"), (int, float)) else []
        card(f"retest-{r.get('item')}-{what}", "retest_due", f"复测：{what}", "到复测时间时，按计划复测，就能看出变化是否超出你的正常波动。",
             nums, None, due=r.get("due_at") or r.get("due"), judgeable=r.get("judgeable"), plan_item=r.get("item"))

    for ref, s in _section(export, ws, "surveillance"):
        ss = s.get("source") or {}
        sources = [{"type": "guideline", "ref": ss.get("title") or ss.get("doi") or ss.get("url"),
                    **({"url": f"https://doi.org/{ss['doi']}" if ss.get("doi") else ss["url"]} if ss.get("doi") or ss.get("url") else {}),
                    **({"quote": ss["quote"]} if ss.get("quote") else {})}] \
            if (ss.get("doi") or ss.get("url") or ss.get("title")) else [{"type": "rule", "ref": ",".join(s.get("rule_ids") or []) or str(s.get("id"))}]
        # a guideline rule carries its guideline label; a retest the plan chose is a plan date (no label, as retest_due)
        lab = L.make_label("guideline_recommendation", sources, guideline_zh=ss.get("title")) if ss.get("title") else None
        nums = [_num(f"{ref}.interval_weeks", s["interval_weeks"], "周", lab["label_zh"] if lab else "方案")] \
            if isinstance(s.get("interval_weeks"), (int, float)) else []
        card(f"surv-{s.get('id') or s.get('marker')}", "surveillance",
             f"{'测一次就够' if s.get('once') else '定期查'}：{s.get('marker')}",
             s.get("note_zh") or ("查一次就够。" if s.get("once") else "按计划定期复查。"), nums, lab, once=s.get("once"), done=s.get("done"))
    return cards


def resolve(ref: str, export: Dict[str, Any], ws: Path) -> Any:
    """The value at a card number's ref: "a[0].b.c" in the export, or "file:<rel path>#a[0].b" in a handover file."""
    if ref.startswith("file:"):
        rel, path = ref[5:].split("#", 1)
        obj: Any = load_json(ws / rel)
    else:
        obj, path = export, ref
    for name, idx in re.findall(r"([^.\[\]]+)|\[(\d+)\]", path):
        obj = obj[int(idx)] if idx else obj[name]
    return obj


def verify(cards: List[Dict[str, Any]], export_path: Path, ws: Path) -> List[str]:
    export = load_json(export_path)
    sha = sha256_file(export_path, limit=None)
    bad = []
    for c in cards:
        if c.get("source", {}).get("sha256") != sha:
            bad.append(f"{c['id']}: not made from this export")
        if not (c["kind"] == "retest_due" or (c["kind"] == "surveillance" and c.get("label") is None)) \
                and L.label_problems(c.get("label")):
            bad.append(f"{c['id']}: label problems {L.label_problems(c.get('label'))}")
        for n in c.get("numbers", []):
            try:
                v = resolve(n["ref"], export, ws)
            except (KeyError, IndexError, TypeError, ValueError, OSError):
                bad.append(f"{c['id']}: {n['ref']} does not resolve")
                continue
            if v != n["value"]:
                bad.append(f"{c['id']}: {n['ref']} is {v!r}, card says {n['value']!r}")
    return bad
