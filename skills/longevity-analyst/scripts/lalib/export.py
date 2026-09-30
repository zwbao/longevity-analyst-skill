"""la-export/1: what a host app (longpi) imports from a finished delivery.

Written by `la.py report` next to report.html, from the same bound files, so it is never newer or older than the
report; the deliver/ folder (and this file) is removed whenever the report stops being current. All text is final:
placeholders are already substituted exactly as the report prints them.
"""
from __future__ import annotations

from datetime import date, timedelta
from pathlib import Path
from typing import Any, Dict, List

from .common import data, load_json, now_iso, sha256_file

SCHEMA = "la-export/1"
# analyst plan category -> longpi-plan/1 category
CATEGORY = {"diet": "diet", "exercise": "exercise", "sleep": "sleep", "supplement": "supplement", "lifestyle": "behavior",
            "test": "other", "referral": "other"}
GROUP = {"genetic_score": "insight", "genetic_finding": "insight", "population_position": "insight", "wearable": "insight",
         "llm_estimate": "organ_ai_estimate"}


def build(st: Dict[str, Any], ws: Path, report_html: Path) -> Dict[str, Any]:
    from .report import _readouts, substitute
    from . import board as B
    from .organs import table
    ro = _readouts(ws)
    m = st["member"]
    today = date.today()
    readouts = []
    for r in ro.values():
        row = {k: r.get(k) for k in ("id", "label_zh", "value", "unit", "kind", "method", "low", "high", "horizon_years") if r.get(k) is not None}
        row["group"] = GROUP.get(r.get("kind"), "method")
        if r.get("provenance_uncertain"):
            row["provenance_uncertain"] = True
        readouts.append(row)
    organs = []
    for o in table(ws, ro):
        if not (o["measured"] or o["indices"] or o["ai_age"] or o["ai_risks"]):
            continue
        organs.append({"organ": o["organ"], "label_zh": o["label_zh"],
                       "measured": [x["id"] for x in o["measured"]], "indices": [x["id"] for x in o["indices"]],
                       "ai_age": o["ai_age"]["id"] if o["ai_age"] else None, "ai_risks": [x["id"] for x in o["ai_risks"]]})
    board = []
    for r in B.rows(ws, st):
        q, f = r["q"], r["f"]
        board.append({"id": q["id"], "title_zh": substitute(q["title_zh"], ro), "hypothesis_zh": substitute(q["hypothesis_zh"], ro),
                      "verdict": (f or {}).get("verdict"), "verdict_zh": B.VERDICTS.get((f or {}).get("verdict"), "未研究" if r.get("skipped") else "未完成"),
                      "confidence": (f or {}).get("confidence"),
                      "summary_zh": substitute(f["summary_zh"], ro) if f else None,
                      "next_step_zh": substitute(f["next_step_zh"], ro) if f else None,
                      "skipped_reason_zh": r.get("skipped")})
    plan_items, retests = [], []
    pl = ws / "work" / "intervene" / "plan.json"
    if pl.exists():
        for i in load_json(pl)["items"]:
            rt = i.get("retest") or {}
            what = substitute(str(rt.get("what") or ""), ro).strip()
            markers = [ro[t]["label_zh"] if t in ro else t for t in i.get("targets", []) if not str(t).startswith("file:")]
            if what and what not in markers:
                markers.append(what)
            plan_items.append({"id": f"la-{i['id'].lower()}", "category": CATEGORY.get(i["category"], "other"),
                               "title": substitute(i["action_zh"], ro)[:60], "detail": substitute(i.get("rationale_zh") or "", ro)[:300],
                               "start": today.isoformat(), "markers": markers[:12], "executor": i.get("executor"),
                               "evidence_grade": i.get("evidence_grade"), "source_item": i["id"]})
            if rt.get("after_weeks"):
                retests.append({"item": i["id"], "what": what, "after_weeks": int(rt["after_weeks"]),
                                "due": (today + timedelta(weeks=int(rt["after_weeks"]))).isoformat()})
    tw = ws / "work" / "twin" / "twin.json"
    return {"schema": SCHEMA, "generated_at": now_iso(), "generation": st.get("generation"),
            "member": {"id": m.get("id"), "age": m.get("age"), "sex": m.get("sex"), "sample_date": m.get("sample_date")},
            "workspace": str(ws),
            "report": {"html": str(report_html), "sha256": sha256_file(report_html, limit=None)},
            "twin": {"path": str(tw), "sha256": sha256_file(tw, limit=None)} if tw.exists() else None,
            "stages": st["stages"], "readouts": readouts, "organs": organs, "board": board,
            "plan": {"title": "深度分析干预方案", "source": "analysis",
                     "note": "来自 longevity-analyst 深度分析；营养师审核后生效。方案只记做什么，不记剂量。", "items": plan_items},
            "retests": retests,
            "boundary_zh": "本结果用于健康管理参考，不是诊断；基因发现需临床确认和遗传咨询；AI 估计不是测量值。"}
