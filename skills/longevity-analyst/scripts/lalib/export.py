"""la-export/1: what a host app (longpi) imports from a finished delivery.

Written by `la.py report` next to report.html, from the same bound files, so it is never newer or older than the
report; the deliver/ folder (and this file) is removed whenever the report stops being current.

Text is plain (no markdown): placeholders are substituted as the report prints them, a PMID is "PMID n", and a
quantity with a unit (5 克, 1500 毫升) becomes "（具体量见报告）", because the host keeps no amounts in a plan and
would otherwise cut the number out of the sentence. The organ table carries the diagnostic-threshold notes that
replace an AI probability. Plan markers are the member's own lab names only, the names the host tracks retests by.
"""
from __future__ import annotations

import re
from datetime import date, timedelta
from pathlib import Path
from typing import Any, Dict, List

from .common import load_json, now_iso, sha256_file

SCHEMA = "la-export/1"
# analyst plan category -> longpi-plan/1 category
CATEGORY = {"diet": "diet", "exercise": "exercise", "sleep": "sleep", "supplement": "supplement", "lifestyle": "behavior",
            "test": "other", "referral": "other"}
GROUP = {"genetic_score": "genetic", "genetic_finding": "genetic", "population_position": "insight", "wearable": "insight",
         "llm_estimate": "organ_ai_estimate"}
GENETIC_METHODS = {"native.apoe"}
QUANTITY = re.compile(r"\d\s*(?:克|千克|公斤|毫克|微克|g|kg|mg|μg|ug|毫升|ml|mL|升|L|杯|片|粒|勺|袋|支|IU|国际单位|千卡|kcal|大卡)(?![A-Za-z])")
MD_LINK = re.compile(r"\[([^\]]*)\]\((https?://[^)]*)\)")
TITLE_MAX, DETAIL_MAX = 60, 300


def plain(text: str, ro: Dict[str, Dict[str, Any]]) -> str:
    """The report's substitution, as plain text a host shows as is."""
    from .report import LITERAL, PMID_PH, substitute
    text = LITERAL.sub(lambda m: "（具体量见报告）" if QUANTITY.search(m.group(1)) else m.group(0), text)
    text = PMID_PH.sub(lambda m: f"PMID {m.group(1)}", text)
    text = substitute(text, ro)
    text = MD_LINK.sub(lambda m: m.group(1), text)
    return re.sub(r"\s+", " ", text).strip()


def _title(text: str) -> str:
    """At most 60 characters, cut at a clause boundary rather than inside a word or a qualifier."""
    if len(text) <= TITLE_MAX:
        return text
    cut = max(text.rfind(p, 0, TITLE_MAX) for p in "，。；、（(")
    return text[:cut].rstrip("，。；、（( ") if cut >= 12 else text[:TITLE_MAX]


def _lab_names(st: Dict[str, Any]) -> Dict[str, str]:
    """folded name -> the member's own lab marker as printed (confirmed rows only)."""
    from .common import skillkit
    from .labnames import usable_rows
    kit = skillkit()
    out: Dict[str, str] = {}
    for r in usable_rows(st):
        out.setdefault(kit.fold_name(r["marker"]), r["marker"])
        if r.get("maps_to"):
            out.setdefault(kit.fold_name(r["maps_to"]), r["marker"])
    return out


def build(st: Dict[str, Any], ws: Path, report_html: Path) -> Dict[str, Any]:
    from .report import _readouts
    from . import board as B
    from .organs import table
    from .common import skillkit
    kit = skillkit()
    ro = _readouts(ws)
    m = st["member"]
    readouts = []
    for r in ro.values():
        if r.get("kind") in ("implausible", "member_answer"):
            continue                                    # the report leaves these out too
        row = {k: r.get(k) for k in ("id", "label_zh", "value", "unit", "kind", "method", "low", "high", "horizon_years",
                                     "confidence", "coverage_pct") if r.get(k) is not None}
        row["group"] = "genetic" if r.get("method") in GENETIC_METHODS else GROUP.get(r.get("kind"), "method")
        if r.get("provenance_uncertain"):
            row["provenance_uncertain"] = True
        readouts.append(row)
    organs = []
    for o in table(ws, ro):
        overrides = [{"disease": x["disease"], "message_zh": x["message_zh"]} for x in o.get("overrides") or []]
        if not (o["measured"] or o["indices"] or o["ai_age"] or o["ai_risks"] or overrides):
            continue
        organs.append({"organ": o["organ"], "label_zh": o["label_zh"],
                       "measured": [x["id"] for x in o["measured"]], "indices": [x["id"] for x in o["indices"]],
                       "ai_age": o["ai_age"]["id"] if o["ai_age"] else None, "ai_risks": [x["id"] for x in o["ai_risks"]],
                       "overrides": overrides})
    board = []
    for r in B.rows(ws, st):
        q, f = r["q"], r["f"]
        board.append({"id": q["id"], "title_zh": plain(q["title_zh"], ro), "hypothesis_zh": plain(q["hypothesis_zh"], ro),
                      "verdict": (f or {}).get("verdict"),
                      "verdict_zh": B.VERDICTS.get((f or {}).get("verdict"), "未研究" if r.get("skipped") else "未完成"),
                      "confidence": (f or {}).get("confidence"),
                      "summary_zh": plain(f["summary_zh"], ro) if f else None,
                      "next_step_zh": plain(f["next_step_zh"], ro) if f else None,
                      "limitations_zh": plain(f["limitations_zh"], ro) if f else None,
                      "public_evidence": list((f or {}).get("public_evidence") or [])[:12],
                      "skipped_reason_zh": r.get("skipped")})
    labs = _lab_names(st)
    plan_items, retests = [], []
    rendered = date.today()
    pl = ws / "work" / "intervene" / "plan.json"
    if pl.exists():
        for i in load_json(pl)["items"]:
            rt = i.get("retest") or {}
            what = plain(str(rt.get("what") or ""), ro)
            markers: List[str] = []
            for t in i.get("targets", []):                 # retests are tracked by the member's own lab names
                name = labs.get(kit.fold_name(ro[t]["label_zh"])) if t in ro else labs.get(kit.fold_name(str(t)))
                if name and name not in markers:
                    markers.append(name)
            if what and labs.get(kit.fold_name(what)) and labs[kit.fold_name(what)] not in markers:
                markers.append(labs[kit.fold_name(what)])
            action = plain(i["action_zh"], ro)
            rationale = plain(i.get("rationale_zh") or "", ro)
            title = _title(action)
            detail = (action + "。" if title != action else "") + rationale
            plan_items.append({"id": f"la-{i['id'].lower()}", "category": CATEGORY.get(i["category"], "other"),
                               "title": title, "detail": detail[:DETAIL_MAX], "markers": markers[:12],
                               "executor": i.get("executor"), "kind": i["category"], "evidence_grade": i.get("evidence_grade"),
                               "source_item": i["id"]})
            if rt.get("after_weeks"):
                retests.append({"item": i["id"], "what": what, "after_weeks": int(rt["after_weeks"]),
                                "due": (rendered + timedelta(weeks=int(rt["after_weeks"]))).isoformat(), "from": rendered.isoformat()})
    tw = ws / "work" / "twin" / "twin.json"
    return {"schema": SCHEMA, "generated_at": now_iso(), "generation": st.get("generation"),
            "member": {"id": m.get("id"), "age": m.get("age"), "sex": m.get("sex"), "sample_date": m.get("sample_date")},
            "workspace": str(ws),
            "report": {"html": str(report_html), "sha256": sha256_file(report_html, limit=None)},
            "twin": {"path": str(tw), "sha256": sha256_file(tw, limit=None)} if tw.exists() else None,
            "stages": st["stages"], "readouts": readouts, "organs": organs, "board": board,
            "plan": {"title": "深度分析干预方案", "source": "analysis",
                     "note": "来自 longevity-analyst 深度分析；营养师审核后生效。方案只记做什么，具体用量见报告。", "items": plan_items},
            "retests": retests,
            "boundary_zh": "本结果用于健康管理参考，不是诊断；基因发现需临床确认和遗传咨询；AI 估计不是测量值。"}
