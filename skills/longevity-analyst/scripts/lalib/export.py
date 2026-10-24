"""la-export/2 (default) and la-export/1: what a host app (longpi) imports from a finished delivery.

Written by `la.py report` next to report.html, from the same bound files, so it is never newer or older than the
report; the deliver/ folder (and this file) is removed whenever the report stops being current.

Text is plain (no markdown): placeholders are substituted as the report prints them, a PMID is "PMID n", and a
quantity with a unit (5 克, 1500 毫升) becomes "（具体量见报告）", because the host keeps no amounts in a plan and
would otherwise cut the number out of the sentence. The organ table carries the diagnostic-threshold notes that
replace an AI probability. Plan markers are the member's own lab names only, the names the host tracks retests by.

/2 adds, from the handover files when they exist (a missing file leaves its field empty): `trends` (twin series,
fixed wording), `concordance` (+ `concordance_transitions`), `surveillance`, `predictions`, `diagnoses`,
`medications`, `identity_alerts`, and `retests[].judgeable / reason / due_at`. /1 is the same document without
them (`la.py export <ws> --schema la-export/1`). Before anything is written, `privacy_problems` refuses names, phone
numbers, any rsid together with a genotype, and the germline fingerprint: those never leave the workspace.
"""
from __future__ import annotations

import re
from datetime import date, timedelta
from pathlib import Path
from typing import Any, Dict, List

from .common import EXIT_INPUT, LAError, load_json, now_iso, sha256_file

SCHEMA = "la-export/2"
SCHEMA_V1 = "la-export/1"
TZ = "+08:00"                       # due dates are wall-clock dates for the member; China time unless told otherwise
V2_FIELDS = ("trends", "concordance", "concordance_note_zh", "concordance_transitions", "surveillance", "predictions",
             "diagnoses", "medications", "identity_alerts")
CONCORDANCE_NOTE_ZH = "基因分值只能解释这一项差异的一部分。"
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
    base = _base_date(m, rendered)
    pl = ws / "work" / "intervene" / "plan.json"
    judge: Dict[str, Dict[str, Any]] = {}
    if pl.exists():
        from .twin import plan_retests
        judge = {r["item"]: r for r in plan_retests(ws, load_json(pl))}
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
                chk = judge.get(i["id"]) or {}
                retests.append({"item": i["id"], "what": what, "after_weeks": int(rt["after_weeks"]),
                                "due": (rendered + timedelta(weeks=int(rt["after_weeks"]))).isoformat(), "from": rendered.isoformat(),
                                "reason": "intervention", "judgeable": bool(chk.get("judgeable")),
                                "due_at": _due_at(base, int(rt["after_weeks"]))})
    retests += _surveillance_retests(ws, base, {r["item"] for r in retests})
    tw = ws / "work" / "twin" / "twin.json"
    doc = {"schema": SCHEMA, "generated_at": now_iso(), "generation": st.get("generation"),
            "member": {"id": m.get("id"), "age": m.get("age"), "sex": m.get("sex"), "sample_date": m.get("sample_date")},
            "workspace": str(ws),
            "report": {"html": str(report_html), "sha256": sha256_file(report_html, limit=None)},
            "twin": {"path": str(tw), "sha256": sha256_file(tw, limit=None)} if tw.exists() else None,
            "stages": st["stages"], "readouts": readouts, "organs": organs, "board": board,
            "plan": {"title": "深度分析干预方案", "source": "analysis",
                     "note": "来自 longevity-analyst 深度分析；营养师审核后生效。方案只记做什么，具体用量见报告。", "items": plan_items},
            "retests": retests,
            "boundary_zh": "本结果用于健康管理参考，不是诊断；基因发现需临床确认和遗传咨询；AI 估计不是测量值。"}
    doc.update(_v2(ws, ro))
    bad = privacy_problems(doc)
    if bad:
        raise LAError("la-export refused (privacy): " + "; ".join(bad[:5]), EXIT_INPUT)
    return doc


def _base_date(m: Dict[str, Any], rendered: date) -> date:
    try:
        return date.fromisoformat(str(m.get("sample_date"))[:10])
    except ValueError:
        return rendered


def _due_at(base: date, weeks: int) -> str:
    return f"{(base + timedelta(weeks=weeks)).isoformat()}T00:00:00{TZ}"


def _items(p: Path) -> List[Dict[str, Any]]:
    """`items` of a handover file; a missing or unreadable file is an absent feature, not an error."""
    try:
        doc = load_json(p) if p.exists() else {}
    except (OSError, ValueError):
        return []
    return list(doc.get("items") or []) if isinstance(doc, dict) else []


def _surveillance_retests(ws: Path, base: date, have: set) -> List[Dict[str, Any]]:
    """Recurring surveillance items (genetic / guideline reasons) join the retests; intervention items are already
    there, and `once` items stay in `surveillance` only."""
    from .twin import retest_check
    from . import series as S
    out = []
    N = None
    for it in _items(ws / "work" / "surveillance" / "plan.json"):
        wk = it.get("interval_weeks")
        if it.get("reason") == "intervention" or it.get("once") or not isinstance(wk, int) or wk <= 0:
            continue
        N = N or S.noise()
        chk = retest_check(N, str(it.get("marker") or ""), [], wk)
        out.append({"item": it.get("id"), "what": it.get("marker"), "after_weeks": wk,
                    "due": (base + timedelta(weeks=wk)).isoformat(), "from": base.isoformat(),
                    "reason": it.get("reason") or "guideline", "judgeable": chk["judgeable"], "due_at": _due_at(base, wk),
                    "executor": it.get("executor")})
    return out


def _v2(ws: Path, ro: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    from .report_trends import trend_items
    sp = ws / "work" / "twin" / "series.json"
    ser = load_json(sp) if sp.exists() else {}
    conc = _items(ws / "work" / "insights" / "concordance.json")
    conc = [{k: v for k, v in c.items() if k != "basis"} for c in conc]
    return {"trends": trend_items(ser) if ser else [],
            "concordance": conc, "concordance_note_zh": CONCORDANCE_NOTE_ZH if conc else None,
            "concordance_transitions": [{k: t.get(k) for k in ("analyte", "from", "to", "from_t", "to_t", "kind", "ask_doctor",
                                                                "real_change", "wording_zh") if k in t}
                                        for t in ser.get("transitions") or []],
            "identity_alerts": [{"t": x.get("t"), "basis": x.get("basis"), "wording_zh": x.get("wording_zh")}
                                for x in ser.get("identity_alerts") or []],
            "surveillance": [{k: v for k, v in it.items()} for it in _items(ws / "work" / "surveillance" / "plan.json")],
            "predictions": _items(ws / "work" / "predict" / "predictions.json"),
            "diagnoses": _items(ws / "work" / "clinical" / "diagnoses.json"),
            "medications": _items(ws / "work" / "clinical" / "medications.json")}


def as_v1(doc: Dict[str, Any]) -> Dict[str, Any]:
    """The same delivery in the la-export/1 shape, for a host that does not read /2 yet."""
    out = {k: v for k, v in doc.items() if k not in V2_FIELDS}
    out["schema"] = SCHEMA_V1
    out["retests"] = [{k: v for k, v in r.items() if k in ("item", "what", "after_weeks", "due", "from")}
                      for r in doc.get("retests") or [] if r.get("reason", "intervention") == "intervention"]
    return out


# ------------------------------------------------------------ privacy (WS-K)
RSID = re.compile(r"\brs\d{3,}\b", re.I)
GENOTYPE = re.compile(r"^(?:[ACGT]{1,2}|[ACGT][/|;:,][ACGT]|[ACGT]{1,2}/[ACGT]{1,2}|DI|II|DD|--)$")
RSID_GT_TEXT = re.compile(r"\brs\d{3,}\s*[:=（(\s]\s*(?:[ACGT]{1,2}\b|[ACGT]\s*[/|;]\s*[ACGT]\b)", re.I)
PHONE = re.compile(r"(?<![0-9A-Za-z])(?:\+?86[-\s]?)?1[3-9]\d[-\s]?\d{4}[-\s]?\d{4}(?![0-9A-Za-z])")
NAME_KEYS = {"name", "full_name", "real_name", "patient_name", "member_name", "姓名", "phone", "mobile", "tel", "手机", "手机号",
             "电话", "genotype_fingerprint", "fingerprint", "calls"}
NAME_TEXT = re.compile(r"(?:姓名|患者姓名|受检者)\s*[:：]\s*\S")


def privacy_problems(obj: Any, path: str = "$") -> List[str]:
    """Names, phone numbers, the germline fingerprint and any rsid + genotype combination never reach a host."""
    out: List[str] = []
    if isinstance(obj, dict):
        rs_vals = [v for v in obj.values() if isinstance(v, str) and RSID.fullmatch(v.strip())]
        gt_vals = [v for k, v in obj.items() if isinstance(v, str) and GENOTYPE.match(v.strip().upper())
                   and str(k).lower() in ("genotype", "gt", "call", "alleles", "diplotype_or_genotype", "value")]
        if rs_vals and gt_vals:
            out.append(f"{path}: an rsid with a genotype")
        for k, v in obj.items():
            kl = str(k).lower()
            if kl in NAME_KEYS:
                out.append(f"{path}.{k}: field not allowed in an export")
                continue
            if RSID.fullmatch(str(k)) and isinstance(v, str) and GENOTYPE.match(v.strip().upper()):
                out.append(f"{path}.{k}: an rsid with a genotype")
                continue
            out += privacy_problems(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += privacy_problems(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        if RSID_GT_TEXT.search(obj):
            out.append(f"{path}: an rsid with a genotype in text")
        if PHONE.search(obj):
            out.append(f"{path}: a phone number")
        if NAME_TEXT.search(obj):
            out.append(f"{path}: a name")
    return out
