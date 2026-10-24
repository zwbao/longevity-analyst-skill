"""Fixed wording for `twin series` results: the report section and the export's `trends`.

Wording is a rule (PLAN WS-A.5, C.2, D.2, D.3): a trend is "超出你个人正常波动的缓慢上升/下降" and the next retest
confirms it; it is never called worse, a disease or a diagnosis. Concurrent interventions are "同期", never the
cause. A marker without a noise model is shown side by side only. The numbers come from series.json, never from
authored prose.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import load_json

UP, DOWN = "上升", "下降"
TREND_ZH = "出现了超出你个人正常波动的缓慢{dir}，下一次复测可以确认。"
CHANGE_ZH = "{to} 和 {frm} 相比的{dir}超出了你个人正常波动。"
PAIR_RETURNED_ZH = "中间有一次变化超出了你个人正常波动，之后又回到了波动范围内。"
PAIR_ZH = "中间有一次变化超出了你个人正常波动。"
WITHIN_ZH = "几次结果的差别都在你个人正常波动之内，还看不出真实变化。"
NOT_JUDGED_ZH = "这一项没有噪声模型，只能并排看几次结果，不能判断是不是真的变化。"
ASK_DOCTOR_ZH = "这个方向对健康不利，而且最近一次已经超出参考范围，建议列入问医生的单子。"
CONCURRENT_ZH = "这段时间你同期在做：{items}。"
CONCURRENT_ITEM_ZH = "{label}（累计 {n} 次）"
SECTION_TITLE = "复测对比与趋势"
INTRO_ZH = ("下面把你几次检测的结果放在一起看：每一项都和你自己的正常波动比（按已发表的个体内生物学变异），"
            "超出波动才说有变化；趋势是指几次结果连起来、超出波动的缓慢变化，它不是诊断。")


def _direction(m: Dict[str, Any]) -> Optional[str]:
    tv = (m.get("trend") or {}).get("verdict")
    if tv in ("rising_beyond_noise", "falling_beyond_noise"):
        return UP if tv == "rising_beyond_noise" else DOWN
    for c in [m.get("baseline_vs_latest") or {}] + list(m.get("pairs") or []):
        if c.get("verdict") in ("increase_beyond_noise", "decrease_beyond_noise"):
            return UP if c["verdict"] == "increase_beyond_noise" else DOWN
    return None


def _range(text: str) -> Optional[Tuple[Optional[float], Optional[float]]]:
    t = str(text or "").replace("～", "-").replace("~", "-").replace("—", "-").replace("–", "-").strip()
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", t)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.fullmatch(r"[<＜≤]\s*=?\s*(\d+(?:\.\d+)?)", t)
    if m:
        return None, float(m.group(1))
    m = re.fullmatch(r"[>＞≥]\s*=?\s*(\d+(?:\.\d+)?)", t)
    if m:
        return float(m.group(1)), None
    return None


def ask_doctor(m: Dict[str, Any]) -> bool:
    """An adverse-direction trend or change whose latest value is outside the member's own reference range."""
    d = _direction(m)
    better = m.get("better")
    if not d or better not in ("lower", "higher", "range") or m.get("summary_verdict") not in ("trend_beyond_noise", "change_beyond_noise"):
        return False
    rng = _range(m.get("latest_ref_range"))
    last = m.get("latest_value")
    if rng is None or last is None:
        return False
    lo, hi = rng
    if d == UP and better in ("lower", "range"):
        return hi is not None and last > hi
    if d == DOWN and better in ("higher", "range"):
        return lo is not None and last < lo
    return False


def concurrent_zh(items: List[Dict[str, Any]]) -> Optional[str]:
    if not items:
        return None
    parts = [CONCURRENT_ITEM_ZH.format(label=x["label_zh"], n=x["cumulative_in_window"])
             if isinstance(x.get("cumulative_in_window"), int) else str(x["label_zh"]) for x in items]
    return CONCURRENT_ZH.format(items="、".join(parts))


def wording(m: Dict[str, Any]) -> str:
    sv = m.get("summary_verdict")
    d = _direction(m)
    if sv == "trend_beyond_noise":
        return TREND_ZH.format(dir=d)
    if sv == "change_beyond_noise":
        b = m.get("baseline_vs_latest") or {}
        if b.get("verdict") in ("increase_beyond_noise", "decrease_beyond_noise"):
            return CHANGE_ZH.format(to=b.get("to"), frm=b.get("from"), dir=d)
        returned = any(c.get("code") == "returned_within_noise" for c in m.get("caveats") or [])
        return PAIR_RETURNED_ZH if returned else PAIR_ZH
    if sv == "within_noise":
        return WITHIN_ZH
    return NOT_JUDGED_ZH


def trend_items(ser: Dict[str, Any]) -> List[Dict[str, Any]]:
    """export `trends`: every judged marker, with its fixed wording."""
    out = []
    for m in ser.get("markers") or []:
        if m.get("summary_verdict") == "not_judged":
            continue
        tr = m.get("trend") or {}
        it = {"marker": m["marker"], "marker_key": m.get("marker_key"), "unit": m.get("unit"),
              "summary_verdict": m["summary_verdict"], "direction": {UP: "up", DOWN: "down"}.get(_direction(m)),
              "wording_zh": wording(m), "observed_at": m.get("observed_at"),
              "points": m.get("points"), "pairs": m.get("pairs"), "baseline_vs_latest": m.get("baseline_vs_latest"),
              "trend": {k: tr.get(k) for k in ("verdict", "why", "n", "span_days", "slope_per_year", "annual_change_pct", "z", "z_crit") if k in tr},
              "ask_doctor": ask_doctor(m), "caveats_zh": [c["zh"] for c in m.get("caveats") or []]}
        if it["ask_doctor"]:
            it["ask_doctor_zh"] = ASK_DOCTOR_ZH
        if m.get("concurrent_exposures"):
            it["concurrent_exposures"] = m["concurrent_exposures"]
            it["concurrent_zh"] = concurrent_zh(m["concurrent_exposures"])
        out.append(it)
    return out


def _fmt(v: Any) -> str:
    return f"{v:.2f}".rstrip("0").rstrip(".") if isinstance(v, float) else str(v)


def sections(ws_root: Path, ctx: Dict[str, Any]) -> List[Tuple[str, str]]:
    sp = Path(ws_root) / "work" / "twin" / "series.json"
    if not sp.exists():
        return []
    ser = load_json(sp)
    from .report import _cell as esc
    L: List[str] = [INTRO_ZH, ""]
    for a in ser.get("identity_alerts") or []:
        L += [f"> **样本核对**：{a['wording_zh']}", ""]
    judged = [m for m in ser.get("markers") or [] if m.get("summary_verdict") != "not_judged"]
    order = {"trend_beyond_noise": 0, "change_beyond_noise": 1, "within_noise": 2}
    for m in sorted(judged, key=lambda x: (order.get(x["summary_verdict"], 3), str(x["marker"]))):
        pts = " → ".join(f"{p['t']} {_fmt(p['value'])}" for p in m.get("points") or [])
        line = f"- **{esc(m['marker'])}**：{wording(m)}（{pts} {m.get('unit') or ''}）".replace(" ）", "）")
        extra = []
        if ask_doctor(m):
            extra.append(ASK_DOCTOR_ZH)
        cz = concurrent_zh(m.get("concurrent_exposures") or [])
        if cz:
            extra.append(esc(cz))
        extra += [c["zh"] + "。" for c in m.get("caveats") or [] if c.get("code") in ("returned_within_noise", "seasonal", "snapshot_excluded_identity")]
        L.append(line + ("".join(extra) if extra else ""))
    rj = [r for r in ser.get("readouts") or [] if (r.get("baseline_vs_latest") or {}).get("verdict") != "not_judged"]
    for r in rj:
        b = r["baseline_vs_latest"]
        txt = (f"{b.get('wording_zh')}（{b['change_abs']} {r.get('unit') or ''}）。{b.get('note_zh')}" if b.get("wording_zh")
               else "两次读数的差别在测量误差之内。")
        L.append(f"- **{esc(r.get('label_zh') or r['readout'])}**：{txt}")
    nj = [m["marker"] for m in ser.get("markers") or [] if m.get("summary_verdict") == "not_judged" and len(m.get("points") or []) > 1]
    nj += [r.get("label_zh") or r["readout"] for r in ser.get("readouts") or [] if (r.get("baseline_vs_latest") or {}).get("verdict") == "not_judged"]
    if nj:
        L += ["", f"以下各项没有噪声模型，只并排展示，不判断好坏：{'、'.join(esc(x) for x in nj)}。"]
    tr = ser.get("transitions") or []
    if tr:
        L += ["", "基因倾向与化验的对照有了更新（基因分值只能解释这一项差异的一部分）："]
        for t in tr:
            if t.get("wording_zh"):
                cz = concurrent_zh(t.get("concurrent_exposures") or [])
                L.append(f"- **{esc(t['analyte'])}**：{t['wording_zh']}{esc(cz) if cz else ''}")
    if len(L) <= 2:
        return []
    return [(SECTION_TITLE, "\n".join(L).strip())]
