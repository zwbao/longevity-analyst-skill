"""la-label/1: the judgement label every prediction, diagnostic conclusion and medication suggestion carries.

The category tag (label_zh) is rendered here from basis_kind, never written by a model. Labelled items live in
three registered files (work/predict/predictions.json, work/clinical/diagnoses.json, work/clinical/medications.json)
and prose cites them only as {{p:<id>}}, {{dx:<id>}}, {{rx:<id>}}: rendering prints the item from its file with its
tag in front and its sources after, so a probability or a dose in the report always comes from a registered item.
`review trace` refuses an item without a label or sources, an AI judgement whose confidence is not low/very_low,
an AI prediction given as a point, and a dose that does not trace to its source (rule text, CPIC text) or to an AI
judgement.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import EXIT_INPUT, LAError, load_json, now_iso, sha256_file, sha_text

SCHEMA = "la-label/1"
BASIS_ZH = {"measurement": "测量", "published_model": "已发表模型", "guideline_criteria": "按《{g}》的标准",
            "guideline_recommendation": "《{g}》推荐", "cpic": "CPIC 指南", "statistical_model": "统计模型估计",
            "ai_estimate": "AI 判断"}
SOURCE_TYPES = ("effects", "biovar", "readout", "lab", "pmid", "doi", "guideline", "cpic", "personal_experiment", "rule", "url")
SOURCE_ZH = {"effects": "效应表", "biovar": "个体内变异表", "readout": "读数", "lab": "化验", "pmid": "PMID", "doi": "DOI",
             "guideline": "指南", "cpic": "CPIC", "personal_experiment": "本人实验", "rule": "规则", "url": "网页"}
AI_LEVELS = ("low", "very_low")
LEVEL_ZH = {"low": "低", "very_low": "很低", "moderate": "中", "high": "高"}
ITEM_PH = re.compile(r"\{\{(p|dx|rx):([A-Za-z0-9_\-]+)\}\}")
FILES = {"p": ("predict", "predictions.json", "la-predictions/1"),
         "dx": ("clinical", "diagnoses.json", "la-diagnoses/1"),
         "rx": ("clinical", "medications.json", "la-medications/1")}
STATE_KEY = {"p": ("predict", None), "dx": ("clinical", "dx"), "rx": ("clinical", "rx")}


def label_zh(basis_kind: str, guideline_zh: Optional[str] = None) -> str:
    if basis_kind not in BASIS_ZH:
        raise LAError(f"basis_kind {basis_kind!r} is not one of {sorted(BASIS_ZH)}", EXIT_INPUT)
    if basis_kind in ("guideline_criteria", "guideline_recommendation"):
        if not guideline_zh:
            raise LAError(f"{basis_kind} needs the guideline name", EXIT_INPUT)
        return BASIS_ZH[basis_kind].format(g=guideline_zh.strip("《》"))
    return BASIS_ZH[basis_kind]


def make_label(basis_kind: str, sources: List[Dict[str, Any]], *, guideline_zh: Optional[str] = None,
               assumptions: Optional[List[str]] = None, confidence: Optional[Dict[str, Any]] = None,
               generated_at: Optional[str] = None, physician_review: str = "none") -> Dict[str, Any]:
    lab: Dict[str, Any] = {"schema": SCHEMA, "basis_kind": basis_kind, "label_zh": label_zh(basis_kind, guideline_zh),
                           "sources": [dict(s) for s in sources]}
    if guideline_zh:
        lab["guideline_zh"] = guideline_zh.strip("《》")
    if assumptions:
        lab["assumptions"] = list(assumptions)
    if confidence is not None:
        lab["confidence"] = dict(confidence)
    lab["generated_at"] = generated_at or now_iso()
    lab["physician_review"] = physician_review
    bad = label_problems(lab)
    if bad:
        raise LAError("label rejected: " + "; ".join(bad), EXIT_INPUT)
    return lab


def label_problems(lab: Any) -> List[str]:
    """Everything that keeps a label from being rendered (the same rules as schema/la-label-1.schema.json)."""
    if not isinstance(lab, dict) or lab.get("schema") != SCHEMA:
        return ["missing label (la-label/1)"]
    out = []
    bk = lab.get("basis_kind")
    if bk not in BASIS_ZH:
        return [f"basis_kind {bk!r} unknown"]
    try:
        want = label_zh(bk, lab.get("guideline_zh"))
        if lab.get("label_zh") != want:
            out.append(f"label_zh {lab.get('label_zh')!r} is not the rendered tag {want!r} (written by the renderer, not a model)")
    except LAError as e:
        out.append(str(e))
    src = lab.get("sources")
    if not isinstance(src, list) or not src:
        out.append("missing sources")
    else:
        for s in src:
            if not isinstance(s, dict) or s.get("type") not in SOURCE_TYPES or not str(s.get("ref") or "").strip():
                out.append(f"source {s!r} needs a type from {SOURCE_TYPES} and a ref")
            elif set(s) - {"type", "ref", "quote", "url"}:
                out.append(f"source fields {sorted(set(s) - {'type', 'ref', 'quote', 'url'})} are not in la-label/1")
    extra = set(lab) - {"schema", "basis_kind", "label_zh", "guideline_zh", "sources", "assumptions", "confidence", "generated_at", "physician_review"}
    if extra:
        out.append(f"fields {sorted(extra)} are not in la-label/1")
    conf = lab.get("confidence")
    if bk == "ai_estimate":
        if not isinstance(conf, dict) or conf.get("level") not in AI_LEVELS:
            out.append("an AI judgement needs confidence level low or very_low")
    if bk == "statistical_model":
        if not isinstance(conf, dict) or not (isinstance(conf.get("interval_80"), list) and len(conf["interval_80"]) == 2):
            out.append("a statistical model estimate needs confidence.interval_80")
        if not lab.get("assumptions"):
            out.append("a statistical model estimate needs its assumptions written out")
    if isinstance(conf, dict) and "level" in conf and conf["level"] not in AI_LEVELS + ("moderate", "high"):
        out.append(f"confidence level {conf['level']!r} unknown")
    if not lab.get("generated_at"):
        out.append("missing generated_at")
    if lab.get("physician_review", "none") not in ("none", "pending", "confirmed"):
        out.append("physician_review must be none, pending or confirmed")
    return out


def render_label(label: Dict[str, Any]) -> str:
    """The tag printed before an item: 〔统计模型估计〕, 〔AI 判断，置信度低〕."""
    tag = label_zh(label["basis_kind"], label.get("guideline_zh"))
    lv = (label.get("confidence") or {}).get("level")
    if label["basis_kind"] == "ai_estimate" and lv:
        tag += f"，置信度{LEVEL_ZH.get(lv, lv)}"
    return f"〔{tag}〕"


def render_sources(label: Dict[str, Any]) -> str:
    parts = []
    for s in label.get("sources") or []:
        t, ref = s.get("type"), str(s.get("ref"))
        if t in ("guideline", "url"):
            parts.append(ref)
        else:
            parts.append(f"{SOURCE_ZH.get(t, t)} {ref}")
    return "来源：" + "；".join(parts)


PERSONAL_ZH = ("这个估计用到了你自己做过的前后对比；它没有对照，同一段时间里其他因素带来的变化也算在里面，"
               "所以把握可能偏高。")


def _about(x: Any) -> str:
    """"约 77%", but "超过 99%" / "不到 1%" without 约 in front."""
    t = _p(x)
    return t if t.startswith(("超过", "不到", "低于")) else f"约 {t}"


def _p(x: Any) -> str:
    if isinstance(x, list):
        return f"{_p(x[0])}–{_p(x[1])}"
    v = 100 * float(x)
    if v > 99:
        return "超过 99%"
    if v < 1:
        return "不到 1%"
    return f"{v:.0f}%"


def _pct(x: float) -> str:
    return f"{x:.1f}%".replace("-", "−")


def item_text(kind: str, it: Dict[str, Any]) -> str:
    """The fixed sentence for one registered item (no tag, no sources)."""
    if kind == "p":
        head = f"{it['intervention_zh']}对{it['marker']}"
        if it["label"]["basis_kind"] == "ai_estimate":
            bits = []
            if it.get("p_detect_favorable") is not None:
                bits.append(f"下次实验能测出有利变化的概率区间 {_p(it['p_detect_favorable'])}")
            if it.get("p_true_favorable") is not None:
                bits.append(f"真实效应有利的概率区间 {_p(it['p_true_favorable'])}")
            return f"{head}：AI 判断，" + "，".join(bits) + "。"
        d = it.get("design") or {}
        lo, hi = it["interval_80_pct"]
        own = any(s.get("type") == "personal_experiment" for s in it["label"].get("sources") or [])
        prior = next((a for a in it["label"].get("assumptions") or [] if a.startswith("只按研究结果")), None)
        body = (f"下次实验能测出有利变化的概率{_about(it['p_detect_favorable'])}，"
                f"真实效应方向有利的概率{_about(it['p_true_favorable'])}，预计变化的 80% 区间 {_pct(lo)} 到 {_pct(hi)}。")
        if own and prior:
            # the trials alone are the answer; the uncontrolled before/after only shifts the expected size, and no
            # "true effect" probability is printed from it (it would read the concurrent change as the cause)
            return (f"{head}：按前后各测 {d.get('repeats_each_side')} 次的设计，{prior}。"
                    f"如果假定你自己那次没有对照的前后对比全部来自这项干预：下次能测出有利变化的概率{_about(it['p_detect_favorable'])}，"
                    f"预计变化的 80% 区间 {_pct(lo)} 到 {_pct(hi)}。{PERSONAL_ZH}")
        return f"{head}：按前后各测 {d.get('repeats_each_side')} 次的设计，{body}" + (PERSONAL_ZH if own else "")
    if kind == "dx":
        g = it["label"].get("guideline_zh")
        c = it["conclusion_zh"]
        if g and c.startswith(f"按《{g}》的标准，"):             # the tag in front already names the guideline
            c = c[len(f"按《{g}》的标准，"):]
        return c + "".join(it.get("caveats_zh") or [])
    txt = it.get("text_zh") or ""
    if it.get("kind") == "cpic" and it.get("label_zh"):
        txt = f"{it['label_zh']}：{txt}"
    if it.get("dose_text_zh"):
        txt += f"（剂量：{it['dose_text_zh']}）"
    return txt + "".join(it.get("caveats_zh") or [])


def render_item(kind: str, it: Dict[str, Any]) -> str:
    """Tag + sentence + sources: the only way a labelled item reaches the member."""
    return f"{render_label(it['label'])}{item_text(kind, it)}（{render_sources(it['label'])}）"


# ------------------------------------------------------------ registry: the three item files
def item_path(ws: Path, kind: str) -> Path:
    sub, name, _ = FILES[kind]
    return ws / "work" / sub / name


def record(st: Dict[str, Any], ws: Path, kind: str, inputs_digest: str) -> Dict[str, Any]:
    """Bind a just-written item file in state: its hash (tamper check) and the digest of its inputs (staleness)."""
    key, sub = STATE_KEY[kind]
    rec = {"path": str(item_path(ws, kind)), "sha256": sha256_file(item_path(ws, kind), limit=None),
           "inputs": inputs_digest, "at": now_iso()}
    if sub:
        box = st.setdefault(key, {})             # (the right-hand side is evaluated first: never index st[key] there)
        box[sub] = dict(box.get(sub) or {}, **rec)
    else:
        st[key] = dict(st.get(key) or {}, **rec)
    return rec


def _state_rec(st: Dict[str, Any], kind: str) -> Optional[Dict[str, Any]]:
    key, sub = STATE_KEY[kind]
    v = st.get(key) or None
    return (v or {}).get(sub) if sub else v


def registry(ws: Path) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """{'p': {id: item}, 'dx': {...}, 'rx': {...}} from the files as they are (checked by `problems`)."""
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for kind in FILES:
        p = item_path(ws, kind)
        try:
            out[kind] = {x["id"]: x for x in load_json(p).get("items", [])} if p.exists() else {}
        except (OSError, ValueError, KeyError, TypeError):
            out[kind] = {}
    return out


class Bound(dict):
    """The readouts dict the report binds placeholders to, carrying the labelled items alongside."""
    labelled: Dict[str, Dict[str, Dict[str, Any]]] = {}


def substitute_items(text: str, reg: Dict[str, Dict[str, Dict[str, Any]]]) -> str:
    def rep(m: re.Match) -> str:
        it = (reg.get(m.group(1)) or {}).get(m.group(2))
        if it is None or label_problems(it.get("label")):
            return f"[缺失条目 {m.group(1)}:{m.group(2)}]"
        return render_item(m.group(1), it)
    return ITEM_PH.sub(rep, text)


def cited(text: str) -> List[tuple]:
    return [(m.group(1), m.group(2)) for m in ITEM_PH.finditer(text)]


def bound_files(ws: Path) -> List[Path]:
    return [item_path(ws, k) for k in FILES] + [ws / "work" / "pgx" / "pgx.json"]


def authored(ws: Path) -> List[tuple]:
    """Prose inside registered items that a model wrote (AI judgements): traced like every other authored text."""
    out = []
    reg = registry(ws)
    for kind, items in reg.items():
        txt = [str(it.get(k) or "") for it in items.values() if (it.get("label") or {}).get("basis_kind") == "ai_estimate"
               for k in ("basis_zh", "condition_zh", "text_zh", "drug_or_class_zh", "intervention_zh")] + \
              [a for it in items.values() if (it.get("label") or {}).get("basis_kind") == "ai_estimate" for a in it["label"].get("assumptions", [])]
        if any(txt):
            out.append((str(item_path(ws, kind)), "\n".join(txt)))
    return out


def problems(st: Dict[str, Any], ws: Path, inputs: Optional[Dict[str, str]] = None) -> List[str]:
    """Why the labelled items may not be rendered: changed or stale files, label problems, untraced numbers."""
    out: List[str] = []
    reg = registry(ws)
    for kind in FILES:
        p = item_path(ws, kind)
        rec = _state_rec(st, kind)
        cmd = {"p": "la.py predict response", "dx": "la.py clinical dx", "rx": "la.py clinical rx"}[kind]
        if p.exists() and (not rec or sha256_file(p, limit=None) != rec.get("sha256")):
            out.append(f"{p.name} does not match what `{cmd}` wrote; run it again")
            continue
        if rec and not p.exists():
            out.append(f"{p.name} is missing; run `{cmd}` again")
            continue
        if rec and inputs is not None and inputs.get(kind) is not None and rec.get("inputs") != inputs[kind]:
            out.append(f"{p.name} was made from inputs that have changed since; run `{cmd}` again")
        for iid, it in reg[kind].items():
            out += [f"{kind}:{iid}: {b}" for b in item_problems(kind, it)]
    rules = _rx_rule_doses()
    pgx_doses = _pgx_doses(ws)
    for iid, it in reg["rx"].items():
        dose = it.get("dose_text_zh")
        if not dose:
            continue
        bk = (it.get("label") or {}).get("basis_kind")
        if bk == "guideline_recommendation" and dose != rules.get(it.get("rule_id")):
            out.append(f"rx:{iid}: dose text is not the rule's verbatim dose text")
        elif bk == "cpic" and dose not in pgx_doses:
            out.append(f"rx:{iid}: dose text is not CPIC's own text in work/pgx/pgx.json")
        elif bk not in ("guideline_recommendation", "cpic", "ai_estimate"):
            out.append(f"rx:{iid}: a dose may come only from a guideline, CPIC or an AI judgement labelled as such")
    return out


def item_problems(kind: str, it: Dict[str, Any]) -> List[str]:
    lab = it.get("label")
    bad = label_problems(lab)
    if bad:
        return bad
    bk = lab["basis_kind"]
    if kind == "p":
        if bk == "ai_estimate":
            for k in ("p_detect_favorable", "p_true_favorable"):
                v = it.get(k)
                if v is not None and not (isinstance(v, list) and len(v) == 2 and 0 <= v[0] <= v[1] <= 1):
                    return [f"an AI prediction gives {k} as a [low, high] interval, never a point"]
            if it.get("p_detect_favorable") is None and it.get("p_true_favorable") is None:
                return ["an AI prediction needs a probability interval"]
            if it.get("expected_change_pct") is not None or it.get("interval_80_pct") is not None:
                return ["an AI prediction gives no point change or model interval"]
        elif bk == "statistical_model":
            if lab["confidence"]["interval_80"] != it.get("interval_80_pct"):
                return ["label interval differs from the item's interval_80_pct"]
            for k in ("p_detect_favorable", "p_true_favorable"):
                if not isinstance(it.get(k), (int, float)) or not 0 <= it[k] <= 1:
                    return [f"{k} must be a probability"]
        else:
            return [f"a prediction is a statistical model estimate or an AI judgement, not {bk}"]
    if kind == "dx" and bk == "ai_estimate":
        pr = it.get("probability_range")
        if not (isinstance(pr, list) and len(pr) == 2 and 0 <= pr[0] <= pr[1] <= 1):
            return ["an AI diagnostic judgement needs probability_range [low, high]"]
    if kind == "rx" and it.get("action") not in ("consider_start", "consider_stop", "consider_switch", "adjust", "avoid", "monitor"):
        return [f"action {it.get('action')!r} unknown"]
    return []


def _rx_rule_doses() -> Dict[str, str]:
    from .common import DATA_DIR
    p = DATA_DIR / "medication_rules.json"
    if not p.exists():
        return {}
    return {r["id"]: r.get("dose_text_zh") or "" for r in load_json(p).get("rules", [])}


def _pgx_doses(ws: Path) -> set:
    p = ws / "work" / "pgx" / "pgx.json"
    if not p.exists():
        return set()
    try:
        d = load_json(p)
        return {str(x.get(k)) for r in d.get("results", []) for x in r.get("drugs", [])
                for k in ("dose_text", "recommendation_text") if x.get(k)}
    except (OSError, ValueError, AttributeError):
        return set()


def digest(obj: Any) -> str:
    return sha_text(json.dumps(obj, ensure_ascii=False, sort_keys=True, default=str))


def check(st: Dict[str, Any], ws: Path) -> List[str]:
    """`problems` with each file's inputs re-derived now (a file made from inputs that have changed is stale)."""
    from . import clinical, predict
    inputs = {"p": predict.inputs_digest(st, ws) if (st.get("predict") or {}).get("sha256") else None,
              "dx": clinical.inputs_digest(st, ws, "dx") if ((st.get("clinical") or {}).get("dx") or {}).get("sha256") else None,
              "rx": clinical.inputs_digest(st, ws, "rx") if ((st.get("clinical") or {}).get("rx") or {}).get("sha256") else None}
    return problems(st, ws, inputs)
