"""Assemble the member report and check that its numbers are bound to results.

The agent writes prose only (summary + per-system analyses + plan text) and
refers to member-specific numbers through {{r:<readout id>}} placeholders.
This module fills tables from readouts.json / plan.json / state.json,
substitutes the placeholders, and `trace` refuses authored text that carries
a number which did not come from a result file.

Approvals are bound to content: `trace` records the sha256 of every file it
checked, `review record` binds the reviewer's file, and `render` refuses if
any of them changed afterwards.
"""
from __future__ import annotations

import html
import json
import re
import subprocess
import unicodedata
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .common import (EXIT_BLOCKED, EXIT_INPUT, LAError, data, load_json, longevity_skills_home,
                     now_iso, sha256_file, write_json)
from .common import sha_text as C_sha
from .integrate import MODES, PLACEHOLDER

KIND_ZH = {"computed": "为你计算", "computed_low_coverage": "为你计算（CpG 覆盖不足，仅供参考）", "implausible": "数值超出合理范围，不采用", "computed_quality_unverified": "为你计算（原始数据缺少质量字段，未核对）",
           "llm_estimate": "AI 估计（非测量、非校准模型）",
           "genetic_score": "遗传倾向（全基因组显著位点风险等位基因计数）", "genetic_finding": "ClinVar 注释",
           "population_position": "人群位置（参照人群见说明）",
           "computed_coverage_unknown": "为你计算（覆盖率未知）", "descriptive": "数据描述", "lookup": "论文名单查表"}
SEX_ZH = {"male": "男", "female": "女"}
BOUNDARY = ("本报告由模型和已发表方法根据你交来的数据计算，是研究性估计，不是诊断，也不是治疗或用药建议。"
            "生物年龄类读数有测量误差，单次结果不足以下结论；是否干预、如何用药由你的医生决定。")

# Authored prose may carry digits ONLY inside placeholders:
#   {{r:<readout id>[|label|value|unit|own words]}}  a member readout (bound to readouts.json)
#   {{n:<literal>}}                                   a non-member literal (IL-6, 每周 3~5 次, GRCh38) - listed for the reviewer
#   {{pmid:<digits>}}                                 a citation, rendered as a PubMed link
# Outside them: no digit of any script, no run of Chinese numerals, no raw URL. Nothing is allow-listed by pattern.
LITERAL = re.compile(r"\{\{n:([^{}]{1,80})\}\}")
PMID_PH = re.compile(r"\{\{pmid:(\d{1,9})\}\}")
CN_NUM = "零一二三四五六七八九十百千万两俩仨〇壹贰叁肆伍陆柒捌玖拾佰仟廿卅卌〡〢〣〤〥〦〧〨〩"
CN_RUN = re.compile(rf"[{CN_NUM}]{{2,}}")
CN_WORDS = {"十分", "一一", "万一", "千万", "三三两两"}   # whole runs that are ordinary words ("十分重要", "千万别"), never values
CN_SINGLE_UNIT = re.compile(rf"[{CN_NUM}]\s*(?:周岁|岁|%|公斤|千克|克|毫克|mmHg|年|个月|月|周|天|小时|分钟|次|步|公里|千米|粒|片|成|倍|分之|分|期|点[{CN_NUM}])|[GＧ]\s*[{CN_NUM}]")
CN_OK_PHRASES = ("一点十分", "一成不变", "一片空白", "十二指肠", "三七粉", "三七", "四分位", "五分位", "十分位", "百岁", "七八分饱",
                 "七分饱", "八分饱", "千千万万", "一日三餐", "三餐", "这一期间", "一期间", "逐一", "唯一", "统一", "单一", "同一", "万一",
                 "五十肩", "早一天", "晚一天", "多一天", "一一", "万万", "三天两头", "百分比", "百分位", "百分点", "百分数", "千克", "千米", "千卡", "千焦")   # ordinary words; blanked before the numeral and dose checks
_EN = r"(?:one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty|thirty|forty|fifty|sixty|seventy|eighty|ninety|hundred)"
EN_NUM = re.compile(rf"\b{_EN}(?:[\s-]+{_EN})*\s*(?:岁|周岁|years?\b|%|percent\b|倍|mmHg)", re.I)


def _blank_ok_phrases(t: str, keep: Tuple[str, ...] = ()) -> str:
    """Blank ordinary words (唯一, 七八分饱 …) - but never one touching another numeral (活到一百岁, 第九十分位, 一三七)."""
    for ph in sorted((x for x in CN_OK_PHRASES if x not in keep), key=len, reverse=True):
        t = re.sub(rf"(?<![{CN_NUM}]){re.escape(ph)}(?![{CN_NUM}])", "□" * len(ph), t)
    return t


ROMAN = re.compile("[\u2160-\u2188]")                               # Ⅰ Ⅱ Ⅲ … ⅻ, checked before NFKC folds them to letters
COMPARE = "大小高低多少快慢早晚老轻差超过满"
CN_OK_UNIT = {"一年", "一分", "十分", "一点", "一天", "一周", "一次", "一步"}   # "一年四季"、"十分重要"、"一点点"、"一次性" - only standalone
ZW = re.compile(r"[\u200b-\u200f\u2060-\u2064\ufeff\u00ad\u034f\u180e]")
URL = re.compile(r"(?:https?|ftp|hxxps?)\s*\\?:\s*/\s*/|www\.|(?<![:\w])//[\w.-]+|"
                 r"\b[a-z0-9-]+(?:\.[a-z0-9-]+)*\.(?:com|cn|net|org|io|gov|edu|info|co|me|app|xyz|top|ly)\b(?:/|\b)", re.I)


def _readouts(ws: Path) -> Dict[str, Dict[str, Any]]:
    """Method readouts plus registered organ AI estimates (their own kind, never mixed into readouts.json)."""
    ro = {r["id"]: r for r in load_json(ws / "work" / "readouts.json")["readouts"]}
    op = ws / "work" / "organs" / "organ_readouts.json"
    if op.exists():
        ro.update({r["id"]: r for r in load_json(op)["readouts"]})
    ip = ws / "work" / "insights" / "insight_readouts.json"
    if ip.exists():
        ro.update({r["id"]: r for r in load_json(ip)["readouts"]})
    return ro


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".") if abs(v) < 1000 else f"{v:.0f}"
    return str(v)


def substitute(text: str, ro: Dict[str, Dict[str, Any]]) -> str:
    def rep(m: re.Match) -> str:
        r = ro.get(m.group(1))
        if r is None:
            return f"[缺失读数 {m.group(1)}]"
        what = (m.group(2) or "value").strip()
        if what not in MODES:
            return what                      # the writer's own wording for the readout's name
        if what == "label":
            return r["label_zh"]
        if what == "unit":
            return r.get("unit", "")
        u = r.get("unit", "")
        if r.get("kind") == "llm_estimate":        # never readable as a measurement, even when cited bare
            if u == "概率":
                return f"{_pct(r['value'])}（AI 估计，区间 {_pct(r.get('low'))}–{_pct(r.get('high'))}，{r.get('horizon_years')} 年）"
            return f"{_fmt(r['value'])} 岁（AI 估计，区间 {_fmt(r.get('low'))}–{_fmt(r.get('high'))} 岁）"
        u = {"a": " 岁", "1": ""}.get(u, f" {u}" if u else "")
        return f"{_fmt(r['value'])}{u}"
    text = PLACEHOLDER.sub(rep, text)
    text = LITERAL.sub(lambda m: m.group(1), text)
    return PMID_PH.sub(lambda m: f"[PMID {m.group(1)}](https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}/)", text)


def normalise(text: str) -> str:
    """What a reader will see: entities decoded, full-width/superscript/circled digits folded, link targets dropped."""
    t = html.unescape(html.unescape(text))
    t = unicodedata.normalize("NFKC", t)
    t = re.sub(r"\[([^\]]*)\]\((https?://pubmed\.ncbi\.nlm\.nih\.gov/\d+/?|https?://doi\.org/[^)]*)\)", r"\1", t)
    t = re.sub(r"PMID\s*:?\s*\d+", " ", t)
    t = re.sub(r"\((https?://[^)]*)\)", " ", t)
    return t


def trace_text(text: str, allowed_values: List[str] = None) -> List[str]:
    """Violations: any digit, numeral run or URL outside the three placeholder kinds."""
    t = PLACEHOLDER.sub(lambda m: f" {m.group(2)} " if m.group(2) and m.group(2).strip() not in MODES else " ", text)
    t = LITERAL.sub(" ", t)
    t = PMID_PH.sub(" ", t)
    t0 = html.unescape(html.unescape(t))
    roman = [m.group(0) for m in ROMAN.finditer(t0)]
    t = unicodedata.normalize("NFKC", t0)
    t = "".join(ch for ch in t if unicodedata.category(ch) not in ("Cf", "Mn", "Me"))   # invisible format / combining marks
    t = re.sub(r"[*_`~]+", "", ZW.sub("", t))                                          # markdown emphasis cannot split a number
    t = _blank_ok_phrases(t)
    t = re.sub(rf"(?<=[{CN_NUM}])\s+(?=[{CN_NUM}])", "", t)        # "八 四", "陆 拾 柒"
    bad = [f"…{r}…  ← 罗马数字（分期、分级也要用占位符或 {{{{n:…}}}}）" for r in roman]
    for m in re.finditer(r"\d+", t):
        s = max(0, m.start() - 10)
        bad.append(f"…{t[s:m.end() + 10].strip()}…  ← `{m.group(0)}`（非会员数字写成 {{{{n:…}}}}，名称里的数字也一样，如 {{{{n:维生素B12}}}}、{{{{n:omega-3}}}}；会员数值用 {{{{r:…}}}}）")
    for m in CN_RUN.finditer(t):
        if m.group(0) not in CN_WORDS and not any(t.startswith(ph, m.start()) or t[max(0, m.start() - 4):m.end() + 4].find(ph) >= 0
                                                 and ph.find(m.group(0)) >= 0 for ph in CN_OK_PHRASES):
            bad.append(f"…{m.group(0)}…  ← 中文数字（写成 {{{{n:…}}}} 或占位符）")
    for m in CN_SINGLE_UNIT.finditer(t):
        standalone = (m.start() == 0 or t[m.start() - 1] not in CN_NUM)
        g = m.group(0).replace(" ", "")
        if m.start() > 0 and t[m.start() - 1] == "第" and g.endswith("步"):
            continue                                                      # "第二步": an ordinal step, not a quantity
        compared = (m.start() > 0 and t[m.start() - 1] in COMPARE) or t[m.end():m.end() + 1] == "半"
        if standalone and not compared and (g in CN_OK_UNIT or t.startswith(CN_OK_PHRASES, m.start())):
            continue                                                      # "这一点十分重要", idioms; not "大一年", "一年半"
        bad.append(f"…{m.group(0)}…  ← 中文数字加单位（写成 {{{{n:…}}}} 或占位符）")
    for m in EN_NUM.finditer(t):
        bad.append(f"…{m.group(0)}…  ← 英文数字（写成占位符或 {{{{n:…}}}}）")
    for m in re.finditer(rf"[{COMPARE}]半\s*(?:年|个月|月|岁)|[{CN_NUM}]\s*(?:来岁|出头|旬)|[{COMPARE}][了得]?[{CN_NUM}]\s*(?:年|岁|个月|周|天)", t):
        bad.append(f"…{m.group(0)}…  ← 数量（写成 {{{{n:…}}}}）")
    for m in URL.finditer(t):
        bad.append("…链接…  ← 文中不放链接；文献用 {{pmid:N}}")
    for lit in literals(text):
        if URL.search(lit):
            bad.append(f"…{{{{n:{lit[:30]}}}}}…  ← 字面量里不能放链接")
    rendered = ZW.sub("", unicodedata.normalize("NFKC", LITERAL.sub(lambda m: m.group(1), PMID_PH.sub(" ", PLACEHOLDER.sub("x", text)))))
    if URL.search(rendered) or re.search(r"\]\s*\(", rendered):
        bad.append("…链接…  ← 拼接后会形成链接；文中不放链接")
    from .evidence import has_dose_outside_literals
    dose_text = _blank_ok_phrases(ZW.sub("", text), keep=("三七粉", "三七"))      # 三七 is itself a supplement
    if has_dose_outside_literals(dose_text):
        bad.append("正文出现剂量或粒数：药物和补充剂剂量只由医生给出；膳食和生活方式用量（盐、饮水、牛奶等）写成 {{n:…}}")
    return bad


APOE_TOPIC = re.compile(r"a[\s\-.]*p[\s\-.]*o[\s\-.]*e\b|apo\s*-?\s*e\d|apolipoprotein\s*-?\s*e|載脂蛋白|载脂蛋白\s*e|ε\s*-?\s*[234二三四]|"
                        r"epsilon\s*-?\s*(?:[234]|two|three|four)|\be\s*-?\s*(?:two|three|four|2|3|4)\b|rs\s*429358|rs\s*7412|"
                        r"ε\s*型?\s*[二三四]\s*号?|[ⅡⅢⅣ]\s*型\s*载脂|four\s+copies", re.I)
# a sentence naming dementia together with genes
_DEM = r"(?:阿尔茨海默|阿尔兹海默|阿兹海默|老年痴呆|痴呆|失智|认知障碍|认知下降|认知衰退|认知功能下降|认知功能减退|记忆力减退|脑退化|(?<![A-Za-z])AD(?![A-Za-z])|alzheimer|dementia)"
APOE_PERSONAL = re.compile(rf"{_DEM}.{{0,14}}(?:基因|等位|携带|易感|遗传|gene|allele|carrier)|"
                           rf"(?:基因|遗传|等位|gene|allele|carrier).{{0,14}}{_DEM}|风险等位基因|易感基因", re.I)


def apoe_mentioned(text: str) -> bool:
    t = _apoe_text(text)
    if APOE_TOPIC.search(t):
        return True
    # declined: no sentence may tie dementia to genes, not even as general knowledge or next to "you chose not to see"
    t = re.sub(r"表观遗传|epigenetic", "□", t, flags=re.I)          # methylation clocks are not genetics
    return any(APOE_PERSONAL.search(sen) for sen in re.split(r"[。！？!?\n；;.]", t))
HOMOGLYPH = str.maketrans({"ᴀ": "A", "ᴘ": "P", "ᴏ": "O", "ᴇ": "E", "ɛ": "ε", "ϵ": "ε", "Ɛ": "ε", "Ε": "E", "Ο": "O", "Ρ": "P", "Α": "A", "А": "A", "Р": "P", "О": "O", "Е": "E",
                           "а": "a", "р": "p", "о": "o", "е": "e", "ο": "o", "ρ": "p", "α": "a"})


def _apoe_text(text: str) -> str:
    return text.translate(HOMOGLYPH) + "\n" + ZW.sub("", unicodedata.normalize("NFKC", text.translate(HOMOGLYPH))).translate(HOMOGLYPH)


def literals(text: str) -> List[str]:
    return [m.group(1) for m in LITERAL.finditer(text)]


def registered_analyses(ws: Path, st: Dict[str, Any] = None) -> List[Path]:
    if st is None:
        sp = ws / "state.json"
        st = load_json(sp) if sp.exists() else {}
    reg = (st.get("integrate") or {}).get("analyses") or {}
    return [Path(v["path"]) for k, v in sorted(reg.items()) if Path(v["path"]).exists()]


VOLATILE = ("events", "updated_at", "review", "report", "review_history", "stages", "review_blocks")


def state_digest(st: Dict[str, Any]) -> str:
    """The whole state except bookkeeping; any change after trace voids trace, review and report."""
    view = {k: v for k, v in st.items() if k not in VOLATILE}
    return C_sha(json.dumps(view, ensure_ascii=False, sort_keys=True, default=str))


def _authored(ws: Path, st: Dict[str, Any] = None) -> List[Tuple[str, str]]:
    out = []
    p = ws / "work" / "report" / "summary.md"
    if p.exists():
        out.append((str(p), p.read_text(encoding="utf-8")))
    for a in registered_analyses(ws, st):
        out.append((str(a), a.read_text(encoding="utf-8")))
    op = ws / "work" / "organs" / "organ_readouts.json"
    if op.exists():
        out.append((str(op), "\n".join(x.get("rationale_zh", "") for x in load_json(op)["readouts"])))   # disease names come from organs.json
    bq = ws / "work" / "insights" / "board" / "questions.json"
    if bq.exists():
        out.append((str(bq), "\n".join(q.get(k, "") for q in load_json(bq)["questions"] for k in ("title_zh", "hypothesis_zh", "why_zh"))))
        from .board import registered_findings
        for fp in sorted(registered_findings(ws, st).values()):
            f = load_json(fp)
            out.append((str(fp), "\n".join(str(f.get(k, "")) for k in ("summary_zh", "next_step_zh", "limitations_zh"))))
    gp = ws / "work" / "insights" / "genotype_phenotype.json"
    if gp.exists() and load_json(gp).get("absent_as_ref"):
        out.append(("absent-as-ref reason", str(load_json(gp)["absent_as_ref"])))       # printed in the report
    pj = ws / "work" / "insights" / "projections.json"
    if pj.exists():
        out.append(("exposure-match reasons", "\n".join(str(x.get("exposure_match") or "") for x in load_json(pj)["projections"])))
    if bq.exists():
        skipped = (((st or {}).get("insights") or {}).get("board") or {}).get("skipped") or {}
        if skipped:                                  # skip reasons are printed in the report too
            out.append(("board skip reasons", "\n".join(str(v.get("reason", "")) for v in skipped.values())))
    pl = ws / "work" / "intervene" / "plan.json"
    if pl.exists():
        plan = load_json(pl)
        fields = []
        for i in plan.get("items", []):
            rt = i.get("retest") or {}
            fields += [i.get("action_zh", ""), i.get("rationale_zh", ""), str(rt.get("what", ""))]
        out.append((str(pl), "\n".join(fields)))
    return out


def bound_files(ws: Path, st: Dict[str, Any] = None) -> List[Path]:
    """Everything the rendered report is made from; any change after trace/review voids them."""
    files = [ws / "work" / "readouts.json", ws / "work" / "report" / "summary.md", ws / "work" / "intervene" / "plan.json",
             ws / "work" / "twin" / "twin.json", ws / "work" / "organs" / "organ_readouts.json"]
    files += sorted((ws / "work" / "organs" / "estimates").glob("*.json"))
    files += [ws / "work" / "insights" / n for n in ("insight_readouts.json", "genotype_phenotype.json", "positions.json", "projections.json")]
    bq = ws / "work" / "insights" / "board" / "questions.json"
    if bq.exists():
        from .board import registered_findings
        files += [bq] + sorted(registered_findings(ws, st).values())
    files += registered_analyses(ws, st)
    return [f for f in files if f.exists()]


def _hashes(ws: Path, st: Dict[str, Any] = None) -> Dict[str, str]:
    h = {str(p): sha256_file(p, limit=None) for p in bound_files(ws, st)}
    h["state"] = state_digest(st if st is not None else load_json(ws / "state.json"))
    return h


def insight_file_problems(ws: Path, st: Dict[str, Any]) -> List[str]:
    """The insight tables are rendered from files `la.py insights` wrote; an edited or hand-made copy is refused."""
    out = []
    files = (st.get("insights") or {}).get("files") or {}
    for name in ("genotype_phenotype.json", "positions.json", "projections.json"):
        fp = ws / "work" / "insights" / name
        if name in files and (not fp.exists() or sha256_file(fp, limit=None) != files[name]):
            out.append(f"{name} does not match what `la.py insights` wrote; run that step again")
        elif name not in files and fp.exists():
            out.append(f"{name} was not written by `la.py insights`; run that step again")
    return out


def trace(st: Dict[str, Any], ws: Path) -> Dict[str, Any]:
    allowed = [str(st["member"].get("age", ""))] if st["member"].get("age") else []
    rp = ws / "work" / "readouts.json"
    want = (st.get("methods") or {}).get("readouts_sha256")
    if not want or sha256_file(rp, limit=None) != want:
        raise LAError("readouts.json does not match what `methods run` wrote; re-run methods", EXIT_BLOCKED)
    op = ws / "work" / "organs" / "organ_readouts.json"
    ow = (st.get("organs") or {}).get("organ_readouts_sha256")
    if op.exists() and (not ow or sha256_file(op, limit=None) != ow):
        raise LAError("organ_readouts.json does not match what `organ register` wrote; register the estimates again", EXIT_BLOCKED)
    ip = ws / "work" / "insights" / "insight_readouts.json"
    iw = (st.get("insights") or {}).get("readouts_sha256")
    if ip.exists() and (not iw or sha256_file(ip, limit=None) != iw):
        raise LAError("insight_readouts.json does not match what `la.py insights` wrote; run it again", EXIT_BLOCKED)
    from .board import integrity
    bad_files = insight_file_problems(ws, st)
    if bad_files:
        raise LAError("; ".join(bad_files), EXIT_BLOCKED)
    bad_board = integrity(ws, st)
    if bad_board:
        raise LAError("; ".join(bad_board), EXIT_BLOCKED)
    ro = _readouts(ws)
    findings = {}
    missing_ids = {}
    lits: Dict[str, List[str]] = {}
    declined = (st["member"].get("answers") or {}).get("genetic_disclosure") == "no"
    for path, text in _authored(ws, st):
        bad = trace_text(text, allowed)
        if declined and apoe_mentioned(LITERAL.sub(lambda m: m.group(1), text)):
            bad.append("会员选择不查看 APOE 结果：正文不得提及 APOE 基因型，也不得把痴呆与基因、遗传写在同一句里")
        if bad:
            findings[path] = bad
        if literals(text):
            lits[path] = literals(text)
        miss = sorted({m.group(1) for m in PLACEHOLDER.finditer(text) if m.group(1) not in ro})
        if miss:
            missing_ids[path] = miss
    cited = sorted({m.group(1) for _, text in _authored(ws, st) for m in PMID_PH.finditer(text)})
    if cited:
        from .evidence import verify_pmids
        live = verify_pmids(cited)           # network failure raises: no unverified citation reaches the report
        dead = [p for p in cited if p not in live]
        if dead:
            findings["citations"] = [f"PMID {p} does not exist or is retracted" for p in dead]
    ok = not findings and not missing_ids
    hashes = _hashes(ws, st) if ok else {}
    trace_id = C_sha(json.dumps(hashes, sort_keys=True))[:16] if ok else None
    from .organs import check_labs, table as organ_table
    organ_checks = [{"organ": r["organ"],
                     "calibrated_or_measured": [f"{x['label_zh']}: {x['value']} {x.get('unit', '')}".strip() for x in r["measured"] + r["indices"]]
                                               + check_labs(st, r["organ"]),
                     "ai_estimates": [f"{x['label_zh']}: {x['value']} [{x.get('low')}–{x.get('high')}]" for x in
                                      ([r["ai_age"]] if r["ai_age"] else []) + r["ai_risks"]]}
                    for r in organ_table(ws, ro) if r["ai_age"] or r["ai_risks"]]
    res = {"ok": ok, "trace_id": trace_id, "unbound_numbers": findings, "unknown_readouts": missing_ids,
           "literals_for_reviewer": lits, "organ_checks_for_reviewer": organ_checks, "pmids_cited": cited, "checked_at": now_iso(), "hashes": hashes}
    write_json(ws / "work" / "review" / "trace.json", res)
    st.setdefault("review", {})["trace"] = {"ok": ok, "at": now_iso(), "hashes": res["hashes"], "trace_id": trace_id}
    st["review"].pop("reviewer", None)           # a new trace needs a new review
    return res


def record_review(st: Dict[str, Any], ws: Path, verdict: str, findings_path: Path) -> Dict[str, Any]:
    tr = (st.get("review") or {}).get("trace") or {}
    if not tr.get("ok"):
        raise LAError("run `la.py review trace` until it passes, then dispatch the reviewer", EXIT_BLOCKED)
    fp = findings_path.resolve()
    if (ws / "work" / "review").resolve() not in fp.parents:
        raise LAError("the reviewer's JSON must be saved under <workspace>/work/review/", EXIT_INPUT)
    try:
        doc = json.loads(fp.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise LAError(f"{findings_path} is not the reviewer's JSON ({e})", EXIT_INPUT)
    if not isinstance(doc, dict) or doc.get("verdict") not in ("pass", "block") \
            or not isinstance(doc.get("findings"), list) or not all(isinstance(x, dict) for x in doc["findings"]):
        raise LAError("reviewer JSON must have verdict pass|block, trace_id and a list of finding objects "
                      "(pass_with_fixes is not a final verdict: apply the fixes, re-trace, review again)", EXIT_INPUT)
    if doc.get("trace_id") != tr.get("trace_id"):
        raise LAError(f"the reviewer's trace_id {doc.get('trace_id')!r} is not the current trace {tr.get('trace_id')!r}; "
                      "the review must be of exactly the traced files", EXIT_INPUT)
    if doc["verdict"] != verdict:
        raise LAError(f"--verdict {verdict} does not match the reviewer's own verdict {doc['verdict']}", EXIT_INPUT)
    blocked = st.setdefault("review_blocks", [])
    if verdict == "pass" and tr["trace_id"] in blocked:
        raise LAError("this exact content was blocked by a reviewer; change it, re-trace and review again", EXIT_INPUT)
    badsev = [i for i, f in enumerate(doc["findings"])
              if unicodedata.normalize("NFKC", str(f.get("severity", ""))).strip().upper() not in ("P0", "P1", "P2")]
    if badsev:
        raise LAError(f"findings {badsev}: every finding needs \"severity\": \"P0\" | \"P1\" | \"P2\" (references/reviewer.md)", EXIT_INPUT)
    worst = {unicodedata.normalize("NFKC", str(f["severity"])).strip().upper() for f in doc["findings"]}
    if verdict == "pass" and worst & {"P0", "P1"}:
        raise LAError("the reviewer listed P0/P1 findings; the verdict cannot be pass", EXIT_INPUT)
    if tr.get("hashes") != _hashes(ws, st):
        raise LAError("files changed after the trace; re-run `la.py review trace` and review again", EXIT_BLOCKED)
    if verdict == "block":
        blocked.append(tr["trace_id"])
    st["review"]["reviewer"] = {"verdict": verdict, "findings": str(fp), "findings_sha256": sha256_file(fp, limit=None),
                                "at": now_iso(), "hashes": tr["hashes"], "trace_id": tr["trace_id"]}
    return st["review"]


def check_bound(st: Dict[str, Any], ws: Path) -> List[str]:
    rv = st.get("review") or {}
    problems = []
    excluded = {f["id"] for f in st["files"] if f.get("excluded")}
    for mod, lst in st["processed"].items():
        if any(set(x.get("source_file_ids") or []) & excluded for x in lst):
            problems.append(f"{mod} data still derives from a file excluded as not this member's")
    if not (rv.get("trace") or {}).get("ok"):
        problems.append("trace has not passed")
    rev = rv.get("reviewer")
    if not rev:
        problems.append("no independent reviewer verdict recorded")
    else:
        if rev["verdict"] == "block":
            problems.append("the reviewer blocked this report")
        if rev.get("hashes") != _hashes(ws, st):
            problems.append("report inputs changed after trace/review")
        if not Path(rev["findings"]).exists() or sha256_file(Path(rev["findings"]), limit=None) != rev.get("findings_sha256"):
            problems.append("the reviewer's findings file changed after it was recorded")
    return problems


def _lib_version() -> str:
    home = longevity_skills_home(required=False)
    if not home:
        return "unknown"
    try:
        v = (home / "VERSION").read_bytes()[:200].decode("utf-8", "replace").strip() if (home / "VERSION").exists() else "?"
    except OSError:
        v = "?"
    try:
        c = subprocess.run(["git", "-c", "core.fsmonitor=false", "-C", str(home), "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        c = ""
    v = re.sub(r"[^0-9A-Za-z.\-]", "", v)[:20]
    c = re.sub(r"[^0-9a-f]", "", c)[:12]
    return f"{v} ({c})" if c else v


def _cell(x: Any) -> str:
    s = html.escape(str(x), quote=False).replace("|", "／").replace("\n", " ")
    return re.sub(r"([\[\]()!*_`#<>])", r"\\\1", s)     # markdown in a cell is shown, never rendered


METHOD_STATUS_ZH = {
    "blocked_platform": "检测平台、组织或数据量纲与该方法的训练数据不一致，按规则不计算",
    "blocked_license": "授权原因，商用模式下不运行",
    "blocked_input_origin": "需要外部算好的年龄值作输入，来源无法核实，商用模式下不运行",
    "needs_answers": "缺少需要本人回答的问题",
    "input_problem": "化验项目缺失或单位无法换算",
    "failed": "运行出错，未得到结果",
    "no_output": "运行完成，但该方法对本人数据不产生数值",
    "declined": "本人选择不在报告中查看此项结果",
    "not_run_excluded": "对应数据文件判定不属于本人，未计算",
    "manual": "输入格式尚未在方法库中声明，本期不自动运行",
}
FILE_STATUS_ZH = {
    "unsupported": "本版本尚不支持该格式", "rejected": "文件内容与其声称的数据类型不符，已拒收（原因见工作区记录）",
    "needs_pipeline": "需要本地生信流程，本期未运行", "needs_judgment": "有待确认的信息",
    "deferred_ask_lab": "检测类型/组织/平台未知，待检测机构说明", "excluded_not_member": "判定不属于本人，已排除",
    "unreadable": "文件损坏或无法读取",
}
PIPE_STATUS_ZH = {"skipped": "未在本机运行（用户决定或本机条件不足）", "void_excluded": "输入文件判定不属于本人，流程作废",
                  "void_stale": "已不在当前计划中", "planned": "已规划，未运行", "stub_ok": "仅做过接线测试，不是结果",
                  "verify_failed": "运行结果未通过校验", "failed_hanging": "运行失败"}


def _prose(text: str) -> str:
    """Authored markdown: raw HTML shown as text; images dropped; only http(s) links kept."""
    text = text.replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"!\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"\[([^\]]*)\]\((?!https?://)[^)]*\)", r"\1", text)
    text = re.sub(r"^\s*\[[^\]]+\]:\s*\S+.*$", "", text, flags=re.M)
    return text


def render(st: Dict[str, Any], ws: Path) -> Dict[str, Any]:
    problems = check_bound(st, ws)
    if problems:
        raise LAError("cannot render: " + "; ".join(problems), EXIT_BLOCKED)
    ro = _readouts(ws)
    sysdef = data("systems.json")["systems"]
    m = st["member"]
    L: List[str] = [f"# 抗衰多组学分析报告 · {_cell(m.get('id'))}", "",
                    f"性别：{SEX_ZH.get(m.get('sex'), '未提供')}　年龄：{m.get('age') or '未提供'} 岁　采样日期：{m.get('sample_date') or '未提供'}　"
                    f"生成时间：{now_iso()[:16]}　模式：{'商用（授权清单待法务确认）' if m.get('mode') == 'commercial' else '内部研究'}", ""]
    ident = st.get("identity") or {}
    excluded = [f for f in st["files"] if f.get("excluded")]
    uncertain = [f for f in st["files"] if f.get("provenance_uncertain") and not f.get("excluded")]
    if uncertain:
        L += [f"> **数据来源提示**：{len(uncertain)} 个文件无法确认属于本人，相关读数在总表中已标注，不作为个人结论的依据。", ""]
    if ident.get("answer") != "consistent" or excluded:
        n_all = len(st["files"])
        what = (f"{len(excluded)}/{n_all} 个文件判定不属于本人，已排除，不参与任何计算" if excluded else
                {"cannot_tell": "无法确认所有文件都属于本人"}.get(ident.get("answer"), "未做数据归属核对"))
        L += [f"> **数据归属提示**：{what}。明细见「本期没有给出结果的项目」。", ""]
    summ = ws / "work" / "report" / "summary.md"
    if summ.exists():
        L += ["## 摘要", "", _prose(substitute(summ.read_text(encoding="utf-8"), ro)).strip(), ""]
    if ro:
        L += ["## 本期读数总表", "", "类型说明：「为你计算」= 用你的数据按论文公式算出；「数据描述」= 质控或计数，不是健康判断。", ""]
        by_sys: Dict[str, List[Dict[str, Any]]] = {}
        for r in ro.values():
            if r.get("kind") in ("llm_estimate", "genetic_score", "genetic_finding", "population_position"):
                continue                              # shown in their own sections (organ table, genes vs labs, position)
            by_sys.setdefault((r.get("systems") or ["overall_aging"])[0], []).append(r)
        for s, rs in sorted(by_sys.items(), key=lambda kv: list(sysdef).index(kv[0]) if kv[0] in sysdef else 99):
            L += [f"### {sysdef.get(s, {}).get('label_zh', s)}", "", "| 读数 | 值 | 单位 | 类型 | 方法 |", "|---|---|---|---|---|"]
            for r in rs:
                cov = f"（CpG 覆盖 {r['coverage_pct']:.0f}%）" if r.get("coverage_pct") is not None else ""
                unc = "（数据来源未确认属于本人）" if r.get("provenance_uncertain") else ""
                unit = {"a": "岁"}.get(r.get("unit", ""), r.get("unit", ""))
                L.append(f"| {_cell(r['label_zh'])} | {_cell(_fmt(r['value']))} | {_cell(unit)} | {KIND_ZH.get(r.get('kind'), r.get('kind'))}{cov}{unc} | {_cell(r.get('method'))} |")
            L.append("")
    L += organ_section(ws, ro)
    L += insight_sections(ws, ro)
    regd = registered_analyses(ws, st)
    if regd:
        L += ["## 分系统解读", ""]
        for a in regd:
            body = _prose(substitute(a.read_text(encoding="utf-8"), ro)).replace("\n## ", "\n#### ")
            L += [f"### {sysdef.get(a.stem, {}).get('label_zh', a.stem)}", "", body.strip(), ""]
    pl = ws / "work" / "intervene" / "plan.json"
    if pl.exists():
        cats = {c["id"]: c["label_zh"] for c in data("intervention_menu.json")["categories"]}
        L += ["## 干预方案（营养师审核后生效）", "", "| # | 类别 | 做什么 | 针对 | 执行人 | 证据 | 复测 |", "|---|---|---|---|---|---|---|"]
        for i in load_json(pl)["items"]:
            tg = "、".join(ro[t]["label_zh"] if t in ro else t for t in i["targets"])
            ev = "；".join(_ev_link(e) for e in i.get("evidence", [])) or "—"
            rt = i.get("retest") or {}
            what = substitute(str(rt.get("what") or ""), ro)
            rtxt = f"{_cell(what)}，{int(rt['after_weeks'])} 周后" if rt.get("after_weeks") else (_cell(what) if what else "—")
            L.append(f"| {_cell(i['id'])} | {cats.get(i['category'], _cell(i['category']))} | {_cell(substitute(i['action_zh'], ro))} | {_cell(tg)} | {_cell(i['executor'])} | {ev} | {rtxt} |")
        L.append("")
    notes = [(it["method"], n) for it in (st.get("methods") or {}).get("items", []) if it.get("status") == "done"
             for n in (it.get("result") or {}).get("notes", [])]
    if m.get("mode") == "commercial" and any(it["method"] == "accelerated-biological-aging-risk" and it.get("status") == "done"
                                             for it in (st.get("methods") or {}).get("items", [])):
        notes.append(("表型年龄", "本报告的「表型年龄」由血液化验计算（Levine 2018 临床公式），与因授权未使用的甲基化版 DNAm PhenoAge 是两个不同方法"))
    if notes:
        L += ["## 方法说明", ""] + [f"- {_cell(m)}：{_cell(n)}" for m, n in notes] + [""]
    L += ["## 本期没有给出结果的项目", "", "| 项目 | 原因 |", "|---|---|"]
    manual = []
    for it in (st.get("methods") or {}).get("items", []):
        if it["status"] == "manual":
            manual.append(it["method"])
        elif it["status"] not in ("done",):
            L.append(f"| {_cell(it['method'])} | {METHOD_STATUS_ZH.get(it['status'], '未计算')} |")
    if manual:
        L.append(f"| 其他 {len(manual)} 个已收录方法 | {METHOD_STATUS_ZH['manual']}：{_cell('、'.join(manual))} |")
    for f in st["files"]:
        if f.get("status") in FILE_STATUS_ZH:
            L.append(f"| 文件 {_cell(f['name'])} | {FILE_STATUS_ZH[f['status']]} |")
    for p in st["pipelines"]:
        if p["status"] != "verified":
            L.append(f"| 流程 {_cell(p['pipeline'])} | {PIPE_STATUS_ZH.get(p['status'], '未运行')} |")
    L.append("")
    tw = ws / "work" / "twin" / "twin.json"
    if tw.exists():
        t = load_json(tw)
        L += ["## 数字孪生与复测计划", "",
              f"本期快照记录了 {len(t['observations'])} 项化验、{len(t['readouts'])} 个读数、{len(t['interventions'])} 条拟执行事项。"
              "下次复测后对比：化验项按个体内生物学变异判断变化是否超出正常波动；没有噪声模型的读数只并列展示，不判好坏。", ""]
        if t["retest_plan"]:
            L += ["| 事项 | 复测什么 | 多久后 |", "|---|---|---|"] + [
                f"| {_cell(r['intervention'])} | {_cell(substitute(str(r.get('what') or ''), ro))} | {int(r['after_weeks']) if r.get('after_weeks') else '—'} 周 |" for r in t["retest_plan"]] + [""]
    L += ["## 数据来源与版本", "", "| 文件 | 类型 | 状态 | 校验 |", "|---|---|---|---|"]
    for f in st["files"]:
        L.append(f"| {_cell(f['name'])} | {_cell(f.get('kind'))} | {_cell(f.get('status'))} | {str(f.get('sha256', ''))[:18]} |")
    L += ["", f"方法库 longevity-skills {_lib_version()}；授权模式 {m.get('mode')}。", ""]
    lic = data("license_policy.json")["modes"][m.get("mode", "commercial")]["excluded"]
    if lic:
        L += ["因商用授权未确认而未使用：" + "；".join(_cell(e["component"]) for e in lic), ""]
    L += ["## 边界", "", BOUNDARY, ""]
    md = "\n".join(L)
    out = ws / "deliver"
    out.mkdir(parents=True, exist_ok=True)
    for old in out.iterdir():                 # only files this render writes are delivered
        if old.is_file():
            old.unlink()
    (out / "report.md").write_text(md, encoding="utf-8")
    (out / "report.html").write_text(_html(md, str(m.get("id"))), encoding="utf-8")
    written = ["report.md", "report.html"]
    if tw.exists():
        (out / "twin.json").write_text(tw.read_text(encoding="utf-8"), encoding="utf-8")
        written.append("twin.json")
    from .export import build as build_export              # what a host app (longpi) imports: same data as the report
    write_json(out / "la-export.json", build_export(st, ws, out / "report.html"))
    written.append("la-export.json")
    st["report"] = {"md": str(out / "report.md"), "html": str(out / "report.html"), "at": now_iso(),
                    "hashes": st["review"]["reviewer"]["hashes"],
                    "deliver_sha256": {n: sha256_file(out / n, limit=None) for n in written}}
    return st["report"]


def _pct(x: Any) -> str:
    return f"{100 * float(x):.1f}%".replace(".0%", "%") if isinstance(x, (int, float)) else "—"


def _qual(x: Dict[str, Any]) -> str:
    """The same quality and provenance qualifiers as the readout table, for a readout shown in the organ table."""
    q = []
    if x.get("kind") not in (None, "computed", "llm_estimate"):
        q.append(KIND_ZH.get(x.get("kind"), str(x.get("kind"))))
    if x.get("provenance_uncertain"):
        q.append("数据来源未确认属于本人")
    return f"（{'；'.join(q)}）" if q else ""


def organ_section(ws: Path, ro: Dict[str, Dict[str, Any]]) -> List[str]:
    from .organs import table
    rows = table(ws, ro)
    if not any(r["measured"] or r["indices"] or r["ai_age"] or r["ai_risks"] or r.get("overrides") for r in rows):
        return []
    L = ["## 器官体检表", "",
         "> **读表说明**：「测量/模型」来自已发表方法对你数据的计算；「公式指数/校准模型」来自已发表公式（其中 TyG、AIP 只作描述）；"
         "**「AI 估计」是大模型依据你的数据和文献推断的区间，不是测量，也不是经过校准的风险模型，误差可能很大，只供参考和提示复查方向。**"
         "单次化验得出的分期或界值提示不是诊断。", "",
         "| 器官 | 器官年龄（测量/模型） | 器官年龄（AI 估计） | 公式指数/校准模型 | 疾病风险（AI 估计） |", "|---|---|---|---|---|"]
    for r in rows:
        meas = "；".join(f"{_cell(x['label_zh'])} {_fmt(x['value'])} 岁{_qual(x)}" for x in r["measured"]) or "本期无测量方法"
        a = r["ai_age"]
        aiage = (f"{_fmt(a['value'])} 岁（{_fmt(a['low'])}–{_fmt(a['high'])}，置信度{'低' if a.get('confidence') == 'low' else '很低'}）{_qual(a)}"
                 if a else "—")
        idx = "；".join(f"{_cell(x['label_zh'])} {_cell(_fmt(x['value']))}{(('' if x['unit'] == '%' else ' ') + _cell(x['unit'])) if x.get('unit') else ''}{_qual(x)}"
                        for x in r["indices"]) or "—"
        risks = "；".join([f"{_cell(x['disease'])}：{_pct(x['value'])}（{_pct(x['low'])}–{_pct(x['high'])}，{x['horizon_years']} 年）{_qual(x)}"
                           for x in r["ai_risks"]] + [f"{_cell(o['disease'])}：{_cell(o['message_zh'])}" for o in r.get("overrides", [])]) or "—"
        L.append(f"| {r['label_zh']} | {meas} | {aiage} | {idx} | {risks} |")
    L.append("")
    notes = [x for r in rows for x in r["ai_risks"] if x.get("rationale_zh")]
    if notes:
        L += ["**AI 估计的依据**", ""]
        for x in notes:
            pm = [e for e in x.get("evidence", []) if e.get("type") == "pubmed" and str(e.get("ref", "")).isdigit()]
            L.append(f"- {_cell(x['disease'])}：{_cell(substitute(x['rationale_zh'], ro))}"
                     + (f"（文献：{'；'.join(_ev_link(e) for e in pm)}）" if pm else ""))
        L.append("")
    return L


def insight_sections(ws: Path, ro: Dict[str, Dict[str, Any]]) -> List[str]:
    from . import board as B
    L: List[str] = []
    gp = ws / "work" / "insights" / "genotype_phenotype.json"
    if gp.exists():
        g = load_json(gp)["analytes"]
        if g:
            L += ["## 基因与化验对照", "",
                  "> 遗传倾向 = 你携带的、已发表全基因组关联研究（GWAS Catalog，p<5×10⁻⁸）中使该指标升高的等位基因个数，"
                  "按 gnomAD 东亚人群频率换算成百分位；它说明「这项指标有多少可能受基因影响」，不是诊断。位点覆盖不足时不计算。", "",
                  "| 指标 | 你的结果 | 与参考区间 | 遗传倾向（东亚百分位） | 位点覆盖 | 单基因致病/可能致病变异（ClinVar） |", "|---|---|---|---|---|---|"]
            for k, a in g.items():
                lab = a["member_lab"]
                flag = {"high": "偏高", "low": "偏低", "in_range": "在范围内", "no_range": "报告未印参考区间"}[lab["flag"]]
                pct = f"第 {_fmt(a['pct_eas'])} 百分位" if a.get("pct_eas") is not None else \
                    ("未计算（公共数据未取到）" if a.get("not_retrieved") else "未计算（在东亚人群有变异的位点太少或覆盖不足）")
                ms = a["monogenic_scan"]
                parts = []
                zmap = {"het": "杂合", "hom": "纯合", "hemi": "半合"}
                for x in ms["pathogenic_or_likely"]:
                    tag = ("，隐性遗传，单个杂合仅为携带者" if x.get("carrier_only") else
                           "，隐性遗传，同一基因两个杂合变异，需验证是否在两条染色体上" if x.get("possible_compound_het") else "")
                    parts.append(f"{_cell(x['gene'])} {_cell(x['title'] or x['variant'])}（{_cell(x['classification'])}，"
                                 f"{zmap.get(x['zygosity'], _cell(x['zygosity']))}{tag}）")
                for x in ms.get("plp_not_unanimous", []):
                    parts.append(f"{_cell(x['gene'])} {_cell(x['title'] or x['variant'])}：ClinVar 分类不一致或缺少审查标准"
                                 f"（{_cell(x['classification'])}），需遗传咨询判断")
                if ms.get("not_scanned"):
                    parts.append("未完整扫描：" + "、".join(_cell(g) for g in sorted(ms["not_scanned"])))
                if ms.get("low_quality_not_counted"):
                    parts.append(f"另有 {len(ms['low_quality_not_counted'])} 个 ClinVar 致病相关变异测序质量不足，未计入，建议验证")
                if ms.get("indels_unresolved"):
                    parts.append(f"{ms['indels_unresolved']} 个插入/缺失变异未能与 ClinVar 比对")
                clean = not (ms["pathogenic_or_likely"] or ms.get("plp_not_unanimous") or ms.get("not_scanned")
                             or ms.get("low_quality_not_counted") or ms.get("indels_unresolved"))
                if clean and ms["genes"]:
                    parts.insert(0, "已扫描，未发现")
                plp_s = "；".join(parts) or "—"
                L.append(f"| {_cell(a['label_zh'])} | {_cell(_fmt(lab['value']))} {_cell(lab['unit'])} | {flag} | {pct} | "
                         f"{a['loci_called']}/{a.get('loci_informative', a['loci_tested'])} | {plp_s} |")
            meta = load_json(gp)
            if meta.get("absent_as_ref"):
                L += ["", f"> 本次 VCF 只列出变异位点；分析时把没列出的位点按「与参考序列相同」处理，依据：{_cell(meta['absent_as_ref'])}。"]
            L.append("")
    pp = ws / "work" / "insights" / "positions.json"
    if pp.exists():
        pos = load_json(pp)
        if pos["labs"] or pos["gmhi"] or pos["ages"]:
            L += ["## 你在同龄人群中的位置", "",
                  "> 化验参照：美国 NHANES 2005–2010 同性别、同年龄段人群（加权百分位）；中国人群分布尚未接入，美国人群的分布可能与中国人不同，只作方向参考。", ""]
            if pos["labs"]:
                L += ["| 指标 | 你的结果 | 参照人群 | 百分位 | 该人群中位数 |", "|---|---|---|---|---|"]
                for x in pos["labs"]:
                    pct = ("低于第 1 百分位" if x["bound"] == "below" else "高于第 99 百分位" if x["bound"] == "above" else
                           f"第 {_fmt(x['pct_low'])}–{_fmt(x['pct_high'])} 百分位（该人群很多人与你同值，多为检测下限）" if x["bound"] == "tie"
                           else f"第 {_fmt(x['pct'])} 百分位")
                    L.append(f"| {_cell(x['label_zh'])} | {_cell(_fmt(x['value']))} {_cell(x['unit'])} | {_cell(x['stratum'])}（{x['n']} 人） | {pct} | {_cell(_fmt(x['median']))} |")
                L.append("")
            if pos["gmhi"]:
                gm = pos["gmhi"]
                parts = [f"全部健康人群（{gm['n_healthy']} 人）第 {_fmt(gm['healthy']['pct'])} 百分位",
                         f"疾病人群（{gm['n_nonhealthy']} 人）第 {_fmt(gm['nonhealthy']['pct'])} 百分位"]
                if gm.get("healthy_east_asia"):
                    parts.insert(1, f"东亚健康人群（{gm['n_healthy_east_asia']} 人）第 {_fmt(gm['healthy_east_asia']['pct'])} 百分位")
                L += [f"**肠道菌群健康指数 GMHI {_fmt(gm['value'])}**：" + "；".join(parts) + "。东亚健康人群的 GMHI 整体偏低，与你对照时以东亚健康人群为主。", ""]
            if len(pos["ages"]) > 1:
                L += ["**几种「年龄」并列**", "", "| 来源 | 数值 | 与实际年龄之差 |", "|---|---|---|"]
                for x in pos["ages"]:
                    rng = f"（{_fmt(x['low'])}–{_fmt(x['high'])}，AI 估计）" if x.get("low") is not None else ""
                    diff = "—" if x.get("minus_age") is None else f"{x['minus_age']:+.1f} 岁"
                    L.append(f"| {_cell(x['label_zh'])} | {_cell(_fmt(x['value']))} 岁{rng} | {diff} |")
                L.append("")
    prj = ws / "work" / "insights" / "projections.json"
    if prj.exists():
        L += ["## 干预的因果推算（孟德尔随机化）", "",
              "> 用公开的孟德尔随机化研究估计（EpiGraphDB）推算「把某项指标降到目标值，你的风险大约会变成多少」。它假定人群中的因果效应也适用于你，只是方向和量级的参考。", "",
              "| 暴露 → 结局 | 你的现值 → 目标 | 基线风险 | 推算后风险（95% 区间） |", "|---|---|---|---|"]
        for x in load_json(prj)["projections"]:
            L.append(f"| {_cell(x['exposure'])} → {_cell(x['outcome'])} | {_cell(_fmt(x['member_value']))} → {_cell(_fmt(x['target']))} {_cell(x['unit'])} | "
                     f"{_pct(x['baseline_risk'])} | {_pct(x['risk_after'])}（{_pct(x['risk_after_ci'][0])}–{_pct(x['risk_after_ci'][1])}） |")
        L.append("")
        for x in load_json(prj)["projections"]:
            if x.get("exposure_match"):
                L.append(f"- 「{_cell(x['exposure'])}」作为{_cell(x['analyte'])}的替代暴露，依据：{_cell(x['exposure_match'])}")
        cav = []
        for x in load_json(prj)["projections"]:
            cav += [c for c in x.get("caveats_zh", []) if c not in cav]
        L += [f"- {_cell(c)}" for c in cav] + ([""] if cav else [])
    rows = B.rows(ws)
    if rows:
        L += ["## 问题看板：为你提出的问题与研究结论", "",
              "> 每个问题由一个独立的 AI 研究员基于你的数据和本次实时检索的公共数据库、文献作出判断。结论只说明「现有证据是否支持」，不是诊断。", "",
              "| 编号 | 问题 | 结论 | 把握 | 下一步 |", "|---|---|---|---|---|"]
        for r in rows:
            q, f = r["q"], r["f"]
            verdict = B.VERDICTS.get((f or {}).get("verdict"), "未研究（" + _cell(r["skipped"]) + "）" if r.get("skipped") else "未完成")
            conf = {"low": "低", "moderate": "中"}.get((f or {}).get("confidence"), "—")
            nxt = _cell(substitute(f.get("next_step_zh", ""), ro)) if f else "—"
            L.append(f"| {q['id']} | {_cell(substitute(q['title_zh'], ro))} | {verdict} | {conf} | {nxt} |")
        L.append("")
        for r in rows:
            q, f = r["q"], r["f"]
            if not f:
                continue
            ev = "；".join(_ev_link({"type": "pubmed", "ref": x[5:]}) if x.startswith("pmid:") else _cell(x) for x in f.get("public_evidence", [])[:8])
            L += [f"**{q['id']} {_cell(substitute(q['title_zh'], ro))}**", "",
                  f"- 假设：{_cell(substitute(q['hypothesis_zh'], ro))}",
                  f"- 结论：{B.VERDICTS[f['verdict']]}。{_cell(substitute(f['summary_zh'], ro))}",
                  f"- 局限：{_cell(substitute(f['limitations_zh'], ro))}"] + ([f"- 依据：{ev}"] if ev else []) + [""]
    return L


def _ev_link(e: Dict[str, Any]) -> str:
    t, r = e.get("type"), str(e.get("ref"))
    if t == "pubmed" and r.isdigit():
        return f"[PMID {r}](https://pubmed.ncbi.nlm.nih.gov/{r}/)"
    if t == "guideline" and re.match(r"^https?://[^\s)?#]+$", r):
        host = re.sub(r"^https?://([^/]+).*", r"\1", r)
        return f"[指南（{_cell(host)}）]({r})"
    return f"{_cell(t)}:{_cell(r)}"


def _html(md: str, title: str) -> str:
    try:
        import markdown  # type: ignore
        body = markdown.markdown(md, extensions=["tables"])
    except ImportError:
        body = f"<pre>{html.escape(md)}</pre>"
    css = """
:root{--bg:#fbfaf7;--fg:#1d1d1f;--muted:#5f6368;--line:#e3e0d8;--accent:#5b4b8a}
@media (prefers-color-scheme:dark){:root{--bg:#141417;--fg:#ececef;--muted:#a1a1aa;--line:#2c2c33;--accent:#b7a6ee}}
body{background:var(--bg);color:var(--fg);font:15px/1.7 -apple-system,"PingFang SC","Noto Sans CJK SC",sans-serif;max-width:920px;margin:0 auto;padding:24px 16px}
h1{font-size:26px;border-bottom:2px solid var(--accent);padding-bottom:8px}h2{margin-top:36px;color:var(--accent)}
table{border-collapse:collapse;width:100%;display:block;overflow-x:auto;font-size:13.5px}th,td{border:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
th{background:color-mix(in srgb,var(--accent) 10%,transparent)}a{color:var(--accent)}blockquote{border-left:3px solid var(--accent);margin:0;padding:4px 12px;color:var(--muted)}
"""
    return (f"<!doctype html><html lang=zh-CN><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
            f"<title>抗衰分析报告 {html.escape(title)}</title><style>{css}</style></head><body>{body}</body></html>")
