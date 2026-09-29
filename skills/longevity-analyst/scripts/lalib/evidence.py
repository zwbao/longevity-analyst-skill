"""Evidence lookup (local library + live PubMed) and intervention-plan checks.

Lookups return rows exactly as stored, each with its id, so a plan can cite
them. PubMed is queried live; a network failure is an error, never a silent
fall back to the local tables.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List

from .common import (EXIT_EXTERNAL, EXIT_INPUT, LAError, data, load_json,
                     longevity_skills_home, now_iso, write_json)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_CN = "零一二三四五六七八九十百千万两半〇壹贰叁肆伍陆柒捌玖拾佰仟廿卅"
_UNITS = (r"(?:mg|mcg|µg|μg|ug|µl|μl|ul|iu|ml|g|caps?(?:ules?)?|tabs?(?:lets?)?|pills?|drops?|softgels?|"
          r"毫克|克|微克|国际单位|单位|粒|片|颗|顆|丸|袋|支|滴|毫升|勺|小勺|大勺|胶囊|包|瓶|茶匙|汤匙|杯|剂|针)")
_CONC = r"(?!\s*/\s*(?:L|dL|mL|l|dl|ml|g|mmol|mol|升|分升|毫升|克)(?![A-Za-z]))"      # 3 mg/L, 50 mg/dL, 30 mg/g: concentrations
DOSE = re.compile(rf"(?:\d+(?:\.\d+)?|[{_CN}]+)\s*(?:k\s*)?(?:个|颗|顆|粒|片)?\s*{_UNITS}{_CONC}", re.I)
DOSE_WORDS = re.compile(rf"\b(?:a|an|one|two|three|four|five|half)\s+{_UNITS}(?![a-z]){_CONC}", re.I)
# Drug and supplement units: never allowed anywhere. Food and lifestyle amounts (克, 毫升, 杯, 包, 勺, 针) are allowed only
# inside an explicit {{n:…}} literal, which the reviewer reads - unless a drug or supplement is named in the same sentence.
_DRUG_UNITS = (r"(?:mg|mcg|µg|μg|ug|iu|caps?(?:ules?)?|tabs?(?:lets?)?|pills?|drops?|softgels?|"
               r"毫克|微克|国际单位|单位|粒|片|丸|胶囊|滴|剂|支|U(?![a-z]))")
DRUG_DOSE = re.compile(rf"(?:\d+(?:\.\d+)?|[{_CN}]+)\s*(?:k\s*)?(?:个|颗|顆|粒|片)?\s*{_DRUG_UNITS}{_CONC}", re.I)
DRUG_WORDS = re.compile(rf"\b(?:a|an|one|two|three|four|five|half)\s+{_DRUG_UNITS}(?![a-z]){_CONC}", re.I)
DRUG_CONTEXT = re.compile(r"(?i)药物|药片|服药|用药|处方|鱼油|他汀|二甲双胍|阿卡波糖|胰岛素|褪黑素|补充剂|辅酶|NMN|NAD\+?|"
                          r"白藜芦醇|雷帕霉素|三七|阿司匹林|司美格鲁肽|钙片|镁剂|锌剂|铁剂|叶酸片|肌酸粉|蛋白粉|胶原蛋白肽|"
                          r"supplement|metformin|statin|insulin|melatonin|rapamycin|aspirin|semaglutide")
# a nutrient name counts only right before the amount ("维生素C 1 g", "肌酸 5 g"), not in a food description
NUTRIENT_BEFORE = re.compile(r"(?i)(?:维生素\s*[A-Z]?\d*|益生菌|叶酸|镁|锌|钙|铁|硒|DHA|EPA|NR|肌酸|牛磺酸|甘氨酸|胶原蛋白|"
                             r"vitamin\s*[A-Z]?\d*|creatine|taurine|glycine)\s*[:：]?\s*$")
NUTRIENT_ANY = re.compile(NUTRIENT_BEFORE.pattern.replace(r"\s*[:：]?\s*$", ""), re.I)
INTAKE_VERB = re.compile(r"补充|服用|冲服|口服|吞服|服下|take|supplement with", re.I)   # not 吃/喝: food is eaten


def _norm(text: str) -> str:
    t = re.sub(r"[\u200b-\u200f\u2060-\u2064\ufeff\u00ad]", "", unicodedata.normalize("NFKC", text or ""))
    return re.sub(r"(?<=[\u4e00-\u9fff])\s+(?=[\u4e00-\u9fff])", "", t)          # "500 毫 克"


def _has_dose(text: str, drug_only: bool = False) -> bool:
    t = _norm(text)
    if drug_only:
        return bool(DRUG_DOSE.search(t) or DRUG_WORDS.search(t))
    return bool(DOSE.search(t) or DOSE_WORDS.search(t))


def has_dose_outside_literals(text: str) -> bool:
    """A dose in the text as it will render. Outside {{n:…}} any dose counts; inside, drug/supplement units count,
    and so does any amount in a sentence that names a drug or supplement. {{r:id|own words}} renders its own words."""
    rtext = re.sub(r"\{\{r:[^{}|]*\|([^{}]*)\}\}", lambda m: m.group(1).split("|")[-1] if m.group(1).split("|")[-1] not in
                   ("label", "value", "unit") else " 1 ", text or "")
    rtext = re.sub(r"\{\{r:[^{}]*\}\}", " 1 ", rtext)
    rest = re.sub(r"\{\{n:[^{}]*\}\}", " ", rtext)
    rendered = re.sub(r"\{\{n:([^{}]*)\}\}", r"\1", rtext)
    if _has_dose(rest) or _has_dose(rendered, drug_only=True):
        return True
    t = _norm(rendered)            # an amount next to a drug/supplement product name, or right after a nutrient name
    for m in list(DOSE.finditer(t)) + list(DOSE_WORDS.finditer(t)):
        near = t[max(0, m.start() - 12):m.end() + 12]
        if DRUG_CONTEXT.search(t[max(0, m.start() - 8):m.end() + 8]) or NUTRIENT_BEFORE.search(t[max(0, m.start() - 12):m.start()]) \
                or (INTAKE_VERB.search(t[max(0, m.start() - 8):m.start()]) and NUTRIENT_ANY.search(near)):
            return True
    return False


def _fold(s: str) -> str:
    return re.sub(r"[\s_\-·]+", "", unicodedata.normalize("NFKC", s).casefold())


def local(terms: List[str]) -> Dict[str, Any]:
    home = longevity_skills_home()
    ft = [_fold(t) for t in terms if t.strip()]
    effects = []
    with open(home / "data" / "effects.jsonl", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            hay = [_fold(x) for x in [r.get("marker", ""), r.get("marker_zh", ""), r.get("intervention", ""),
                                       r.get("intervention_zh", ""), *r.get("keywords", [])]]
            if any(t in h or h in t for t in ft for h in hay if h):
                effects.append(r)
    claims = []
    cpath = home / "skills" / "longevity-evidence" / "data" / "claims.jsonl"
    with open(cpath, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            hay = [_fold(x) for x in [r.get("entity", ""), *r.get("aliases", [])]]
            if any(t == h for t in ft for h in hay if h) and r.get("direction") not in (None, "named"):
                claims.append(r)
    return {"terms": terms, "effects": effects, "claims": claims[:200],
            "note": "effects = intervention→marker trial effects with quotes; claims = paper statements (human/animal/cell). Claims with direction 'named' or null are excluded (no effect direction)."}


def _get(url: str, retries: int = 3) -> bytes:
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (i + 1))
    raise LAError(f"PubMed unreachable ({last}); report the gap, do not substitute remembered citations", EXIT_EXTERNAL)


def pubmed(ws_root: Path, query: str, retmax: int = 10) -> Dict[str, Any]:
    q = urllib.parse.urlencode({"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json", "sort": "relevance"})
    ids = json.loads(_get(EUTILS + "esearch.fcgi?" + q))["esearchresult"].get("idlist", [])
    items = []
    if ids:
        s = json.loads(_get(EUTILS + "esummary.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "json"})))
        for pmid in ids:
            d = s["result"].get(pmid, {})
            items.append({"pmid": pmid, "title": d.get("title"), "journal": d.get("fulljournalname"),
                          "pubdate": d.get("pubdate"), "pubtype": d.get("pubtype", []),
                          "retracted": any("retract" in t.lower() for t in d.get("pubtype", [])),
                          "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        abstracts = _abstracts(ids)
        for it in items:
            it["abstract"] = abstracts.get(it["pmid"], "")
    res = {"query": query, "retrieved_at": now_iso(), "count": len(items), "items": items}
    key = hashlib.sha1(query.encode()).hexdigest()[:10]
    write_json(ws_root / "work" / "evidence" / f"pubmed_{key}.json", res)
    return res


def _abstracts(ids: List[str]) -> Dict[str, str]:
    """PubMed abstracts (plain text, first 2500 characters) so an estimator can read baseline figures, not only titles.
    A failure here leaves abstracts empty; it never blocks the search result."""
    try:
        import xml.etree.ElementTree as ET
        x = _get(EUTILS + "efetch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "xml"}))
        out = {}
        for art in ET.fromstring(x).iter("PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""
            parts = []
            for a in art.iter("AbstractText"):
                label = a.get("Label")
                text = "".join(a.itertext()).strip()
                parts.append(f"{label}: {text}" if label else text)
            out[pmid] = " ".join(parts)[:2500]
        return out
    except Exception:  # noqa: BLE001
        return {}


def _cached_pmids(ws_root: Path) -> set:
    out = set()
    for p in (ws_root / "work" / "evidence").glob("pubmed_*.json"):
        out |= {i["pmid"] for i in load_json(p)["items"]}
    return out


def verify_pmids(pmids: List[str]) -> Dict[str, Any]:
    """Live check that each PMID exists (a cache file alone can be forged)."""
    if not pmids:
        return {}
    s = json.loads(_get(EUTILS + "esummary.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "json"})))
    res = s.get("result", {})
    ok = {}
    for p in pmids:
        d = res.get(p)
        if not d or d.get("error"):
            continue
        types = [t.lower() for t in d.get("pubtype", [])]
        if any("retract" in t for t in types):
            continue                                   # retracted articles are not evidence
        ok[p] = d
    return ok


def fetch_url(ws_root: Path, url: str) -> Dict[str, Any]:
    """Record that a guideline URL was actually retrieved in this run (status + title)."""
    if not re.match(r"^https?://[^\s?#]+$", url):
        raise LAError("only plain http(s) URLs without a query string or fragment (the link is shown to the member; "
                      "use the guideline page's own address)", EXIT_INPUT)
    path = urllib.parse.unquote(re.sub(r"^https?://[^/]+", "", url)).replace("-", " ").replace("_", " ")
    if _has_dose(path) or DOSE_WORDS.search(path) or (DRUG_CONTEXT.search(path) and re.search(r"\d", path)) or re.search(r"(?i)(?<![a-z])(?:iu|mg|mcg|ug|caps?|capsules?|tablets?|pills?|softgels?|dose|dosage)(?![a-z])", path):
        raise LAError("the URL path names a dose; a guideline link is shown to the member and must not carry one", EXIT_INPUT)
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "longevity-analyst"}), timeout=20) as r:
            body = r.read(200000).decode("utf-8", errors="replace")
            status = r.status
    except Exception as e:  # noqa: BLE001
        raise LAError(f"could not retrieve {url}: {e}", EXIT_EXTERNAL)
    m = re.search(r"<title[^>]*>(.*?)</title>", body, re.S | re.I)
    rec = {"url": url, "status": status, "title": (m.group(1).strip()[:200] if m else ""), "retrieved_at": now_iso()}
    key = hashlib.sha1(url.encode()).hexdigest()[:10]
    write_json(ws_root / "work" / "evidence" / f"url_{key}.json", rec)
    return rec


def _cached_urls(ws_root: Path) -> set:
    return {load_json(p)["url"] for p in (ws_root / "work" / "evidence").glob("url_*.json")}


def check_plan(st: Dict[str, Any], ws_root: Path, plan_path: Path) -> Dict[str, Any]:
    plan = load_json(plan_path)
    menu = data("intervention_menu.json")
    cats = {c["id"]: c for c in menu["categories"]}
    ro_ids = {r["id"] for r in load_json(ws_root / "work" / "readouts.json")["readouts"]}
    op = ws_root / "work" / "organs" / "organ_readouts.json"
    if op.exists():
        ro_ids |= {r["id"] for r in load_json(op)["readouts"]}
    lab_names = {l["marker"] for l in st["labs"]}
    home = longevity_skills_home()
    eff_ids = {json.loads(l)["id"] for l in open(home / "data" / "effects.jsonl", encoding="utf-8")}
    claim_ids = {json.loads(l)["id"] for l in open(home / "skills" / "longevity-evidence" / "data" / "claims.jsonl", encoding="utf-8")}
    pmids = _cached_pmids(ws_root)
    urls = _cached_urls(ws_root)
    cited_urls = sorted({str(e.get("ref")) for it in (plan.get("items") or []) for e in it.get("evidence", []) if e.get("type") == "guideline"})
    live_urls = set()
    for u in cited_urls:
        if u in urls:
            try:
                fetch_url(ws_root, u)          # re-fetched now: a hand-written cache entry is not evidence
                live_urls.add(u)
            except LAError:
                pass
    urls = live_urls
    cited = sorted({str(e.get("ref")) for it in (plan.get("items") or []) for e in it.get("evidence", []) if e.get("type") == "pubmed"})
    live = verify_pmids([p for p in cited if p.isdigit()]) if cited else {}
    problems: List[str] = []
    items = plan.get("items")
    if not isinstance(items, list) or not items:
        raise LAError("plan.json needs a non-empty `items` list", EXIT_INPUT)
    item_keys = {"id", "category", "action_zh", "rationale_zh", "targets", "evidence", "evidence_grade", "executor", "retest"}
    for it in items:
        if not isinstance(it, dict):
            raise LAError("every plan item must be an object", EXIT_INPUT)
        iid = it.get("id", "?")
        if set(it) - item_keys:
            problems.append(f"{iid}: unknown fields {sorted(set(it) - item_keys)} (see references/plan-schema.md)")
        if it.get("retest") is not None and (not isinstance(it["retest"], dict) or set(it["retest"]) - {"what", "after_weeks"}):
            problems.append(f"{iid}: retest has only `what` and `after_weeks`")
        for e in it.get("evidence") or []:
            if not isinstance(e, dict) or set(e) != {"type", "ref"}:
                problems.append(f"{iid}: each evidence entry is exactly {{type, ref}}")
        for k in ("id", "category", "action_zh", "targets", "evidence", "executor"):
            if k not in it:
                problems.append(f"{iid}: missing field {k}")
        cat = cats.get(it.get("category"))
        if not cat:
            problems.append(f"{iid}: category {it.get('category')!r} is not on the intervention menu")
        elif it.get("executor") not in cat["executors"]:
            problems.append(f"{iid}: executor {it.get('executor')!r} not allowed for {cat['id']} (allowed {cat['executors']})")
        for fld in ("action_zh", "rationale_zh"):
            if has_dose_outside_literals(it.get(fld, "")) or (fld == "action_zh" and has_dose_outside_literals(str((it.get("retest") or {}).get("what", "")))):
                problems.append(f"{iid}: {fld} names a dose or pill count; doses are for the physician, remove it")
        file_ids = {f["id"] for f in st["files"]}
        for t in it.get("targets", []):
            if t.startswith("file:") and t[5:] in file_ids and it.get("category") in ("test", "referral"):
                continue            # e.g. "ask the lab what assay F005 is"
            if t not in ro_ids and t not in lab_names:
                problems.append(f"{iid}: target {t!r} is neither a readout id, a lab marker of this member, nor file:<id> on a test/referral item")
        ev = it.get("evidence", [])
        if not ev and it.get("category") not in ("test", "referral"):
            problems.append(f"{iid}: no evidence cited")
        for e in ev:
            typ, ref = e.get("type"), str(e.get("ref", ""))
            ok = (typ == "effects" and ref in eff_ids) or (typ == "claims" and ref in claim_ids) or \
                 (typ == "pubmed" and ref in pmids and ref in live) or (typ == "guideline" and ref in urls)
            if not ok:
                problems.append(f"{iid}: evidence {typ}:{ref} was not retrieved and verified in this run (effects/claims id, a PMID from `la.py evidence pubmed` that PubMed confirms, or a URL fetched with `la.py evidence fetch`)")
        rt = it.get("retest") or {}
        if rt.get("after_weeks") is not None and (isinstance(rt["after_weeks"], bool) or not isinstance(rt["after_weeks"], int)):
            problems.append(f"{iid}: retest.after_weeks must be a JSON integer")
        if it.get("evidence_grade") not in (None, "human_rct", "human_cohort", "animal_or_cell", "guideline", "not_applicable"):
            problems.append(f"{iid}: evidence_grade must be human_rct|human_cohort|animal_or_cell|guideline|not_applicable")
        if it.get("category") not in ("test", "referral") and it.get("evidence_grade") == "not_applicable":
            problems.append(f"{iid}: only test/referral items may have evidence_grade not_applicable")
        if not re.fullmatch(r"I[1-9][0-9]?", str(it.get("id", "")), re.ASCII):
            problems.append(f"{iid!r}: item id must be I1, I2, … (the digits are item numbers only)")
        if it.get("category") not in ("test", "referral"):
            if not rt.get("what") or not rt.get("after_weeks"):
                problems.append(f"{iid}: retest needs `what` and `after_weeks`")
        if rt.get("after_weeks") is not None:
            try:
                wk = int(rt["after_weeks"])
                if not 1 <= wk <= 260:
                    raise ValueError
            except (TypeError, ValueError):
                problems.append(f"{iid}: retest.after_weeks must be a whole number of weeks between 1 and 260")
    if problems:
        raise LAError("plan rejected:\n- " + "\n- ".join(problems), EXIT_INPUT)
    dst = ws_root / "work" / "intervene" / "plan.json"
    write_json(dst, plan)
    st["intervene"] = {"plan": str(dst), "items": len(items), "registered_at": now_iso()}
    return st["intervene"]
