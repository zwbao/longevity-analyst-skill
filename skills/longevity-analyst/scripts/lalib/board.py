"""The question board: the member's own open questions, each investigated by a researcher subagent.

The agent proposes the questions and the researchers judge them (references/researcher.md). The harness checks form
and binding only: question ids, that every basis item exists for this member, that every public citation is a
record retrieved in this workspace, and the allowed verdicts and confidence levels.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Set

from . import pubdata
from .common import EXIT_INPUT, LAError, load_json, now_iso, sha256_file, write_json

VERDICTS = {"supported": "证据支持", "not_supported": "证据不支持", "insufficient": "证据不足"}
CONFIDENCE = ("low", "moderate")
QID = re.compile(r"^Q(10|[1-9])$")


def member_ids(st: Dict[str, Any], ws: Path) -> Set[str]:
    ids: Set[str] = {"member.age", "member.sex"} | {f"member.answers.{k}" for k in (st["member"].get("answers") or {})}
    for p in (ws / "work" / "readouts.json", ws / "work" / "organs" / "organ_readouts.json", ws / "work" / "insights" / "insight_readouts.json"):
        if p.exists():
            ids |= {r["id"] for r in load_json(p)["readouts"]}
    from .labnames import usable_rows
    ids |= {l["marker"] for l in usable_rows(st)}                 # a refused or unconfirmed row is not the member's data
    gp = ws / "work" / "insights" / "genotype_phenotype.json"
    if gp.exists():
        for key, a in load_json(gp)["analytes"].items():
            ids.add(f"gen:{key}")
            ids |= {f"gt:{l['rsid']}" for l in a["loci"] if l.get("member_risk_alleles") is not None}
            ids |= {f"gt:{v['hgvs']}" for v in a["monogenic_scan"]["pathogenic_or_likely"]}
    ids |= {f"file:{f['id']}" for f in st["files"] if not f.get("excluded")}
    return ids


def public_refs(ws: Path) -> Set[str]:
    refs = set(pubdata.known_refs(ws))
    for p in (ws / "work" / "evidence").glob("pubmed_*.json"):
        refs |= {f"pmid:{i['pmid']}" for i in load_json(p)["items"]}
    return refs


def register_questions(st: Dict[str, Any], ws: Path, path: Path) -> Dict[str, Any]:
    try:
        doc = load_json(path)
    except (OSError, ValueError) as e:
        raise LAError(f"{path} is not JSON ({e})", EXIT_INPUT)
    qs = doc.get("questions") if isinstance(doc, dict) else None
    if not isinstance(qs, list) or not 3 <= len(qs) <= 10:
        raise LAError("questions JSON needs 3-10 `questions`", EXIT_INPUT)
    known = member_ids(st, ws)
    problems, seen = [], set()
    for i, q in enumerate(qs):
        w = f"questions[{i}]"
        if not isinstance(q, dict) or set(q) - {"id", "title_zh", "hypothesis_zh", "basis", "why_zh"}:
            problems.append(f"{w}: fields are id, title_zh, hypothesis_zh, basis, why_zh")
            continue
        if not QID.match(str(q.get("id", ""))) or q["id"] in seen:
            problems.append(f"{w}: id must be Q1..Q10, unique")
        seen.add(q.get("id"))
        for k in ("title_zh", "hypothesis_zh", "why_zh"):
            if not isinstance(q.get(k), str) or not q[k].strip():
                problems.append(f"{w}: {k} is empty")
        basis = q.get("basis")
        if not isinstance(basis, list) or not basis:
            problems.append(f"{w}: basis lists this member's data the question starts from")
        else:
            bad = [b for b in basis if b not in known]
            if bad:
                problems.append(f"{w}: basis {bad} is not this member's data (readout ids, lab markers, member.answers.*, gen:<analyte>, gt:<rsid>, file:<id>)")
    if problems:
        raise LAError("questions rejected:\n- " + "\n- ".join(problems), EXIT_INPUT)
    dst = ws / "work" / "insights" / "board" / "questions.json"
    write_json(dst, {"questions": qs, "registered_at": now_iso()})
    b = st.setdefault("insights", {}).setdefault("board", {})
    b.update(questions=[q["id"] for q in qs], questions_sha256=sha256_file(dst, limit=None), findings={}, skipped={})
    return {"questions": [q["id"] for q in qs]}


def register_finding(st: Dict[str, Any], ws: Path, qid: str) -> Dict[str, Any]:
    b = (st.get("insights") or {}).get("board") or {}
    if qid not in (b.get("questions") or []):
        raise LAError(f"{qid} is not a registered question", EXIT_INPUT)
    path = ws / "work" / "insights" / "board" / f"{qid}.json"
    if not path.exists():
        raise LAError(f"{path} not written", EXIT_INPUT)
    try:
        f = load_json(path)
    except ValueError as e:
        raise LAError(f"{path} is not JSON ({e})", EXIT_INPUT)
    allowed = {"id", "verdict", "confidence", "member_evidence", "public_evidence", "summary_zh", "next_step_zh", "limitations_zh"}
    problems: List[str] = []
    if not isinstance(f, dict):
        raise LAError(f"{path} must be an object", EXIT_INPUT)
    if set(f) - allowed:
        problems.append(f"unknown fields {sorted(set(f) - allowed)}")
    if f.get("id") != qid:
        problems.append("`id` does not match")
    if f.get("verdict") not in VERDICTS:
        problems.append(f"verdict must be one of {sorted(VERDICTS)}")
    if f.get("confidence") not in CONFIDENCE:
        problems.append(f"confidence must be one of {CONFIDENCE} (one person's data never carries high confidence)")
    known = member_ids(st, ws)
    me = f.get("member_evidence")
    if not isinstance(me, list) or not me:
        problems.append("member_evidence lists the member's data the verdict rests on")
    else:
        bad = [x for x in me if x not in known]
        if bad:
            problems.append(f"member_evidence {bad} is not this member's data")
    pe = f.get("public_evidence")
    pubs = public_refs(ws)
    if not isinstance(pe, list):
        problems.append("public_evidence must be a list of retrieved refs (gwas:, clinvar:, mr:, proj:, myvariant:, pmid:)")
    else:
        bad = [x for x in pe if x not in pubs]
        if bad:
            problems.append(f"public_evidence {bad[:5]} were not retrieved in this workspace")
        if f.get("verdict") in ("supported", "not_supported") and not pe:
            problems.append("a supported / not_supported verdict needs at least one retrieved public record")
        pm = sorted({x[5:] for x in pe if isinstance(x, str) and x.startswith("pmid:") and x in pubs})
        if pm and not problems:                      # a cache file can be edited; PMIDs are checked live at registration
            from .evidence import verify_pmids
            live = verify_pmids(pm)
            gone = [p for p in pm if p not in live]
            if gone:
                problems.append(f"PMIDs {gone[:5]} do not resolve on PubMed or are retracted")
    for k in ("summary_zh", "next_step_zh", "limitations_zh"):
        if not isinstance(f.get(k), str) or not f[k].strip():
            problems.append(f"{k} is empty")
    if isinstance(f.get("summary_zh"), str):
        from .report import trace_text
        tt = trace_text(f["summary_zh"] + "\n" + str(f.get("next_step_zh", "")) + "\n" + str(f.get("limitations_zh", "")))
        if tt:
            problems.append("the text breaks the number rules: " + "; ".join(tt[:3]))
    if problems:
        b.get("findings", {}).pop(qid, None)
        raise LAError(f"{qid} finding rejected:\n- " + "\n- ".join(problems), EXIT_INPUT)
    b.setdefault("findings", {})[qid] = {"path": str(path), "sha256": sha256_file(path, limit=None), "at": now_iso()}
    (b.get("skipped") or {}).pop(qid, None)
    return b["findings"][qid]


def skip(st: Dict[str, Any], qid: str, reason: str) -> Dict[str, Any]:
    b = (st.get("insights") or {}).get("board") or {}
    if qid not in (b.get("questions") or []):
        raise LAError(f"{qid} is not a registered question", EXIT_INPUT)
    b.setdefault("skipped", {})[qid] = {"reason": reason, "at": now_iso()}
    (b.get("findings") or {}).pop(qid, None)
    return {"skipped": qid}


def complete(st: Dict[str, Any]) -> bool:
    b = (st.get("insights") or {}).get("board") or {}
    qs = set(b.get("questions") or [])
    return bool(qs) and qs <= set(b.get("findings") or {}) | set(b.get("skipped") or {})


def integrity(ws: Path, st: Dict[str, Any]) -> List[str]:
    """Questions and findings are bound to the bytes the harness accepted; an edit after registration voids them."""
    b = (st.get("insights") or {}).get("board") or {}
    out = []
    qp = ws / "work" / "insights" / "board" / "questions.json"
    if b.get("questions_sha256") and (not qp.exists() or sha256_file(qp, limit=None) != b["questions_sha256"]):
        out.append("board/questions.json changed after `board questions`; register the questions again")
    for q, v in (b.get("findings") or {}).items():
        p = Path(v["path"])
        if not p.exists() or sha256_file(p, limit=None) != v["sha256"]:
            out.append(f"board/{q}.json changed after `board finding --id {q}`; register it again")
    return out


def registered_findings(ws: Path, st: Dict[str, Any] = None) -> Dict[str, Path]:
    """Only findings the harness accepted count; a rejected or skipped file on disk is never rendered or traced."""
    if st is None:
        st = load_json(ws / "state.json")
    fnd = (((st.get("insights") or {}).get("board") or {}).get("findings")) or {}
    return {q: Path(v["path"]) for q, v in fnd.items() if Path(v["path"]).exists()}


def rows(ws: Path, st: Dict[str, Any] = None) -> List[Dict[str, Any]]:
    qp = ws / "work" / "insights" / "board" / "questions.json"
    if not qp.exists():
        return []
    if st is None:
        st = load_json(ws / "state.json")
    reg = registered_findings(ws, st)
    skipped = ((st.get("insights") or {}).get("board") or {}).get("skipped") or {}
    out = []
    for q in load_json(qp)["questions"]:
        out.append({"q": q, "f": load_json(reg[q["id"]]) if q["id"] in reg else None,
                    "skipped": (skipped.get(q["id"]) or {}).get("reason")})
    return out
