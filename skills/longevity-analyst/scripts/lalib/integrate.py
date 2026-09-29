"""Group readouts by body system and check the analysts' write-ups.

`bundle` is a pure projection of readouts.json + labs + method notes.
`register` only checks structure: the write-up exists, has the required
sections, and every {{r:<id>}} it cites is a real readout. Whether the
reasoning is sound is the reviewer subagent's job, not this script's.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List

from .labnames import refused_rows, usable_rows
from .common import EXIT_INPUT, LAError, data, load_json, now_iso, write_json

PLACEHOLDER = re.compile(r"\{\{r:([A-Za-z0-9_.\-]+)(?:\|([^{}]*))?\}\}")
MODES = ("label", "value", "unit")
REQUIRED_SECTIONS = ["## 结论", "## 证据", "## 不确定性", "## 建议复测"]


def readouts(ws_root: Path) -> List[Dict[str, Any]]:
    p = ws_root / "work" / "readouts.json"
    if not p.exists():
        raise LAError("no readouts yet; run `la.py methods run`", EXIT_INPUT)
    return load_json(p)["readouts"]


def bundle(st: Dict[str, Any], ws_root: Path) -> Dict[str, Any]:
    ro = readouts(ws_root)
    sysdef = data("systems.json")["systems"]
    notes = {it["method"]: (it.get("result") or {}).get("notes", []) for it in (st.get("methods") or {}).get("items", [])}
    not_run = [{"method": it["method"], "status": it["status"],
                "reasons": it.get("reasons") or [str(x.get("message_zh") or x) for x in ((it.get("result") or {}).get("problems") or [])][:5]
                or ([str((it.get("result") or {}).get("error"))[:300]] if (it.get("result") or {}).get("error") else [])}
               for it in (st.get("methods") or {}).get("items", []) if it["status"] not in ("done",)]
    out_dir = ws_root / "work" / "integrate" / "bundles"
    made = []
    by_sys: Dict[str, List[Dict[str, Any]]] = {}
    for r in ro:
        for s in r.get("systems") or ["overall_aging"]:
            by_sys.setdefault(s, []).append(r)
    for s, rs in sorted(by_sys.items()):
        b = {
            "system": s, "label_zh": sysdef.get(s, {}).get("label_zh", s),
            "member": {k: st["member"].get(k) for k in ("id", "age", "sex", "mode", "sample_date", "answers")},
            "readouts": rs,
            "method_notes": {m: notes.get(m, []) for m in sorted({r["method"] for r in rs})},
            "labs": [{k: v for k, v in l.items() if k not in ("confirm", "_row_key")} for l in usable_rows(st)],
            "lab_rows_not_used": [{"marker": l["marker"], "why": (l.get("confirm") or {}).get("why", "not confirmed")} for l in refused_rows(st)],
            "not_run_methods": not_run,
            "identity": st.get("identity") or {"answer": "not_checked"},
            "files_not_analysed": [{"id": f["id"], "name": f["name"], "status": f.get("status"), "reason": f.get("reason") or (f.get("excluded") or {}).get("reason")}
                                   for f in st["files"] if f.get("status") in ("deferred_ask_lab", "excluded_not_member", "rejected", "unsupported", "needs_pipeline")],
            "contract": "references/system-analyst.md",
            "write_to": str(ws_root / "work" / "integrate" / "analyses" / f"{s}.md"),
        }
        write_json(out_dir / f"{s}.json", b)
        made.append({"system": s, "readouts": len(rs), "bundle": str(out_dir / f"{s}.json")})
    st.setdefault("integrate", {})["bundles"] = made
    return {"bundles": made}


def check_placeholders(text: str, ids: set) -> List[str]:
    return sorted({m.group(1) for m in PLACEHOLDER.finditer(text) if m.group(1) not in ids})


def register(st: Dict[str, Any], ws_root: Path, system: str) -> Dict[str, Any]:
    made = {b["system"] for b in (st.get("integrate") or {}).get("bundles", [])}
    if system not in made:
        raise LAError(f"{system!r} is not one of this run's bundles {sorted(made)}", EXIT_INPUT)
    path = ws_root / "work" / "integrate" / "analyses" / f"{system}.md"
    if not path.exists():
        raise LAError(f"{path} not written", EXIT_INPUT)
    text = path.read_text(encoding="utf-8")
    missing = [h for h in REQUIRED_SECTIONS if h not in text]
    ids = {r["id"] for r in readouts(ws_root)}
    bad = check_placeholders(text, ids)
    cited = sorted({m.group(1) for m in PLACEHOLDER.finditer(text)})
    problems = []
    if missing:
        problems.append(f"missing sections: {missing}")
    if bad:
        problems.append(f"cites readouts that do not exist: {bad}")
    if not cited:
        problems.append("cites no readout; every claim about this member must point at {{r:<id>}}")
    if problems:
        raise LAError(f"{system}: " + "; ".join(problems), EXIT_INPUT)
    reg = st.setdefault("integrate", {}).setdefault("analyses", {})
    reg[system] = {"path": str(path), "cited": cited, "registered_at": now_iso()}
    return reg[system]
