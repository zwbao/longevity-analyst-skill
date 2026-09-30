"""Read a member's labs and daily wearable values from Mirobody through the member's personal MCP (read only), into CSV
files in a data folder before intake. Intake treats them like any other file: lab rows are still confirmed by the
agent (`labs candidates` / `labs confirm`), watch days are still mapped by the agent (`insights wearable`).

- Labs are written one file per checkup date (mirobody_labs_<date>.csv), so rows from different years are never
  mixed into one checkup; Mirobody's raw view carries no printed reference range.
- A result Mirobody cut (row cap, render limit) is asked again in smaller pieces; one that is still cut stops the
  pull. Nothing is written until everything was read.
- The MCP URL carries the member's secret: it is read from a file or the environment, never printed or stored.
Mirobody has no file-download tool, so omics raw files (VCF, methylation, stool, proteomics) still come from folders.
"""
from __future__ import annotations

import csv
import json
import os
import re
import urllib.request
from datetime import date as _date, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import EXIT_EXTERNAL, EXIT_INPUT, LAError, now_iso

PROTOCOL = "2025-03-26"
DAILY_MIN_ROWS = 14        # a daily series: at least this many readings, on at least half of the days it spans
CUT = re.compile(r"^(?:…|\.\.\.)\s*cut at")
PAIR = re.compile(r"(?:^|, )([a-z_][a-z0-9_]*)=")
OUR_FILES = re.compile(r"^mirobody_(labs(_[0-9a-z-]+)?|wearable_daily)\.csv$")


def _scrub(text: Any, url: str) -> str:
    return str(text).replace(url, "<mirobody-mcp>")


class Mcp:
    def __init__(self, url: str, timeout: int = 90):
        if not re.match(r"^https?://", url or ""):
            raise LAError("the Mirobody MCP URL must start with http(s)://", EXIT_INPUT)
        self.url, self.timeout, self.session, self.n = url, timeout, None, 0

    def _post(self, payload: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        headers = {"Content-Type": "application/json", "Accept": "application/json, text/event-stream"}
        if self.session:
            headers["Mcp-Session-Id"] = self.session
        req = urllib.request.Request(self.url, data=json.dumps(payload).encode(), headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                self.session = r.headers.get("Mcp-Session-Id") or self.session
                body = r.read().decode("utf-8")
                ctype = r.headers.get("Content-Type", "")
        except Exception as e:  # noqa: BLE001 - the URL holds the member's secret: never echoed
            raise LAError(f"Mirobody MCP unreachable ({type(e).__name__}: {_scrub(e, self.url)[:120]})", EXIT_EXTERNAL) from None
        if not body.strip():
            return None
        msgs: List[Any] = []
        if "text/event-stream" in ctype:
            for event in re.split(r"\r?\n\r?\n", body):
                data = "\n".join(l[5:].lstrip() for l in event.splitlines() if l.startswith("data:"))
                if data:
                    try:
                        msgs.append(json.loads(data))
                    except ValueError:
                        continue
        else:
            try:
                msgs.append(json.loads(body))
            except ValueError:
                raise LAError(f"Mirobody MCP answered something that is not JSON ({ctype or 'no content type'})", EXIT_EXTERNAL) from None
        want = payload.get("id")
        for m in msgs:                                  # the response to this request, not a notification beside it
            if isinstance(m, dict) and (want is None or m.get("id") == want):
                return m
        return None if want is not None else (msgs[-1] if msgs and isinstance(msgs[-1], dict) else None)

    def _rpc(self, method: str, params: Dict[str, Any]) -> Dict[str, Any]:
        self.n += 1
        res = self._post({"jsonrpc": "2.0", "id": self.n, "method": method, "params": params}) or {}
        if res.get("error"):
            raise LAError(f"Mirobody MCP {method}: {_scrub((res['error'] or {}).get('message'), self.url)[:200]}", EXIT_EXTERNAL)
        return res.get("result") or {}

    def open(self) -> List[str]:
        self._rpc("initialize", {"protocolVersion": PROTOCOL, "capabilities": {},
                                 "clientInfo": {"name": "longevity-analyst", "version": "0.7"}})
        self._post({"jsonrpc": "2.0", "method": "notifications/initialized"})
        return [t["name"] for t in self._rpc("tools/list", {}).get("tools", [])]

    def call(self, tool: str, args: Dict[str, Any]) -> Tuple[str, bool]:
        """The table text and whether Mirobody says it was cut."""
        res = self._rpc("tools/call", {"name": tool, "arguments": args})
        text = "\n".join(c.get("text", "") for c in res.get("content", []) if c.get("type") == "text")
        if res.get("isError"):
            raise LAError(f"Mirobody {tool} refused: {_scrub(text, self.url)[:300]}", EXIT_EXTERNAL)
        cut = False
        try:                                              # {"result": "<table>", "status": ..., "truncated": ...}
            env = json.loads(text)
        except ValueError:
            env = None
        if isinstance(env, dict) and "result" in env:
            if env.get("status") not in (None, "ok"):
                raise LAError(f"Mirobody {tool}: {_scrub(env.get('result'), self.url)[:300]}", EXIT_EXTERNAL)
            cut = bool(env.get("truncated"))
            text = str(env["result"])
        if any(CUT.match(l) for l in text.splitlines()) or re.search(r"(?m)^\(.*\btruncated\b", text):
            cut = True
        return text, cut


def parse_table(text: str) -> List[Dict[str, str]]:
    """Mirobody's compact table: an optional '(constants: k=v, …)' line, a header 'a|b|c', rows, then a blank line and
    meta lines. A value is not escaped, so extra cells belong to the last column. A table whose every column is
    constant has no header: the constants are its one row. '… cut at N characters' is never a row."""
    constants: Dict[str, str] = {}
    header: Optional[List[str]] = None
    rows: List[Dict[str, str]] = []
    for line in text.splitlines():
        if header is None:
            if line.startswith("(constants:"):
                body = line[len("(constants:"):].strip()
                body = body[:-1] if body.endswith(")") else body
                starts = [(m.start(), m.group(1), m.end()) for m in PAIR.finditer(body)]
                for i, (_, key, v) in enumerate(starts):
                    end = starts[i + 1][0] if i + 1 < len(starts) else len(body)
                    constants[key] = body[v:end]
                continue
            if not line.strip():
                if constants:
                    break                                # constants, then the blank line: a one-row table
                continue
            if line.startswith("(") or CUT.match(line) or line.startswith(("notes:", "next:")):
                continue
            header = line.split("|")
            continue
        if not line.strip() or CUT.match(line):
            break                                        # the table ends; meta lines follow
        cells = line.split("|")
        if len(cells) > len(header):
            cells = cells[:len(header) - 1] + ["|".join(cells[len(header) - 1:])]
        rows.append({**constants, **dict(zip(header, cells))})
    if header is None and constants:
        rows.append(constants)
    return rows


def _pick(row: Dict[str, str], *keys: str) -> str:
    for k in keys:
        if row.get(k) not in (None, ""):
            return row[k]
    return ""


def _read_all(mcp: Mcp, names: List[str], extra: Dict[str, Any], batch: int) -> List[Dict[str, str]]:
    """Every row for these indicators; a cut result is asked again in halves down to one indicator."""
    out: List[Dict[str, str]] = []
    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        text, cut = mcp.call("query_health_indicators", {"indicators": chunk, **extra})
        if cut:
            if len(chunk) == 1:
                raise LAError(f"Mirobody cut the readings of {chunk[0]!r} even alone; narrow the window (--days) and pull again",
                              EXIT_EXTERNAL)
            out += _read_all(mcp, chunk, extra, max(1, len(chunk) // 2))
            continue
        for r in parse_table(text):
            name = _pick(r, "indicator", "name")
            if not name and len(chunk) == 1:          # a one-indicator table carries no indicator column
                r["indicator"] = name = chunk[0]
            if name in chunk:                         # an unknown name makes Mirobody answer with the catalogue instead
                out.append(r)
    return out


def _write(path: Path, fields: List[str], rows: List[Dict[str, str]]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def pull(url: str, out_dir: Path, days: int = 120) -> Dict[str, Any]:
    """Write mirobody_labs_<date>.csv (one per checkup date) and mirobody_wearable_daily.csv into the folder."""
    mcp = Mcp(url)
    tools = mcp.open()
    if "query_health_indicators" not in tools:
        raise LAError("this Mirobody account has no health data yet (query_health_indicators is not offered)", EXIT_INPUT)
    text, cut = mcp.call("query_health_indicators", {})
    if cut:
        raise LAError("Mirobody cut the list of indicators; the pull would miss some", EXIT_EXTERNAL)
    daily, labs = [], []
    for r in parse_table(text):
        n = _pick(r, "indicator", "name")
        if not n:
            continue
        try:
            span = (_date.fromisoformat(r.get("last_date", "")[:10]) - _date.fromisoformat(r.get("first_date", "")[:10])).days + 1
            count = int(r.get("count") or 0)
        except ValueError:
            span, count = 0, 0
        (daily if count >= DAILY_MIN_ROWS and count * 2 >= span else labs).append(n)
    # read everything first; nothing is written if any read fails
    by_date: Dict[str, List[Dict[str, str]]] = {}
    for r in _read_all(mcp, labs, {"view": "raw"}, 10):
        code = _pick(r, "code")
        d = _pick(r, "date", "time", "start_time")[:10]
        by_date.setdefault(d if re.fullmatch(r"\d{4}-\d{2}-\d{2}", d) else "undated", []).append({
            "marker": _pick(r, "indicator", "name"), "value": _pick(r, "value"), "unit": _pick(r, "unit"),
            "ref_range": _pick(r, "reference_range", "ref_range", "range"), "date": d,
            "loinc": code if re.fullmatch(r"\d{1,7}-\d", code or "") else "", "source": "mirobody"})
    day_rows: Dict[str, Dict[str, str]] = {}
    units: Dict[str, str] = {}
    today = _date.today()
    for ind in daily:                                    # one indicator and 30 days per call
        for k in range(0, days, 30):
            a, b = today - timedelta(days=min(days, k + 30) - 1), today - timedelta(days=k)
            for r in _read_all(mcp, [ind], {"view": "day", "start": a.isoformat(), "end": b.isoformat()}, 1):
                d = _pick(r, "period", "date", "day", "time", "start_time")[:10]
                if d:
                    day_rows.setdefault(d, {"date": d})[ind] = _pick(r, "avg", "value", "mean")
                    units.setdefault(ind, _pick(r, "unit"))
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.iterdir():                        # only files this command writes are replaced
        if OUR_FILES.match(old.name):
            old.unlink()
    lab_files = []
    for d, rows in sorted(by_date.items()):
        p = out_dir / f"mirobody_labs_{d}.csv"
        _write(p, ["marker", "value", "unit", "ref_range", "date", "loinc", "source"], rows)
        lab_files.append({"file": str(p), "date": d, "rows": len(rows)})
    wp = out_dir / "mirobody_wearable_daily.csv"
    if day_rows:
        cols = {ind: f"{ind} ({units[ind]})" if units.get(ind) else ind for ind in daily}
        _write(wp, ["date"] + list(cols.values()),
               [{"date": r["date"], **{cols[k]: v for k, v in r.items() if k in cols}}
                for r in (day_rows[d] for d in sorted(day_rows, reverse=True))])
    return {"at": now_iso(), "indicators": len(labs) + len(daily), "lab_files": lab_files,
            "lab_rows": sum(f["rows"] for f in lab_files), "wearable_days": len(day_rows),
            "wearable_file": str(wp) if day_rows else None, "wearable_indicators": daily,
            "note": "One lab file per checkup date: confirm the rows of the checkup this report is about; earlier dates are "
                    "history. Mirobody's raw view carries no printed reference range. The watch file is mapped with "
                    "`la.py insights wearable`, not transcribed."}


def url_from(arg: Optional[str]) -> str:
    """--mcp-url-file <path> or the env var LONGPI_MCP_URL; a literal URL on the command line would land in logs."""
    if arg:
        p = Path(arg).expanduser()
        if not p.exists():
            raise LAError(f"{p} not found", EXIT_INPUT)
        url = p.read_text(encoding="utf-8").strip()
        if not url:
            raise LAError(f"{p} is empty (the URL copy is cleared after an import; start a new run)", EXIT_INPUT)
        return url
    env = os.environ.get("LONGPI_MCP_URL", "").strip()
    if not env:
        raise LAError("give --mcp-url-file <file holding the member's MCP URL> or set LONGPI_MCP_URL", EXIT_INPUT)
    return env
