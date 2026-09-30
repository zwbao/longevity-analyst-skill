"""Read a member's confirmed labs and daily wearable values from Mirobody through the member's personal MCP (read only),
into CSV files in the member's data folder, before intake. Intake then treats them like any other file: lab rows are
still transcribed/confirmed by the agent (`labs candidates` / `labs confirm`), watch days still mapped by the agent.

The MCP URL carries the member's secret: it is read from a file or the environment, never printed or stored in state.
Mirobody has no file-download tool, so omics raw files (VCF, methylation, stool, proteomics) still come from the folder.
"""
from __future__ import annotations

import csv
import json
import os
import re
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .common import EXIT_EXTERNAL, EXIT_INPUT, LAError, now_iso

PROTOCOL = "2025-03-26"
DAILY_MIN_ROWS = 14        # an indicator with at least this many readings over at least this many days is a daily series


class Mcp:
    def __init__(self, url: str, timeout: int = 90):
        if not re.match(r"^https?://", url or ""):
            raise LAError("the Mirobody MCP URL must start with http(s)://", EXIT_INPUT)
        self.url, self.timeout, self.session, self.n, self.truncated = url, timeout, None, 0, False

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
        except Exception as e:  # noqa: BLE001 - the URL is never echoed: it holds the member's secret
            raise LAError(f"Mirobody MCP unreachable ({type(e).__name__}: {str(e).split('/mcp/')[0][:120]})", EXIT_EXTERNAL)
        if not body.strip():
            return None
        if "text/event-stream" in ctype:
            data = [l[5:].strip() for l in body.splitlines() if l.startswith("data:")]
            body = data[-1] if data else "{}"
        return json.loads(body)

    def _rpc(self, method: str, params: Dict[str, Any]) -> Dict[str, Any]:
        self.n += 1
        res = self._post({"jsonrpc": "2.0", "id": self.n, "method": method, "params": params}) or {}
        if res.get("error"):
            raise LAError(f"Mirobody MCP {method}: {str(res['error'].get('message'))[:200]}", EXIT_EXTERNAL)
        return res.get("result") or {}

    def open(self) -> List[str]:
        self._rpc("initialize", {"protocolVersion": PROTOCOL, "capabilities": {},
                                 "clientInfo": {"name": "longevity-analyst", "version": "0.7"}})
        self._post({"jsonrpc": "2.0", "method": "notifications/initialized"})
        return [t["name"] for t in self._rpc("tools/list", {}).get("tools", [])]

    def call(self, tool: str, args: Dict[str, Any]) -> str:
        res = self._rpc("tools/call", {"name": tool, "arguments": args})
        text = "\n".join(c.get("text", "") for c in res.get("content", []) if c.get("type") == "text")
        if res.get("isError"):
            raise LAError(f"Mirobody {tool} refused: {text[:300]}", EXIT_EXTERNAL)
        try:                                              # {"result": "<table>", "status": ..., "truncated": ...}
            env = json.loads(text)
        except ValueError:
            return text
        if isinstance(env, dict) and "result" in env:
            if env.get("status") not in (None, "ok"):
                raise LAError(f"Mirobody {tool}: {str(env.get('result'))[:300]}", EXIT_EXTERNAL)
            if env.get("truncated"):
                self.truncated = True
            return str(env["result"])
        return text


def parse_table(text: str) -> List[Dict[str, str]]:
    """Mirobody's compact table: an optional '(constants: k=v, …)' line, a header 'a|b|c', rows with '\\|' escaped,
    then a blank line and meta lines. Constants are merged into every row."""
    constants: Dict[str, str] = {}
    header: Optional[List[str]] = None
    rows: List[Dict[str, str]] = []
    for line in text.splitlines():
        if header is None:
            if line.startswith("(constants:"):
                for part in re.findall(r"(\w[\w ]*)=([^,)]*)", line[len("(constants:"):]):
                    constants[part[0].strip()] = part[1].strip()
                continue
            if not line.strip() or line.startswith("("):
                continue
            header = line.split("|")
            continue
        if not line.strip():
            break                                        # the table ends; meta lines follow
        cells = [c.replace("\\|", "|").replace("\\n", "\n").replace("\\\\", "\\")
                 for c in re.split(r"(?<!\\)\|", line)]
        rows.append({**constants, **dict(zip(header, cells))})
    if header is None and constants:
        rows.append(constants)
    return rows


def _pick(row: Dict[str, str], *keys: str) -> str:
    for k in keys:
        if row.get(k) not in (None, ""):
            return row[k]
    return ""


def pull(url: str, out_dir: Path, days: int = 120) -> Dict[str, Any]:
    """Write mirobody_labs.csv and mirobody_wearable_daily.csv into the member's data folder."""
    mcp = Mcp(url)
    tools = mcp.open()
    if "query_health_indicators" not in tools:
        raise LAError("this Mirobody account has no health data yet (query_health_indicators is not offered)", EXIT_INPUT)
    catalog = parse_table(mcp.call("query_health_indicators", {}))
    from datetime import date as _d
    daily, labs = [], []
    for r in catalog:
        n = _pick(r, "indicator", "name")
        if not n:
            continue
        try:
            span = (_d.fromisoformat(r.get("last_date", "")[:10]) - _d.fromisoformat(r.get("first_date", "")[:10])).days
            count = int(r.get("count") or 0)
        except ValueError:
            span, count = 0, 0
        (daily if count >= DAILY_MIN_ROWS and span >= DAILY_MIN_ROWS else labs).append(n)
    names, wearable = labs + daily, daily
    out_dir.mkdir(parents=True, exist_ok=True)
    lab_rows: List[Dict[str, str]] = []
    for i in range(0, len(labs), 40):
        for r in parse_table(mcp.call("query_health_indicators", {"indicators": labs[i:i + 40], "view": "raw"})):
            code = _pick(r, "code")
            lab_rows.append({"marker": _pick(r, "indicator", "name"), "value": _pick(r, "value"), "unit": _pick(r, "unit"),
                             "ref_range": _pick(r, "reference_range", "ref_range", "range"),
                             "date": _pick(r, "date", "time", "start_time")[:10],
                             "loinc": code if re.fullmatch(r"\d{1,7}-\d", code or "") else "", "source": "mirobody"})
    lp = out_dir / "mirobody_labs.csv"
    with open(lp, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["marker", "value", "unit", "ref_range", "date", "loinc", "source"])
        w.writeheader()
        w.writerows(lab_rows)
    day_rows: Dict[str, Dict[str, str]] = {}
    if wearable:                     # one indicator and 30 days per call: a longer table is cut at Mirobody's render limit
        from datetime import date, timedelta
        today = date.today()
        for ind in wearable[:30]:
            for k in range(0, days, 30):
                a, b = today - timedelta(days=min(days, k + 30) - 1), today - timedelta(days=k)
                text = mcp.call("query_health_indicators", {"indicators": [ind], "view": "day", "start": a.isoformat(), "end": b.isoformat()})
                if "… cut at" in text:
                    mcp.truncated = True
                for r in parse_table(text):
                    d = _pick(r, "period", "date", "day", "time", "start_time")[:10]
                    if d:
                        day_rows.setdefault(d, {"date": d})[ind] = _pick(r, "avg", "value", "mean")
    wp = out_dir / "mirobody_wearable_daily.csv"
    if day_rows:
        cols = ["date"] + wearable
        with open(wp, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for d in sorted(day_rows, reverse=True):
                w.writerow(day_rows[d])
    return {"at": now_iso(), "indicators": len(names), "lab_rows": len(lab_rows), "lab_file": str(lp),
            "wearable_days": len(day_rows), "wearable_file": str(wp) if day_rows else None,
            "truncated": mcp.truncated,
            "note": "rows are Mirobody's; intake still asks you to confirm them like any lab file. Mirobody's raw view "
                    "carries no printed reference range, so ranges come from the checkup PDF if it is also in the folder"}


def url_from(arg: Optional[str]) -> str:
    """--mcp-url-file <path> or the env var LONGPI_MCP_URL; a literal URL on the command line would land in logs."""
    if arg:
        p = Path(arg).expanduser()
        if not p.exists():
            raise LAError(f"{p} not found", EXIT_INPUT)
        return p.read_text(encoding="utf-8").strip()
    env = os.environ.get("LONGPI_MCP_URL", "").strip()
    if not env:
        raise LAError("give --mcp-url-file <file holding the member's MCP URL> or set LONGPI_MCP_URL", EXIT_INPUT)
    return env
