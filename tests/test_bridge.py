"""longpi bridge: Mirobody pull (against a fake MCP server) and the la-export/1 shape."""
from __future__ import annotations

import csv
import json
import sys
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import common as C  # noqa: E402
from lalib import mirobody  # noqa: E402

CATALOG = "indicator|system|code|count|first_date|last_date\n" \
          "低密度脂蛋白胆固醇 LDL-C|mirobody-device|x|2|2025-09-12|2026-09-10\n" \
          "restingHeartRate|mirobody-device|restingHeartRate|90|2026-06-12|2026-09-09\n\n(meta)"
RAW = "(constants: system=mirobody-device)\nindicator|name|time|value|unit|code\n" \
      "低密度脂蛋白胆固醇 LDL-C|低密度脂蛋白胆固醇 LDL-C|2026-09-10 08:00:00|3.8|mmol/L|13457-7\n" \
      "低密度脂蛋白胆固醇 LDL-C|低密度脂蛋白胆固醇 LDL-C|2025-09-12 08:00:00|3.5|mmol/L|13457-7\n\n(window=…)"


def _server(tools=("query_health_indicators",), fail=False):
    calls = []

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_POST(self):
            req = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            calls.append(req)
            if "id" not in req:
                self.send_response(202)
                self.end_headers()
                return
            m = req["method"]
            if m == "initialize":
                res = {"protocolVersion": "2025-03-26"}
            elif m == "tools/list":
                res = {"tools": [{"name": t} for t in tools]}
            else:
                args = req["params"]["arguments"]
                if fail:
                    res = {"isError": True, "content": [{"type": "text", "text": "unknown parameter aggregate"}]}
                else:
                    if not args:
                        table = CATALOG
                    elif args.get("view") == "raw":
                        table = RAW
                    else:
                        table = "indicator|period|avg|unit\nrestingHeartRate|2026-09-09|67.0|bpm\nrestingHeartRate|2026-09-08|71.0|bpm\n\n(m)"
                    res = {"content": [{"type": "text", "text": json.dumps({"result": table, "status": "ok", "truncated": False})}]}
            body = json.dumps({"jsonrpc": "2.0", "id": req["id"], "result": res}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Mcp-Session-Id", "s1")
            self.end_headers()
            self.wfile.write(b"event: message\ndata: " + body + b"\n\n")

    srv = HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_port}/mcp/SECRET123", calls


def test_parse_table_constants_and_escapes():
    rows = mirobody.parse_table("(constants: system=x, n=1)\na|b\n1|two\\|pipes\n\n(meta)\nignored|row")
    assert rows == [{"system": "x", "n": "1", "a": "1", "b": "two|pipes"}]


def test_pull_writes_labs_and_daily_and_never_echoes_secret(tmp_path):
    srv, url, calls = _server()
    try:
        res = mirobody.pull(url, tmp_path, days=30)
    finally:
        srv.shutdown()
    labs = list(csv.DictReader(open(tmp_path / "mirobody_labs.csv", encoding="utf-8")))
    assert [r["value"] for r in labs] == ["3.8", "3.5"] and labs[0]["loinc"] == "13457-7" and labs[0]["date"] == "2026-09-10"
    days = list(csv.DictReader(open(tmp_path / "mirobody_wearable_daily.csv", encoding="utf-8")))
    assert days[0] == {"date": "2026-09-09", "restingHeartRate": "67.0"}
    assert "SECRET123" not in json.dumps(res)
    sent = [c["params"]["arguments"] for c in calls if c.get("method") == "tools/call"]
    assert all(set(a) <= {"keywords", "indicators", "start", "end", "view"} for a in sent)     # Mirobody 1.5.3 schema only


def test_pull_errors_are_clean(tmp_path):
    srv, url, _ = _server(tools=("convert_unit",))
    try:
        with pytest.raises(C.LAError, match="no health data"):
            mirobody.pull(url, tmp_path)
    finally:
        srv.shutdown()
    with pytest.raises(C.LAError) as e:
        mirobody.pull("http://127.0.0.1:9/mcp/SECRET123", tmp_path)
    assert "SECRET123" not in str(e.value)


def test_url_only_from_file_or_env(tmp_path, monkeypatch):
    monkeypatch.delenv("LONGPI_MCP_URL", raising=False)
    with pytest.raises(C.LAError, match="mcp-url-file"):
        mirobody.url_from(None)
    f = tmp_path / "u.txt"
    f.write_text("http://h/mcp/x\n")
    assert mirobody.url_from(str(f)) == "http://h/mcp/x"
