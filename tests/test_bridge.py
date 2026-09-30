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
          "脂蛋白(a) Lp(a)|mirobody-device|x|1|2026-09-10|2026-09-10\n" \
          "空腹血糖|mirobody-device|x|20|2026-07-01|2026-09-10\n" \
          "restingHeartRate|mirobody-device|restingHeartRate|90|2026-06-12|2026-09-09\n\n(meta)"
ROWS = {"低密度脂蛋白胆固醇 LDL-C": [("2026-09-10 08:00:00", "3.8"), ("2025-09-12 08:00:00", "3.5")],
        "脂蛋白(a) Lp(a)": [("2026-09-10 08:00:00", "38")],
        "空腹血糖": [(f"2026-{7 + i // 10:02d}-{1 + i % 10:02d} 08:00:00", "5.9") for i in range(20)]}


def _raw(names):
    rows = [(n, t, v) for n in names for t, v in ROWS.get(n, [])]
    if len(rows) == 1:                                   # every column constant: Mirobody prints no header at all
        n, t, v = rows[0]
        return f"(constants: indicator={n}, name={n}, time={t}, value={v}, unit=mg/dL, system=mirobody-device, code=x)\n\n(window=…)\nnotes: …"
    body = "\n".join(f"{n}|{n}|{t}|{v}|u|13457-7" for n, t, v in rows)
    return "(constants: system=mirobody-device)\nindicator|name|time|value|unit|code\n" + body + "\n\n(window=…)"


def _server(tools=("query_health_indicators",), fail=False, cut_over=2):
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
                truncated = False
                if fail:
                    res = {"isError": True, "content": [{"type": "text", "text": "unknown parameter aggregate"}]}
                    table = None
                elif not args:
                    table = CATALOG
                elif args.get("view") == "raw":
                    names = args["indicators"]
                    truncated = len(names) > cut_over               # like the 250-row cap on a big batch
                    table = _raw(names) if not truncated else _raw(names[:1]) + "\n… cut at 8000 characters"
                else:                                               # day view, one indicator: no indicator column
                    table = "(constants: n=1, unit=/min, system=mirobody-device, code=restingHeartRate)\nperiod|avg|min|max\n" \
                            "2026-09-09|67.0|67|67\n2026-09-08|71.0|71|71\n\n(window=…)"
                if table is not None:
                    res = {"content": [{"type": "text", "text": json.dumps({"result": table, "status": "ok", "truncated": truncated})}]}
            body = json.dumps({"jsonrpc": "2.0", "id": req["id"], "result": res}).encode()
            note = json.dumps({"jsonrpc": "2.0", "method": "notifications/message", "params": {}}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Mcp-Session-Id", "s1")
            self.end_headers()
            self.wfile.write(b"event: message\ndata: " + note + b"\n\nevent: message\ndata: " + body + b"\n\n")

    srv = HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_port}/mcp/SECRET123", calls


def test_parse_table_real_shapes():
    assert mirobody.parse_table("(constants: system=x, n=1)\na|b\n1|two|pipes\n\n(meta)\nignored|row") == \
        [{"system": "x", "n": "1", "a": "1", "b": "two|pipes"}]           # values are not escaped: extra cells join the last
    one = "(constants: indicator=脂蛋白(a) Lp(a), value=38, unit=mg/dL)\n\n(window=…)\nnotes: no data means …"
    assert mirobody.parse_table(one) == [{"indicator": "脂蛋白(a) Lp(a)", "value": "38", "unit": "mg/dL"}]
    assert mirobody.parse_table("a|b\n1|2\n… cut at 8000 characters\n") == [{"a": "1", "b": "2"}]
    assert mirobody.parse_table("(no rows)\n\n(window=…)") == []


def test_pull_writes_per_date_labs_daily_series_and_never_echoes_secret(tmp_path):
    (tmp_path / "mirobody_labs.csv").write_text("stale")               # an older pull's file is replaced
    (tmp_path / "member_own.csv").write_text("keep")
    srv, url, calls = _server()
    try:
        res = mirobody.pull(url, tmp_path, days=30)
    finally:
        srv.shutdown()
    names = sorted(p.name for p in tmp_path.iterdir())
    assert "mirobody_labs.csv" not in names and "member_own.csv" in names
    new = list(csv.DictReader(open(tmp_path / "mirobody_labs_2026-09-10.csv", encoding="utf-8")))
    old = list(csv.DictReader(open(tmp_path / "mirobody_labs_2025-09-12.csv", encoding="utf-8")))
    assert {r["marker"]: r["value"] for r in new}["脂蛋白(a) Lp(a)"] == "38"   # a one-row table is not lost
    assert [r["value"] for r in old] == ["3.5"]                             # another year is another file
    assert sum(f["rows"] for f in res["lab_files"]) == 23                    # 2 + 1 + 20 glucose (a repeated lab, not a watch)
    days = list(csv.DictReader(open(tmp_path / "mirobody_wearable_daily.csv", encoding="utf-8")))
    assert days[0] == {"date": "2026-09-09", "restingHeartRate (/min)": "67.0"}
    assert res["wearable_indicators"] == ["restingHeartRate"]
    assert "SECRET123" not in json.dumps(res)
    sent = [c["params"]["arguments"] for c in calls if c.get("method") == "tools/call"]
    assert all(set(a) <= {"keywords", "indicators", "start", "end", "view"} for a in sent)     # Mirobody 1.5.3 schema only
    assert any(len(a.get("indicators", [])) == 1 and a.get("view") == "raw" for a in sent)      # a cut batch was split


def test_pull_stops_when_one_indicator_is_still_cut(tmp_path):
    srv, url, _ = _server(cut_over=0)
    try:
        with pytest.raises(C.LAError, match="even alone"):
            mirobody.pull(url, tmp_path)
    finally:
        srv.shutdown()
    assert not list(tmp_path.glob("mirobody_*.csv"))                         # nothing written after a failed read


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


def test_export_text_is_plain_and_keeps_no_amounts():
    from lalib import export
    ro = {"native.egfr": {"id": "native.egfr", "label_zh": "eGFR", "value": 88.0, "unit": "mL/min/1.73m²"},
          "organ.kidney.risk.1": {"id": "organ.kidney.risk.1", "label_zh": "慢性肾病", "value": 0.08, "unit": "概率",
                                  "kind": "llm_estimate", "low": 0.04, "high": 0.16, "horizon_years": 10}}
    t = export.plain("每日食盐不超过 {{n:5 克}}，饮水 {{n:1500 毫升}}，坚持 {{n:12 周}}，参见 {{pmid:34554658}}；eGFR {{r:native.egfr}}", ro)
    assert "5 克" not in t and "1500" not in t and t.count("（具体量见报告）") == 2
    assert "12 周" in t and "PMID 34554658" in t and "](" not in t and "88 mL/min" in t
    risk = export.plain("风险 {{r:organ.kidney.risk.1}}", ro)
    assert "AI 估计" in risk
    long = "每天快走三十分钟以上并逐步增加到每周五次，同时减少久坐时间，每坐一小时起身活动五分钟，晚饭后散步二十分钟，周末安排一次较长距离的户外徒步活动"
    title = export._title(long)
    assert len(title) <= 60 and long.startswith(title) and title[-1] not in "，、"
