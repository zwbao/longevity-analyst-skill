"""Wearable daily summaries: the agent maps the file's columns (a judgment about what each column is); the harness
computes 30- and 90-day means and the 90-day trend per 30 days. Descriptive only."""
from __future__ import annotations

import csv
import json
import statistics
from datetime import date
from pathlib import Path
from typing import Any, Dict, List

from .common import EXIT_INPUT, LAError, now_iso

FIELDS = {"steps": ("日均步数", "步"), "rhr": ("静息心率", "次/分"), "hrv": ("心率变异性 HRV", "ms"),
          "sleep_h": ("睡眠时长", "小时"), "deep_h": ("深睡时长", "小时"), "mvpa_min": ("中高强度运动", "分钟/天"),
          "spo2_min": ("夜间最低血氧", "%")}


def summarize(st: Dict[str, Any], ws: Path, file_id: str, mapping: Dict[str, str]) -> Dict[str, Any]:
    f = next((x for x in st["files"] if x["id"] == file_id and not x.get("excluded")), None)
    if not f:
        raise LAError(f"no usable file {file_id}", EXIT_INPUT)
    if "date" not in mapping or not set(mapping) - {"date"} <= set(FIELDS):
        raise LAError(f"--map needs 'date' plus any of {sorted(FIELDS)} -> column names", EXIT_INPUT)
    with open(f["path"], encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    missing = [c for c in mapping.values() if rows and c not in rows[0]]
    if missing:
        raise LAError(f"columns {missing} are not in {f['name']}", EXIT_INPUT)
    series: Dict[str, List] = {k: [] for k in mapping if k != "date"}
    for r in rows:
        try:
            d = date.fromisoformat(r[mapping["date"]].strip()[:10])
        except ValueError:
            continue
        for k in series:
            try:
                series[k].append((d, float(r[mapping[k]])))
            except (ValueError, TypeError):
                pass
    readouts = []
    for k, pts in series.items():
        if len(pts) < 14:
            continue
        pts.sort()
        last = pts[-1][0]
        m30 = [v for d, v in pts if (last - d).days < 30]
        m90 = [v for d, v in pts if (last - d).days < 90]
        xs = [(d - pts[0][0]).days for d, _ in pts]
        ys = [v for _, v in pts]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        sxx = sum((x - mx) ** 2 for x in xs)
        slope30 = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx * 30) if sxx else 0.0
        label, unit = FIELDS[k]
        base = {"kind": "descriptive", "method": "native.wearable", "source_file": file_id, "systems": [], "days": len(pts)}
        readouts += [{"id": f"wear.{k}.mean30", "label_zh": f"{label}（近 30 天均值）", "value": round(statistics.mean(m30), 2), "unit": unit, **base},
                     {"id": f"wear.{k}.mean90", "label_zh": f"{label}（近 90 天均值）", "value": round(statistics.mean(m90), 2), "unit": unit, **base},
                     {"id": f"wear.{k}.trend30", "label_zh": f"{label}（每 30 天变化趋势）", "value": round(slope30, 3), "unit": f"{unit}/30天", **base}]
    if not readouts:
        raise LAError("fewer than 14 dated values in every mapped column", EXIT_INPUT)
    return {"file": file_id, "mapping": mapping, "readouts": readouts, "at": now_iso()}
