"""Minimal synthetic workspaces for the clinical lane's unit tests (synthetic data, not a real member)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))
from lalib import common as C  # noqa: E402


def make_ws(tmp: Path, labs=(), age=49, sex="male", sample_date="2026-09-01", answers=None, readouts=()):
    ws = tmp / "ws"
    (ws / "work").mkdir(parents=True, exist_ok=True)
    st = C.new_state({"id": "synthetic", "age": age, "sex": sex, "sample_date": sample_date, "answers": answers or {}}, tmp)
    for row in labs:
        marker, value, unit = row[:3]
        r = {"marker": marker, "value": str(value), "unit": unit, "source_file": "F001", "maps_to": marker}
        if len(row) > 3:
            r["collected_at"] = row[3]
        st["labs"].append(r)
    C.write_json(ws / "work" / "readouts.json", {"readouts": list(readouts)})
    C.write_json(ws / "state.json", st)
    return ws, st


def write_plan(ws: Path, items):
    C.write_json(ws / "work" / "intervene" / "plan.json", {"items": items})


SALT_ITEM = {"id": "I1", "category": "diet", "action_zh": "减盐", "rationale_zh": "收缩压偏高。", "targets": ["收缩压"],
             "evidence": [{"type": "effects", "ref": "salt-reduction-sbp"}], "evidence_grade": "human_rct",
             "executor": "nutritionist", "retest": {"what": "收缩压", "after_weeks": 6}}


def series(ws: Path, markers):
    C.write_json(ws / "work" / "twin" / "series.json", {"schema": "la-series/1", "markers": markers})
