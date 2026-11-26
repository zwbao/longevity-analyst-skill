"""Minimal stand-in for the engine's Library: only the root that holds data/reference/ and a provenance stamp."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


class LKError(Exception):
    def __init__(self, message: str, code: str = "error", exit_code: int = 1):
        super().__init__(message)
        self.code = code
        self.exit_code = exit_code


class Library:
    def __init__(self, root: Path):
        self.root = Path(root)

    def stamp(self) -> Dict[str, Any]:
        return {"provided_by": "longevity-analyst age-evidence"}
