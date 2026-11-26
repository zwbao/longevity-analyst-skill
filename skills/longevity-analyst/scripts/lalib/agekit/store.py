"""Minimal stand-in for the engine's store: panel text comes straight from the payload."""
from typing import Optional


def resolve(payload: dict) -> Optional[str]:
    return payload.get("panel_text")
