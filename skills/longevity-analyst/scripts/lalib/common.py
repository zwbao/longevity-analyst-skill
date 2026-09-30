"""Shared state handling for the longevity-analyst harness.

The workspace is one directory per member. `state.json` is the only state
file; every command reads it, changes it through a typed function, and
writes it back atomically. The LLM never edits state.json by hand.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

SCHEMA = "la-state/1"
SKILL_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = SKILL_DIR / "data"

STAGES = [
    "intake",
    "preflight",
    "pipelines",
    "methods",
    "integrate",
    "organs",
    "insights",
    "intervene",
    "twin",
    "review",
    "report",
]

# Exit codes shared by every command.
EXIT_OK = 0
EXIT_USAGE = 2
EXIT_INPUT = 3          # input did not pass checks
EXIT_BLOCKED = 4        # a gate refused (license, platform, preflight red)
EXIT_PENDING = 5        # judgments are still pending
EXIT_EXTERNAL = 6       # external tool / network failed
EXIT_BUSY = 7           # another command holds the workspace lock


class LAError(Exception):
    def __init__(self, msg: str, code: int = EXIT_INPUT):
        super().__init__(msg)
        self.code = code


def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def load_json(path: Path) -> Any:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=".tmp-", suffix=".json")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def data(name: str) -> Any:
    return load_json(DATA_DIR / name)


def sha256_file(path: Path, limit: Optional[int] = 4 * 1024 ** 3) -> str:
    """Full hash for files up to `limit` bytes; head+tail+size hash above it."""
    h = hashlib.sha256()
    size = path.stat().st_size
    with open(path, "rb") as fh:
        if limit is None or size <= limit:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
            return h.hexdigest()
        h.update(fh.read(8 << 20))
        fh.seek(-8 << 20, os.SEEK_END)
        h.update(fh.read())
        h.update(str(size).encode())
        return "partial:" + h.hexdigest()


def longevity_skills_home(required: bool = True) -> Optional[Path]:
    """Location of the longevity-skills method library.

    Taken only from LONGEVITY_SKILLS_HOME; never guessed from a private path.
    """
    env = os.environ.get("LONGEVITY_SKILLS_HOME")
    if env and (Path(env) / "catalog.json").exists():
        return Path(env)
    if required:
        raise LAError(
            "LONGEVITY_SKILLS_HOME is not set or has no catalog.json. "
            "Clone https://github.com/zwbao/longevity-skills and export "
            "LONGEVITY_SKILLS_HOME=<that directory>.",
            EXIT_INPUT,
        )
    return None


_KIT = None


def skillkit():
    """The longevity-skills skillkit module (name matching, unit handling)."""
    global _KIT
    if _KIT is None:
        import importlib.util
        home = longevity_skills_home()
        spec = importlib.util.spec_from_file_location("la_skillkit", home / "tools" / "skillkit" / "skillkit.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["la_skillkit"] = mod
        spec.loader.exec_module(mod)  # type: ignore[union-attr]
        _KIT = mod
    return _KIT


# What each stage's outputs depend on. Changing an upstream stage voids everything after it.
DOWNSTREAM = {s: STAGES[i + 1:] for i, s in enumerate(STAGES)}


def invalidate_after(st: Dict[str, Any], stage: str, why: str) -> List[str]:
    """Mark every stage after `stage` as todo and drop approvals that were bound to old content."""
    voided = []
    for s in DOWNSTREAM[stage]:
        if st["stages"].get(s) not in (None, "todo"):
            voided.append(s)
        st["stages"][s] = "todo"
    if "methods" in DOWNSTREAM[stage] and st.get("methods"):
        st["methods"] = None
    if "integrate" in DOWNSTREAM[stage]:
        st.pop("integrate", None)
    if "organs" in DOWNSTREAM[stage]:
        st.pop("organs", None)
    if "insights" in DOWNSTREAM[stage]:
        st.pop("insights", None)
    if "intervene" in DOWNSTREAM[stage]:
        st.pop("intervene", None)
    if "twin" in DOWNSTREAM[stage]:
        st.pop("twin", None)
    if "review" in DOWNSTREAM[stage]:
        st.pop("review", None)
    st.pop("report", None)
    if voided:
        st.setdefault("events", []).append({"t": now_iso(), "event": "invalidated", "stages": voided, "why": why})
    return voided


def require(st: Dict[str, Any], *stages: str) -> None:
    """Refuse to run a stage whose inputs are not finished (stage order is enforced, not advised)."""
    missing = [s for s in stages if st["stages"].get(s) != "done"]
    if missing:
        raise LAError(f"finish {', '.join(missing)} first (see `la.py observe`)", EXIT_BLOCKED)


def sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha_paths(paths: List[Path]) -> Dict[str, str]:
    out = {}
    for p in paths:
        out[str(p)] = sha256_file(p, limit=None) if p.exists() else "missing"
    return out


def validate_member(member: Dict[str, Any]) -> None:
    age = member.get("age")
    if age is not None and not (0 < int(age) < 125):
        raise LAError(f"age {age} is not plausible", EXIT_USAGE)
    if member.get("sex") not in (None, "male", "female"):
        raise LAError(f"sex must be male/female (m/f/男/女), got {member.get('sex')!r}", EXIT_USAGE)
    sd = member.get("sample_date")
    if sd:
        import datetime as _d
        try:
            _d.date.fromisoformat(sd)
        except ValueError:
            raise LAError(f"sample date must be YYYY-MM-DD, got {sd!r}", EXIT_USAGE)


OPEN_WORKSPACES: List["Workspace"] = []


class Workspace:
    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.state_path = self.root / "state.json"
        self._lock_fh = None
        OPEN_WORKSPACES.append(self)

    def _generation(self) -> Optional[str]:
        try:
            return load_json(self.state_path).get("generation")
        except (OSError, ValueError):
            return None

    def lock(self, wait: float = 60.0) -> None:
        """One command at a time per workspace; parallel analysts queue here instead of failing.

        The lock file sits next to the workspace (not inside it), so `init --force` renaming the workspace cannot
        hand a queued command a different member's state; a generation id catches a re-init that happened while
        the command waited.
        """
        import fcntl
        import time as _t
        self.root.parent.mkdir(parents=True, exist_ok=True)
        before = self._generation()
        try:
            self._lock_fh = open(self.root.parent / f".{self.root.name}.la-lock", "w")
        except PermissionError:
            raise LAError(f"cannot create the workspace lock in {self.root.parent} (not writable); "
                          "put the workspace under a writable folder", EXIT_USAGE)
        deadline = _t.time() + wait
        while True:
            try:
                fcntl.flock(self._lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if _t.time() >= deadline:
                    raise LAError("another la.py command is using this workspace; try again when it finishes", EXIT_BUSY)
                _t.sleep(0.2)
        if before is not None and self._generation() != before:
            self.unlock()
            raise LAError("the workspace was re-initialised while this command waited; re-run it deliberately", EXIT_BLOCKED)

    def unlock(self) -> None:
        import fcntl
        if self._lock_fh:
            fcntl.flock(self._lock_fh, fcntl.LOCK_UN)
            self._lock_fh.close()
            self._lock_fh = None

    # ---- paths -------------------------------------------------------
    def p(self, *parts: str) -> Path:
        return self.root.joinpath(*parts)

    # ---- state -------------------------------------------------------
    def exists(self) -> bool:
        return self.state_path.exists()

    def load(self) -> Dict[str, Any]:
        if not self.exists():
            raise LAError(f"no workspace at {self.root}; run `la.py init` first", EXIT_USAGE)
        st = load_json(self.state_path)
        if st.get("schema") != SCHEMA:
            raise LAError(f"state schema {st.get('schema')} != {SCHEMA}", EXIT_INPUT)
        return st

    def save(self, st: Dict[str, Any]) -> None:
        st["updated_at"] = now_iso()
        d = self.root / "deliver"
        if st.get("stages", {}).get("report") != "done" and d.exists():
            # a report that is no longer current is removed: it may hold data later excluded or declined
            import shutil
            shutil.rmtree(d)
        orp = self.root / "work" / "organs" / "organ_readouts.json"
        if orp.exists() and not (st.get("organs") or {}).get("organ_readouts_sha256"):
            orp.unlink()                                # derived from registrations that were voided upstream
        irp = self.root / "work" / "insights" / "insight_readouts.json"
        if irp.exists() and not (st.get("insights") or {}).get("readouts_sha256"):
            irp.unlink()                                # derived from insight runs voided upstream
        if "insights" not in st and (self.root / "work" / "insights").exists():
            # Voided upstream: nothing here is registered any more, so none of it reaches the report. The folder is
            # kept aside (researcher findings cost real work) as insights.previous; the agent may bring files back
            # and register them again, where every check runs again.
            import shutil
            prev = self.root / "work" / "insights.previous"
            if prev.exists():
                shutil.rmtree(prev)
            (self.root / "work" / "insights").rename(prev)
        if "organs" not in st and (self.root / "work" / "organs").exists():
            import shutil                               # upstream changed: bundles and estimates describe old data
            shutil.rmtree(self.root / "work" / "organs")
        from . import pubdata                            # evidence caches written by the harness in this command
        mine = {k: v for k, v in pubdata.WRITTEN.items() if k.startswith(str(self.root))}
        if mine:
            st.setdefault("evidence_ledger", {}).update(mine)
        write_json(self.state_path, st)

    def log(self, st: Dict[str, Any], event: str, **kw: Any) -> None:
        st.setdefault("events", []).append({"t": now_iso(), "event": event, **kw})


def new_state(member: Dict[str, Any], raw_dir: Path) -> Dict[str, Any]:
    import uuid
    return {
        "schema": SCHEMA,
        "generation": uuid.uuid4().hex,
        "created_at": now_iso(),
        "updated_at": now_iso(),
        "member": member,
        "raw_dir": str(raw_dir.resolve()),
        "files": [],
        "labs": [],
        "processed": {},
        "preflight": None,
        "pipelines": [],
        "methods": None,
        "judgments": [],
        "stages": {s: "todo" for s in STAGES},
        "events": [],
    }


def pending_judgments(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [j for j in st.get("judgments", []) if j.get("status") == "pending"]


def add_judgment(st: Dict[str, Any], kind: str, target: str, question: str,
                 options: Optional[List[str]] = None, hints: Optional[List[str]] = None) -> None:
    for j in st.setdefault("judgments", []):
        if j["kind"] == kind and j["target"] == target and j["status"] == "pending":
            return
    st["judgments"].append({
        "id": f"J{len(st['judgments']) + 1:03d}",
        "kind": kind,
        "target": target,
        "question": question,
        "options": options or [],
        "hints": hints or [],
        "status": "pending",
    })


def resolve_judgment(st: Dict[str, Any], kind: str, target: str, answer: str, reason: str) -> None:
    """Resolve the pending judgment (kind, target). Raises if there is none to resolve."""
    if not reason or not reason.strip():
        raise LAError("a judgment needs --reason: say what you read that decided it", EXIT_USAGE)
    for j in st.get("judgments", []):
        if j["kind"] == kind and j["target"] == target and j["status"] == "pending":
            if j["options"] and answer not in j["options"]:
                raise LAError(f"answer {answer!r} not in options {j['options']}", EXIT_USAGE)
            j.update(status="resolved", answer=answer, reason=reason, resolved_at=now_iso())
            return
    raise LAError(f"no pending {kind} judgment for {target}", EXIT_USAGE)


def emit(obj: Any) -> None:
    json.dump(obj, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")


def file_by_id(st: Dict[str, Any], fid: str) -> Dict[str, Any]:
    for f in st["files"]:
        if f["id"] == fid:
            return f
    raise LAError(f"no file with id {fid}", EXIT_USAGE)
