#!/usr/bin/env python3
"""longevity-analyst harness CLI.

Every state change goes through one of these typed commands; the agent never
edits <workspace>/state.json by hand. Run `la.py observe <ws>` at the start
of every step: it re-renders what is done, what is pending and what to do next.
Stage order is enforced, and changing an earlier stage voids the later ones.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lalib import board, causal, genomics, reference, wearable
from lalib import common as C  # noqa: E402
from lalib import evidence, integrate, intake, methods, organs, pipelines, preflight, report, twin  # noqa: E402

SEX = {"m": "male", "male": "male", "男": "male", "f": "female", "female": "female", "女": "female"}


def _open(a, lock: bool = True):
    ws = C.Workspace(Path(a.workspace))
    if lock:
        ws.lock()
    return ws, ws.load()


# ---------------------------------------------------------------- setup
def cmd_init(a):
    ws = C.Workspace(Path(a.workspace))
    raw = Path(a.raw_dir)
    if not raw.is_dir():
        raise C.LAError(f"raw dir {raw} not found", C.EXIT_USAGE)
    if not re.fullmatch(r"[A-Za-z0-9_\-]{1,40}", a.member_id or ""):
        raise C.LAError("--member-id must be 1-40 letters, digits, _ or -", C.EXIT_USAGE)
    if ws.root.exists() and any(ws.root.iterdir()):
        if not a.force:
            raise C.LAError(f"{ws.root} is not empty (use --force to start over; the old workspace is moved aside, not reused)", C.EXIT_USAGE)
        ws.lock(wait=0)                                    # nobody may be using the old workspace
        try:
            old = C.load_json(ws.state_path) if ws.state_path.exists() else {}
        except ValueError:
            old = {}
        live = [p["run_id"] for p in old.get("pipelines", []) if pipelines.run_alive(p)]
        if live:
            raise C.LAError(f"pipelines still running in this workspace ({live}); `pipeline stop` them first", C.EXIT_BLOCKED)
        if not ws.state_path.exists() or not old:
            import subprocess as _sp
            ps = _sp.run(["ps", "-axo", "command"], capture_output=True, text=True).stdout
            if str(ws.root) in ps:
                raise C.LAError("a process is still using files in this workspace; stop it first", C.EXIT_BLOCKED)
        backup = ws.root.with_name(f"{ws.root.name}.old-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}-{os.getpid()}")
        n = 1
        while backup.exists():
            backup = backup.with_name(f"{backup.name}-{n}")
            n += 1
        ws.root.rename(backup)
        ws.unlock()
        print(f"moved the previous workspace to {backup}", file=sys.stderr)
    lib = a.longevity_skills or os.environ.get("LONGEVITY_SKILLS_HOME")
    if not lib or not (Path(lib) / "catalog.json").exists():
        raise C.LAError("pass --longevity-skills <path to a longevity-skills clone> (it must contain catalog.json)", C.EXIT_USAGE)
    sex = SEX.get((a.sex or "").lower(), a.sex)
    member = {"id": a.member_id, "age": a.age, "sex": sex, "sample_date": a.sample_date, "mode": a.mode, "answers": {}}
    C.validate_member(member)
    ws.root.mkdir(parents=True, exist_ok=True)
    ws.lock()
    st = C.new_state(member, raw)
    st["workspace"] = str(ws.root)
    st["longevity_skills_home"] = str(Path(lib).resolve())
    ws.log(st, "init", mode=a.mode)
    ws.save(st)
    C.emit({"workspace": str(ws.root), "member": member})


def cmd_member_set(a):
    ws, st = _open(a)
    reserved = set(C.data("license_policy.json")["reserved_answer_keys"])
    for kv in a.pairs:
        if "=" not in kv:
            raise C.LAError(f"expected key=value, got {kv}", C.EXIT_USAGE)
        k, v = kv.split("=", 1)
        if k in reserved:
            raise C.LAError(f"{k} is chosen by the harness, not by answers", C.EXIT_USAGE)
        if k == "id":
            raise C.LAError("the member id is fixed at init; start a new workspace for another person", C.EXIT_USAGE)
        if k == "age":
            try:
                st["member"]["age"] = int(v)
            except ValueError:
                raise C.LAError(f"age must be a whole number, got {v!r}", C.EXIT_USAGE)
        elif k in ("sex", "sample_date"):
            st["member"][k] = SEX.get(v.lower(), v) if k == "sex" else v
        else:
            if k == "genetic_disclosure" and v not in ("yes", "no"):
                raise C.LAError("genetic_disclosure must be yes or no", C.EXIT_USAGE)
            st["member"].setdefault("answers", {})[k] = v
    C.validate_member(st["member"])
    known = {"genetic_disclosure"} | {i["key"] for sk in C.load_json(C.longevity_skills_home() / "catalog.json")["skills"]
                                      for i in (sk.get("inputs") or []) if i.get("from") in ("argument", "profile")}
    unknown = [kv.split("=", 1)[0] for kv in a.pairs if kv.split("=", 1)[0] not in known | {"age", "sex", "sample_date"}]
    C.invalidate_after(st, "pipelines", f"member facts changed: {a.pairs}")      # answers feed methods, not intake
    ws.log(st, "member_set", pairs=a.pairs, source=a.source)
    ws.save(st)
    C.emit({**st["member"], **({"warning_unknown_keys": unknown,
                                "note": "no method reads these keys; check the key names `methods plan` asks for"} if unknown else {})})


# ---------------------------------------------------------------- intake
def _intake_done(st) -> None:
    from lalib import labnames
    if C.pending_judgments(st):
        st["stages"]["intake"] = "pending_judgments"
    elif labnames.pending(st):
        st["stages"]["intake"] = "pending_lab_confirmations"
    else:
        st["stages"]["intake"] = "done"


def cmd_intake(a):
    ws, st = _open(a)
    res = intake.run_intake(st, Path(st["raw_dir"]))
    _intake_done(st)
    if res["added"]:
        C.invalidate_after(st, "intake", "new files")
    ws.log(st, "intake", **res)
    ws.save(st)
    C.emit({**res, "pending_judgments": C.pending_judgments(st),
            "pending_lab_confirmations": __import__("lalib.labnames", fromlist=["pending"]).pending(st),
            "files": [{k: f.get(k) for k in ("id", "name", "kind", "status", "reason")} for f in st["files"]]})


def cmd_assign(a):
    ws, st = _open(a)
    if a.what == "identity":
        if a.value not in ("consistent", "inconsistent", "cannot_tell"):
            raise C.LAError("--what identity takes --value consistent | inconsistent | cannot_tell", C.EXIT_USAGE)
        excl = [x.strip() for x in (a.files or "").split(",") if x.strip()]
        if a.value == "inconsistent" and not excl and not st.get("identity"):
            raise C.LAError("--value inconsistent needs --files F00x,F00y: the files that are not this member's", C.EXIT_USAGE)
        if a.value == "cannot_tell" and not [x for x in (a.uncertain or "").split(",") if x.strip()]:
            raise C.LAError("--value cannot_tell needs --uncertain F00x,...: the files whose ownership cannot be confirmed", C.EXIT_USAGE)
        if any(j["kind"] == "identity" and j["status"] == "pending" for j in st["judgments"]):
            C.resolve_judgment(st, "identity", "member", a.value, a.reason)
        prev = st.get("identity")
        if prev:
            st.setdefault("identity_history", []).append(prev)      # a revision keeps the earlier verdict on record
        unc = [x.strip() for x in (a.uncertain or "").split(",") if x.strip()]
        if excl:
            intake.exclude_files(st, excl, a.reason)
        if unc:
            intake.mark_uncertain(st, unc, a.reason)
        st["identity"] = {"answer": a.value, "reason": a.reason,
                          "excluded_files": sorted(set(excl) | set((prev or {}).get("excluded_files", []))),
                          "uncertain_files": sorted((set(unc) | set((prev or {}).get("uncertain_files", [])))
                                                    - set(excl) - set((prev or {}).get("excluded_files", []))), "at": C.now_iso()}
        if a.value == "consistent" and st["identity"]["excluded_files"]:
            st["identity"]["answer"] = "inconsistent"      # exclusions are never undone by a later "consistent"
    elif a.what == "transcribe":
        if a.value == "transcribed":
            raise C.LAError("transcribed rows are added with `la.py labs add` (it resolves this judgment)", C.EXIT_USAGE)
        f = C.file_by_id(st, a.file)
        C.resolve_judgment(st, "transcribe", a.file, a.value, a.reason)
        f["status"] = "ready" if a.value == "no_lab_values" else "unreadable"
        f["reason"] = a.reason
    elif a.what == "primary":
        intake.choose_primary(st, a.file, a.value, a.reason)
    else:
        intake.assign(st, a.file, a.what, a.value, a.reason)
    _intake_done(st)
    C.invalidate_after(st, "intake", f"judgment {a.what}={a.value}")
    ws.log(st, "assign", file=a.file, what=a.what, value=a.value, reason=a.reason)
    ws.save(st)
    C.emit({"what": a.what, "value": a.value, "pending": len(C.pending_judgments(st))})


def _check_map_target(key: str) -> None:
    """A lab row may only be mapped to a routine-lab input of a verified method (never a protein or CpG input) or to an
    organ-index analyte name in data/organs.json."""
    kit = C.skillkit()
    for a in C.data("organs.json")["lab_aliases"].values():
        if kit.fold_name(key) in {kit.fold_name(n) for n in a["names"]}:
            return
    cat = C.load_json(C.longevity_skills_home() / "catalog.json")
    for s in cat["skills"]:
        if s.get("inputs_status") != "verified" or "蛋白、代谢与器官年龄" in s.get("domains", []):
            continue
        for i in s.get("inputs") or []:
            if i.get("from") == "measurements" and i["key"] == key and not i["key"].startswith("cg"):
                return
    raise C.LAError(f"{key!r} is not a routine-lab input key of any verified method", C.EXIT_USAGE)


def cmd_labs(a):
    if a.action == "candidates":
        ws, st = _open(a, lock=False)
        from lalib import labnames
        C.emit({"pending": labnames.pending(st),
                "how": "every numeric row needs an answer. yes = the row is from THIS visit (the sample date) and is what its name "
                       "says, measured the standard way (serum/plasma for blood chemistry, venous blood for blood counts, fasting "
                       "where the candidate says so); no = an earlier or later visit, urine, post-meal, OGTT, random or capillary "
                       "values, another analyte, or you cannot tell. "
                       "`candidates` is what the harness would match the row to (empty: the harness does not recognise the name; "
                       "answer from what the row really is). "
                       "Answer with `la.py labs confirm <ws> --answers <json> --reason <why>`: "
                       "{\"answers\": [{\"row_key\": \"…\", \"answer\": \"yes|no\", \"why\": \"…\"}]}"})
        return
    if not (a.reason or "").strip():
        raise C.LAError(f"labs {a.action} needs --reason", C.EXIT_USAGE)
    ws, st = _open(a)
    if a.action == "confirm":
        from lalib import labnames
        if not a.answers:
            raise C.LAError("labs confirm needs --answers <json>", C.EXIT_USAGE)
        try:
            doc = C.load_json(Path(a.answers))
        except (OSError, ValueError) as e:
            raise C.LAError(f"{a.answers} is not JSON ({e})", C.EXIT_INPUT)
        ans = doc.get("answers") if isinstance(doc, dict) else None
        if not isinstance(ans, list) or not ans:
            raise C.LAError("answers JSON needs a non-empty `answers` list", C.EXIT_INPUT)
        by_key = labnames.keyed(st)
        dup = [k for k in {x.get("row_key") for x in ans if isinstance(x, dict)} if sum(1 for x in ans if isinstance(x, dict) and x.get("row_key") == k) > 1]
        if dup:
            raise C.LAError(f"row_key answered twice: {dup}", C.EXIT_INPUT)
        problems = []
        for i, x in enumerate(ans):
            if not isinstance(x, dict) or x.get("row_key") not in by_key or x.get("answer") not in ("yes", "no") \
                    or not isinstance(x.get("why"), str) or not x["why"].strip():
                problems.append(f"answers[{i}]: needs a current row_key (from `labs candidates`), answer yes|no and why")
        if problems:
            raise C.LAError("\n".join(problems), C.EXIT_INPUT)
        for x in ans:
            by_key[x["row_key"]]["confirm"] = {"row_key": x["row_key"], "answer": x["answer"], "why": x["why"], "at": C.now_iso()}
        res = {"confirmed": sum(x["answer"] == "yes" for x in ans), "refused": sum(x["answer"] == "no" for x in ans),
               "still_pending": len(labnames.pending(st))}
    elif a.action == "add":
        if not a.csv or not a.file:
            raise C.LAError("labs add needs --csv and --file", C.EXIT_USAGE)
        res = {"added": intake.add_labs(st, Path(a.csv), a.file, a.reason), "labs_total": len(st["labs"])}
    elif a.action == "map":
        if not a.marker or not a.to:
            raise C.LAError("labs map needs --marker and --to", C.EXIT_USAGE)
        _check_map_target(a.to)
        hit = [r for r in st["labs"] if r["marker"] == a.marker]
        if a.row:
            hit = hit[a.row - 1:a.row] if 0 < a.row <= len(hit) else []
        if not hit:
            raise C.LAError(f"no lab row named {a.marker!r}" + (f" number {a.row}" if a.row else ""), C.EXIT_USAGE)
        refused = [r for r in hit if (r.get("confirm") or {}).get("answer") == "no"]
        if refused:
            raise C.LAError(f"{len(refused)} of these rows were answered `no` in labs confirm; pick the row with --row, "
                            "or answer it again with labs confirm first", C.EXIT_USAGE)
        for r in hit:
            r["maps_to"] = a.to
            r["map_reason"] = a.reason
        res = {"mapped": len(hit), "marker": a.marker, "to": a.to}
    else:
        if not a.marker:
            raise C.LAError(f"labs {a.action} needs --marker", C.EXIT_USAGE)
        res = {a.action: intake.fix_lab(st, a.marker, a.action, a.value, a.unit, a.reason, row=a.row), "marker": a.marker}
    _intake_done(st)
    C.invalidate_after(st, "intake", f"labs {a.action}")
    ws.log(st, f"labs_{a.action}", reason=a.reason)
    ws.save(st)
    C.emit(res)


# ---------------------------------------------------------------- pipelines
def cmd_preflight(a):
    ws, st = _open(a)
    C.require(st, "intake")
    pf = preflight.run_preflight(st, ws.p("work"), check_network=not a.no_network)
    st["preflight"] = pf
    st["stages"]["preflight"] = "done"
    if not [r for r in pf["pipelines"] if not r["optional"]] and \
            all(p["status"] in ("verified", "skipped", "void_excluded", "void_stale") for p in st["pipelines"] if not p.get("optional")):
        st["stages"]["pipelines"] = "done"          # nothing required to run is a finished stage, not a stuck one
    md = preflight.render_md(pf)
    ws.p("work").mkdir(parents=True, exist_ok=True)
    ws.p("work", "preflight.md").write_text(md, encoding="utf-8")
    ws.log(st, "preflight", overall=pf["overall"])
    ws.save(st)
    print(md)


def cmd_pipeline(a):
    ws, st = _open(a)
    C.require(st, "intake", "preflight")
    if a.action == "plan":
        res = pipelines.plan(st, ws.root, pgs_id=a.pgs_id)
        res = [{k: p[k] for k in ("run_id", "pipeline", "revision", "tier", "optional", "estimate", "blocked", "cmd")} for p in res]
    elif a.action == "run":
        if not a.run_id:
            raise C.LAError("--run-id is required", C.EXIT_USAGE)
        res = pipelines.launch(st, a.run_id, a.confirmed or "", stub=a.stub, force_yellow=a.accept_yellow, reopen=a.reopen)
    elif a.action == "status":
        res = [pipelines.status(st, r["run_id"]) for r in st["pipelines"] if not a.run_id or r["run_id"] == a.run_id]
    elif a.action == "stop":
        res = pipelines.stop(st, a.run_id)
    elif a.action == "skip":
        res = pipelines.skip(st, a.run_id, a.reason or "")
    else:
        res = pipelines.verify(st, a.run_id)
        if res.get("ok"):
            C.invalidate_after(st, "pipelines", f"verified {a.run_id}")
    closed = ("verified", "skipped", "void_excluded", "void_stale")
    required = [p for p in st["pipelines"] if not p.get("optional")]
    optional = [p for p in st["pipelines"] if p.get("optional")]
    if all(p["status"] in closed for p in required) and all(p["status"] in closed + ("planned",) for p in optional):
        st["stages"]["pipelines"] = "done"
    else:
        st["stages"]["pipelines"] = "todo"
    ws.log(st, f"pipeline_{a.action}", run_id=a.run_id)
    ws.save(st)
    C.emit(res)


# ---------------------------------------------------------------- analysis
def cmd_methods(a):
    ws, st = _open(a)
    C.require(st, "intake", "preflight", "pipelines")
    if a.action == "plan":
        res = methods.plan(st, ws.root)
        C.invalidate_after(st, "methods", "methods re-planned")
        st["stages"]["methods"] = "planned"
        ws.save(st)
        C.emit({"counts": res["counts"], "clearance": res.get("clearance"), "unmatched_lab_names": res["unmatched_lab_names"],
                "withheld_lab_rows": res.get("withheld_lab_rows") or None,
                "items": [{k: i.get(k) for k in ("method", "tier", "status", "reasons", "needs", "problems", "optional_unanswered") if i.get(k)} for i in res["items"]]})
        return
    if not st.get("methods"):
        raise C.LAError("run `la.py methods plan` first", C.EXIT_BLOCKED)
    res = methods.run(st, ws.root, None, checkpoint=lambda: ws.save(st))
    C.invalidate_after(st, "methods", "methods ran")
    st["stages"]["methods"] = "done"
    ws.log(st, "methods_run")
    ws.save(st)
    C.emit(res)


def cmd_integrate(a):
    ws, st = _open(a)
    C.require(st, "methods")
    if a.action == "bundle":
        res = integrate.bundle(st, ws.root)
        if not res["bundles"]:
            st["stages"]["integrate"] = "done"
    else:
        if not a.system:
            raise C.LAError("--system is required", C.EXIT_USAGE)
        res = integrate.register(st, ws.root, a.system)
        C.invalidate_after(st, "organs", f"analysis {a.system} registered")   # organ bundles do not read the analyses
        made = {b["system"] for b in st.get("integrate", {}).get("bundles", [])}
        if made <= set(st["integrate"].get("analyses", {})):
            st["stages"]["integrate"] = "done"
    ws.log(st, f"integrate_{a.action}", system=a.system)
    ws.save(st)
    C.emit(res)


def cmd_evidence(a):
    ws, st = _open(a, lock=False)
    if a.action == "local":
        C.emit(evidence.local([t.strip() for t in a.terms.split(",") if t.strip()]))
    elif a.action == "pubmed":
        C.emit(evidence.pubmed(ws.root, a.query, a.max))
    else:
        C.emit(evidence.fetch_url(ws.root, a.url))


def cmd_organ(a):
    ws, st = _open(a)
    C.require(st, "integrate")
    spec = C.data("organs.json")["organs"]
    if a.action == "bundle":
        res = organs.bundle(st, ws.root)
    elif a.action == "register":
        if not a.organ:
            raise C.LAError("--organ is required", C.EXIT_USAGE)
        try:
            res = organs.register(st, ws.root, a.organ, pmid_check=evidence.verify_pmids)
        except C.LAError:
            if st["stages"].get("organs") == "done":
                st["stages"]["organs"] = "todo"
            C.invalidate_after(st, "organs", f"organ {a.organ} registration rejected")
            ws.save(st)                               # the voided earlier estimate is recorded, not only on disk
            raise
        C.invalidate_after(st, "organs", f"organ {a.organ} registered")
    else:
        if not a.organ or not (a.reason or "").strip():
            raise C.LAError("organ skip needs --organ and --reason", C.EXIT_USAGE)
        res = organs.skip(st, ws.root, a.organ, a.reason)
        C.invalidate_after(st, "organs", f"organ {a.organ} skipped")
    done = set((st.get("organs") or {}).get("estimates", {})) | set((st.get("organs") or {}).get("skipped", {}))
    if (st.get("organs") or {}).get("bundles") and set(spec) <= done:
        st["stages"]["organs"] = "done"
    elif st["stages"].get("organs") == "done":             # an organ was voided: the table is no longer complete
        st["stages"]["organs"] = "todo"
        C.invalidate_after(st, "organs", f"organ {a.action} left organs incomplete")
    ws.log(st, f"organ_{a.action}", organ=a.organ)
    ws.save(st)
    C.emit(res)



def _insight_readouts(st, ws, group, rows):
    """Replace one group's rows (genomics / position / wearable) in work/insights/insight_readouts.json."""
    p = ws.root / "work" / "insights" / "insight_readouts.json"
    old = C.load_json(p)["readouts"] if p.exists() else []
    keep = [r for r in old if r.get("group") != group]
    new = [{**r, "group": group} for r in rows]
    C.write_json(p, {"generated_at": C.now_iso(), "readouts": keep + new})
    ins = st.setdefault("insights", {})
    ins["readouts_sha256"] = C.sha256_file(p, limit=None)
    ins.setdefault("done", {})[group] = C.now_iso()


def _insights_stage(st):
    ins = st.get("insights") or {}
    d = ins.get("done") or {}
    genomics_ok = "genomics" in d or "genomics_skipped" in ins
    if "position" in d and genomics_ok and board.complete(st):
        st["stages"]["insights"] = "done"
    elif st["stages"].get("insights") == "done":
        st["stages"]["insights"] = "todo"
        C.invalidate_after(st, "insights", "insights incomplete")


def cmd_insights(a):
    ws, st = _open(a)
    C.require(st, "organs")
    if a.action == "genomics":
        res = genomics.explain(st, ws.root, absent_as_ref=a.absent_as_ref)
        if a.absent_as_ref:                          # a recorded judgment, shown in the report next to the scores
            st.setdefault("insights", {})["absent_as_ref"] = {"reason": a.absent_as_ref, "at": C.now_iso()}
        else:
            (st.get("insights") or {}).pop("absent_as_ref", None)
        _insight_readouts(st, ws, "genomics", res.pop("readouts"))
    elif a.action == "skip-genomics":
        if not (a.reason or "").strip():
            raise C.LAError("skip-genomics needs --reason (e.g. no VCF, member declined genetic results)", C.EXIT_USAGE)
        st.setdefault("insights", {})["genomics_skipped"] = {"reason": a.reason, "at": C.now_iso()}
        res = {"skipped": "genomics"}
    elif a.action == "position":
        res = reference.position(st, ws.root)
        _insight_readouts(st, ws, "position", res.pop("readouts"))
    elif a.action == "wearable":
        if not a.file or not a.map:
            raise C.LAError("wearable needs --file and --map '{\"date\": \"日期\", \"steps\": \"步数\", ...}'", C.EXIT_USAGE)
        try:
            mapping = json.loads(a.map)
        except ValueError as e:
            raise C.LAError(f"--map is not JSON ({e})", C.EXIT_USAGE)
        res = wearable.summarize(st, ws.root, a.file, mapping)
        _insight_readouts(st, ws, "wearable", res.pop("readouts"))
    elif a.action == "mr":
        if not a.exposure or not a.outcome:
            raise C.LAError("mr needs --exposure and --outcome trait names (e.g. 'LDL cholesterol', 'Coronary heart disease')", C.EXIT_USAGE)
        from lalib import pubdata
        res = {"records": pubdata.mr(ws.root, a.exposure, a.outcome)}
    else:  # project
        if not (a.mr_ref and a.analyte and a.target is not None and a.baseline):
            raise C.LAError("project needs --mr-ref --analyte --target --baseline", C.EXIT_USAGE)
        res = causal.project(st, ws.root, a.mr_ref, a.analyte, a.target, a.baseline, exposure_match=a.exposure_match)
    written = {"genomics": "genotype_phenotype.json", "position": "positions.json", "project": "projections.json"}.get(a.action)
    if written:                                      # bound to its bytes: trace refuses an edited copy
        fp = ws.root / "work" / "insights" / written
        st.setdefault("insights", {}).setdefault("files", {})[written] = C.sha256_file(fp, limit=None)
    C.invalidate_after(st, "insights", f"insights {a.action}")
    _insights_stage(st)
    ws.log(st, f"insights_{a.action}")
    ws.save(st)
    C.emit(res)


def cmd_board(a):
    ws, st = _open(a)
    C.require(st, "organs")
    if a.action == "questions":
        if not a.file:
            raise C.LAError("board questions needs --file <questions.json>", C.EXIT_USAGE)
        res = board.register_questions(st, ws.root, Path(a.file))
    elif a.action == "finding":
        if not a.id:
            raise C.LAError("board finding needs --id Q<n>", C.EXIT_USAGE)
        try:
            res = board.register_finding(st, ws.root, a.id)
        except C.LAError:
            C.invalidate_after(st, "insights", f"finding {a.id} rejected")
            _insights_stage(st)
            ws.save(st)
            raise
    else:
        if not a.id or not (a.reason or "").strip():
            raise C.LAError("board skip needs --id and --reason", C.EXIT_USAGE)
        res = board.skip(st, a.id, a.reason)
    C.invalidate_after(st, "insights", f"board {a.action}")
    _insights_stage(st)
    ws.log(st, f"board_{a.action}", id=a.id)
    ws.save(st)
    C.emit(res)

def cmd_intervene(a):
    ws, st = _open(a)
    C.require(st, "integrate", "organs", "insights")
    try:
        res = evidence.check_plan(st, ws.root, Path(a.plan))
    except C.LAError:
        if st.get("intervene"):                  # a rejected new plan must not leave the old one in force
            C.invalidate_after(st, "organs", "plan registration rejected")
            ws.save(st)
        raise
    C.invalidate_after(st, "intervene", "plan registered")
    st["stages"]["intervene"] = "done"
    ws.log(st, "intervene_register", items=res["items"])
    ws.save(st)
    C.emit(res)


def cmd_twin(a):
    if a.action == "compare":
        if not a.prev or not a.cur:
            raise C.LAError("twin compare needs --prev and --cur", C.EXIT_USAGE)
        C.emit(twin.compare(Path(a.prev), Path(a.cur)))
        return
    ws, st = _open(a)
    C.require(st, "intervene")
    res = twin.build(st, ws.root)
    C.invalidate_after(st, "twin", "twin rebuilt")
    st["stages"]["twin"] = "done"
    ws.log(st, "twin_build")
    ws.save(st)
    C.emit(res)


def cmd_review(a):
    ws, st = _open(a)
    C.require(st, "twin")
    if a.action == "trace":
        res = report.trace(st, ws.root)
        st["stages"]["review"] = "todo"
        st["stages"]["report"] = "todo"
        ws.save(st)
        C.emit(res)
        return C.EXIT_OK if res["ok"] else C.EXIT_BLOCKED
    if not a.verdict or not a.findings:
        raise C.LAError("review record needs --verdict and --findings", C.EXIT_USAGE)
    res = report.record_review(st, ws.root, a.verdict, Path(a.findings))
    st["stages"]["review"] = "done" if a.verdict == "pass" else "todo"
    ws.log(st, "review_record", verdict=a.verdict)
    ws.save(st)
    C.emit(res)


def cmd_report(a):
    ws, st = _open(a)
    C.require(st, "intake", "preflight", "pipelines", "methods", "integrate", "organs", "insights", "intervene", "twin", "review")
    res = report.render(st, ws.root)
    st["stages"]["report"] = "done"
    ws.log(st, "report")
    ws.save(st)
    C.emit(res)


def cmd_mirobody(a):
    from lalib import mirobody
    res = mirobody.pull(mirobody.url_from(a.mcp_url_file), Path(a.folder).expanduser(), days=a.days)
    C.emit(res)


def cmd_export(a):
    ws, st = _open(a, lock=False)
    C.require(st, "report")
    p = ws.root / "deliver" / "la-export.json"
    want = ((st.get("report") or {}).get("deliver_sha256") or {}).get("la-export.json")
    if not p.exists() or C.sha256_file(p, limit=None) != want:
        raise C.LAError("deliver/la-export.json is missing or changed since `la.py report`; run la.py report again", C.EXIT_BLOCKED)
    d = C.load_json(p)
    C.emit({"path": str(p), "schema": d["schema"], "report": d["report"]["html"], "readouts": len(d["readouts"]),
            "organs": len(d["organs"]), "board": len(d["board"]), "plan_items": len(d["plan"]["items"]), "retests": len(d["retests"])})


NEXT = {
    "intake": "resolve pending judgments (la.py assign / la.py labs add), then confirm lab rows (la.py labs candidates / labs confirm), see workflows/01-intake.md",
    "preflight": "la.py preflight <ws>; show the result to the user (workflows/02-preflight-pipelines.md)",
    "pipelines": "plan, get the user's consent, run and verify, or skip with the reason (workflows/02-preflight-pipelines.md)",
    "methods": "la.py methods plan / run (workflows/03-methods.md)",
    "integrate": "bundle, dispatch one analyst per system, register (workflows/04-integrate.md)",
    "organs": "organ bundle, one estimator per organ, organ register (workflows/04b-organs.md)",
    "insights": "insights genomics / position / wearable, then the question board with one researcher per question (workflows/04c-insights.md)",
    "intervene": "evidence lookups, write plan.json, register (workflows/05-intervene.md)",
    "twin": "la.py twin build (workflows/06-twin-report.md)",
    "review": "summary.md, trace, independent reviewer, record (workflows/06-twin-report.md)",
    "report": "la.py report (workflows/06-twin-report.md)",
}


def cmd_observe(a):
    ws, st = _open(a, lock=False)
    files = {}
    for f in st["files"]:
        files.setdefault(f.get("status"), []).append(f"{f['id']} {f['name']} [{f.get('kind')}]")
    nxt = next((s for s in C.STAGES if st["stages"][s] != "done"), None)
    C.emit({
        "member": st["member"],
        "stages": st["stages"],
        "identity": st.get("identity"),
        "files_by_status": files,
        "labs": len(st["labs"]),
        "processed": {k: len(v) for k, v in st["processed"].items()},
        "pending_judgments": C.pending_judgments(st),
        "pending_lab_confirmations": len(__import__("lalib.labnames", fromlist=["pending"]).pending(st)),
        "pipelines": [{k: p.get(k) for k in ("run_id", "pipeline", "tier", "status", "optional")} for p in st["pipelines"]],
        "methods": (st.get("methods") or {}).get("counts"),
        "review": {k: (v.get("ok") if k == "trace" else v.get("verdict")) for k, v in (st.get("review") or {}).items()},
        "next": f"{nxt}: {NEXT[nxt]}" if nxt else "all stages done; deliver/report.html is ready",
    })


def cmd_validate(a):
    ws, st = _open(a, lock=False)
    problems = []
    ids = set()
    excluded = {f["id"] for f in st["files"] if f.get("excluded")}
    for f in st["files"]:
        if f["id"] in ids:
            problems.append(f"duplicate file id {f['id']}")
        ids.add(f["id"])
        p = Path(f["path"])
        if not p.exists():
            problems.append(f"{f['id']} missing on disk: {p}")
        elif p.stat().st_size != f.get("size", -1):
            problems.append(f"{f['id']} changed size since intake (start a new workspace)")
    for mod, lst in st["processed"].items():
        for x in lst:
            if not Path(x["path"]).exists():
                problems.append(f"processed {mod} points at missing {x['path']}")
            if set(x.get("source_file_ids") or []) & excluded:
                problems.append(f"processed {mod} still derives from an excluded file")
    for j in st["judgments"]:
        if j["status"] == "resolved" and not j.get("reason"):
            problems.append(f"judgment {j['id']} resolved without a reason")
    rp = ws.p("work", "readouts.json")
    if rp.exists():
        rid = [r["id"] for r in C.load_json(rp)["readouts"]]
        if len(rid) != len(set(rid)):
            problems.append("duplicate readout ids")
    if st["stages"]["report"] == "done":
        problems += [f"report: {p}" for p in report.check_bound(st, ws.root)]
        if (st.get("report") or {}).get("hashes") != report._hashes(ws.root, st):
            problems.append("report: inputs changed after the report was rendered")
        for name, sha in ((st.get("report") or {}).get("deliver_sha256") or {}).items():
            p = ws.p("deliver", name)
            if not p.exists() or C.sha256_file(p, limit=None) != sha:
                problems.append(f"deliver/{name} was changed after rendering")
    elif ws.p("deliver").exists():
        problems.append("deliver/ exists but the report stage is not done (stale deliverable)")
    for lst in st["processed"].values():
        for x in lst:
            if x.get("sha256") and Path(x["path"]).exists() and C.sha256_file(Path(x["path"]), limit=None) != x["sha256"]:
                problems.append(f"derived file changed after intake: {x['path']}")
    for f in st["files"]:
        p = Path(f["path"])
        if not f.get("excluded") and f.get("sha256") not in (None, "unreadable") and p.exists() and C.sha256_file(p) != f["sha256"]:
            problems.append(f"{f['id']} content changed since intake")
    C.emit({"ok": not problems, "problems": problems})
    return C.EXIT_OK if not problems else C.EXIT_INPUT


def build_parser():
    p = argparse.ArgumentParser(prog="la.py", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="create a member workspace")
    s.add_argument("raw_dir")
    s.add_argument("workspace")
    s.add_argument("--member-id", required=True)
    s.add_argument("--age", type=int)
    s.add_argument("--sex")
    s.add_argument("--sample-date")
    s.add_argument("--mode", choices=["commercial", "research"], default="commercial")
    s.add_argument("--longevity-skills", help="path to the longevity-skills clone (stored in the workspace)")
    s.add_argument("--force", action="store_true", help="move an existing workspace aside and start fresh")
    s.set_defaults(fn=cmd_init)

    s = sub.add_parser("member", help="record facts the user told you (age, sex, questionnaire answers)")
    s.add_argument("workspace")
    s.add_argument("pairs", nargs="+", help="key=value")
    s.add_argument("--source", required=True, help="where the value came from (quote the user)")
    s.set_defaults(fn=cmd_member_set)

    for name, fn in (("intake", cmd_intake), ("observe", cmd_observe), ("validate", cmd_validate), ("report", cmd_report)):
        s = sub.add_parser(name)
        s.add_argument("workspace")
        s.set_defaults(fn=fn)

    s = sub.add_parser("assign", help="record a judgment")
    s.add_argument("workspace")
    s.add_argument("--file", default="member", help="file id; `member` for identity; a modality name for primary")
    s.add_argument("--what", required=True, choices=["modality", "platform", "value_scale", "tissue", "assembly",
                                                     "sample_column", "transcribe", "identity", "primary"])
    s.add_argument("--value", required=True)
    s.add_argument("--reason", required=True)
    s.add_argument("--files", help="with --what identity: comma-separated ids of files that are not this member's (excluded)")
    s.add_argument("--uncertain", help="with --what identity: ids of files whose ownership cannot be confirmed (analysed, labelled)")
    s.set_defaults(fn=cmd_assign)

    s = sub.add_parser("labs", help="add transcribed rows, map a name to a method input, fix or drop a row")
    s.add_argument("action", choices=["add", "map", "fix", "drop", "candidates", "confirm"])
    s.add_argument("workspace")
    s.add_argument("--answers", help="for confirm: JSON {answers: [{row_key, answer: yes|no, why}]}")
    s.add_argument("--csv")
    s.add_argument("--file")
    s.add_argument("--marker")
    s.add_argument("--to")
    s.add_argument("--value")
    s.add_argument("--unit")
    s.add_argument("--row", type=int, help="for fix/drop: which of several rows with this name (1-based, in file order)")
    s.add_argument("--reason")
    s.set_defaults(fn=cmd_labs)

    s = sub.add_parser("preflight")
    s.add_argument("workspace")
    s.add_argument("--no-network", action="store_true")
    s.set_defaults(fn=cmd_preflight)

    s = sub.add_parser("pipeline")
    s.add_argument("action", choices=["plan", "run", "status", "stop", "verify", "skip"])
    s.add_argument("workspace")
    s.add_argument("--run-id")
    s.add_argument("--confirmed", help="the user's own words agreeing to run (quoted)")
    s.add_argument("--accept-yellow", action="store_true")
    s.add_argument("--reopen", action="store_true", help="run a previously skipped pipeline after the user changed their mind")
    s.add_argument("--stub", action="store_true", help="developer wiring check (needs LA_DEV=1); never a member result")
    s.add_argument("--reason", help="for skip: why (quote the user if they declined)")
    s.add_argument("--pgs-id", help="for plan: PGS Catalog ids the user confirmed, e.g. PGS000018,PGS000337")
    s.set_defaults(fn=cmd_pipeline)

    s = sub.add_parser("methods")
    s.add_argument("action", choices=["plan", "run"])
    s.add_argument("workspace")
    s.set_defaults(fn=cmd_methods)

    s = sub.add_parser("integrate")
    s.add_argument("action", choices=["bundle", "register"])
    s.add_argument("workspace")
    s.add_argument("--system")
    s.set_defaults(fn=cmd_integrate)

    s = sub.add_parser("evidence")
    s.add_argument("action", choices=["local", "pubmed", "fetch"])
    s.add_argument("workspace")
    s.add_argument("--terms", default="")
    s.add_argument("--query", default="")
    s.add_argument("--url", default="")
    s.add_argument("--max", type=int, default=10)
    s.set_defaults(fn=cmd_evidence)

    s = sub.add_parser("organ", help="organ checkup table: bundles, AI estimates, skips")
    s.add_argument("action", choices=["bundle", "register", "skip"])
    s.add_argument("workspace")
    s.add_argument("--organ")
    s.add_argument("--reason")
    s.set_defaults(fn=cmd_organ)

    s = sub.add_parser("insights", help="genotype-phenotype, population position, wearables, MR projections")
    s.add_argument("action", choices=["genomics", "skip-genomics", "position", "wearable", "mr", "project"])
    s.add_argument("workspace")
    s.add_argument("--reason")
    s.add_argument("--file")
    s.add_argument("--map", help="wearable: JSON mapping date/steps/rhr/hrv/sleep_h/deep_h/mvpa_min/spo2_min -> column names")
    s.add_argument("--exposure")
    s.add_argument("--outcome")
    s.add_argument("--mr-ref")
    s.add_argument("--analyte")
    s.add_argument("--target", type=float)
    s.add_argument("--baseline")
    s.add_argument("--exposure-match", help="project: why an MR exposure with another name measures the same quantity as --analyte")
    s.add_argument("--absent-as-ref", help="genomics: judgment that sites absent from a variant-only VCF are reference, with the reason")
    s.set_defaults(fn=cmd_insights)

    s = sub.add_parser("board", help="question board: register questions, researcher findings, skips")
    s.add_argument("action", choices=["questions", "finding", "skip"])
    s.add_argument("workspace")
    s.add_argument("--file")
    s.add_argument("--id")
    s.add_argument("--reason")
    s.set_defaults(fn=cmd_board)

    s = sub.add_parser("mirobody", help="before init: pull the member's confirmed labs and daily wearable values from Mirobody into the data folder")
    s.add_argument("action", choices=["pull"])
    s.add_argument("folder", help="the member's data folder (the one you will pass to init)")
    s.add_argument("--mcp-url-file", help="file holding the member's personal MCP URL (or set LONGPI_MCP_URL)")
    s.add_argument("--days", type=int, default=120)
    s.set_defaults(fn=cmd_mirobody)

    s = sub.add_parser("export", help="the la-export/1 file a host app (longpi) imports; written by `report`")
    s.add_argument("workspace")
    s.set_defaults(fn=cmd_export)

    s = sub.add_parser("intervene")
    s.add_argument("action", choices=["register"])
    s.add_argument("workspace")
    s.add_argument("--plan", required=True)
    s.set_defaults(fn=cmd_intervene)

    s = sub.add_parser("twin")
    s.add_argument("action", choices=["build", "compare"])
    s.add_argument("workspace", nargs="?")
    s.add_argument("--prev")
    s.add_argument("--cur")
    s.set_defaults(fn=cmd_twin)

    s = sub.add_parser("review")
    s.add_argument("action", choices=["trace", "record"])
    s.add_argument("workspace")
    s.add_argument("--verdict", choices=["pass", "block"])
    s.add_argument("--findings")
    s.set_defaults(fn=cmd_review)
    return p


def _adopt_workspace_env(a) -> None:
    """Shell state does not persist between an agent's tool calls; the workspace remembers the library path."""
    wsp = getattr(a, "workspace", None)
    if wsp and not os.environ.get("LONGEVITY_SKILLS_HOME"):
        sp = Path(wsp) / "state.json"
        if sp.exists():
            lib = C.load_json(sp).get("longevity_skills_home")
            if lib:
                os.environ["LONGEVITY_SKILLS_HOME"] = lib


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        _adopt_workspace_env(a)
        rc = a.fn(a)
        return rc or C.EXIT_OK
    except C.LAError as e:
        print(f"error: {e}", file=sys.stderr)
        return e.code
    finally:
        for w in C.OPEN_WORKSPACES:
            w.unlock()
        C.OPEN_WORKSPACES.clear()


if __name__ == "__main__":
    sys.exit(main())
