"""Plan, launch, watch and verify local pipeline runs.

The command line and parameter file come from data/pipelines.json templates
(pinned revision, fixed params). The agent chooses nothing here except which
planned run to start, and only after the user confirmed the preflight.
"""
from __future__ import annotations

import csv
import glob
import hashlib
import json
import os
import platform
import re
import unicodedata
import signal
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import (EXIT_BLOCKED, EXIT_INPUT, LAError, add_judgment, data, file_by_id,
                     load_json, now_iso, write_json)


def _run_dir(ws_root: Path, run_id: str) -> Path:
    return ws_root / "work" / "pipelines" / run_id


def pipeline_cache(repo: str, revision: str) -> Path:
    base = Path(os.environ.get("LA_PIPELINE_CACHE", Path.home() / ".cache" / "longevity-analyst" / "pipelines"))
    return base / repo.replace("/", "__") / revision


def fetch(repo: str, revision: str) -> Dict[str, Any]:
    """Clone the pinned revision with the system git (honours the user's git/ssh config, unlike Nextflow's JGit)."""
    dst = pipeline_cache(repo, revision)
    if (dst / "main.nf").exists():
        commit = subprocess.run(["git", "-C", str(dst), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        return {"path": str(dst), "commit": commit, "cached": True}
    dst.parent.mkdir(parents=True, exist_ok=True)
    p = subprocess.run(["git", "clone", "--quiet", "--depth", "1", "--branch", revision,
                        f"https://github.com/{repo}.git", str(dst)], capture_output=True, text=True, timeout=600)
    if p.returncode != 0:
        raise LAError(f"could not fetch {repo}@{revision}: {p.stderr.strip()[:400]}", EXIT_BLOCKED)
    commit = subprocess.run(["git", "-C", str(dst), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {"path": str(dst), "commit": commit, "cached": False}


def _samplesheet(pid: str, st: Dict[str, Any], file_ids: List[str], out: Path) -> Optional[Path]:
    files = [file_by_id(st, i) if not i.startswith("path:") else {"id": i, "path": i[5:], "kind": "variants_vcf", "name": Path(i[5:]).name} for i in file_ids]
    member = st["member"]["id"]
    sheet = out / "samplesheet.csv"
    if pid == "sarek":
        r1 = [f for f in files if f["kind"] == "reads_fastq" and f.get("mate") == "1"]
        rows = []
        for i, f in enumerate(sorted(r1, key=lambda x: x["name"])):
            mate = next((g for g in files if g.get("pair_key") == f.get("pair_key") and g.get("mate") == "2"), None)
            rows.append({"patient": member, "sample": member, "lane": f"L{i + 1}", "fastq_1": f["path"], "fastq_2": mate["path"] if mate else ""})
        crams = [f for f in files if f.get("kind") == "alignment"]
        if crams and not rows:
            with open(sheet, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["patient", "sample", "cram", "crai"])
                w.writeheader()
                for c in crams:
                    w.writerow({"patient": member, "sample": member, "cram": c["path"], "crai": c["path"] + ".crai"})
            return sheet
        cols = ["patient", "sample", "lane", "fastq_1", "fastq_2"]
    elif pid == "methylseq":
        r1 = [f for f in files if f.get("mate") in ("1", None)]
        rows = []
        for f in r1:
            mate = next((g for g in files if g.get("pair_key") == f.get("pair_key") and g.get("mate") == "2"), None)
            rows.append({"sample": member, "fastq_1": f["path"], "fastq_2": mate["path"] if mate else "", "genome": ""})
        cols = ["sample", "fastq_1", "fastq_2", "genome"]
    elif pid == "taxprofiler":
        r1 = [f for f in files if f.get("mate") in ("1", None)]
        rows = []
        for i, f in enumerate(r1):
            mate = next((g for g in files if g.get("pair_key") == f.get("pair_key") and g.get("mate") == "2"), None)
            rows.append({"sample": member, "run_accession": f"run{i + 1}", "instrument_platform": "ILLUMINA",
                         "fastq_1": f["path"], "fastq_2": mate["path"] if mate else "", "fasta": ""})
        cols = ["sample", "run_accession", "instrument_platform", "fastq_1", "fastq_2", "fasta"]
        db = os.environ.get("LA_REF_METAPHLAN_DB", "")
        with open(out / "databases.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["tool", "db_name", "db_params", "db_type", "db_path"])
            w.writerow(["metaphlan", "mpa", "", "short", db])
    elif pid == "pgsc_calc":
        v = files[0]
        rows = [{"sampleset": re.sub(r"[^A-Za-z0-9]", "", member) or "member",
                 "path_prefix": re.sub(r"\.g?\.?vcf\.gz$", "", v["path"]), "chrom": "", "format": "vcf"}]
        cols = ["sampleset", "path_prefix", "chrom", "format"]
    else:
        return None
    with open(sheet, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        w.writerows(rows)
    return sheet


def _paths_for(st: Dict[str, Any], ids: List[str]) -> List[str]:
    out = []
    for i in ids:
        if i.startswith("path:"):
            out.append(i[5:])
        else:
            out.append(file_by_id(st, i)["path"])
    return out


def plan(st: Dict[str, Any], ws_root: Path, pgs_id: Optional[str] = None) -> List[Dict[str, Any]]:
    pf = st.get("preflight")
    if not pf:
        raise LAError("run `la.py preflight` first", EXIT_INPUT)
    cfg = data("pipelines.json")["pipelines"]
    mode = st["member"].get("mode", "commercial")
    plans = []
    for r in pf["pipelines"]:
        pid = r["pipeline"]
        p = cfg[pid]
        run_id = f"{pid}-{hashlib.sha1(','.join(r['file_ids']).encode()).hexdigest()[:8]}"
        out = _run_dir(ws_root, run_id)
        out.mkdir(parents=True, exist_ok=True)
        params = dict(p["params"])
        params["outdir"] = str(out / "results")
        if pid == "sarek" and os.environ.get("LA_REF_GRCH38"):
            params["igenomes_base"] = os.environ["LA_REF_GRCH38"]
        if pid == "methylseq" and os.environ.get("LA_REF_BISMARK"):
            params["bismark_index"] = os.environ["LA_REF_BISMARK"]
        if pid == "sarek" and any(not i.startswith("path:") and file_by_id(st, i).get("kind") == "alignment" for i in r["file_ids"]):
            params["step"] = "variant_calling"          # start from the delivered CRAM/BAM
        if pid == "methylseq" and any(not i.startswith("path:") and file_by_id(st, i).get("modality") == "em_seq" for i in r["file_ids"]):
            params["em_seq"] = True
        if pid == "taxprofiler":
            params["perform_runmerging"] = True
        if pid == "pgsc_calc":
            if pgs_id and not re.fullmatch(r"PGS\d{6}(,PGS\d{6})*", pgs_id):
                raise LAError("--pgs-id must look like PGS000018,PGS000337", EXIT_INPUT)
            if pgs_id:
                params["pgs_id"] = pgs_id
            if os.environ.get("LA_REF_PGSC_PANEL"):
                params["run_ancestry"] = os.environ["LA_REF_PGSC_PANEL"]
        sheet = _samplesheet(pid, st, r["file_ids"], out)
        if sheet:
            params["input"] = str(sheet)
        if pid == "taxprofiler":
            params["databases"] = str(out / "databases.csv")
        write_json(out / "params.json", params)
        if p["engine"] == "nextflow":
            profile = "docker,arm" if platform.machine() in ("arm64", "aarch64") and p.get("has_arm_profile") else "docker"
            cmd = ["nextflow", "run", p["repo"], "-r", p["revision"], "-profile", profile,
                   "-params-file", str(out / "params.json"), "-work-dir", str(out / "work"),
                   "-with-trace", str(out / "trace.txt"), "-resume"]
        elif p["engine"] == "java":
            jar = os.environ.get("LA_PHARMCAT_JAR", "pharmcat-3.4.0-all.jar")
            vcf = _paths_for(st, r["file_ids"])[0]
            cmd = ["java", "-jar", jar, "-vcf", vcf, "-o", str(out / "results"), "-reporterJson"]
        else:
            idats = sorted({file_by_id(st, i)["path"].rsplit("_", 1)[0] for i in r["file_ids"]})
            cmd = ["Rscript", str(Path(__file__).resolve().parents[1] / "idat_to_beta.R"), str(out / "results" / "betas.csv"), *idats]
        blocked = None
        if mode == "commercial" and not r.get("commercial_ok", True):
            blocked = "license: not cleared for commercial use"
        if pid == "pgsc_calc" and not params.get("pgs_id"):
            blocked = "needs the user-confirmed PGS Catalog ids: `pipeline plan --pgs-id PGS000018,...`"
        if sheet and sum(1 for _ in open(sheet)) < 2:
            blocked = "no usable read pairs for this pipeline (sarek needs paired-end FASTQ or a CRAM)"
        if pid == "taxprofiler" and not os.environ.get("LA_REF_METAPHLAN_DB"):
            blocked = "set LA_REF_METAPHLAN_DB to a downloaded MetaPhlAn database first (taxprofiler cannot run without it)"
        sources = set()
        for i in r["file_ids"]:
            if i.startswith("path:"):
                for lst in st["processed"].values():
                    for x in lst:
                        if x["path"] == i[5:]:
                            sources |= set(x.get("source_file_ids") or [])
            else:
                sources.add(i)
        tissues = {file_by_id(st, i).get("tissue") for i in r["file_ids"] if not i.startswith("path:")} - {None}
        plans.append({"run_id": run_id, "pipeline": pid, "revision": p["revision"], "tier": r["tier"],
                      "tissue": tissues.pop() if len(tissues) == 1 else None, "source_file_ids": sorted(sources),
                      "optional": r["optional"], "file_ids": r["file_ids"], "cmd": cmd, "dir": str(out),
                      "estimate": r["estimate"], "blocked": blocked, "status": "planned"})
    excluded = {f["id"] for f in st["files"] if f.get("excluded")}
    current = {pl["run_id"] for pl in plans}
    for x in st["pipelines"]:
        if x["run_id"] not in current and x["status"] == "planned":
            x.update(status="void_stale", blocked="no longer in the current plan")
        if set(x.get("source_file_ids") or x.get("file_ids") or []) & excluded:
            x.update(status="void_excluded", blocked="input files were excluded as not this member's")
    existing = {x["run_id"]: x for x in st["pipelines"]}
    for pl in plans:
        old = existing.get(pl["run_id"])
        if old and old["status"] not in ("planned",):
            continue                     # keep history of runs that started, were skipped or verified
        existing[pl["run_id"]] = pl
    st["pipelines"] = list(existing.values())
    return plans


def launch(st: Dict[str, Any], run_id: str, confirmed_by_user: str, stub: bool = False, force_yellow: bool = False,
           reopen: bool = False) -> Dict[str, Any]:
    run = next((r for r in st["pipelines"] if r["run_id"] == run_id), None)
    if not run:
        raise LAError(f"no planned run {run_id}; run `la.py pipeline plan`", EXIT_INPUT)
    if run.get("blocked") or run["status"].startswith("void"):
        raise LAError(f"{run_id} blocked: {run.get('blocked') or run['status']}", EXIT_BLOCKED)
    excluded = {f["id"] for f in st["files"] if f.get("excluded")}
    if set(run.get("source_file_ids") or run.get("file_ids") or []) & excluded:
        raise LAError(f"{run_id} uses files excluded as not this member's", EXIT_BLOCKED)
    if run["status"] == "skipped" and not reopen:
        raise LAError(f"{run_id} was skipped ({run.get('skip_reason')}); pass --reopen with the user's new words to run it", EXIT_BLOCKED)
    words = re.sub(r"[\u200b-\u200f\u2060-\u2064\ufeff\u00ad\s]+", " ", unicodedata.normalize("NFKC", confirmed_by_user or "")).strip()
    meaningful = re.findall(r"[A-Za-z\u4e00-\u9fff]", words)
    # A floor, not the judgment: reading the user's reply as consent is the agent's job. Words that are empty,
    # attributed to the agent, or carry a refusal or non-answer are never accepted as consent.
    refusal = (re.search(r"(?:不|别|勿|没|未|拒绝?)\s*(?:要|用|想|再|必|需要?)?\s*(?:同意|允许|确认|跑|运行|执行|开始|启动|做|回复|答复|授权)", words)
               or re.search(r"不(?:要|行|可以|准|许|好)(?![问担怕急])|算了|取消|以后再说|改天|等等|等一下|等下|等我|再想想|考虑|稍后|不了|不用|不必|先不|再说吧|^再说|明天|不太想|还是别|没兴趣|不确定|随便|你决定|你看着办|[?？]",
                         re.sub(r"不用再问|不用担心|不用确认了?|不必再确认|不好意思", "", words))
               or re.search(r"(?i)\b(?:no|nope|nah|not|stop|cancel|refuse|decline|never|later|wait)\b|\bn/?a\b|\bnone\b|"
                            r"\b(?:don'?t|do not)\s+(?:run|go|start|proceed|do|launch)|hold on|tomorrow|please don'?t|not sure|whatever|up to you", words)
               or re.fullmatch(r"(?i)[\W_]*(?:n|否|不|无|没有|拒绝)[\W_]*", words))
    single_yes = re.fullmatch(r"(?i)[\W_]*(?:好|行|嗯|可|是|对|跑|y)[\W_]*", words)
    affirm = re.search(r"(?i)好|行|可以|同意|跑吧|跑一下|帮我跑|开始|确认|没问题|嗯|是的|^对|(?<![a-z])(?:ok|okay|yes|yep|sure|go|proceed|run it)(?![a-z])|^y$", words)
    # questions, conditions and postponements go back to the agent to clarify; they are never consent
    unclear = re.search(r"吗|呢|什么|怎么|为什么|解释|告诉我|但是?|只跑|除了?|别给|下周|下次|下个月|过几天|回头|暂停|不需要|不是现在|晚点|等会|先放|没说|谁让|"
                        r"(?i:hold off|not now|explain|what is|why)", words)
    if ((len(meaningful) < 2 and not single_yes) or not affirm or unclear or re.fullmatch(r"[xX.\-_ ]+", words) or refusal
            or re.search(r"(?i)^[\W_]*(agent|assistant|claude|system|ai|助手|代理)", words) or re.search(r"(?i)(ai|agent|助手|代理)\s*(确认|confirm)", words)):
        raise LAError("--confirmed needs the user's own words of consent (quote them, e.g. \"好的，跑吧\"); "
                      "show the preflight first", EXIT_BLOCKED)
    if stub and os.environ.get("LA_DEV") != "1":
        raise LAError("--stub is a developer wiring check (set LA_DEV=1); it is never part of a member analysis", EXIT_BLOCKED)
    live = next((p for p in (st.get("preflight") or {}).get("pipelines", [])
                 if p["pipeline"] == run["pipeline"] and p["file_ids"] == run["file_ids"]), None)
    if not live:
        raise LAError(f"{run_id} is not in the latest preflight; re-run preflight and plan", EXIT_BLOCKED)
    run["tier"] = live["tier"]                    # the latest preflight decides, not the tier at planning time
    if run["tier"] == "red" and not stub:
        raise LAError(f"{run_id} preflight is red; fix the listed problems and re-run preflight", EXIT_BLOCKED)
    if run["tier"] == "yellow" and not (force_yellow or stub):
        raise LAError(f"{run_id} preflight is yellow; pass --accept-yellow after the user accepted the listed costs", EXIT_BLOCKED)
    if run_alive(run):
        raise LAError(f"{run_id} still has a live process (pid {run['pid']}, status {run['status']}); `pipeline status` / `pipeline stop` first", EXIT_BLOCKED)
    cmd = list(run["cmd"])
    if cmd[0] == "nextflow":
        cfg = data("pipelines.json")["pipelines"][run["pipeline"]]
        got = fetch(cfg["repo"], cfg["revision"])
        i = cmd.index(cfg["repo"])
        cmd[i] = got["path"]
        del cmd[i + 1:i + 3]  # drop "-r <revision>": the local clone is already at that tag
        run["pipeline_commit"] = got["commit"]
    if stub:
        if cmd[0] != "nextflow":
            raise LAError("--stub only applies to nextflow pipelines", EXIT_INPUT)
        i = cmd.index("-profile")
        cmd[i + 1] = "test," + cmd[i + 1]   # stubs still call `tool --version` inside the containers
        cmd.append("-stub-run")
        # The stub run checks our samplesheet and command wiring; reference-database params come from
        # the pipeline's own test profile because no real database is downloaded for a stub.
        params = load_json(Path(run["dir"]) / "params.json")
        params["outdir"] = str(Path(run["dir"]) / "results_stub")     # stub outputs never land where verify looks
        params.pop("input", None)                                      # the test profile's own inputs, never member data
        for k in ("genome", "igenomes_base"):
            params.pop(k, None)
        if "databases" in params:
            # stub processes never open the database; a local empty directory satisfies schema validation offline
            dummy = Path(run["dir"]) / "stub_db"
            dummy.mkdir(exist_ok=True)
            with open(Path(run["dir"]) / "databases.stub.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["tool", "db_name", "db_params", "db_type", "db_path"])
                w.writerow(["metaphlan", "stub", "", "short", str(dummy)])
            params["databases"] = str(Path(run["dir"]) / "databases.stub.csv")
        write_json(Path(run["dir"]) / "params.stub.json", params)
        cmd[cmd.index("-params-file") + 1] = str(Path(run["dir"]) / "params.stub.json")
    out = Path(run["dir"])
    attempt = int(run.get("attempt", 0)) + 1
    for name in ("trace.txt", "run.log"):   # keep earlier attempts, start each launch clean
        if (out / name).exists():
            (out / name).rename(out / f"{name}.attempt{attempt - 1}")
    run["attempt"] = attempt
    log = open(out / "run.log", "wb")
    (out / "exitcode").unlink(missing_ok=True)
    wrapped = ["/bin/sh", "-c", '"$@"; echo $? > exitcode', "la-run", *cmd]
    proc = subprocess.Popen(wrapped, cwd=out, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    run.update(status="running", pid=proc.pid, pid_started=_start_time(proc.pid), started_at=now_iso(), stub=stub,
               confirmed_by_user=confirmed_by_user, launched_cmd=cmd)
    run.setdefault("attempts", []).append({"attempt": attempt, "stub": stub, "started_at": run["started_at"],
                                            "confirmed_by_user": confirmed_by_user})
    return run


def _start_time(pid: int) -> Optional[str]:
    """Process start time in a fixed locale and time zone, so the value is comparable across shells."""
    try:
        env = {**os.environ, "LC_ALL": "C", "LANG": "C", "TZ": "UTC"}
        out = subprocess.run(["ps", "-o", "lstart=", "-p", str(pid)], capture_output=True, text=True, timeout=5, env=env).stdout.strip()
        return out or None
    except (OSError, subprocess.TimeoutExpired):
        return None


def _alive(pid: int, started: Optional[str] = None) -> bool:
    """Alive and, when we know it, still the process we launched (pids are reused)."""
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    if started is not None and _start_time(pid) != started:
        return False
    return True


def run_alive(run: Dict[str, Any]) -> bool:
    return bool(run.get("pid")) and _alive(run["pid"], run.get("pid_started"))


def status(st: Dict[str, Any], run_id: str) -> Dict[str, Any]:
    run = next((r for r in st["pipelines"] if r["run_id"] == run_id), None)
    if not run:
        raise LAError(f"no run {run_id}", EXIT_INPUT)
    out = Path(run["dir"])
    counts: Dict[str, int] = {}
    trace = out / "trace.txt"
    if trace.exists():
        with open(trace) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                counts[row.get("status", "?")] = counts.get(row.get("status", "?"), 0) + 1
    tail = ""
    if (out / "run.log").exists():
        tail = "\n".join((out / "run.log").read_text(errors="replace").splitlines()[-15:])
    hint = None
    if run["status"] == "running" and run.get("pid") and not run_alive(run):
        run["status"] = "exited"
        run["ended_at"] = now_iso()
    elif run["status"] == "running" and ("Pipeline failed" in tail or "Execution cancelled" in tail):
        # Nextflow can linger after a failure while it waits for container tasks; say so instead of "running".
        run["status"] = "failed_hanging"
        hint = "the pipeline failed but the Nextflow process is still alive; run `pipeline stop`, read run.log, fix, then relaunch (-resume reuses finished tasks)"
    return {"run_id": run_id, "status": run["status"], "tasks": counts, "log_tail": tail, "hint": hint}


def stop(st: Dict[str, Any], run_id: str) -> Dict[str, Any]:
    import time
    run = next((r for r in st["pipelines"] if r["run_id"] == run_id), None)
    if not run:
        raise LAError(f"no run {run_id}", EXIT_INPUT)
    if run_alive(run):
        try:
            os.killpg(os.getpgid(run["pid"]), signal.SIGTERM)
            for _ in range(20):
                if not run_alive(run):
                    break
                time.sleep(1)
            if run_alive(run):
                os.killpg(os.getpgid(run["pid"]), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
    run["status"] = "stopped"
    run["stopped_at"] = now_iso()
    return run


def verify(st: Dict[str, Any], run_id: str) -> Dict[str, Any]:
    run = next((r for r in st["pipelines"] if r["run_id"] == run_id), None)
    if not run:
        raise LAError(f"no run {run_id}", EXIT_INPUT)
    cfg = data("pipelines.json")["pipelines"][run["pipeline"]]
    res_dir = Path(run["dir"]) / "results"
    found: Dict[str, List[str]] = {}
    import datetime as _dt
    started = _dt.datetime.fromisoformat(run["started_at"]).timestamp() if run.get("started_at") else 0
    for modality, pattern in cfg["outputs"].items():
        found[modality] = sorted(p for p in glob.glob(str(res_dir / pattern)) if os.path.getmtime(p) >= started - 1)
    missing = [m for m, v in found.items() if not v]
    st_trace = status(st, run_id)["tasks"]
    failed = st_trace.get("FAILED", 0)
    if run.get("status") == "running":
        raise LAError(f"{run_id} is still running; verify after it exits", EXIT_BLOCKED)
    excluded = {f["id"] for f in st["files"] if f.get("excluded")}
    if run["status"].startswith("void") or set(run.get("source_file_ids") or run.get("file_ids") or []) & excluded:
        raise LAError(f"{run_id} used files excluded as not this member's; its outputs are never registered", EXIT_BLOCKED)
    if not run.get("attempts"):
        raise LAError(f"{run_id} was never launched; nothing to verify", EXIT_BLOCKED)
    log = (Path(run["dir"]) / "run.log").read_text(errors="replace") if (Path(run["dir"]) / "run.log").exists() else ""
    ec_file = Path(run["dir"]) / "exitcode"
    exit_ok = ec_file.exists() and ec_file.read_text().strip() == "0"
    err = re.search(r"ERROR\s*~|^Error(?: in|:)|^Execution halted|Traceback \(most recent call last\)", log, re.M)
    succeeded = exit_ok and not err and ("Pipeline completed successfully" in log if cfg["engine"] == "nextflow" else True)
    empty = [p for v in found.values() for p in v if os.path.getsize(p) == 0 or not _gz_ok(p)]
    ok = not missing and failed == 0 and not run.get("stub") and succeeded and not empty
    prov = {"origin": "pipeline", "pipeline": cfg["repo"], "revision": run["revision"], "commit": run.get("pipeline_commit"), "run_id": run_id,
            "params_file": str(Path(run["dir"]) / "params.json"), "verified_at": now_iso(),
            "trace_tasks": st_trace, "stub": run.get("stub", False)}
    if ok:
        for modality, paths in found.items():
            key = {"variants": "variants", "cpg_report": "methylation_cov", "gut_profile": "gut_profile",
                   "prs": "prs", "pgx": "pgx", "methylation": "methylation"}[modality]
            for pth in paths:
                lst = st["processed"].setdefault(key, [])
                if not any(x["path"] == pth for x in lst):
                    unc = any(f.get("provenance_uncertain") for f in st["files"] if f["id"] in set(run.get("source_file_ids") or []))
                    lst.append({"file_id": None, "source_file_ids": list(run.get("source_file_ids") or run["file_ids"]), "path": pth,
                                **({"provenance_uncertain": True} if unc else {}),
                                "provenance": {**prov, **({"tissue": run.get("tissue")} if run.get("tissue") else {}),
                                               **({"assembly": "GRCh38", "gvcf": True} if key == "variants" else {})}})
                    if len(lst) > 1:
                        add_judgment(st, "primary", key, f"{len(lst)} {key} entries (delivered and pipeline output): which one is this visit's?",
                                     [x.get("file_id") or x["path"] for x in lst])
        run["status"] = "verified"
        for fid in run["file_ids"]:
            if fid.startswith("path:"):
                continue
            f = file_by_id(st, fid)
            if f.get("status") == "needs_pipeline":
                f["status"] = "processed"
    elif run.get("stub"):
        log = (Path(run["dir"]) / "run.log")
        text = log.read_text(errors="replace") if log.exists() else ""
        run["status"] = "stub_ok" if failed == 0 and st_trace and not re.search(r"ERROR\s*~", text) else "stub_failed"
    else:
        run["status"] = "verify_failed"
    return {"run_id": run_id, "ok": ok, "outputs": found, "missing": missing, "failed_tasks": failed,
            "pipeline_reported_success": succeeded, "empty_outputs": empty,
            "stub": run.get("stub", False), "status": run["status"]}


def _gz_ok(path: str) -> bool:
    if not path.endswith(".gz"):
        return True
    import gzip as _g
    try:
        with _g.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def skip(st: Dict[str, Any], run_id: str, reason: str) -> Dict[str, Any]:
    """Record that a planned run will not happen (user declined, red preflight, lab will deliver the file)."""
    if not reason.strip():
        raise LAError("--reason is required; quote the user if they declined", EXIT_INPUT)
    run = next((r for r in st["pipelines"] if r["run_id"] == run_id), None)
    if not run:
        raise LAError(f"no run {run_id}", EXIT_INPUT)
    if run_alive(run):
        raise LAError(f"{run_id} is running; stop it first", EXIT_BLOCKED)
    run.update(status="skipped", skip_reason=reason, skipped_at=now_iso())
    return run
