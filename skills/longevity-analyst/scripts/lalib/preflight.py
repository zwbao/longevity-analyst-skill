"""Can this machine run the planned pipelines, and how long will it take?

Every check actually executes the tool (`java -version`, `docker info`,
`nextflow -version`, `Rscript -e library(sesame)`) instead of trusting
`which`: macOS ships /usr/bin/java as a stub that prints "Unable to locate a
Java Runtime", and a docker client can exist without a running daemon.
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

from .common import data, now_iso

GREEN, YELLOW, RED = "green", "yellow", "red"


def _run(cmd: List[str], timeout: int = 25, env: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        return {"rc": p.returncode, "out": (p.stdout + p.stderr).strip()[:2000]}
    except FileNotFoundError:
        return {"rc": 127, "out": "not found"}
    except subprocess.TimeoutExpired:
        return {"rc": 124, "out": f"timed out after {timeout}s"}


def mem_gb() -> Optional[float]:
    try:
        if platform.system() == "Darwin":
            return int(_run(["sysctl", "-n", "hw.memsize"])["out"]) / 1024 ** 3
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal"):
                    return int(line.split()[1]) / 1024 ** 2
    except (ValueError, OSError):
        return None
    return None


def java_status() -> Dict[str, Any]:
    r = _run(["java", "-version"])
    ok = r["rc"] == 0 and "version" in r["out"]
    res: Dict[str, Any] = {"ok": ok, "detail": r["out"].splitlines()[0] if r["out"] else ""}
    if not ok:
        cands = sorted(Path("/opt/homebrew/opt").glob("openjdk*/bin/java")) + sorted(Path("/usr/local/opt").glob("openjdk*/bin/java"))
        if cands:
            home = cands[-1].parents[1] / "libexec" / "openjdk.jdk" / "Contents" / "Home"
            res["fix"] = f"a JDK is installed but not on PATH: export JAVA_HOME={home if home.exists() else cands[-1].parents[1]}"
        else:
            res["fix"] = "install a JDK 17+: `brew install openjdk@17` (macOS) or `apt install openjdk-17-jre` (Linux)"
    return res


def docker_status() -> Dict[str, Any]:
    if not shutil.which("docker"):
        return {"ok": False, "detail": "docker not installed", "fix": "install Docker Desktop / docker engine, or use conda/singularity"}
    r = _run(["docker", "info", "--format", "{{.ServerVersion}}|{{.NCPU}}|{{.MemTotal}}|{{.Architecture}}"])
    if r["rc"] != 0 or "|" not in r["out"]:
        return {"ok": False, "detail": "docker client present but the daemon is not running",
                "fix": "start Docker Desktop (`open -a Docker`) and wait until `docker info` succeeds"}
    ver, ncpu, mem, arch = r["out"].splitlines()[-1].split("|")
    return {"ok": True, "detail": f"daemon {ver}", "ncpu": int(ncpu), "mem_gb": int(mem) / 1024 ** 3, "arch": arch}


def nextflow_status() -> Dict[str, Any]:
    r = _run(["nextflow", "-version"], timeout=60)
    if r["rc"] != 0:
        return {"ok": False, "detail": r["out"][:200], "fix": "install Nextflow: `brew install nextflow` or `curl -s https://get.nextflow.io | bash`"}
    ver = next((l.strip() for l in r["out"].splitlines() if "version" in l), "")
    return {"ok": True, "detail": ver}


def sesame_status() -> Dict[str, Any]:
    r = _run(["Rscript", "-e", "suppressMessages(library(sesame)); cat('ok')"], timeout=120)
    if r["rc"] == 0 and r["out"].endswith("ok"):
        return {"ok": True, "detail": "sesame loads"}
    return {"ok": False, "detail": (r["out"].splitlines() or ["Rscript failed"])[-1][:200],
            "fix": "in R: install.packages('BiocManager'); BiocManager::install('sesame'); sesame::sesameDataCache()"}


def network_status(urls: List[str]) -> Dict[str, Any]:
    import urllib.error
    res = {}
    for u in urls:
        try:
            with urllib.request.urlopen(urllib.request.Request(u, method="HEAD"), timeout=8) as r:
                res[u] = True
        except urllib.error.HTTPError:
            res[u] = True          # the server answered (401/403/405): reachable
        except Exception:  # noqa: BLE001 - timeouts, resets, DNS: unreachable
            res[u] = False
    return {"ok": all(res.values()), "detail": res}


PROBE_IMAGE = "community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690"  # ~20 MB, used by taxprofiler 2.0.1


def registry_status() -> Dict[str, Any]:
    """Pull one small nf-core container through the docker daemon.

    The daemon has its own network route (proxy settings), so a host-side HTTP check can pass while every
    container pull fails; only a real pull tells the truth.
    """
    r = _run(["docker", "pull", "--quiet", "--platform", "linux/amd64", PROBE_IMAGE], timeout=180)
    if r["rc"] == 0:
        return {"ok": True, "detail": f"docker pulled {PROBE_IMAGE.split('/')[-1]}"}
    return {"ok": False, "detail": "docker cannot pull from community.wave.seqera.io: " + r["out"].splitlines()[-1][:160],
            "fix": "give the Docker daemon a working route to community.wave.seqera.io / quay.io (Docker Desktop → Settings → Resources → Proxies), or pre-pull the images on another network"}


def git_status(pipes: List[Dict[str, Any]], check_network: bool) -> Dict[str, Any]:
    """Pipelines are fetched with the system git, so the user's git/ssh setup must reach GitHub."""
    from .pipelines import pipeline_cache
    missing = [p for p in pipes if not (pipeline_cache(p["repo"], p["revision"]) / "main.nf").exists()]
    if not missing:
        return {"ok": True, "detail": "all pinned pipelines already cached"}
    if not check_network:
        return {"ok": True, "detail": f"{len(missing)} pipeline(s) to fetch; network not checked"}
    bad = []
    for p in missing:
        r = _run(["git", "ls-remote", f"https://github.com/{p['repo']}.git", f"refs/tags/{p['revision']}"], timeout=40)
        if r["rc"] != 0 or not r["out"]:
            bad.append(f"{p['repo']}@{p['revision']}: {r['out'][:120]}")
    if bad:
        return {"ok": False, "detail": "; ".join(bad), "fix": "make `git ls-remote https://github.com/nf-core/sarek.git` work (network/proxy/ssh keys)"}
    return {"ok": True, "detail": f"{len(missing)} pipeline(s) reachable, fetched at launch"}


def _input_gb(st: Dict[str, Any], file_ids: List[str]) -> float:
    by = {f["id"]: f for f in st["files"]}
    total = 0
    for i in file_ids:
        if i in by:
            total += by[i]["size"]
        elif i.startswith("path:") and os.path.exists(i[5:]):
            total += os.path.getsize(i[5:])
    return total / 1024 ** 3


def plan_needed(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Which pipelines the current files call for (structural: kind + assigned modality)."""
    cfg = data("pipelines.json")["pipelines"]
    need: List[Dict[str, Any]] = []
    groups: Dict[str, List[str]] = {}
    for f in st["files"]:
        if f.get("status") != "needs_pipeline":
            continue
        mod = f.get("modality") or ("methylation_idat" if f["kind"] == "methylation_idat" else None)
        for pid, p in cfg.items():
            if mod in p["for_modalities"] and f["kind"] in p["inputs"]:
                groups.setdefault(pid, []).append(f["id"])
    for pid, ids in groups.items():
        need.append({"pipeline": pid, "file_ids": sorted(ids), "optional": False})
    if st["processed"].get("variants"):
        ids = [x["file_id"] or f"path:{x['path']}" for x in st["processed"]["variants"]][:1]
        for pid in ("pharmcat", "pgsc_calc"):
            need.append({"pipeline": pid, "file_ids": ids, "optional": True})
    return need


def run_preflight(st: Dict[str, Any], workdir: Path, check_network: bool = True) -> Dict[str, Any]:
    cfg = data("pipelines.json")
    pipes = cfg["pipelines"]
    arch = platform.machine()
    cpus = os.cpu_count() or 1
    mem = mem_gb()
    workdir.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(workdir).free / 1024 ** 3
    needed = plan_needed(st)
    engines = {pipes[n["pipeline"]]["engine"] for n in needed}
    env: Dict[str, Any] = {"arch": arch, "cpus": cpus, "mem_gb": round(mem, 1) if mem else None,
                           "free_disk_gb": round(free_gb, 1), "workdir": str(workdir)}
    tools: Dict[str, Any] = {}
    env["docker"] = docker_status()
    if engines & {"nextflow"}:
        tools["nextflow"] = nextflow_status()
        tools["docker"] = docker_status()
    if engines & {"java"}:
        tools["java"] = java_status()
    if engines & {"rscript"}:
        tools["sesame"] = sesame_status()
    if check_network and needed:
        tools["network"] = network_status(["https://github.com", "https://quay.io", "https://ngi-igenomes.s3.amazonaws.com"])
    if engines & {"nextflow"} and check_network and tools.get("docker", {}).get("ok"):
        tools["registry"] = registry_status()
    if engines & {"nextflow"}:
        tools["git"] = git_status([pipes[n["pipeline"]] for n in needed if pipes[n["pipeline"]]["engine"] == "nextflow"], check_network)

    results = []
    for n in needed:
        p = pipes[n["pipeline"]]
        reasons: List[str] = []
        fixes: List[str] = []
        tier = GREEN
        in_gb = _input_gb(st, n["file_ids"])
        eff_cpu = cpus
        if p["engine"] == "nextflow":
            for t in ("nextflow", "docker", "git", "registry"):
                if t in tools and not tools[t]["ok"]:
                    tier = RED
                    reasons.append(f"{t}: {tools[t]['detail']}")
                    fixes.append(tools[t]["fix"])
            if tools["docker"].get("ok"):
                eff_cpu = min(cpus, tools["docker"]["ncpu"])
                if tools["docker"]["mem_gb"] < p["mem_gb_min"]:
                    tier = RED
                    reasons.append(f"docker VM has {tools['docker']['mem_gb']:.1f} GB RAM, pipeline needs ≥{p['mem_gb_min']} GB")
                    fixes.append(f"raise Docker Desktop memory to ≥{p['mem_gb_recommended']} GB (Settings → Resources)")
        elif p["engine"] == "java" and not tools["java"]["ok"]:
            tier = RED
            reasons.append(f"java: {tools['java']['detail']}")
            fixes.append(tools["java"]["fix"])
        elif p["engine"] == "rscript" and not tools["sesame"]["ok"]:
            tier = RED
            reasons.append(f"sesame: {tools['sesame']['detail']}")
            fixes.append(tools["sesame"]["fix"])

        if mem is not None and mem < p["mem_gb_min"]:
            tier = RED
            reasons.append(f"machine RAM {mem:.0f} GB < minimum {p['mem_gb_min']} GB")
        elif mem is not None and mem < p["mem_gb_recommended"] and tier != RED:
            tier = YELLOW
            reasons.append(f"RAM {mem:.0f} GB below recommended {p['mem_gb_recommended']} GB (slower, may retry)")

        refs_missing = [r for r in p["references"] if not (os.environ.get(r.get("env", "")) and Path(os.environ[r["env"]]).exists())]
        ref_gb = sum(r["size_gb"] for r in refs_missing)
        work_gb = in_gb * p["work_disk_factor"] + 1
        disk_need = work_gb + ref_gb
        if disk_need > free_gb:
            tier = RED
            reasons.append(f"needs ~{disk_need:.0f} GB free, {free_gb:.0f} GB available")
            fixes.append("free disk space or point the workspace at a larger volume")
        elif disk_need > 0.8 * free_gb and tier != RED:
            tier = YELLOW
            reasons.append(f"uses ~{disk_need:.0f} of {free_gb:.0f} GB free disk")
        if refs_missing:
            if tier == GREEN:
                tier = YELLOW
            reasons.append("references to download first: " + ", ".join(f"{r['label']} (~{r['size_gb']} GB)" for r in refs_missing))
            if check_network and tools.get("network") and not tools["network"]["ok"] and tier != RED:
                tier = RED
                reasons.append("network checks failed, references cannot be downloaded")

        b = p["benchmark"]
        cpu_h = b.get("cpu_hours_fixed", 0) + b.get("cpu_hours_per_input_gb", 0) * in_gb
        factor = 1.0
        if p.get("containers_amd64_only") and arch in ("arm64", "aarch64"):
            factor = cfg["arch_emulation_factor"]["value"]
            if tier == GREEN:
                tier = YELLOW
            reasons.append(f"{arch} host running amd64 containers under emulation (×{factor} planning factor, not measured here)")
        use_cpu = max(1, min(eff_cpu, p["cpus_recommended"] * 2))
        hours = cpu_h / (use_cpu * 0.7) * factor
        if hours > 24 and tier == GREEN:
            tier = YELLOW
            reasons.append(f"estimated {hours:.0f} h wall time")
        results.append({
            "pipeline": n["pipeline"], "label_zh": p["label_zh"], "optional": n["optional"],
            "revision": p["revision"], "file_ids": n["file_ids"], "input_gb": round(in_gb, 3),
            "tier": tier, "reasons": reasons, "fixes": sorted(set(fixes)),
            "estimate": {"cpu_hours": round(cpu_h, 2), "wall_hours": round(hours, 2), "cores_used": use_cpu,
                          "disk_gb": round(disk_need, 1), "basis": b["basis"], "is_estimate": b.get("estimate", False)},
            "commercial_ok": p.get("commercial_ok", True),
            "note": p.get("note") or p.get("needs_user_choice"),
        })
    # Files whose assay is unknown still get a what-if estimate, so the user hears how long each possibility takes.
    deferred = []
    for f in st["files"]:
        if f.get("status") in ("deferred_ask_lab", "excluded_not_member") and f.get("kind") in ("reads_fastq", "alignment"):
            gb = f["size"] / 1024 ** 3
            opts = []
            for pid, p in pipes.items():
                if f["kind"] in p["inputs"] and p["engine"] == "nextflow" and "variants" not in p["for_modalities"]:
                    b = p["benchmark"]
                    cpu_h = b.get("cpu_hours_fixed", 0) + b.get("cpu_hours_per_input_gb", 0) * gb
                    opts.append({"if_modality": p["for_modalities"], "pipeline": pid, "cpu_hours": round(cpu_h, 3),
                                 "reference_download_gb": sum(r["size_gb"] for r in p["references"]), "is_estimate": b.get("estimate", False)})
            deferred.append({"file": f["id"], "name": f["name"], "size_gb": round(gb, 4), "status": f.get("status"), "what_if": opts})
    overall = RED if any(r["tier"] == RED and not r["optional"] for r in results) else (
        YELLOW if any(r["tier"] != GREEN for r in results) else GREEN)
    return {"checked_at": now_iso(), "env": env, "tools": tools, "pipelines": results, "deferred_files": deferred,
            "overall": overall if results else "none_needed"}


def render_md(pf: Dict[str, Any]) -> str:
    icon = {GREEN: "🟢", YELLOW: "🟡", RED: "🔴"}
    e = pf["env"]
    lines = [f"# 运行前检查（preflight）  {pf['checked_at']}", "",
             f"本机：{e['arch']}，{e['cpus']} 核，内存 {e['mem_gb']} GB，工作目录可用磁盘 {e['free_disk_gb']} GB。", ""]
    if pf.get("deferred_files"):
        lines += ["检测类型未知、暂不运行的测序文件（按可能的检测类型给出估算，真实耗时取决于确认后的类型）：", "",
                  "| 文件 | 大小 | 若是 | 流程 | CPU·h（估算） | 需先下载参考 |", "|---|---|---|---|---|---|"]
        for d in pf["deferred_files"]:
            for o in d["what_if"]:
                lines.append(f"| {d['name']} | {d['size_gb']} GB | {'/'.join(o['if_modality'])} | {o['pipeline']} | {o['cpu_hours']} | 约 {o['reference_download_gb']} GB |")
        lines.append("")
    if not pf["pipelines"]:
        lines.append("当前数据不需要跑任何生信流程（都已是可直接分析的中间格式）。")
        return "\n".join(lines)
    lines += ["| 流程 | 结论 | 预计耗时 | 磁盘 | 原因 |", "|---|---|---|---|---|"]
    for r in pf["pipelines"]:
        est = r["estimate"]
        lines.append(f"| {r['label_zh']} ({r['pipeline']} {r['revision']}){'（可选）' if r['optional'] else ''} | {icon[r['tier']]} "
                     f"| 约 {est['wall_hours']} 小时{'（估算）' if est['is_estimate'] else ''} | 约 {est['disk_gb']} GB | {'；'.join(r['reasons']) or '—'} |")
    fixes = sorted({f for r in pf["pipelines"] for f in r["fixes"]})
    if fixes:
        lines += ["", "要先解决：", *[f"- {f}" for f in fixes]]
    notes = [f"- {r['pipeline']}: {r['note']}" for r in pf["pipelines"] if r.get("note")]
    if notes:
        lines += ["", "注意：", *notes]
    lines += ["", "耗时 = 基准 CPU·h ÷（可用核数 × 0.7）× 架构系数。基准多为估算，跑完后按实测回填。"]
    return "\n".join(lines)
