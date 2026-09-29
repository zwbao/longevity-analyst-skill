"""Decide which methods can run on this member's data, then run them.

Planning is structural: a method is runnable when its declared inputs are
present (matched with the longevity-skills skillkit itself), its data
platform and value scale pass the platform gate, and the licence mode
allows it. The agent resolves only what needs judgment (a lab name the
aliases do not cover, a questionnaire answer the user must give).
"""
from __future__ import annotations

import csv
import gzip
import io
import os
import re
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from . import native
from .common import (EXIT_BLOCKED, EXIT_INPUT, LAError, data, load_json, longevity_skills_home,
                     now_iso, sha256_file, skillkit, write_json)

CPG_KEY = re.compile(r"^cg\d{8}$")
AUTO_TIERS = {"A", "tool"}
COVERAGE_ROW = re.compile(r"`([a-z0-9_]+)`[）)]?\s*\|[^|]*\|[^|]*\|\s*([\d.]+)%（缺\s*(\d+)/(\d+)）")


def _name_for_kit(name: str) -> str:
    from .labnames import name_for_kit
    return name_for_kit(name)


def _lab_rows(st: Dict[str, Any], blocked_markers: Optional[List[str]] = None,
              withheld: Optional[List[Dict[str, str]]] = None) -> List[Dict[str, str]]:
    """Lab rows the methods may match: a row that names a known analyte reaches a method only after the agent
    confirmed it (`la.py labs confirm`) or mapped it (`la.py labs map`). Qualitative results (阴性, +) never do."""
    from .labnames import usable_rows
    ok = {id(r) for r in usable_rows(st)}
    kit_fold = lambda s: re.sub(r"[\s_\-·]+", "", s.casefold())
    bad = [kit_fold(b) for b in (blocked_markers or [])]
    out = []
    for r in st["labs"]:
        name = r.get("maps_to") or _name_for_kit(r["marker"])
        if any(b in kit_fold(r["marker"]) or b in kit_fold(name) for b in bad):
            continue                     # a licensed clock value typed in as a "lab" does not enter commercial methods
        if not re.search(r"\d", str(r.get("value", ""))):
            continue                     # 阴性 / 阳性 / + : not a number any method can use
        if id(r) not in ok:              # every numeric row needs the agent's yes (or a map); nothing is matched unconfirmed
            if withheld is not None:
                withheld.append({"marker": r["marker"], "why": "answered no" if (r.get("confirm") or {}).get("answer") == "no" else "not confirmed"})
            continue
        from .labnames import clean_value
        out.append({"marker": name, "value": clean_value(r["value"]), "unit": r.get("unit", ""), "source_file": r.get("source_file")})
    return out


def pick(st: Dict[str, Any], modality: str) -> Optional[Dict[str, Any]]:
    """The one processed entry for a modality: the only one, or the one marked primary."""
    lst = st["processed"].get(modality) or []
    if len(lst) == 1:
        return lst[0]
    prim = [x for x in lst if x.get("primary")]
    return prim[0] if len(prim) == 1 else None


def _beta_lookup(path: Path, keys: List[str]) -> Dict[str, float]:
    want = set(keys)
    out: Dict[str, float] = {}
    fh = io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8-sig")
    with fh:
        for line in fh:
            parts = line.strip().split(",")
            if parts and parts[0] in want and len(parts) > 1:
                try:
                    out[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return out


def beta_sanity(path: Path, home: Path) -> Optional[str]:
    """Real blood β values correlate strongly with the whole-blood median profile shipped with epiage.

    Detection p-values, 1-β, constants, M-values squashed into [0,1] or another tissue do not. Returns a reason to
    refuse, or None.
    """
    import statistics
    ref_path = home / "skills" / "epiage" / "data" / "sesame_450k_median.csv"
    if not ref_path.exists():
        return "reference blood profile (epiage/data/sesame_450k_median.csv) not found; cannot verify β values"
    ref = {}
    for line in open(ref_path, encoding="utf-8"):
        p = line.strip().split(",")
        try:
            ref[p[0]] = float(p[1])
        except (ValueError, IndexError):
            pass
    got = _beta_lookup(path, list(ref))
    keys = [k for k in got if k in ref]
    if len(keys) < 1000:
        return f"only {len(keys)} of the {len(ref)} reference CpGs are present; too few to verify β values or run clocks"
    xs = [got[k] for k in keys]
    ys = [ref[k] for k in keys]
    if statistics.pstdev(xs) < 1e-6:
        return "all β values are the same number"
    r = statistics.correlation(xs, ys)
    if r < 0.9:
        return (f"correlation with the whole-blood reference profile is {r:.2f} (real blood arrays give 0.95-0.99); "
                "these are not whole-blood β values (p-values, 1-β, another tissue or a wrong column?)")
    slope, icpt = statistics.linear_regression(ys, xs)
    mad = statistics.mean(abs(x - y) for x, y in zip(xs, ys))
    # real whole-blood arrays across 7 GEO series / 3 preprocessing styles: r 0.955-0.987, slope 0.80-1.02,
    # intercept -0.016..0.067, MAD 0.040-0.116. The limits below only catch gross rescaling (x0.5, beta^2, min-max
    # squashed M-values); small global shifts are indistinguishable from preprocessing and are a known limitation.
    if not (0.7 <= slope <= 1.15) or abs(icpt) > 0.09 or mad > 0.14:
        return (f"β values are shifted or rescaled against the whole-blood reference (slope {slope:.2f}, intercept "
                f"{icpt:.3f}, mean |Δ| {mad:.3f}; real arrays: about 0.9 / 0.01 / 0.05). Clocks would be biased; ask the "
                "lab for unmodified SeSAMe/minfi β values")
    return None


def _member_value(member: Dict[str, Any], key: str) -> Optional[str]:
    if key in ("age", "sex"):
        v = member.get(key)
        return None if v in (None, "") else str(v)
    return (member.get("answers") or {}).get(key)


def _prot_gate(name: str, prot: Optional[Dict[str, Any]], gate: Dict[str, Any]) -> Optional[str]:
    rule = gate["proteomics"]["skills"].get(name, gate["proteomics"]["default"])
    plat = prot["provenance"].get("platform") if prot else None
    scale = prot["provenance"].get("value_scale") if prot else None
    if plat not in rule["allowed_platforms"]:
        return f"member proteomics platform {plat} not in {rule['allowed_platforms'] or 'any declared platform'}: {rule['reason']}"
    need = rule.get("required_value_scale")
    if need and scale != need:
        return f"platform matches but values are {scale}, method needs {need}: {gate['proteomics'].get('note_value_scale', '')}"
    return None


def plan(st: Dict[str, Any], ws_root: Path) -> Dict[str, Any]:
    home = longevity_skills_home()
    kit = skillkit()
    cat = load_json(home / "catalog.json")
    gate = data("platform_requirements.json")["modalities"]
    pol = data("license_policy.json")
    mode = st["member"].get("mode", "commercial")
    lic = pol["modes"][mode]
    member = st["member"]
    from .labnames import pending as labs_pending
    withheld: List[Dict[str, str]] = []
    labs = _lab_rows(st, lic.get("blocked_lab_markers"), withheld)
    lab_pending = labs_pending(st)
    confirm_reason = (f"{len(lab_pending)} lab rows are not yet confirmed: `la.py labs candidates <ws>`, then "
                      "`la.py labs confirm <ws> --answers <json>` (workflows/01-intake.md)")
    uncertain_files = {f["id"] for f in st["files"] if f.get("provenance_uncertain")}
    labs_uncertain = any(r.get("source_file") in uncertain_files for r in labs)
    methyl = pick(st, "methylation")
    prot = pick(st, "proteomics")
    gut = pick(st, "gut_profile")
    variants = pick(st, "variants")
    items: List[Dict[str, Any]] = []
    unmatched_names: set = set()
    tissue_ok = bool(methyl) and methyl["provenance"].get("tissue") in gate["methylation"]["default"]["allowed_tissues"]
    beta_problem = beta_sanity(Path(methyl["path"]), home) if (methyl and tissue_ok) else None
    if beta_problem:
        tissue_ok = False
        items.append({"method": "(all methylation methods)", "tier": "-", "status": "blocked_platform", "reasons": [beta_problem]})

    for s in cat["skills"]:
        name = s["name"]
        if not s.get("has_script") or s["tier"] not in AUTO_TIERS:
            continue
        intents_data = {d for i in cat["intents"] if i["id"] in s.get("intents", []) for d in i["data"]}
        inputs = s.get("inputs") or []
        meas = [i for i in inputs if i.get("from") == "measurements"]
        is_cpg = any(CPG_KEY.match(i["key"]) for i in meas)
        uses_prot = "proteomics" in intents_data or "蛋白、代谢与器官年龄" in s.get("domains", [])
        takes_clock_menu = any(i["key"] == "clocks" for i in inputs)
        item: Dict[str, Any] = {"method": name, "tier": s["tier"], "domains": s.get("domains", []),
                                "blurb_zh": s.get("blurb_zh", ""), "status": None, "reasons": []}

        if name == "epiage":
            if not methyl:
                continue
            if not tissue_ok:
                why = beta_problem or f"tissue {methyl['provenance'].get('tissue')}: {gate['methylation']['default']['reason']}"
                item.update(status="blocked_platform", reasons=[why])
            elif not member.get("age") or not member.get("sex"):
                item.update(status="needs_answers", reasons=["age and sex are required"], needs=["age", "sex"])
            else:
                item.update(status="ready", run={"kind": "epiage", "betas": methyl["path"], "clocks": lic["epiage_clocks"]},
                            uncertain=bool(methyl.get("provenance_uncertain")))
            items.append(item)
            continue
        if takes_clock_menu:
            # pyaging and similar tools pick clocks from an open menu; licence and tissue checks cannot be enforced there.
            if methyl or prot:
                item.update(status="blocked_license", reasons=["takes a free clock menu (may include GrimAge/PC clocks); not auto-run"])
                items.append(item)
            continue

        takes_age_value = any(i.get("unit") == "a" and i["key"] != "age" for i in meas)
        if takes_age_value and mode == "commercial":
            if labs:
                item.update(status="blocked_input_origin", reasons=["takes a pre-computed biological age typed in as a lab value; "
                                                                    "its origin (possibly a licensed clock) cannot be checked"])
                items.append(item)
            continue
        have = {"routine_labs": bool(labs), "methylation": bool(methyl), "proteomics": bool(prot), "genotype": bool(variants)}
        if s.get("inputs_status") != "verified":
            if any(have.get(d) for d in intents_data):
                why = _prot_gate(name, prot, gate) if (uses_prot and prot) else None
                if why:
                    item.update(status="blocked_platform", reasons=[why])
                else:
                    item.update(status="manual", reasons=["inputs not declared in skill.json (inputs_status != verified); not auto-run"])
                items.append(item)
            continue

        if not meas:
            req_args = [i["key"] for i in inputs if i.get("from") == "argument" and i.get("required")]
            if not req_args or any(_member_value(member, k) is None for k in req_args):
                continue
        if uses_prot:
            if not prot:
                continue
            why = _prot_gate(name, prot, gate)
            if why:
                item.update(status="blocked_platform", reasons=[why])
                items.append(item)
                continue

        manifest = load_json(home / "skills" / name / "skill.json")
        if is_cpg:
            if not methyl:
                continue
            if not tissue_ok:
                item.update(status="blocked_platform", reasons=[beta_problem or "methylation tissue is not whole blood"])
                items.append(item)
                continue
            betas = _beta_lookup(Path(methyl["path"]), [i["key"] for i in meas if CPG_KEY.match(i["key"])])
            rows = [{"marker": k, "value": str(v), "unit": ""} for k, v in betas.items()] + labs
        else:
            rows = labs
        if meas and lab_pending and not is_cpg and "routine_labs" in intents_data:
            item.update(status="needs_answers", reasons=[confirm_reason])
            items.append(item)
            continue
        if meas:
            col = kit.collect_measurements(rows, manifest)
            missing = [p.key for p in col.problems if p.kind == "missing"]
            other = [p.as_dict() for p in col.problems if p.kind != "missing"]
            req = [i["key"] for i in meas if i.get("required")]
            have_any = any(k in col.values for k in req) or other
            if (missing and not have_any) or (not col.values and not other):
                lab_missing = [k for k in missing if any(i["key"] == k and i.get("loinc") for i in meas)]
                if labs and not is_cpg and "routine_labs" in intents_data and lab_missing and len(lab_missing) * 2 >= len(missing):
                    unmatched_names.update(col.unmatched)
                    item.update(status="input_problem",
                                reasons=[f"no lab row matched this method's inputs ({', '.join(missing)}); check `unmatched_lab_names` and map real equivalents with `la.py labs map`"],
                                unmatched_rows=col.unmatched[:30])
                    items.append(item)
                continue
            if missing or other:
                unmatched_names.update(col.unmatched)
                item.update(status="input_problem", reasons=[f"missing: {', '.join(missing)}"] if missing else [],
                            problems=other, unmatched_rows=col.unmatched[:30])
                items.append(item)
                continue
            item["matched"] = {k: col.sources[k] for k in col.values}
        args_missing = []
        arg_vals: Dict[str, str] = {}
        for i in inputs:
            if i.get("from") not in ("argument", "profile") or i["key"] in pol["reserved_answer_keys"]:
                continue
            v = _member_value(member, i["key"])
            if v is None:
                if i.get("required"):
                    args_missing.append(i["key"])
            else:
                arg_vals[i.get("flag") or f"--{i['key']}"] = v
        if args_missing:
            item.update(status="needs_answers", needs=args_missing,
                        reasons=["ask the user: " + ", ".join(f"{i['key']} ({i.get('label_zh', '')})" for i in inputs if i["key"] in args_missing)])
            items.append(item)
            continue
        optional_unanswered = [f"{i['key']} ({i.get('label_zh', '')}{'；' + i['note_zh'] if i.get('note_zh') else ''})"
                               for i in inputs if i.get("from") == "argument" and not i.get("required")
                               and _member_value(member, i["key"]) is None and i["key"] not in pol["reserved_answer_keys"]]
        item.update(status="ready", run={"kind": "skill", "rows": rows if meas else [], "args": arg_vals, "cpg": is_cpg},
                    optional_unanswered=optional_unanswered,
                    uncertain=labs_uncertain or bool(is_cpg and methyl and methyl.get("provenance_uncertain")))
        items.append(item)

    # native methods
    if lab_pending:
        items.append({"method": "native.organ_indices", "tier": "native", "status": "needs_answers", "reasons": [confirm_reason]})
    elif labs:
        items.append({"method": "native.organ_indices", "tier": "native", "status": "ready", "reasons": [],
                      "uncertain": labs_uncertain, "run": {"kind": "native", "fn": "organ_indices"}})
    if variants:
        consent = (member.get("answers") or {}).get("genetic_disclosure")
        if consent not in ("yes", "no"):
            items.append({"method": "native.apoe", "tier": "native", "status": "needs_answers", "needs": ["genetic_disclosure"],
                          "reasons": ["ask whether the member wants APOE (Alzheimer-associated) results in the report: "
                                      "`la.py member <ws> genetic_disclosure=yes|no --source \"<quote>\"`"]})
        elif consent == "no":
            items.append({"method": "native.apoe", "tier": "native", "status": "declined", "reasons": ["member declined APOE results"]})
        else:
            items.append({"method": "native.apoe", "tier": "native", "status": "ready", "reasons": [], "uncertain": bool(variants.get("provenance_uncertain")),
                          "run": {"kind": "native", "fn": "apoe", "path": variants["path"], "assembly": variants["provenance"].get("assembly")}})
        items.append({"method": "native.vcf_summary", "tier": "native", "status": "ready", "reasons": [],
                      "run": {"kind": "native", "fn": "vcf_summary", "path": variants["path"]}})
    if gut:
        items.append({"method": "native.gut_diversity", "tier": "native", "status": "ready", "reasons": [], "uncertain": bool(gut.get("provenance_uncertain")),
                      "run": {"kind": "native", "fn": "gut_diversity", "path": gut["path"]}})
        ver = (gut["provenance"].get("mpa_version") or "").lower()
        body = Path(gut["path"]).read_text(encoding="utf-8", errors="replace")
        if "|t__" in body or "SGB" in body:
            ver = "mpa_v3_or_later (SGB/strain rows found)"
        if ver.startswith("mpa_v2"):
            items.append({"method": "native.gmhi", "tier": "native", "status": "ready", "reasons": [], "uncertain": bool(gut.get("provenance_uncertain")),
                          "run": {"kind": "native", "fn": "gmhi", "path": gut["path"]}})
        else:
            items.append({"method": "native.gmhi", "tier": "native", "status": "blocked_platform",
                          "reasons": [f"GMHI species names are MetaPhlAn2; this profile is {ver or 'of unknown MetaPhlAn version'}; renamed species would be missed and can flip the sign"]})
    if prot:
        items.append({"method": "native.proteomics_describe", "tier": "native", "status": "ready", "reasons": [],
                      "run": {"kind": "native", "fn": "proteomics_describe", "path": prot["path"], "platform": prot["provenance"].get("platform")}})
    excluded_kinds = {f.get("kind") for f in st["files"] if f.get("excluded")}
    for kind, mod, names in (("variants_vcf", "variants", ["native.apoe", "native.vcf_summary"]),
                             ("metaphlan_profile", "gut_profile", ["native.gut_diversity", "native.gmhi"]),
                             ("protein_matrix", "proteomics", ["native.proteomics_describe"]),
                             ("methylation_beta", "methylation", ["epiage"])):
        if kind in excluded_kinds and not pick(st, mod):
            for n in names:
                if not any(i["method"] == n for i in items):
                    items.append({"method": n, "tier": "-", "status": "not_run_excluded",
                                  "reasons": ["the only input file of this kind was excluded as not this member's"]})
    for mod in ("methylation", "proteomics", "gut_profile", "variants"):
        if len(st["processed"].get(mod) or []) > 1 and pick(st, mod) is None:
            items.append({"method": f"(all {mod} methods)", "tier": "-", "status": "needs_answers",
                          "reasons": [f"several {mod} files and none chosen as this visit's: `la.py assign --what primary --file {mod} --value <file id>`"]})

    counts: Dict[str, int] = {}
    for it in items:
        counts[it["status"]] = counts.get(it["status"], 0) + 1
    res = {"planned_at": now_iso(), "mode": mode, "clearance": lic.get("clearance"), "counts": counts, "items": items,
           "unmatched_lab_names": sorted(unmatched_names), "withheld_lab_rows": withheld, "license_excluded": lic["excluded"]}
    st["methods"] = res
    return res


def _python() -> str:
    return os.environ.get("LA_PYTHON", sys.executable)


def _epiage_coverage(report_md: Path) -> Dict[str, Dict[str, Any]]:
    cov: Dict[str, Dict[str, Any]] = {}
    if not report_md.exists():
        return cov
    for line in report_md.read_text(encoding="utf-8").splitlines():
        m = COVERAGE_ROW.search(line)
        if m:
            key, pct, miss, n = m.group(1), float(m.group(2)), int(m.group(3)), int(m.group(4))
            cov[f"dnam_{key}"] = {"coverage_pct": pct, "missing_cpgs": miss, "model_cpgs": n}
    return cov


def _run_skill(home: Path, it: Dict[str, Any], member: Dict[str, Any], out: Path) -> Dict[str, Any]:
    name = it["method"]
    manifest = load_json(home / "skills" / name / "skill.json")
    entry = manifest["entry"]
    script = home / "skills" / name / entry["script"]
    cmd = [_python(), str(script)]
    run = it["run"]
    out.mkdir(parents=True, exist_ok=True)
    if run["kind"] == "epiage":
        cmd += ["--betas", run["betas"], "--clocks", ",".join(run["clocks"])]
    if run.get("rows"):
        mpath = out / "measurements.csv"
        hdr = entry.get("measurements_header") or ["marker", "value", "unit"]
        with open(mpath, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(hdr)
            for r in run["rows"]:
                w.writerow([r["marker"], r["value"], r["unit"]])
        cmd += [entry.get("measurements_flag", "--measurements"), str(mpath)]
    flags_done = set()
    for flag, v in (run.get("args") or {}).items():
        cmd += [flag, v]
        flags_done.add(flag)
    if entry.get("age_flag") and member.get("age") and entry["age_flag"] not in flags_done:
        cmd += [entry["age_flag"], str(member["age"])]
    if entry.get("sex_flag") and member.get("sex") and entry["sex_flag"] not in flags_done:
        cmd += [entry["sex_flag"], str(member["sex"])]
    cmd += [entry.get("out_flag", "--out"), str(out)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    (out / "stderr.txt").write_text(p.stderr[-20000:])
    res: Dict[str, Any] = {"exit": p.returncode, "cmd": cmd}
    rj = out / "result.json"
    if p.returncode == 0 and rj.exists():
        res["readouts"] = _readouts_from(name, manifest, load_json(rj))
        if name == "epiage":
            cov = _epiage_coverage(out / "report.md")
            allowed = {f"dnam_{c}" for c in run["clocks"]}
            kept = []
            for ro in res["readouts"]:
                key = ro["id"].split(".", 1)[1]
                if key not in allowed:
                    continue                      # never surface a clock the licence mode excludes
                c = cov.get(key)
                if c:
                    ro.update(c)
                    if c["coverage_pct"] < 50:
                        continue                 # mostly imputed: not a measurement of this person, not reported
                    if c["coverage_pct"] < 90:
                        ro["kind"] = "computed_low_coverage"
                else:
                    ro["kind"] = "computed_coverage_unknown"
                if ro.get("unit") == "a" and isinstance(ro["value"], (int, float)) and not 0 <= ro["value"] <= 125:
                    raise LAError(f"{ro['id']} = {ro['value']} years is outside any human range; the β input is not what "
                                  "it claims to be. Nothing from this run is reported.", EXIT_INPUT)
                kept.append(ro)
            res["readouts"] = kept
    elif (out / "problems.json").exists():
        res["problems"] = load_json(out / "problems.json")
    res["report_md"] = str(out / "report.md") if (out / "report.md").exists() else None
    return res


def _readouts_from(name: str, manifest: Dict[str, Any], result: Dict[str, Any]) -> List[Dict[str, Any]]:
    outs = {o["key"]: o for o in manifest.get("outputs", [])}
    vals = result.get("outputs", result) if isinstance(result.get("outputs", result), dict) else result
    ro = []
    for key, spec in outs.items():
        v = vals.get(key)
        if isinstance(v, dict):
            v = v.get("value")
        if v is None:
            continue
        ro.append({"id": f"{name}.{key}", "label_zh": spec.get("label_zh", key), "value": v, "unit": spec.get("unit", ""),
                   "kind": "computed", "note_zh": spec.get("note_zh"), "method": name, "tier": manifest.get("tier")})
    return ro


def run(st: Dict[str, Any], ws_root: Path, only: Optional[List[str]] = None,
        checkpoint: Optional[Callable[[], None]] = None) -> Dict[str, Any]:
    if not st.get("methods"):
        raise LAError("run `la.py methods plan` first", EXIT_INPUT)
    home = longevity_skills_home()
    member = st["member"]
    done = []
    sysmap = data("systems.json")
    changed = []
    for f in st["files"]:
        p = Path(f["path"])
        if f.get("excluded") or f.get("status") not in ("ready", "processed"):
            continue
        if not p.exists() or p.stat().st_size != f.get("size") or sha256_file(p) != f.get("sha256"):
            changed.append(f["name"])
    for lst in st["processed"].values():
        for x in lst:
            if x.get("sha256") and (not Path(x["path"]).exists() or sha256_file(Path(x["path"]), limit=None) != x["sha256"]):
                changed.append(x["path"])
    if changed:
        raise LAError(f"input files changed since intake: {changed[:5]}; start a new workspace", EXIT_BLOCKED)
    order = sorted(st["methods"]["items"], key=lambda i: 0 if i["method"] == "epiage" else 1)
    methyl_bad = None
    for it in order:
        if it["status"] not in ("ready", "failed") or (only and it["method"] not in only):
            continue
        uses_methyl = it["method"] == "epiage" or (it.get("run") or {}).get("cpg")
        if uses_methyl and methyl_bad:
            it.update(status="failed", result={"exit": 3, "error": methyl_bad}, reasons=[methyl_bad])
            continue
        out = ws_root / "work" / "methods" / it["method"]
        try:
            if it["run"]["kind"] == "native":
                r = it["run"]
                fn = getattr(native, r["fn"], None)
                if r["fn"] == "organ_indices":
                    from . import organs as _organs
                    res = _organs.compute_indices(st)
                elif r["fn"] == "apoe":
                    res = fn(Path(r["path"]), r.get("assembly"))
                elif r["fn"] == "proteomics_describe":
                    res = fn(Path(r["path"]), r.get("platform"))
                else:
                    res = fn(Path(r["path"]))
                out.mkdir(parents=True, exist_ok=True)
                write_json(out / "result.json", res)
                for ro in res.get("readouts", []):
                    ro.update(method=it["method"], tier="native")
                it["result"] = {"exit": 0, "readouts": res.get("readouts", []), "notes": res.get("notes", [])}
            else:
                it["result"] = _run_skill(home, it, member, out)
        except LAError as e:
            it["result"] = {"exit": e.code, "error": str(e)}
        except Exception as e:  # noqa: BLE001 - one bad input must not lose the whole run
            it["result"] = {"exit": 1, "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-2000:]}
        if it["method"] == "epiage" and "outside any human range" in str(it["result"].get("error", "")):
            methyl_bad = it["result"]["error"]      # the same β feeds every CpG method: none of them is reported
        ok = it["result"].get("exit") == 0
        it["status"] = "failed" if not ok else ("done" if it["result"].get("readouts") else "no_output")
        if it["status"] == "no_output":
            notes = it["result"].get("notes") or []
            it["reasons"] = notes or ["ran without error but produced no value for this member; see " + str(out / "report.md")]
        it["ran_at"] = now_iso()
        doms = it.get("domains") or []
        systems = sorted({s for d in doms for s in sysmap["domain_to_systems"].get(d, [])} |
                         set(sysmap["native_to_systems"].get(it["method"], [])))
        for ro in it["result"].get("readouts", []) or []:
            ro["systems"] = ro.get("systems") or systems
            ro["source"] = str(out / "result.json")
            if it.get("uncertain"):
                ro["provenance_uncertain"] = True
        done.append({"method": it["method"], "status": it["status"], "readouts": len(it["result"].get("readouts") or [])})
        if checkpoint:
            try:
                checkpoint()
            except OSError as e:           # disk full: stop cleanly, results so far stay in memory for the caller
                raise LAError(f"could not save progress ({e}); free space and re-run `methods run`", EXIT_INPUT)
    readouts: List[Dict[str, Any]] = []
    for it in st["methods"]["items"]:
        if it.get("status") == "done":
            readouts.extend(it["result"].get("readouts") or [])
    rp = ws_root / "work" / "readouts.json"
    write_json(rp, {"generated_at": now_iso(), "member": member.get("id"), "mode": st["methods"].get("mode"), "readouts": readouts})
    st["methods"]["readouts_sha256"] = sha256_file(rp, limit=None)      # trace refuses a readouts file edited after this
    return {"ran": done, "readouts": len(readouts)}
