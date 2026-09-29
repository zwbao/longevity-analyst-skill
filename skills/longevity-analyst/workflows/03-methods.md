# 03 · Methods

## Use When
Intake is resolved and every pipeline the user agreed to has been verified (or declined).

## Workflow
1. `$LA methods plan <ws>`. Each candidate method gets one status:
   - `ready`: inputs present, gates passed.
   - Recording an answer with `la.py member` voids the method plan and everything after it. Plan again; this is
     expected, not an error.
   - `native.apoe` asks `genetic_disclosure=yes|no` first. APOE ε4 is an Alzheimer-associated finding, and the member
     decides whether it appears in the report.
   - `needs_answers`: ask the user the listed questions (for example China-PAR asks for antihypertensive use, smoking,
     diabetes, north/south). Record each with `$LA member <ws> key=value --source "<quote>"` and plan again.
   - `input_problem`: a lab row did not match or has a wrong unit. Look at `unmatched_lab_names` and `problems`.
     If a row clearly is the missing input under another name, map it:
     `$LA labs map <ws> --marker "超敏C反应蛋白(hs-CRP)" --to crp_mg_dl --reason "same analyte; unit mg/L stated"`.
     Never change a value or invent a unit.
   - `ready` items may list `optional_unanswered` inputs. Some are needed for part of the population (China-PAR's male
     equation needs `urban` and `family_history`). Ask them if the user can answer; otherwise the method may
     stop with a problems file, which is reported as not computed.
   - `blocked_platform`: the data platform does not match the method's training assay. Final. Report it.
   - `blocked_license`: the method takes an open clock menu (pyaging) or is excluded in this licence mode. Final.
   - `manual`: the skill's inputs are not declared. It is not auto-run and is listed as not computed.
2. `$LA methods run <ws>`. It runs every ready method and writes `work/readouts.json`. Each readout carries a
   stable id, value, unit, kind (`computed` / `descriptive`), method and body systems.
3. Failed or `no_output` methods: read their `work/methods/<name>/report.md` or `stderr.txt` and note why.
   Do not rerun them with altered inputs to force a number.

## Output Format
1. Counts by status.
2. The readouts that matter most for this member (by id).
3. What could not run and why.

## Guardrails
- Ask the user every `optional_unanswered` answer before `methods run`: some methods need them in practice (China-PAR's
  male equation fails without `urban` and `family_history`), and a method that fails reports only its problems file.
- `methods plan` lists `withheld_lab_rows`: rows you answered `no` in `labs confirm` (or left unanswered). They reach
  no method. Qualitative results (阴性, +) never reach a method.
- Commercial mode runs only the provisional clock list in `data/license_policy.json` (GrimAge, DunedinPACE, DNAm
  PhenoAge and derivatives, skin & blood and DNAmTL are excluded). That list has not been cleared by counsel, and
  the report says so. `--mode research` is for internal R&D only.
- An input file changed after intake (raw or derived) stops `methods run`; start a new workspace.
- A clock age outside any human range marks the readout `implausible`. Several such values stop the run, because the
  input was not what it claimed to be.
- Clock readouts carry their CpG coverage; below 90% they are marked `computed_low_coverage` and shown as
  "仅供参考". Do not build a finding on a low-coverage clock.
- Native modules: APOE genotype (needs both rs429358 and rs7412 called), gut diversity, GMHI (Gupta 2020, reproduces
  the published GMHI.R output), proteomics QC description (no protein age).
