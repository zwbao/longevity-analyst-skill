# 02 · Preflight and local pipelines

## Use When
Some files need a pipeline (`observe` shows files `needs_pipeline`, or a VCF invites the optional PGx/PRS runs),
or the user asks whether their machine can run the analysis.

## Inputs Needed
- Required: a workspace past intake.
- Optional: reference paths already on disk, exported as `LA_REF_GRCH38`, `LA_REF_METAPHLAN_DB`, `LA_REF_BISMARK`,
  `LA_REF_PGSC_PANEL`, `LA_PHARMCAT_JAR`, `LA_SESAME_CACHE`.

## Workflow
1. `$LA preflight <ws>` writes `work/preflight.md`. **Show that table to the user as it is.** It lists, for each pipeline:
   green / yellow / red, the estimated wall time, the disk it needs and the reasons. Say that the times are estimates.
2. Red: tell the user what blocks it and the listed fix. Offer the alternative: ask the lab for the intermediate file
   (VCF/CRAM, β values, MetaPhlAn table), or run on a server. Do not launch.
   A run that will not happen (red, or the user said no) is closed with
   `$LA pipeline skip <ws> --run-id <id> --reason "<why, quoting the user if they declined>"`. The report lists it.
3. Yellow or green: ask plainly, e.g. "Run sarek locally? Estimated 19 h, 310 GB disk. Yes/No".
   Optional pipelines (PharmCAT, pgsc_calc) are asked separately. pgsc_calc also needs the user to confirm which
   PGS Catalog scores to compute; propose them, wait, then plan with `--pgs-id PGS000018,...`.
4. `$LA pipeline plan <ws>`, then after a yes:
   `$LA pipeline run <ws> --run-id <id> --confirmed "<the user's words>" [--accept-yellow]`.
   The run is detached; tell the user it runs in the background and how to check it.
5. Check with `$LA pipeline status <ws> [--run-id <id>]` (task counts + log tail). An interrupted run is resumed
   by launching it again (`-resume` is built in).
6. When it has exited: `$LA pipeline verify <ws> --run-id <id>`. Only a verified run registers its outputs
   (with pipeline revision and params in the provenance) for the method stage.
7. `--stub` (developers only, needs `LA_DEV=1`) runs `nextflow -stub-run -profile test,docker`: a wiring check
   with no real compute. Its outputs go to `results_stub/` and can never be registered as member results.
8. A skipped run is only started again with `--reopen` and the user's new words.
9. Stage order is enforced: `methods` refuses to start until every required pipeline is verified or skipped. When
   nothing needs to run, `preflight` closes the pipelines stage by itself.

## Output Format
1. The preflight table.
2. The user's decision per pipeline, quoted.
3. Run ids with their status; verified outputs.

## Guardrails
- Never add parameters, change the revision or swap the container profile; the plan is fixed by `data/pipelines.json`.
- A pipeline that did not run is not an error in the report. It is listed with its reason.
