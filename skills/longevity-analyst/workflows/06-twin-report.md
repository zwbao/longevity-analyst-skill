# 06 · Twin, trace, review, report

## Use When
Methods have run and the analyses and plan are registered (or deliberately skipped). This stage produces the
deliverables. It is also used later to compare a new visit with an old snapshot.

## Workflow
1. `$LA twin build <ws> [--exposures <la-exposures/1 file>]` writes `work/twin/twin.json` (`la-twin/2`): lab
   observations with their own sampling date (`t`, `t_source` row/file/workspace), readouts, data files with hashes,
   pipeline revisions, interventions (status and since from the exposures file), the retest plan with `judgeable`,
   and, when the genomics steps wrote them, `concordance`, `genotype_fingerprint` and `surveillance`. The exposures
   file comes from the host (LongPi site: adoption status only; coach member file: also weekly check-ins); a field
   outside the schema rejects it. A rebuild without `--exposures` keeps the ones registered before.
1b. A member with earlier visits: `$LA twin series --snapshots <old1/work/twin/twin.json> <old2/...> --ws <ws>` adds
   this workspace's twin and writes `work/twin/series.json` (report and export read it; run it again after any
   `twin build`). Per marker: adjacent pairs, first vs latest, and a trend test (≥ 3 points over ≥ 180 days);
   `summary_verdict` is `trend_beyond_noise`, `change_beyond_noise`, `within_noise` or `not_judged`. It also lists
   interventions under way at the same time (`concurrent_exposures`), changes of the gene-lab class (`transitions`)
   and `identity_alerts` (a visit whose germline genotype disagrees is left out of every comparison). The report
   section uses fixed wording; never rephrase a trend as worse, a disease or a cause. Snapshots of another member
   are refused. `la-twin/1` snapshots work with a warning (their observation dates are the snapshot date).
1c. Labelled items read the twin and the series, so run them now, in this order: `$LA predict response <ws>`
   (`--repeats N` for the experiment design, e.g. 14 home readings each side), `$LA clinical dx <ws> [--guideline-set cn|ada]`,
   `$LA clinical rx <ws>`. Trace refuses predictions or clinical items made before the latest exposures or series.
2. Write `work/report/summary.md` (Chinese, member-facing, at most about 300 characters of prose plus a short list):
   what stands out, where modalities agree or conflict, what is not known yet. Digits only inside placeholders:
   `{{n:…}}` for non-member numbers, `{{pmid:N}}` for citations, and member numbers only as `{{r:<id>}}` (prints value **and unit**, so do not repeat the unit); `{{r:<id>|label}}` prints the readout's
   name; `{{r:<id>|your own wording}}` prints your wording for its name (digits inside it are still traced).
3. `$LA review trace <ws>`. It fails on any digit, Chinese numeral run or URL outside placeholders, records a hash of
   every file the report is built from plus a digest of the state, and prints a `trace_id` and the list of `{{n:…}}`
   literals for the reviewer. Lab values have no placeholder: describe their direction ("高于参考范围"). The report
   never prints the free-text reasons you record in state; it prints fixed Chinese status wording instead. Dates,
   schedules and list numbers are not exceptions: write `{{n:2026 年 9 月}}`, `{{n:12 周后}}`, or use markdown list
   bullets. Each cited `{{pmid:N}}` is checked live against PubMed (retracted or unknown PMIDs fail), and doses in
   prose fail as in the plan. Fix the text rather than weakening it.
4. Dispatch an **independent reviewer subagent** with [`../references/reviewer.md`](../references/reviewer.md) and the
   workspace path. It returns a JSON verdict; save each round to a new file (`work/review/reviewer-1.json`,
   `reviewer-2.json`, …) and record it:
   `$LA review record <ws> --verdict pass|block --findings work/review/reviewer-N.json`.
   Every printed plan field (`action_zh`, `rationale_zh`, `retest.what`) follows the same digit rule as prose.
   The reviewer's JSON must carry the current `trace_id`, sit under `work/review/`, and its verdict is `pass` or
   `block` (any P0/P1 means block). `--verdict` must equal it.
   Any edit after the trace (summary, analyses, plan, twin, readouts) voids the trace and the review: re-run
   trace, send the changed files to the reviewer again, record again. `report` refuses otherwise.
   `block` means fix, re-run trace and send the new trace to the reviewer; the same trace can never be recorded as
   `pass` after a `block`.
5. `$LA report <ws>` → `deliver/report.html`, `deliver/report.md`, `deliver/twin.json` and `deliver/la-export.json`
   (what longpi imports; `$LA export <ws>` checks it and prints its counts).
6. `$LA validate <ws>` must print `ok: true`.

## Comparing with a later visit
Prefer `twin series` (above): it sees slow drifts that every single step hides. `$LA twin compare --prev
<old>/work/twin/twin.json --cur <new>/work/twin/twin.json` remains for two snapshots and is the two-point case of
series (same output as before; units declared convertible are converted). `deliver/twin.json` has no germline
calls, so use the `work/` copies when the fingerprint check matters. Lab changes are judged against
the reference change value from published within-person biological variation (log-normal markers use the
asymmetric form). Readouts without a noise model are shown side by side as `not_judged`. Never call them
better or worse.

## Labs-only retest (fast route)
When `init --reuse-germline-from <earlier ws>` reused the germline results and the new folder holds only labs and
wearables, `observe` prints `fast_route`: skip pipelines with the reason, run `twin series --snapshots ... --ws <ws>`
early (before `twin build` it reads the confirmed lab rows and prints, without writing), and redo organ estimates
and board questions only for systems whose markers have `summary_verdict` `trend_beyond_noise` or
`change_beyond_noise`. That is the definition of "a system with change": series verdicts beyond noise, not values
outside a reference range. Skip the other organs with the reason. `observe` also prints each stage's duration
(`timing`) for cost calibration.

## Output Format
1. Paths of `deliver/report.html`, `deliver/report.md` and `deliver/twin.json`.
2. The trace result and the reviewer verdict with its P0/P1 findings and how each was fixed.
3. What the member report says was not computed, in one line each.

## Guardrails
- After any edit: re-register what changed (integrate / organ), then `intervene register`, `twin build`, `review trace`,
  in that order; `observe` shows what went back to todo.
- Save each reviewer reply verbatim as `work/review/reviewer-N.json`; never abridge it.
- A `pass` with P2 findings may be recorded as is. If you fix P2s anyway, that is a new trace and a new review.
- Citing an organ AI estimate in the summary (`{{r:organ.<organ>.…}}`) renders it with "AI 估计" and its range; do not
  restate it as a plain number in your own words.
- Never render with a failing trace, and never render over a `block` verdict.
- The reviewer is a separate subagent with no stake in the report. Do not review your own writing and call it
  independent.
