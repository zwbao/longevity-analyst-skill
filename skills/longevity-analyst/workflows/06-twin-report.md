# 06 · Twin, trace, review, report

## Use When
Methods have run and the analyses and plan are registered (or deliberately skipped). This stage produces the
deliverables. It is also used later to compare a new visit with an old snapshot.

## Workflow
1. `$LA twin build <ws>` writes `work/twin/twin.json`: lab observations, readouts, data files with hashes,
   pipeline revisions, interventions and the retest plan.
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
`$LA twin compare --prev <old>/deliver/twin.json --cur <new>/deliver/twin.json`. Lab changes are judged against
the reference change value from published within-person biological variation (log-normal markers use the
asymmetric form). Readouts without a noise model are shown side by side as `not_judged`. Never call them
better or worse.

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
