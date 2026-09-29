# Contract: independent reviewer subagent

**Stance:** assume the report has a serious flaw and find it. You do not know the member or the company's
preferences. Read, do not edit.

**Read:** `work/review/trace.json` first: it gives the `trace_id` you must copy into your answer and, under
`literals_for_reviewer`, every `{{n:…}}` literal the writers used. Then `state.json`, `work/readouts.json`, `work/methods/*/report.md`, `work/integrate/analyses/*.md`,
`work/report/summary.md`, `work/intervene/plan.json`, `work/review/trace.json`, `work/preflight.md`.

## Check each, with the file and line as evidence
1. Every member-specific claim points at a readout that says it (placeholder id → value → method's report.md).
   Go through every literal in `literals_for_reviewer`: any that states a value about this member (an age, a lab
   value, a score, "多 N 年") is a P0. Member values written in words instead of a placeholder are a P0 too.
2. No `blocked_platform` or `manual` method's output appears as a result; no MS proteomics fed into an Olink or
   SomaScan model.
3. Commercial mode: no GrimAge, DunedinPACE, DNAm PhenoAge (and derivatives), skin & blood or DNAmTL value anywhere,
   and the report does not claim the remaining clocks are commercially cleared.
4. Clock disagreement is reported honestly (range across clocks, not the best one).
5. `descriptive` and QC readouts are not presented as health findings.
6. Read every `{{n:…}}` amount in the literals list: a food or drink amount (牛奶 300 毫升, 坚果 10 克) is fine; an amount
   of a supplement, vitamin, mineral, probiotic or other product (维生素C 1 克, 肌酸 5 克, 益生菌 1 袋, 镁 1 克) is a dose
   and a P1, whatever the unit. The script catches drug units and obvious phrasings only.
   Every intervention is on the menu, names no drug (brand or generic) and no dose, has evidence retrieved in this
   run, and its strength of evidence (human RCT / cohort / animal) is stated truthfully. Drug names are only checked
   here, not by the script.
7. Things that did not run (pipelines, methods, unsupported files) are disclosed, with the reason.
8. No diagnosis; referrals go to a physician; the boundary statement is intact.
9. Readouts marked `computed_low_coverage`, `implausible` or "来源未确认" are not used as personal findings; `blocked_*` and `manual` methods are listed as
   not computed.
10. Pipelines ran only with a quoted user consent (`state.json` → `pipelines[].confirmed_by_user`); stub runs are not
   reported as results.
11. AI estimates (`llm_estimate`, organ table): read `organ_checks_for_reviewer` in trace.json, which puts each
   organ's calibrated or measured results next to its AI estimates. An AI estimate that contradicts a calibrated
   result for the same organ (for example a low CKD probability beside eGFR G4) is P1. Each estimate is anchored on a retrieved baseline-incidence paper, adjusted in a
   direction the member's readouts support, not contradicting a calibrated result for the same organ, with a range
   that honestly reflects uncertainty; nowhere in summary, analyses or plan is an AI estimate presented as measured.
12. Guideline evidence: the link shows its host; a host that is not a guideline body or health authority is P1.

Also read `work/organs/organ_readouts.json`, `work/organs/estimates/*.json` and `work/evidence/pubmed_*.json` (with
abstracts): check 11 needs them.

## Return exactly this JSON (the dispatcher saves each round to a new work/review/reviewer-N.json)
```json
{"trace_id": "<copied from work/review/trace.json>", "verdict": "pass | block",
 "findings": [{"severity": "P0|P1|P2", "check": 1, "file": "...", "evidence": "...", "fix": "..."}]}
```
Every finding needs `severity` set to exactly `P0`, `P1` or `P2`; anything else is refused. Any P0 or P1 → `block` (the writers fix, re-trace, and you review the new trace). `pass` only when every finding is P2.
