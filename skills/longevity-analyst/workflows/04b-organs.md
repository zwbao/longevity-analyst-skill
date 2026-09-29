# 04b · Organ checkup table

## Use When
The system analyses are registered and the member (or the plan) needs a per-organ view: organ age and disease risk.

## Inputs Needed
- Required: a workspace past `integrate`.
- The deterministic organ indices (eGFR with KDIGO stage, FIB-4 with tier, TyG, AIP) were already computed in
  `methods run` from lab rows; missing inputs are listed in the method's notes.

## Workflow
1. `$LA organ bundle <ws>` writes one bundle per organ to `work/organs/bundles/<organ>.json` with the measured organ
   ages, the calibrated indices, all readouts and labs, and the diseases to estimate.
2. Dispatch **one estimator subagent per organ, in parallel**, with [`../references/organ-estimator.md`](../references/organ-estimator.md),
   the bundle path, its `bundle_sha256` (from the `organ bundle` output) and the output path, plus: "You do not know
   the member or the company's preferences."
3. `$LA organ register <ws> --organ <organ>` for each returned file. Rejections name what to fix; send them back to
   the same estimator. A rejected re-registration also voids that organ's earlier estimate. Re-running
   `organ bundle` after an upstream change voids estimates and skips made from the old bundle.
4. An organ that should not be estimated (the member declined, no relevant signal and no baseline literature):
   `$LA organ skip <ws> --organ <organ> --reason "<why>"`.
5. When every organ is registered or skipped, the stage is done and the report gets an "器官体检表" section.

## Output Format
1. Per organ: measured/model age, AI age range, indices/calibrated models, AI disease probability ranges with horizon.
2. The PMIDs each estimate rests on.
3. Organs skipped and why.

## Guardrails
- `organ register` applies the report's number rules to `rationale_zh` (a `{{r:…}}` takes a readout id, never a value;
  the member's own age or lab values are never `{{n:…}}`; 十年 is `{{n:十年}}`). Fix and re-register the same organ.
- A new registration voids the plan, twin and review; after any change re-register in order: organ → intervene →
  twin build → review trace. An upstream change (labs, methods, member answers) deletes `work/organs/` entirely.
- If the user does not want to see a disease (for example "不想看老年痴呆"), tell the estimator to leave that disease
  out and say so in the report; for APOE use `member genetic_disclosure=no`.
- PubMed may answer 429 when nine estimators search at once: wait and retry the search; never cite an unretrieved PMID.
- The organ indices use only lab rows confirmed in intake (`labs confirm`). Several confirmed rows with different
  values for one analyte stop that index until the extra rows are fixed or dropped.
- A single eGFR below 60 already meets the stage-3 GFR threshold: do not estimate "慢性肾脏病 3 期及以上" then (the
  harness refuses it and the table says a repeat test is needed).
- AI estimates are a separate kind (`llm_estimate`), always shown in their own column with the note that they are
  neither measurements nor calibrated models. Never copy an AI estimate into a summary as if it were measured.
- A calibrated result (China-PAR, eGFR stage, FIB-4 tier) always outranks an AI estimate for the same organ.
- Plans may target AI estimates, but the reviewer checks that the plan does not overstate them.
