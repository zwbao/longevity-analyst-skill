# 05 · Intervention plan

## Use When
The system analyses are registered.

## Workflow
1. Read `data/intervention_menu.json`, or the nutrition team's own menu if they provided one. The plan may use only
   its categories and each category's allowed executors.
2. For each finding worth acting on, retrieve evidence:
   - `$LA evidence local <ws> --terms "DunedinPACE,甘油三酯,Akkermansia"` → trial effects (`effects`, with quotes)
     and paper claims (`claims`, with human/animal/cell labels).
   - `$LA evidence pubmed <ws> --query "<intervention> AND <marker> AND randomized controlled trial[pt]"` (live).
   - A guideline page you want to cite: `$LA evidence fetch <ws> --url https://...` first (records that it was retrieved).
   If PubMed fails, stop citing new PubMed evidence and say so. Do not fill in from memory. At registration every
   cited PMID is re-checked live against PubMed; a PMID that PubMed does not return is rejected.
3. A `test` or `referral` item may target a file as `file:F005` (e.g. "ask the lab which assay produced F005").
   Write `work/intervene/plan.json` following [`../references/plan-schema.md`](../references/plan-schema.md). Each item
   has a target (readout id, organ estimate id such as `organ.liver.risk.1`, or lab marker), evidence refs retrieved in this run, an executor and a retest.
   Human RCT evidence ranks above cohort, and cohort ranks above animal or cell. An item backed only by animal
   evidence is phrased as "under study", or left out.
4. `$LA intervene register <ws> --plan work/intervene/plan.json`. It rejects unknown categories, disallowed executors,
   doses in the text, targets that are not this member's data, and evidence not retrieved in this run.

## Output Format
1. The registered plan: item id, category, action, target, executor, evidence grade, retest.
2. Findings deliberately left without an intervention, and why (weak evidence, needs a physician).
3. Evidence searches that returned nothing or failed.

## Guardrails
- An item that targets an AI estimate (`organ.*`) says in `rationale_zh` that the target is an estimate; the plan never
  treats it as a measured risk.
- No drugs and no doses. The script rejects doses and pill counts in `action_zh` and `rationale_zh` (digits or
  Chinese numerals); naming a drug is checked by the independent reviewer. A finding that needs medical follow-up (APOE ε4/ε4, ACMG-reportable variant,
  high China-PAR risk) becomes a `referral` or `test` item for the physician.
- Fewer, well-supported items beat a long list. The nutritionist reviews the plan before it reaches the member.
