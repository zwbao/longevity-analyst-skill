# Contract: researcher subagent (question board)

**Scope (fixed):** one question about one member. Inputs: the question (id, title, hypothesis, basis), and read-only
access to `work/insights/*.json`, `work/readouts.json`, `work/organs/organ_readouts.json`, `state.json` (labs,
answers). Output: exactly one file, `work/insights/board/Q<n>.json`. You may run these retrieval commands (their
results are cached and become citable): `la.py evidence pubmed <ws> --query "…"`, `la.py insights mr <ws> --exposure …
--outcome …`. You may not edit any other file or register anything.

**Stance:** you do not know the member or the company's preferences. You are a skeptical physician-scientist: look
for the evidence that would make the hypothesis wrong as hard as for the evidence that supports it.

## How to investigate
1. State what would have to be true for the hypothesis to hold, in terms of this member's data and public evidence.
2. Check the member's data: lab values and their printed ranges, genotypes (`gt:<rsid>` / `gt:<hgvs>` in
   `genotype_phenotype.json`), genetic percentiles and their coverage, population positions, readouts, answers.
3. Check public evidence retrieved in this run: GWAS records (`gwas:…`), ClinVar (`clinvar:…`), MR (`mr:…`) and
   projections (`proj:…`), PubMed (`pmid:<n>` from `la.py evidence pubmed`; read the abstract in the cache).
4. Decide: `supported` (the member's data and retrieved evidence agree with the hypothesis), `not_supported` (they
   contradict it), `insufficient` (the data cannot decide). One person's data never gives more than `moderate`
   confidence; a common-variant score alone is at most `low`.

## Return file (JSON, exactly these fields)
```json
{"id": "Q1", "verdict": "supported | not_supported | insufficient", "confidence": "low | moderate",
 "member_evidence": ["<readout id | lab marker | member.answers.<k> | gen:<analyte> | gt:<rsid or hgvs>>"],
 "public_evidence": ["<gwas:… | clinvar:… | mr:… | proj:… | myvariant:… | pmid:<n>>"],
 "summary_zh": "<2–4 sentences: what the member's data shows, what the public evidence says, the verdict>",
 "next_step_zh": "<one concrete next step: a test, a referral, a retest, or none>",
 "limitations_zh": "<what could make this wrong>"}
```
Text follows the report's number rules: member values only as `{{r:<id>}}`, other numbers as `{{n:…}}`, no bare
digits, no doses, no drug advice. `supported` and `not_supported` need at least one retrieved public record.

## Return to the dispatcher (exactly these fields)
`{"id": "Q<n>", "file": "…", "verdict": "…", "public_evidence_count": <int>}`

## Citing the member's own numbers
- A number the member told us (pack years, a parent's age at a heart attack, their age): `{{r:member.answers.<key>}}`
  or `{{r:member.age}}`, where `<key>` is the answer's key in `la.py member`. Never `{{n:…}}` for these.
- Lab values have no placeholder: say the direction against the printed reference range (偏高 / 在范围内).
- `{{r:…}}` takes a readout id (letters, digits, `._-`), never a lab name or a value; anything else stays as raw
  braces and is refused at registration and at trace.
