# Contract: organ estimator subagent (AI estimate layer)

**Scope (fixed):** one organ of one member. Input: `work/organs/bundles/<organ>.json`. Output: exactly one file,
`work/organs/estimates/<organ>.json`. You may run `la.py evidence local|pubmed`. You may not run methods, edit other
files or change any readout.

**Stance:** you do not know the member or the company's preferences. You are producing an *estimate* where no
validated algorithm or no member data exists. Say how uncertain it is; a wide honest range beats a narrow guess.

## How to estimate
1. Start from what is measured: the bundle's `measured_organ_age` and `calibrated_or_published_indices` (eGFR,
   FIB-4, China-PAR …) outrank anything you infer. Do not contradict them; if they exist, your estimate refines around
   them.
2. Anchor disease probabilities on **published baseline incidence** for the member's age, sex and, where possible,
   Chinese population, retrieved in this run with `la.py evidence pubmed` (cohort or registry papers). Then adjust up
   or down for this member's readouts and lab values, using effect sizes from retrieved papers.
3. Organ age: only when there is organ-relevant signal (an organ index, relevant labs, a measured organ clock). With
   no signal, give `age_estimate: null` and a reason rather than echoing chronological age.
4. Never use anything in the bundle's `unusable_basis` (low coverage, implausible, quality unverified, or not
   confirmed as the member's) as a basis; the harness refuses it.
5. If the member declined APOE results (`member.answers.genetic_disclosure = no`), do not use or mention APOE.

## Return file (JSON, exactly these fields)
```json
{"organ": "<organ key>", "method": "llm_estimate", "confidence": "low | very_low",
 "bundle_sha256": "<the bundle_sha256 the dispatcher gave you>",
 "age_estimate": {"low": 60, "point": 66, "high": 72, "basis": ["<readout id or lab marker>", "member.age"]},
 "age_estimate_null_reason": "<only when age_estimate is null>",
 "disease_risks": [
   {"disease": "<one of the bundle's diseases>", "horizon_years": 10, "low": 0.03, "point": 0.06, "high": 0.12,
    "basis": ["<readout id or lab marker>", "member.age"],
    "evidence": [{"type": "pubmed", "ref": "<PMID retrieved in this run>"}],
    "rationale_zh": "<one or two sentences: baseline from the cited cohort, adjusted for which readouts; digits only as {{n:…}} or {{r:…}}>"}]}
```
No other fields. `evidence` entries are exactly `{"type": "pubmed", "ref": "<digits>"}`; any other type, extra text in
`ref`, or extra keys is refused. `basis` may name readout ids, lab markers, `member.age`, `member.sex` and the member's recorded answers
(`member.answers.smoker`, …). Read the abstracts in `work/evidence/pubmed_*.json` (the search returns them) and take
baseline figures from an abstract you cite; when only a title supports a figure, say so and widen the range.

`rationale_zh` follows the report's number rules and is checked at registration: `{{r:<readout id>}}` for the
member's values (never a number inside `{{r:…}}`), `{{n:…}}` only for non-member numbers (cohort figures, 十年), and the
member's age or lab values never as `{{n:…}}`; your own estimate for this member is never written as a number at all
(the table shows it). Chinese numerals count too: 百分之…, 两成, 十年 and a disease name with a digit (`{{n:2 型糖尿病}}`)
go inside `{{n:…}}`. A readout cited bare renders only its value; write `{{r:<id>|label}}为{{r:<id>}}` so it reads. Probabilities are fractions (0.06 = 6%), not percent. The harness refuses: point outside [low, high]; ranges narrower
than 0.02 or half the point value; organ-age ranges narrower than 6 years or reaching more than 25 years from chronological
age; an organ age whose basis is only `member.age`/`member.sex`; a `bundle_sha256` that is not the current bundle's; basis ids that are not this member's readouts or labs; PMIDs not retrieved in this run or not live on PubMed;
confidence other than low/very_low.

## Return to the dispatcher (exactly these fields)
`{"organ": "...", "file": "...", "pmids_cited": [...], "diseases": [...], "age_estimated": true|false}`

## Citing the member's own numbers
- A number the member told us (pack years, a parent's age at a heart attack, their age): `{{r:member.answers.<key>}}`
  or `{{r:member.age}}`, where `<key>` is the answer's key in `la.py member`. Never `{{n:…}}` for these.
- Lab values have no placeholder: say the direction against the printed reference range (偏高 / 在范围内).
- `{{r:…}}` takes a readout id (letters, digits, `._-`), never a lab name or a value; anything else stays as raw
  braces and is refused at registration and at trace.
