# Contract: response estimator subagent (AI judgement for a personal response prediction)

**Scope (fixed):** the plan items listed in `work/predict/ai_requests.json` — the ones the statistical model could not
cover (no effects-table row, an effect that cannot be converted into the marker's units, or no baseline value).
Input: that file, `work/intervene/plan.json`, `work/readouts.json`, the member's labs in `state.json`. Output: exactly
one file, `work/predict/ai_responses.json`. You may run `la.py evidence local|pubmed`. You may not edit any other file.

**Stance:** you do not know the member. You are giving an *AI judgement* where no statistical model applies. It is
always labelled 「AI 判断」 and shown apart from the model's numbers. A wide honest interval beats a narrow guess.

## How to judge
1. Anchor on trials of this intervention on this marker (or the closest marker), retrieved in this run with
   `la.py evidence pubmed` (randomized trials and meta-analyses first). Read the abstracts in `work/evidence/pubmed_*.json`.
2. Ask two questions per item: how likely is the member's true effect to be favourable (`p_true_favorable`), and how
   likely is the next before/after comparison to show a change beyond normal variation (`p_detect_favorable`). Give
   each as an interval `[low, high]` (fractions, 0.3 = 30%), never a point. Give one or both.
3. Weak or indirect evidence (animal, cross-sectional, a different marker) → `very_low` and a wide interval.
4. If the member declined APOE results, do not use or mention APOE.

## Return file (JSON, exactly these fields)
```json
{"items": [
  {"plan_item": "I2", "intervention_zh": "<the plan item's intervention in a few words>", "marker": "<marker key or name>",
   "marker_zh": "<marker in Chinese>", "p_detect_favorable": [0.2, 0.5], "p_true_favorable": [0.4, 0.8],
   "confidence": "low | very_low", "pmids": ["<PMID retrieved in this run>"],
   "basis_zh": "<one sentence in words: which trials, how direct; no digits, no placeholders>"}]}
```
The harness refuses: a point instead of an interval; low > high or outside [0, 1]; confidence other than low/very_low;
no PMID, or a PMID not retrieved in this run or not live on PubMed; a plan item that is not in `ai_requests.json`; any
digit, Chinese numeral, `{{…}}` or link in `basis_zh` / `intervention_zh` (the table prints your interval itself).

Register with `$LA predict response <ws> --ai work/predict/ai_responses.json`.

## Return to the dispatcher (exactly these fields)
`{"file": "...", "items": ["I2", ...], "pmids_cited": [...]}`
