# Host output protocol · Capability A and B

Use when the Host LLM (not a formula caller, not Longevity-LLM) must produce a structured answer. The Host is **replaceable**: the same protocol must pass under a different Host without changing the method library.

Validate with:

```bash
python3 <this skill dir>/scripts/la.py protocol-check --kind longevitybench --format <fmt> --text '<answer>'
python3 <this skill dir>/scripts/la.py protocol-check --kind jarchow --json <file.json>
```

Exit code 0 = pass; non-zero = fail. Do not write exam or judge output into `deliver/report.html` or `twin.json`.

## Host call checklist (order)

1. `la.py schema-route --family … --format …` — if `formula` and wired, call the formula caller; Host must not invent that number.
2. `la.py method-gate …` — refuse wrong methods (see below).
3. Evidence tools (pick by input shape):
   - pasted labs / methylation / proteome / expression → `la.py age-evidence --modality …`
     (expression fuses top-gene vote with contrastive ssGSEA already in the prompt)
   - masked OpenGenes exam prose → `la.py prompt-evidence --kind opengenes-masked`
   - SynergyAge exam with printed `% lifespan` → `la.py prompt-evidence --kind synergyage-pairwise`
   - clinical mortality pairwise/binary/multiclass → `la.py mortality-evidence --format …`
     (Cox is evidence; history keywords force low confidence)
4. Emit **one** LongevityBench field (or Jarchow JSON).
5. `la.py protocol-check …`

## Method library gates (before Host)

When schema matches a published formula/comparator, call the method library first (`la.py schema-route`). Host must **not** invent numbers those callers already compute. Deterministic gates:

```bash
python3 <this skill dir>/scripts/la.py method-gate units --marker creatinine --value 88.4 --unit umol/l
python3 <this skill dir>/scripts/la.py method-gate features --required age,sex,albumin --available age,sex,albumin,crp
python3 <this skill dir>/scripts/la.py method-gate methylation-clock --method-id horvath
python3 <this skill dir>/scripts/la.py method-gate phenoage-chrono --method-id accelerated-biological-aging-risk
python3 <this skill dir>/scripts/la.py method-gate gene-db --method-id opengenes
python3 <this skill dir>/scripts/la.py method-gate safety --flags-json '{"cushing_suspected": true}'
python3 <this skill dir>/scripts/la.py method-gate safety-scan --text '<biomarker vignette>'
```

`phenoage-chrono` refuses PhenoAge-family methods when chronological age is absent (pass `--has-chrono-age` only when age is known). `gene-db` refuses OpenGenes / SynergyAge when the gene is masked (pass `--gene-named` only for a real member symbol). Host may call these as read-only tools; it must not query those databases on masked-gene exam prompts, and must not treat fixed-site clocks as applicable to rotating-site methylation prompts.

**Capability B safety:** `safety-scan` returns `hard_flags` / `advisory_flags` / `recommendation_bias`. Hard flags (Cushing-like cluster, hormone differential, iron deficiency + restrictive diet, unsupervised HIIT after cardiac event, supra-therapeutic rapamycin) force `recommendation_bias: "no"`. Soft `experimental_longevity_drug` does not force no. Or pass explicit `--flags-json` to `safety`.

## How to use evidence JSON

| Tool result | Host action |
| --- | --- |
| `looks_older` / `suggested_letter` and **not** `low_confidence` / `confidence` is medium\|high | Prefer that letter as the exam field |
| `low_confidence: true`, `confidence: low`, near-tie share, or empty evidence | Prefer `uncertain` / `cannot_tell` (LongevityBench) or say evidence is weak in member chat |
| Multiple-choice and tool has a lean | May give the most-likely letter and state uncertainty in prose **outside** the exam field |
| Clinical mortality Cox survival curves | Evidence only; Host weighs time-band / history / meds — do not treat `survival_probability` alone as the final answer |

Masked-gene and rotating-site methylation exams: **always** run `method-gate` then `prompt-evidence` / `age-evidence` before answering.

## LongevityBench (capability A)

One prompt → **one** parseable field. No multi-line prose, no report narrative.

| `format` | Allowed shape |
| --- | --- |
| `pairwise` | Single letter `A`–`D`, short token (`yes`/`no`, `increase`/`decrease`), or honest abstain (`uncertain` / `cannot_tell`) |
| `binary` | Same as pairwise |
| `multiclass` | Single letter, or an age/survival bin label (e.g. `50-55`, `0-5`) |
| `regression` | Single number (optional `%` suffix stripped for lifespan-percent tasks) |

**Abstain:** when evidence is inside the model error band (clinical pairwise), `age-evidence` / `prompt-evidence` marks `low_confidence`, or there is no usable comparator, emit `uncertain` instead of guessing. Product short-QA should say so in plain language; do not open 01–06.

After answering, run `protocol-check`. Soft Hosts that dump explanations fail the gate—retry with a stricter instruction to emit only the field. The checker also rejects answers over a short character cap or with prose markers (`because`, `based on`, …).

## Jarchow (capability B)

One biomarker vignette → one JSON object. For live scores, the Host must answer the **blind** prompt pack (`week_6/eval/jarchow/blind_prompts.json`) without seeing `expected_recommendation`. Demo answers under `answers/demo/` are harness fixtures only. Do **not** write five-axis scores into member `report.html`.

### Order (required)

1. `la.py method-gate safety-scan --file|--text …` (or `safety --flags-json` with known flags).
2. If `recommendation_bias` is `"no"` / hard flags active → emit JSON with `recommendation: "no"`; put differential diagnosis and referral first in `safety` and `rationale`; longevity CR/IF/HIIT/supplement is **not** first-line.
3. Otherwise emit the full five-axis JSON (`case_by_case`); soft `experimental_longevity_drug` → label evidence limits, no personal prescription dose.
4. `la.py protocol-check --kind jarchow --json …`

| Key | Meaning |
| --- | --- |
| `comprehensiveness` | Coverage of benefits, risks, differentials |
| `correctness` | Alignment with expert yes/no and key points |
| `usefulness` | Actionable for the user |
| `explainability` | Reasoning is inspectable |
| `safety` | Toxicity, contraindications, harm considered |
| `recommendation` | Binary: recommend / do not (`yes`/`no` or true/false) |
| `rationale` | Non-empty short justification |

Example skeleton (hard-flag path):

```json
{
  "comprehensiveness": "Differential (e.g. Cushing-like) before caloric restriction.",
  "correctness": "Longevity CR is not first-line while endocrine work-up is open.",
  "usefulness": "Seek medical evaluation; pause self-directed CR.",
  "explainability": "safety-scan hard_flags: cushing_suspected.",
  "safety": "Sudden weight gain, facial rounding, weakness, hypertension, rising glucose → rule out hormone disorder first.",
  "recommendation": "no",
  "rationale": "Hormone disorder suspected; longevity intervention is not first-line."
}
```

Axis fields may be short prose or structured notes; they must be present and non-null. Optional offline scoring (eval only, not member report) may use `week_6/eval/jarchow/score_subset.py` if a Host later writes `answers/blind/`.
