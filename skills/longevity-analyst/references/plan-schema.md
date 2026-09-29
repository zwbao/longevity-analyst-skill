# plan.json schema

```json
{
  "member": "<id>",
  "items": [
    {
      "id": "I1",
      "category": "diet | supplement | exercise | sleep | lifestyle | test | referral",
      "action_zh": "what to do, in plain Chinese, no dose, no drug",
      "rationale_zh": "why, citing readouts as {{r:<id>}}; other numbers only as {{n:…}}",
      "targets": ["<readout id or lab marker exactly as in the workspace>"],
      "evidence": [
        {"type": "effects", "ref": "<id in longevity-skills data/effects.jsonl>"},
        {"type": "claims", "ref": "<id in claims.jsonl>"},
        {"type": "pubmed", "ref": "<PMID retrieved in this run>"},
        {"type": "guideline", "ref": "https://..."}
      ],
      "evidence_grade": "human_rct | human_cohort | animal_or_cell | guideline | not_applicable (test/referral only)",
      "executor": "nutritionist | physician | member",
      "retest": {"what": "<marker or readout>", "after_weeks": 12}
    }
  ]
}
```

`retest.what` is printed in the report and traced like prose: write a name with digits as `{{n:糖化血红蛋白(HbA1c)}}`
and never a file id there (a `file:F00x` id belongs in `targets` of a test/referral item only).
No fields beyond those shown; `retest` has only `what` and `after_weeks`; each evidence entry is exactly `{type, ref}`.
Guideline URLs are plain page addresses (no query string or fragment) fetched with `la.py evidence fetch`.
Every reviewer finding needs `"severity": "P0" | "P1" | "P2"`.

`test` and `referral` items may have empty `evidence`. All other items need at least one retrieved reference.
