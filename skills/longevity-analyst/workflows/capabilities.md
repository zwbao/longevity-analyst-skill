# Capabilities · product path (not exam routing)

LongPi analyst work is organized by **three transferable capabilities**. Benchmarks are mirrors that score those capabilities through adapters under `week_6/eval/`. Exam task ids (`LB-xxxx`, `JARCHOW-n`) never define the product path.

## Capability map

| Cap | What the analyst does for a person | Primary mirror (eval only) | This round |
| --- | --- | --- | --- |
| **A · Verifiable table reading** | Single-modality table or vignette → parseable age / who-older / bin / 10-year survival | LongevityBench | Accept |
| **B · Bounded intervention advice** | Biomarker vignette → structured advice with safety gates | Jarchow five axes (fixed 10-item subset) | Accept |
| **C · Intervention response reading** | Same-person pre/post → change vs noise; clocks: prefer reliable panel, multi-clock agree, short-term ≠ lifespan | Eval: capability_c rules (epigenetic probe); discovery mining is appendix only | **Narrative thin accept** (coach retest-response); no full-trial library in default path |

## How the analyst chooses a path (schema, not task id)

1. Inspect **input schema** (family of features + question format).
2. If a published formula / comparator in the method library matches → call it (`path: formula`).
3. Else → Host LLM answers under [`host-protocol.md`](host-protocol.md), optionally calling read-only tools (`path: host`).
4. If neither applies → report uncovered (`path: uncovered`).
5. Clinical pairwise formula may return `uncertain` when the predicted-age gap is inside ~2× regression MAE (model error band)—do not force A/B.
6. Transcriptome top-gene / GSEA vignettes: method library has **no** personal transcriptomic-age formula (tier B/C papers only look up tables). Prefer Host `uncertain` over a confident guess.

```bash
# Product-facing: family + format only (no LongevityBench id required)
python3 <this skill dir>/scripts/la.py schema-route --family clinical_age --format pairwise
```

LongevityBench `LB-xxxx` ids appear only in [`bench.md`](bench.md) as an **eval adapter** hint for scoring tables. Do not download Longevity-LLM as the default answer path; `--control-model` is optional对照.

## What must not enter the member report

- LongevityBench / Jarchow exam scores
- Longevity-LLM control answers
- Task ids as if they were clinical findings
- Intervention-response heatmaps, board ranks, or short-term clock deltas presented as proven lifespan gain

Member-facing ages still come from fixed-feature clocks in workflow 03. Capability C may narrate pre/post change vs noise; it must not overwrite that age with a short clock delta.

## Capability C (thin)

When the Host or coach narrates methylation / clock retests: follow `week_6/eval/capability_c/sehgal_rules.md` (preferred clocks, dual consistency, noise band, short-term ≠ lifespan). Do **not** add a default dependency on any full multi-trial clock-response platform.
