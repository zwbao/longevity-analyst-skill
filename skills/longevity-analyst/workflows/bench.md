# Bench · LongevityBench route

## Use when

The input is one LongevityBench prompt, not a member's folder. Do not start at workflow 01. Do not write `deliver/report.html` or `twin.json`. Exam scores stay out of the member report.

## Command

```bash
python3 <this skill dir>/scripts/la.py bench --prompt <file> [--lb-id LB-0034] [--follow-up <text>] [--run-model] [--model-path <dir>]
```

`--run-model` loads the local Longevity-LLM named in the route. Without it, a model route prints the path and stays `callable: false`.

## Published best (Table S6)

Parity uses the better of ML and Best LLM. Pairwise = concordance, binary/multiclass = balanced accuracy, regression = MAE.

| Task | Id | Call | Published best | Local score |
| --- | --- | --- | --- | --- |
| Clinical age, pairwise | LB-0034 | `clinical_hgb` | 0.81 | 0.811 |
| Clinical age, multiclass | LB-0030 | `clinical_elasticnet` | 0.43 | 0.452 |
| Clinical age, regression | LB-0038 | `clinical_elasticnet` | 11.55 | 11.97 |
| Clinical mortality, binary | LB-0042 | `clinical_cox` | 0.81 | 0.807 |
| Clinical mortality, pairwise | LB-0050 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.684 | scoring via `run_step3_scores.sh` |
| Clinical mortality, multiclass | LB-0046 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.389 | same |
| Methylation, pairwise | LB-0006 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.868 | stride-5 interim then full |
| Methylation, regression | LB-0010 | `L-Qwen3.5-9B` (MLX 4-bit) | 12.96 | same |
| Methylation, multiclass | LB-0002 | `L-Qwen3.5-9B` (MLX 4-bit) | 9B Table S3 (task-best grok-4 0.203) | same |
| Transcriptome, pairwise | LB-0018 / LB-0220 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.676 | 0.624 (n=900) |
| Transcriptome, multiclass | LB-0014 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.223 | scoring next |
| Plasma protein, pairwise | LB-0299 | `protein_clock` | 0.75 | 0.752 |
| Plasma protein, multiclass | LB-0074 | `protein_clock` | 0.56 | 0.636 |
| Plasma protein, regression | LB-0082 | `protein_elasticnet` | 3.94 | 4.37 |
| OpenGenes, masked name | LB-0086 | `L-LFM2-1.2B` | 0.780 | 0.765 |
| SynergyAge, pairwise | LB-0122 | `L-Qwen3.5-9B` (MLX 4-bit) | 0.846 | not scored locally yet (S3: 9B 0.731) |
| SynergyAge, regression | LB-0126 | `L-Qwen3.5-9B` (MLX 4-bit) | 11.96 | not scored locally yet (S3: 9B 12.38) |

Local Cox control on mortality pairwise / multiclass: concordance 0.631 and balanced_acc 0.329. Those are not the parity targets.

## Open models

| Paper name | Repository | Local path |
| --- | --- | --- |
| L-Qwen3.5-9B | `insilicomedicine/longevity-llm` | `week_6/sample/data/hf/longevity-llm` (+ MLX 4-bit at `longevity-llm-mlx-4bit`) |
| L-LFM2-1.2B | `LiquidAI/LFM2-1.2B-Longevity` | `week_6/sample/data/hf/LFM2-1.2B-Longevity` |

## Protein

- Pairwise / multiclass (`protein_clock`): `fit_protein.py --goeminne <Conventional_coefs.csv>` loads OrganAge Conventional fold-0, then affine-calibrates on Allen train (eval kitGuids dropped). Weights: `data/protein_clock.json`.
- Regression (`protein_elasticnet`): `fit_protein.py <olink_wide.csv>` fits ElasticNetCV on ~392 consensus proteins, Allen Imm-of-Aging after dropping 42 eval kitGuids (train n=178). **Not** the paper's official Olink train split. Weights: `data/protein_elasticnet.json`.
- Details: `week_6/sample/data/olink/SOURCE.txt`. Aggregate scores: `week_6/sample/data/scores_formula.json`.

## Clinical features

Prompt union mask (LB-0034/0030/0038/0042) is `data/clinical_feature_mask.json` (29 labs = `clinical_age._FEATURES`). HGB / elasticnet / Cox train on NHANES 2003–2017 only; 1999/2001 are eval.

## Output

One field is added only after `callable` is true. A member-facing age still comes from the fixed-feature clocks in workflow 03.
