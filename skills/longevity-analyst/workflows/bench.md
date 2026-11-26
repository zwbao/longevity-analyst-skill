# Bench · Capability A eval adapter (LongevityBench)

> **非产品主路径。** 产品短问答由 **longevity-coach** 在对话里答（见教练 `references/short-qa.md`）；会员整包数据走 01–06。本文件只服务离线 LongevityBench 刷分 / 对照。

## Use when

The input is one LongevityBench prompt for **batch scoring**, not a member's folder and not a live coach chat. Do not start at workflow 01. Do not write `deliver/report.html` or `twin.json`. Exam scores stay out of the member report.

For schema policy without an exam id, use [`capabilities.md`](capabilities.md) and:

```bash
python3 <this skill dir>/scripts/la.py schema-route --family clinical_age --format pairwise
```

This workflow is an **evaluation adapter** that maps a LongevityBench prompt onto the same schema policy. It is not a product path that downloads Longevity-LLM. Default routes are `formula` | `host` | `uncovered` (see `la.py bench` JSON `path`). The `_TASKS` / `LB-xxxx` table in `bench_route.py` is **eval metadata** for scoring tables—not the Skill capability definition.

## Command

```bash
python3 <this skill dir>/scripts/la.py bench --prompt <file> [--lb-id LB-0034] [--follow-up <text>]
```

- Formula tasks fill `answer` when weights are present.
- Host tasks return `path: host` and `host_protocol: <format>`; the **Host LLM** produces the exam field. Validate with [`host-protocol.md`](host-protocol.md) / `la.py protocol-check`.
- Optional对照 only:

```bash
python3 <this skill dir>/scripts/la.py bench --prompt <file> --lb-id LB-0006 --control-model [--model-path <dir>]
```

`--control-model` writes `control_answer` from a paper Longevity-LLM. It never becomes the default `answer`. Bulk control scoring lives under `week_6/sample/` (see `week_6/eval/README.md`).

`--lb-id` is an eval-table hint. Routing is by family/format and whether a formula can run, not by "task id → download model".

## Default path table

| Task | Id | Default path | Target / note | Published reference |
| --- | --- | --- | --- | --- |
| Clinical age, pairwise | LB-0034 | formula | `clinical_hgb` | 0.81 |
| Clinical age, multiclass | LB-0030 | formula | `clinical_elasticnet` | 0.43 |
| Clinical age, regression | LB-0038 | formula | `clinical_elasticnet` | 11.55 MAE |
| Clinical mortality, binary | LB-0042 | formula | `clinical_cox` | 0.81 |
| Clinical mortality, pairwise | LB-0050 | Host+tools | host; control: L-Qwen3.5-9B | 0.684 |
| Clinical mortality, multiclass | LB-0046 | Host+tools | host; control: L-Qwen3.5-9B | 0.389 |
| Methylation, pairwise | LB-0006 | Host+tools | host; fixed clocks in do_not_call | 0.868 |
| Methylation, regression | LB-0010 | Host+tools | host | 12.96 MAE |
| Methylation, multiclass | LB-0002 | Host+tools | host | see paper_best |
| Transcriptome, pairwise | LB-0018 / LB-0220 | Host+tools | host | 0.676 |
| Transcriptome, multiclass | LB-0014 | Host+tools | host | 0.223 |
| Plasma protein, pairwise | LB-0299 | formula | `protein_clock` | 0.75 |
| Plasma protein, multiclass | LB-0074 | formula | `protein_clock` | 0.56 |
| Plasma protein, regression | LB-0082 | formula | `protein_elasticnet` | 3.94 MAE |
| OpenGenes, masked name | LB-0086 | Host+tools | host; control: L-LFM2-1.2B | 0.780 |
| SynergyAge, pairwise | LB-0122 | Host+tools | host; do not query SynergyAge | 0.846 |
| SynergyAge, regression | LB-0126 | Host+tools | host; do not query SynergyAge | 11.96 MAE |

Published reference numbers are对照 baselines, not Skill pass criteria. Missing formula weights → `path: uncovered`.

## Formula assets

- Clinical: `data/clinical_models.joblib`, `data/clinical_cox.json` (NHANES 2003–2017 train; 1999/2001 eval).
- Protein clock / elasticnet: `data/protein_clock.json`, `data/protein_elasticnet.json` (see `week_6/sample/data/olink/SOURCE.txt`).
- Aggregate formula scores: `week_6/sample/data/scores_formula.json`.

## Output

One exam field is added only after a formula caller sets `callable: true`, or after the Host supplies a field that passes `protocol-check`. A member-facing age still comes from the fixed-feature clocks in workflow 03.
