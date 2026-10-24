---
name: longevity-analyst
description: >
  Turns one person's longevity / anti-aging test data (WGS or VCF, DNA
  methylation, gut metagenome, plasma proteomics, lab reports and PDFs) into a
  multi-omics report, a digital-twin snapshot and an evidence-cited
  intervention plan, running local nf-core pipelines after a preflight check.
  Use when a user hands over a folder of health or omics test results and asks
  for biological age, an aging or longevity analysis, a digital twin, or a
  personalised intervention plan; also when a new retest must be compared with
  an earlier snapshot. Triggers on 抗衰分析, 长寿检测, 生物学年龄, 甲基化年龄,
  数字孪生, 多组学报告, 干预方案, 复测对比, longevity report, biological age,
  epigenetic clock, digital twin. Do NOT use for cancer treatment decisions,
  diagnosing disease, or reading a single lab value. When the longevity-coach
  skill is installed it is the entry point for personal conversations and calls
  this skill for the analysis; inside LongPi, start only after run_deep_analysis.
license: MIT
category: health-omics
metadata:
  author: CancerDAO / zwbao
  version: "0.8.0"
  tags: longevity aging multiomics digital-twin epigenetic-clock nf-core bioinformatics
---

# longevity-analyst

A member's raw test folder goes in; `deliver/report.html`, `deliver/report.md`
and `deliver/twin.json` come out. The work runs through one harness,
`scripts/la.py`. The harness owns state, file sniffing, pipeline commands,
method dispatch, gates and number binding. You own the judgment calls: what
a file is, what a PDF says, what the results mean for this person, and which
interventions to propose.

## Setup (once per machine)

```bash
git clone https://github.com/zwbao/longevity-skills   # method library, 172 published methods
pip install pandas numpy markdown                      # epiage clocks and the HTML report need them
```

`$LA` in the workflows stands for `python3 <this skill dir>/scripts/la.py`. Shell variables do not survive
between your tool calls, and zsh does not split a multi-word variable, so write the full command every time (or
define a shell function `la() { python3 <dir>/scripts/la.py "$@"; }` inside the same call). `init --longevity-skills <clone>` stores the
library path in the workspace, so later commands find it without any environment variable.

## Running inside dsh / longpi
- In longpi the AI (you, in a health chat) decides to run: `run_deep_analysis` with the reason, when the snapshot line
  深度分析 says new data is waiting, or when the member asks. It gives you the data folder, the workspace, the
  member id, the member's age/sex and the method library (`~/longpi/longevity-skills`). Use exactly those; the member
  id stays the same across runs so longpi can compare a retest with `twin compare`. When the report is done, call
  longpi's `import_analysis` yourself and read the plan back for the member to adopt. Items whose executor is the
  physician, and supplements, tests and referrals, go into longpi's doctor brief, not the member's plan.
- When it says the member's Mirobody is connected, first `$LA mirobody pull <data folder> --mcp-url-file <file it gives>`
  (never paste the URL into a command: it is the member's secret). Labs land as one CSV per checkup date and watch
  days as one CSV; they go through intake like any other file (workflows/01-intake.md, step 0).
- `$LA report` also writes `deliver/la-export.json` (schema `la-export/1`), which `import_analysis` reads. The plan
  becomes a longpi plan only after the member hears it read back and confirms.

## Routing Rules

1. New data folder → start at workflow 01 and walk 01 → 02 → 03 → 04 → 04b → 04c → 05 → 06 in order. Never skip
   a stage: `$LA observe <ws>` names the next one.
2. Re-analysis after a retest → run 01–06 on the new folder in a new workspace (`init --reuse-germline-from <old ws>`
   when the germline data are the same; pass the host's exposures with `twin build --exposures`), then
   `$LA twin series --snapshots <every earlier work/twin/twin.json> --ws <new ws>` before review (workflow 06).
   `twin compare` is only the two-snapshot case.
3. The user only asks "can my machine run this?" → workflow 02.
4. Before every step run `$LA observe <ws>` and act on its `next` line.
   It is the durable state; your memory of earlier turns is not.
5. Open only the workflow for the current stage.

## Workflow Index

| Stage | Workflow |
|---|---|
| Set up the workspace, inventory files, resolve what each file is, transcribe PDFs | [`workflows/01-intake.md`](workflows/01-intake.md) |
| Check the machine, tell the user the cost, run pipelines only after consent | [`workflows/02-preflight-pipelines.md`](workflows/02-preflight-pipelines.md) |
| Run the published methods and native modules through the platform and licence gates | [`workflows/03-methods.md`](workflows/03-methods.md) |
| One analyst subagent per body system, cross-modality synthesis | [`workflows/04-integrate.md`](workflows/04-integrate.md) |
| Organ checkup table: organ ages, organ indices, AI estimates of organ age and disease probability | [`workflows/04b-organs.md`](workflows/04b-organs.md) |
| Insights: genes vs labs, population position, MR projections, the question board (one researcher per question) | [`workflows/04c-insights.md`](workflows/04c-insights.md) |
| Evidence retrieval and the intervention plan | [`workflows/05-intervene.md`](workflows/05-intervene.md) |
| Twin snapshot, number trace, independent review, render | [`workflows/06-twin-report.md`](workflows/06-twin-report.md) |

## Evidence Contract

- **Sources of record:** the member's files; `work/readouts.json` (every number
  a method computed, with its source file); the longevity-skills library
  (`catalog.json`, `data/effects.jsonl`, `claims.jsonl`,
  `biological_variation.json`); live PubMed via `$LA evidence pubmed`; live GWAS Catalog, ClinVar/gnomAD
  (myvariant.info), Ensembl and EpiGraphDB MR via `$LA insights`; population references in
  `data/ref_population.json` (NHANES 2005–2010, GMHI 4347 metagenomes).
- **Fallback:** a missing input is named and reported as not computed. A
  pipeline that cannot run is listed with its preflight reason. PubMed down
  means "evidence not retrieved", not remembered citations.
- **Never fabricate:** measured lab values, computed clock values, PMIDs, DOIs,
  effect sizes, guideline thresholds. In prose, cite member numbers only as
  `{{r:<readout id>}}`; `$LA review trace` refuses any other number.
- **Doses:** a dose comes only from the verbatim text of a registered source (a rule in
  `data/medication_rules.json`, CPIC text in `work/pgx/pgx.json`) or is an AI judgement labelled as such. It reaches
  prose only inside a registered item cited as `{{rx:<id>}}`; a dose written anywhere else is refused.
- **Labelled items:** predictions (`$LA predict response`), diagnostic conclusions (`$LA clinical dx`) and medication
  suggestions (`$LA clinical rx`) are registered items, each with a `la-label/1` label: its category (测量, 已发表模型,
  按《指南》的标准, 《指南》推荐, CPIC 指南, 统计模型估计, AI 判断, rendered from `basis_kind`, never written by you) and
  its sources. Prose cites them as `{{p:<id>}}`, `{{dx:<id>}}`, `{{rx:<id>}}`; rendering prints the tag, the sentence
  and the sources. AI judgements carry confidence low/very_low and are shown in their own column.
- **AI estimates (allowed in this skill):** where no validated algorithm or no
  member data exists, organ ages and disease probabilities may be estimated by
  an LLM estimator subagent. Every estimate is a range with its basis (readout
  ids / lab markers / retrieved PMIDs), carries confidence low or very_low, is
  registered as its own readout kind `llm_estimate`, and is shown in its own
  column labelled "AI 估计". It is never presented as a measurement or a
  calibrated model output.

## Hard lines

| Line | Why |
|---|---|
| No pipeline starts without the preflight shown and the user's words quoted in `--confirmed` | a WGS run can take a day and hundreds of GB |
| Never pass MS proteomics into Olink/SomaScan clocks; never override `blocked_platform` | raw MS intensities produce "brain age 140" with no error |
| `--mode commercial` (default) runs only the provisional clock list (no GrimAge, DunedinPACE, DNAm PhenoAge, skin & blood, DNAmTL); counsel has not cleared it yet | licences are listed in `data/license_policy.json` |
| A multi-sample file is never read by position; the member's column is a judgment | another person's column gives that person's results |
| Any edit after `review trace` voids the trace and the review | the report is built only from the files that were checked |
| Diagnoses, medication suggestions and predictions only as registered, labelled items (`clinical dx`, `clinical rx`, `predict response`); an AI judgement is always labelled AI 判断; an emergency means "拨打 120" first | every such sentence must show where it came from and what kind of judgement it is |
| Change workspace state only through `la.py` commands; result files are read-only for you | the trace and the twin depend on them |

## Subagents

Workflow 04 dispatches one analyst per body system, workflow 04b one estimator per organ, workflow 04c one
researcher per board question, and workflow 06 one independent reviewer. Each has a
fixed scope and returns exactly the fields named in its contract: [`references/system-analyst.md`](references/system-analyst.md)
returns `{system, file, readouts_cited, pubmed_cited, open_questions}`, [`references/organ-estimator.md`](references/organ-estimator.md)
returns `{organ, file, pmids_cited, diseases, age_estimated}`, [`references/researcher.md`](references/researcher.md)
returns `{id, file, verdict, public_evidence_count}`, [`references/clinical-estimator.md`](references/clinical-estimator.md)
returns `{files, diagnoses, medications, pmids_cited}`, [`references/response-estimator.md`](references/response-estimator.md)
returns `{file, items, pmids_cited}`, and [`references/reviewer.md`](references/reviewer.md)
returns `{verdict, findings[]}`. Free prose instead of those fields counts as a failed dispatch; send it back.

## Output Requirements

1. `deliver/report.html` + `report.md`: summary, readout table, per-system
   analyses, intervention plan, what was not computed and why, twin and retest
   plan, data provenance, boundary statement.
2. `deliver/twin.json`: the member's state for the next comparison.
3. In chat: where the files are, which stages ran, what did not run and why,
   and what the user must decide next (for example, whether to run WGS locally).
