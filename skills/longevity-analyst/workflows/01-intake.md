# 01 · Intake

## Use When
A new data folder arrives. Output: a workspace whose `observe` shows no pending judgments.

## Inputs Needed
- Required: the raw folder; the member's age and sex (ask if unknown; never infer from a report without saying so).
- Optional: sample date, lab delivery notes, questionnaire answers.

## Workflow
1. `$LA init <raw_dir> <ws> --member-id <id> --age <n> --sex <m|f> --longevity-skills <clone> [--sample-date YYYY-MM-DD] [--mode research]`.
   `--mode commercial` is the default and is right for paid member work.
2. `$LA intake <ws>`. It sniffs every file and lists `pending_judgments`.
3. Resolve each judgment with evidence you actually read, and write that evidence in `--reason`:
   - `modality` (FASTQ/BAM/CRAM): which assay made the reads. Look for a lab delivery note, an order form or the user's words.
     If nothing says, **ask the user**. A file name like `R1` says nothing about the assay.
     If the user cannot say either, answer `unknown`: the file is parked as `deferred_ask_lab` and listed in the
     report, and nothing runs on it. Never guess an assay to keep the pipeline moving.
     Paired mates (R1/R2) are resolved together: answering for R1 closes R2's judgment too.
     `$LA assign <ws> --file F005 --what modality --value metagenome_shotgun --reason "lab note: 粪便宏基因组"`
   - `platform` (protein matrix): Olink, SomaScan or mass spec. The file header, the lab name, or the user decides it.
     This choice gates which methods may run, so never pick Olink to "get more results".
   - `tissue` (methylation β, IDAT, WGBS reads): whole blood, saliva and so on. Look for the lab's delivery note, the
     order form or the user's words; if none says, ask the user. `unknown` stops every methylation clock, the most
     costly single answer in intake, so say that in the question.
   - `value_scale` (protein matrix): raw intensity, log2 intensity, NPX, RFU or cohort z-scores. Published clocks
     need a specific scale; the platform alone is not enough.
   - `sample_column` (any table or VCF with several samples, or a methylation export whose p-value columns were set
     aside): which column is this member. Use the lab's sample sheet or the user's words. Never pick the first column
     by default: another person's column gives that person's results.
   - `assembly` (VCF without a build in its header): GRCh38 or GRCh37.
   - `primary` (two files of the same kind, e.g. two methylation exports): which one is this visit's measurement:
     `$LA assign <ws> --what primary --file methylation --value F004 --reason "..."`.
   - `transcribe` (PDF/images): follow [`../references/pdf-transcription.md`](../references/pdf-transcription.md),
     then `$LA labs add <ws> --csv <rows.csv> --file F001 --reason "..."`.
     If the document holds no lab values: `$LA assign ... --what transcribe --value no_lab_values --reason "..."`.
4. `identity` (always asked once): do all files belong to this member? Compare the VCF sample ids and header lines
   shown in the hints, and the name/sex/age on transcribed documents, with the member record. A public reference
   sample id (e.g. HG00096), a "synthetic" header or a mismatched sex means `inconsistent`; say which files.
   `$LA assign <ws> --what identity --value consistent|inconsistent|cannot_tell --reason "..." [--files F002,F004] [--uncertain F003]`.
   A file with no identity information at all (a CSV, a MetaPhlAn table, a β matrix) is `consistent` unless something
   contradicts it. A file that declares itself synthetic, demo or a public reference sample (in a header, title or
   footer) is excluded, whatever the user believes; say so in the report. `--uncertain` files contribute no personal
   finding anywhere (analysts, estimators and the twin comparison skip them), so use it only when that is the point.
   `--files` excludes files that are not this member's. `--uncertain` keeps a file whose ownership cannot be confirmed in
   the analysis, but every readout from it is labelled "来源未确认" and must not carry a personal conclusion.
   With `inconsistent`, list the files that are not this member's in `--files`: they are excluded from every
   method, from the twin baseline and from the readout table (still listed in the report as excluded).
   The report prints a warning unless the answer is `consistent`. Decide identity before methods run; a later
   revision (`assign --what identity` again) is allowed, keeps the earlier verdict on record, and voids every later
   stage so nothing computed from an excluded file survives.
5. A lab value read wrongly (e.g. `6.8↑` where the arrow is the out-of-range flag) is corrected with
   `$LA labs fix <ws> --marker "葡萄糖" --value 6.8 --reason "..."`, or dropped with `$LA labs drop`. The old
   value is kept in the row's history.
6. Confirm the lab rows. `$LA labs candidates <ws>` lists every numeric row with its value, unit, reference range and
   the analyte the harness would match it to (empty: no method uses it; answer it anyway) (for example "空腹血糖（静脉血浆/血清；不是餐后、OGTT、随机或尿糖）").
   For each row decide from the row itself, its reference range and the source document: `yes` only if it is from
   this visit (the sample date; an earlier year's table is `no`, it belongs to that year's workspace) and is that
   analyte in serum/plasma (fasting may be inferred from the document's fasting reference range; say so in `why`). Physical measurements
(blood pressure, waist) are `yes` when they were measured at this visit (venous blood for blood counts) and fasting where the candidate says so; `no` for urine,
   post-meal, OGTT, random or capillary values, a different analyte with a similar abbreviation (甲状腺球蛋白(TG)),
   or anything you cannot tell. Write `{"answers": [{"row_key": "…", "answer": "yes|no", "why": "…"}]}` and run
   `$LA labs confirm <ws> --answers <json> --reason "<what you checked>"`. Intake is not done until every numeric
   row is answered; editing a row later voids its answer. Rows answered `no` reach no method, index, estimator or twin. A row that names the right analyte in words the harness
   does not know is mapped instead: `$LA labs map <ws> --marker "<row>" --to <method input key or 空腹血糖/肌酐/…> --reason "…"`.
7. Record the user's own facts with `$LA member <ws> key=value --source "<quote>"` (also `sample_date=YYYY-MM-DD` when
   the date comes from a document; name that document in `--source`)
   (for example smoker=no, when `methods plan` asks for it later).
8. `$LA observe <ws>` → `pending_judgments` must be empty before going on.

## Output Format
1. One line per file: what it is and its status.
2. The judgments you made and the evidence behind each.
3. Anything unsupported (it stays listed in the report as not analysed).

A file whose content contradicts what it is said to be comes back `rejected` with the reason: duplicate probe ids,
a "β" column shaped like detection p-values, abundances that do not sum to a percentage, or a declared protein scale
the numbers do not fit. Tell the user and ask the lab. Do not re-answer the judgment until one fits.

## Guardrails
- A `rejected` methylation file (M-values or percentages) is re-exported as β by the lab, never rescaled by you.
- Do not transcribe values you cannot read clearly. Leave them out and say so.
