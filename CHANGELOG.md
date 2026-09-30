# Changelog

## [0.6.0] - 2026-09-30

### Added

- Insights stage: genes vs labs, population position, MR projections, question board
- insights genomics: live GWAS Catalog lead loci x member genotype (gVCF-aware), East Asian risk-allele percentile, ClinVar scan of monogenic genes via myvariant/Ensembl
- insights position: NHANES 2005-2010 weighted lab percentiles (27 analytes incl. fasting LDL/TG/glucose), GMHI vs 4347 metagenomes incl. East Asian healthy, multi-age comparison
- insights mr/project: EpiGraphDB MR estimates projected onto member value, target and baseline risk
- question board: 3-10 member-specific questions, researcher subagent per question, verdict bound to retrieved records (references/researcher.md)
- insights wearable: 30/90-day means and trends from a mapped daily export

### Fixed (v0.6 adversarial review: 8 P0, 13 P1)
- Genotype reads: every VCF row at a position is read (split multi-allelic records); the effect allele is matched
  exactly; a reference-base effect allele at a multi-allelic site is scored as 2 minus every alt copy; no tabix binary
  falls back to one streaming pass that keeps only the wanted positions and gene regions.
- Genetic percentile uses only loci that vary in East Asians (risk-allele frequency 1-99%, at least 8); a
  variant-only VCF gives no score unless the agent records `--absent-as-ref "<reason>"`, which the report prints.
- Monogenic scan: only alleles the member carries (0/2 annotates allele 2 only), indels annotated through genomic
  HGVS, calls must pass FILTER/GQ/DP and ≥3 alt reads (failures listed separately), ClinVar must be unanimously P/LP
  with review criteria (a "no assertion criteria" record no longer counts), zygosity + inheritance mode per gene;
  a het in a recessive gene is a carrier finding and is not counted; a source failure marks the gene "not scanned"
  and suppresses the count readout instead of writing 0.
- Board: questions and findings are bound to their registered hashes (trace blocks an edit afterwards); skip reasons
  are traced; cited PMIDs are verified live at registration; basis ids come only from confirmed lab rows.
- MR projection: baseline must be a percent risk; exposure must be a known name for the analyte or carry
  `--exposure-match "<reason>"`; lifetime-exposure and outcome-mismatch caveats printed; projection ref includes the
  baseline; an analyte without a reference gives a clear error.
- Population position: values tied at a detection limit report the percentile span; HDL now includes NHANES 2005-06
  (column LBDHDD), raising n in each stratum.
- Network: a host that failed is not retried per item in the same run; Ensembl gene regions are cached.

## [0.5.0] - 2026-09-29

### Added

- Organ checkup table: per-organ age and disease risk
- Published organ indices from labs: eGFR CKD-EPI 2021 + KDIGO G, FIB-4 with age-adjusted tiers, TyG, AIP
- AI estimate layer (llm_estimate): estimator subagent per organ, ranges, confidence low/very_low, basis ids, ≥1 retrieved live PMID; own report column 'AI 估计'
- New stage 'organs' between integrate and intervene; organ skip/register; stale organ readouts removed when upstream changes
- AI estimates cited in prose always render with 'AI 估计' and range
- eGFR refuses unknown sex instead of defaulting to the male equation

### Fixed (fifth independent re-attack: 3 P0, 9 P1, 16 P2)
- Lab rows (fifth and sixth re-attacks): name rules alone either let 尿肌酐/餐后血糖 through or withheld 肌酐(酶法)/
  葡萄糖(空腹), so row identity is now a judgment: `labs candidates` lists rows that name a known analyte and
  `labs confirm` records the agent's yes/no (serum/plasma, fasting where needed); only confirmed or mapped rows reach
  methods and organ indices, and editing a row voids its answer. Duplicate confirmed rows with different values,
  non-finite, implausible or mis-unit values stop the index; eGFR/FIB-4 only for adults; FIB-4 tier on the unrounded
  value, worded as a screening cut-off. A met diagnostic threshold blocks the matching AI probability: eGFR < 60 (CKD
  stage ≥ 3), fasting glucose ≥ 7.0 or HbA1c ≥ 6.5 (type 2 diabetes), FIB-4 > 2.67 (advanced fibrosis).
- Organ estimates: evidence is PubMed ids only (no links, doses or notes via `ref`); strict schema; bound to the
  current bundle hash; a rejected re-registration voids the old one; bases that are low-quality or not confirmed as
  the member's are refused; organ age needs organ-relevant data; AI estimates no longer appear in the readout table;
  the organ table carries quality/provenance labels; trace lists calibrated results beside AI estimates for review.
- Number binding: Chinese numerals containing 十一/一万 etc., numerals with measure words (次/天/个月/倍/成/分之),
  runs split by spaces or zero-width characters, uppercase numerals in doses, bare domains and disguised links are
  refused. APOE-declined check covers epsilon/ɛ/Cyrillic/Greek look-alikes and colloquial Chinese names.
- Consent refuses refusals and non-answers (negation tied to a run/consent verb, so "可以，不用再问了" still counts);
  re-bundling that voids an organ moves the organs stage back to todo; re-registering an analysis or a rejected plan
  no longer discards organ estimates; markdown emphasis, invisible format characters and Roman numerals cannot hide a
  number; drug and supplement doses stay refused everywhere while food and lifestyle amounts (盐, 饮水, 牛奶) may be
  written as `{{n:…}}`; common words (逐一, 唯一, 这一期间, 七八分饱) no longer trip the numeral check; a dose in a guideline
  URL path is refused; license reasons in the member report are Chinese; (seventh re-attack) candidates use the same
  name rewriting as the method matcher (肌酐 Cr, 空腹血糖 GLU), `{{r:id|own words}}` cannot carry a dose unit, 单位/g/克
  count as a dose next to a drug or supplement name, concentrations (3 mg/L) are not doses, overrides fire on any
  confirmed value past the threshold (HbA1c in mmol/mol, 126 mg/dL = 7.0, unrounded FIB-4), exemption phrases never
  touch another numeral, APOE-declined refuses any sentence tying dementia to genes; (eighth re-attack) a drug name
  counts only within a few characters of an amount (diet advice with 维生素/叶酸/因素 in the sentence passes), lab flags
  (8.4↑, 8.4 H) are stripped, overrides read blank-unit and duplicate rows, consent needs an affirmative word and
  refuses questions and conditions, 百分比/千克 are not numerals, 表观遗传 is not genetics; (ninth re-attack) more flag
  forms are stripped (8.4H, 8.4(↑), 8.4 偏高, 7.2%↑) and an organ with an unreadable confirmed threshold row cannot be
  estimated until the row is fixed; a lab row from another visit is answered no; nutrient names count as dose context
  only right before an amount (food descriptions with 维生素/益生菌 pass); (tenth re-attack) 补充/服用 + a nutrient counts
  as a dose, a blank-unit threshold row is readable only if plausible in the canonical unit, and the reviewer checks every
  supplement amount in `{{n:…}}` (the script is a floor for that judgment); (final E2E feedback) twin.json keeps the
  provenance-unconfirmed flag and twin compare does not judge such values, `organ register` applies the number rules
  to `rationale_zh` (no value inside `{{r:…}}`, no member age as a literal), an identity revision removes excluded files
  from the uncertain list, a failed method's reasons reach the analysts, the transcription header asks for ref_range,
  and the workflows now state the re-registration order, the APOE sentence rule, `retest.what` rules and when to ask
  optional answers; plan items and retest accept only declared fields (twin.json no longer
  copies extra retest fields); guideline URLs must be plain addresses and show their host; reviewer severity must be
  P0/P1/P2; identity values are fixed; MS header words are matched on word boundaries and in file/column names;
  pipeline logs with R/Python errors fail verification; a non-UTF-8 library VERSION no longer crashes the report.

## [0.4.1] - 2026-09-29
### Fixed (fourth independent re-attack of v0.4.0: 0 P0, 9 P1, 15 P2)
- β shape check adds slope/intercept/mean-deviation limits (refuses ×0.5, +0.05, β², compressed M-values); an
  implausible clock age stops every method fed by the same β, not only epiage.
- Process identity uses start time read in a fixed locale and time zone; stub runs never read member inputs.
- A trace once blocked by a reviewer can never be recorded as pass; `{{pmid:…}}` verified live at trace;
  `{{n:…}}` cannot hold links or hide doses; ftp:// and //host links refused; plan item ids are I1…I99.
- APOE-declined check covers spelling variants; GQ=-1 no longer crashes; stale deliverables are deleted.
- Rejected plan registration voids the previously registered plan; validate re-hashes raw files; β-check failures
  and typed-in-age blocks are reported with their true reason; GMHI version check reads the whole profile;
  build inference cross-checks chr1 and chr2; 25(OH)D comparisons across seasons carry a caveat.

## [0.4.0] - 2026-09-29
### Changed (third independent re-attack of v0.3.0: 1 P0 + 20 P1 new, 36 partial)
- Number binding redesigned: authored prose may contain digits only inside `{{r:…}}` (member readouts), `{{n:…}}`
  (explicit non-member literals, listed in trace.json for the reviewer) and `{{pmid:…}}`; no digit, Chinese-numeral
  run or URL elsewhere. The regex allow-list (and its 42% false-positive rate and bypass channels) is gone.
- Reviewer JSON must carry the current `trace_id` (binds the review to exactly the traced content); final verdicts
  are pass|block; severity parsed across key variants; any P0/P1 forbids pass; a trace once blocked can never pass.
- `{{pmid:…}}` citations in prose are verified live at trace time; literals may not carry links or hide doses; pid
  identity uses a locale- and time-zone-independent start time; APOE-declined text checks cover spelling variants.
- The state digest covers the whole state except bookkeeping; the report prints fixed Chinese status wording and
  never free-text reasons (no trace deadlock, no reason channel).
- Methylation β must correlate (r ≥ 0.8, ≥ 1000 CpGs) with the whole-blood reference profile shipped with epiage;
  clocks below 50% CpG coverage are not reported; any clock age outside 0–125 stops epiage and every other method fed
  by the same β values. The β check also refuses shifted or rescaled values (slope/intercept/mean deviation against
  the reference, calibrated on three public whole-blood 450K samples).
- verify refuses runs whose inputs were excluded (also while running); excluding files stops their live runs;
  pid identity includes start time; gzip outputs are checked; only "Pipeline completed successfully" counts.
- Commercial mode: methods fed a typed-in age are not run; APOE record GQ/AD/ploidy checks; declined APOE text is
  refused by trace; lock file lives beside the workspace with a generation check; lock timeout exit code 7;
  stale deliverables moved to .stale/ on every save; EPICv2 replicates that disagree are refused; rejected answers
  can be corrected; `labs fix --row`; dose check covers word forms; `tests/test_cli_chain.py` runs the whole chain.

## [0.3.0] - 2026-09-29
### Fixed (independent re-attack of v0.2.0: 102 original P0/P1 re-run; 5 new P0 and 10 new P1 found and closed)
- Duplicate probe ids refused; β columns shaped like detection p-values refused; clock ages outside human range
  marked implausible and stop the run; `ch.` probes kept in derived β files.
- APOE reference calls from gVCF blocks need GQ ≥ 20 and depth ≥ 10; RefCall/LowQual/low-GQ calls are not called.
- Excluded inputs void their pipeline runs (plan, launch and verify refuse them); verify needs a launched, successful
  attempt with non-empty outputs; launch uses the latest preflight tier.
- `init --force` takes the old workspace's lock and refuses while its pipelines run; unique backup names; the member
  id is fixed at init; locks queue for up to 120 s instead of failing.
- Trace/review/render bind a digest of the state fields the report prints; free-text reasons are traced; only
  registered analyses are used; reviewer JSON must live in work/review, be newer than the trace and not be reused;
  severity parsing normalised; stale deliverables are moved aside on invalidation; validate checks deliver hashes.
- Readouts and derived files are hash-bound between methods and trace; raw files re-checked before methods run.
- Trace: identifier/value tokens (LDL160, GrimAge58), traditional Chinese numerals, heart-rate and unit-bearing
  years, month offsets are caught; markdown images and non-http links are neutralised.
- Plan check: NFKC-normalised dose detection with capsule/tablet/spoon units; retracted PMIDs rejected; integer weeks.
- APOE results need the member's `genetic_disclosure` answer; `identity --uncertain` labels readouts from files whose
  ownership cannot be confirmed; Weidner and Ying clocks removed from the commercial list; licensed clock values typed
  as lab rows are dropped in commercial mode.
- Note: v0.3.0 claims about reviewer reuse, clock-age stopping and state digest coverage were only partly true;
  corrected in 0.4.0.
- Intake: Olink long tables, `sp|…|` UniProt ids, NA-tolerant sample columns, declared protein scale checked against
  the values, header synonyms (`检验项目名称`, `Test Name`), unknown tables routed to transcription, MetaPhlAn version
  read from content, abundance sums checked, content-identical copies of excluded files excluded.

## [0.2.0] - 2026-09-29
### Changed (after three adversarial reviews: 31 P0 / 70 P1 / 61 P2, and two E2E rounds)
- Trace, reviewer verdict and render are bound to file hashes; any edit after trace voids them. Stage order enforced;
  upstream changes invalidate downstream stages. `init --force` moves the old workspace aside.
- Multi-sample tables and VCFs ask which column is the member; detection p-value columns are set aside; every table is
  rewritten into a normalised single-sample file. Protein matrices need platform and value scale; VCFs without a build
  ask for it; several files of one kind ask which is this visit's.
- Identity verdict can exclude files; nothing derived from them survives. Reviewer JSON must match the recorded verdict.
- Commercial clock list narrowed (no GrimAge, DunedinPACE, DNAm PhenoAge and derivatives, skin & blood, DNAmTL) and
  marked pending legal review; pyaging not auto-run.
- APOE merges split multi-allelic records and reads gVCF reference blocks; sarek now uses DeepVariant (gVCF).
  GMHI only on MetaPhlAn2 profiles. epiage coverage parsed per clock; low coverage flagged.
- Trace normalises text (NFKC, HTML entities, link text, Chinese numerals, number+unit tokens); every rendered field
  is traced or escaped. Plan check: doses in Chinese numerals and pill counts, live PMID verification, fetched URLs.
- Intake: GBK/UTF-16/xlsx, title rows, header synonyms, symlink dedupe, hidden files skipped, bounded reads, corrupt
  files marked unreadable. `labs fix|drop`. Workspace lock. Per-method exception isolation with checkpoints.
- Pipelines: stub needs LA_DEV=1 and writes to results_stub/; skipped runs need --reopen; pgsc_calc needs --pgs-id;
  outputs carry source file ids; verify ignores outputs older than the current attempt.

## [0.1.0] - 2026-09-29
### Added
- Router skill `longevity-analyst` with six stage workflows and the `la.py` harness.
- Intake sniffing for FASTQ/BAM/CRAM/VCF/IDAT/β tables/MetaPhlAn/protein matrices/lab tables/PDFs, with pending judgments.
- Preflight that actually executes java/docker/nextflow/git/sesame and estimates wall time and disk per pipeline.
- Pinned local pipelines (sarek 3.10.0, methylseq 4.2.0, taxprofiler 2.0.1, pgsc_calc v2.3.0, PharmCAT v3.4.0, SeSAMe);
  fetched with the system git; launch needs the user's quoted consent.
- Method dispatch over longevity-skills with platform gate (MS proteomics refused by affinity clocks) and licence modes.
- Native methods: APOE genotype, VCF summary, gut diversity, GMHI (validated against GMHI_2020 output), proteomics QC.
- System bundles + analyst/reviewer contracts, evidence lookup (local + live PubMed), plan checks, twin build/compare
  with RCV, number trace, report rendering.
