# longevity-analyst-skill

This skill turns one member's longevity test data into three deliverables:

- a multi-omics report;
- a digital-twin snapshot (`twin.json`) that can be compared at the next retest;
- an intervention plan that cites evidence and that the nutritionist reviews.

It accepts WGS (FASTQ, CRAM or VCF), DNA methylation (β-value tables or IDAT), gut metagenome data (FASTQ or a
MetaPhlAn table), plasma proteomics matrices, and lab reports as tables, PDFs or photos. It runs nf-core pipelines
locally. Before each run it checks whether the machine can handle the job and estimates run time and disk use. The
run starts only after the user agrees.

中文说明见 [README.zh-CN.md](README.zh-CN.md).

## Install

```bash
npx skills add zwbao/longevity-analyst-skill
git clone https://github.com/zwbao/longevity-skills
export LONGEVITY_SKILLS_HOME=$PWD/longevity-skills
pip install pandas numpy markdown
# optional, for local pipelines: nextflow, docker (≥16–64 GB RAM for the VM), JDK 17+, R + sesame
```

Then give your agent a data folder: *"我 67 岁，女，检测数据在 ~/data/me，帮我做一份抗衰分析报告和干预方案。"*

## What is computed, and what is not

| Data | Computed in this version | Not computed, and why |
|---|---|---|
| Blood methylation β | epiage clocks with per-clock CpG coverage; eFRS. Commercial mode runs a provisional list (Horvath 2013, Hannum, Lin, Vidal-Bralo, Zhang, EpiTOC1, StocH, StocZ) that has not yet been cleared by counsel | GrimAge, DunedinPACE, DNAm PhenoAge and derivatives, skin & blood, DNAmTL, Weidner, Ying in commercial mode (licensing); non-blood tissue |
| Routine labs (table, xlsx, PDF, photo) | PhenoAge (9 markers), China-PAR 10-year ASCVD risk, other verified methods | methods whose inputs are not declared (listed as `manual`) |
| VCF / gVCF | APOE genotype (merges split multi-allelic records, reads gVCF reference blocks); optional PharmCAT and pgsc_calc (user picks PGS ids) | ACMG secondary findings (not in v0.2); APOE from a variant-only VCF when a site is absent |
| Gut MetaPhlAn profile | richness, Shannon, Gini-Simpson; GMHI on MetaPhlAn2 profiles (matches the published GMHI.R output to within 0.001) | GMHI on MetaPhlAn3/4 (species renamed); HUMAnN |
| Proteomics matrix | QC description only | every protein clock unless platform **and** value scale match its training assay |
| Organ checkup table (9 organs) | Measured organ ages when an organ clock ran; published indices from labs: eGFR (CKD-EPI 2021, race-free) with KDIGO G category, FIB-4 with age-adjusted tiers, TyG, AIP; China-PAR for the heart. **AI estimates** (organ-age range and disease-probability ranges) by an estimator subagent anchored on PubMed papers retrieved in the run, shown in their own "AI 估计" column | Calibrated organ-disease probabilities: no validated model exists for most organs. Most organ ages are AI estimates, because measured organ clocks need Olink/SomaScan proteomics |
| Insights (v0.6) | Genes vs labs: for each confirmed lab, the member's genotype at the genome-wide significant loci retrieved live from GWAS Catalog, a risk-allele count placed in the gnomAD East Asian distribution, and a ClinVar scan of the member's own variants in monogenic genes. Population position: lab percentiles in the NHANES 2005–2010 same-sex, same-age-band population; GMHI against 4,347 public metagenomes (healthy, East Asian healthy, non-healthy); the several ages side by side. Wearable 30/90-day summaries. MR projections: published EpiGraphDB estimates applied to the member's own value, target and baseline risk. The question board: 5–10 member-specific questions, one researcher subagent each, verdict bound to retrieved records | Chinese population percentiles (CHARLS/CKB need registration); genome-wide polygenic scores (need a full WGS/gVCF, not a subset VCF); a researcher's verdict is judgment on retrieved evidence, not a validated model |
| FASTQ / CRAM / IDAT | sarek (DeepVariant gVCF), taxprofiler, methylseq, SeSAMe (pinned revisions), after preflight and the user's quoted consent | nothing runs without consent; stub runs are developer checks and never results |

## Known limitations (v0.6)

- Genes vs labs uses an unweighted count of trait-raising alleles at clumped genome-wide significant loci (GWAS
  Catalog effect sizes are in mixed units, so they are not summed). It is a tendency with its locus coverage, not a
  polygenic score from a validated PGS, and needs the loci to be present in the VCF (a gVCF or full VCF).
- Lab percentiles come from the US NHANES population; glucose and triglycerides also have fasting-subsample tables.
  No Chinese population table is bundled.
- An MR projection assumes the published population causal estimate (exposure GWAS in SD units) applies to the member.
  MR reflects lifelong exposure, so the projected benefit of lowering a value in adulthood is an upper bound.
- The monogenic scan covers the genes listed per analyte in `data/trait_map.json` (Ensembl span ± 2 kb), small variants
  with a canonical SPDI in ClinVar (structural variants and haplotype records are out of scope), and is complete for a
  gene only when the VCF covers its region. Phase is not known, so two het P/LP variants in a recessive gene are
  "possible compound heterozygous". First retrieval of ClinVar for all genes takes minutes; it is cached 30 days.
  It is a screen, not a clinical genetic test.
- Board citations: records must be retrieved in the workspace by the harness (cache files are checked against a hash
  ledger in state.json); PMIDs are also checked live. An agent that edits state.json by hand can defeat any of these
  bindings; the harness assumes changes go through `la.py`.

- AI estimates were benchmarked against real outcomes on public cohorts (NHANES 2005–2008 with 2019 mortality
  linkage; the Framingham teaching dataset): Claude, one person per call, matched a cross-validated logistic
  regression in discrimination (10-year death C 0.855 vs 0.847; incident CVD 0.733 vs 0.741; hypertension 0.747 vs
  0.755) and beat age + sex alone, but over-predicted incidence (O/E 0.61 for CVD, 0.74 for stroke). Dementia, CKD,
  diabetes, fracture, sarcopenia and organ age have not been validated, and no Chinese cohort has been tested yet.

- Organ AI estimates are not calibrated. The harness checks their form (ranges no narrower than 6 years or ±0.02 and
  half the point value, confidence low/very_low, basis ids that exist, ≥1 PubMed article retrieved in the run and live);
  whether the numbers are reasonable is judged by the reviewer. Disease probabilities for the same member can differ
  between two estimator runs.
- eGFR and FIB-4 come from a single measurement; a CKD or fibrosis diagnosis needs repeat testing. The report says so.
- Which lab row is which analyte is a judgment: the harness proposes candidates by name, and every numeric row reaches a
  method, an organ index, an estimator or the twin only after the agent confirmed it (`labs confirm`: serum/plasma,
  fasting where needed) or mapped it. The confirmation is as good as the agent's reading of the source document.
  Rows carry no collection date, so two confirmed rows with different values stop the index instead of picking one.
- "PMID retrieved in this run" is checked against the workspace's PubMed cache files, which a misbehaving agent could
  write by hand; that the PMID exists and is not retracted is checked live. Pipeline verification likewise trusts the
  run directory's `exitcode` and `run.log`.
- Where a diagnostic threshold is already met (eGFR < 60, fasting glucose ≥ 7.0, HbA1c ≥ 6.5, FIB-4 > 2.67) no AI
  probability is given for that disease. Other contradictions between an AI estimate and calibrated results (for example
  an organ age far below chronological age beside eGFR G4) are judged by the reviewer from `organ_checks_for_reviewer`.
- Drug and supplement doses are refused by unit (mg, IU, 粒, 片 …) everywhere, and by any amount written within a few
  characters of a known drug or supplement name. A drug name the list does not know, written with a food unit (0.1 g),
  is left to the reviewer; the skill's own rule is never to name a drug or a dose. Telling a supplement amount from a
  food amount is a judgment: the script catches drug units and obvious phrasings (补充/服用 + a nutrient, 维生素C 1 g), and
  the reviewer is instructed to treat every other supplement amount in `{{n:…}}` as a P1 (in the last re-attack the script
  alone caught 11 of 20 supplement-dose phrasings and let 18 of 20 food lines through).
- Some real lab names are not known to the method library (红细胞平均体积 without MCV, ALKP, TRIG); those rows reach no
  method until mapped with `labs map`, and the report lists the method as missing an input.
- Guideline links are checked for being fetched live and for not carrying a dose in the URL; whether the host is a real
  guideline body is the reviewer's call (no domain allow-list).

- The commercial clock list is provisional and has not been reviewed by counsel. The report says so.
- PharmCAT and pgsc_calc outputs are registered after a verified run but not yet summarised into readouts; PharmCAT's
  VCF preprocessor is not called by the harness.
- Drug names in the plan are checked by the independent reviewer, not by the script; the script's dose check is
  best-effort (numbers, Chinese numerals and common dose words) and the reviewer is the backstop.
- Whether a `{{n:…}}` literal is really a member value is judged by the reviewer, not the script. Member values spelled
  out in English words or Roman numerals ("sixty-seven") are also left to the reviewer; the script catches digits,
  Chinese numerals and dose words only.
- The β shape check (correlation, slope, intercept and mean deviation against the epiage whole-blood reference) was
  calibrated on three public whole-blood 450K samples; arrays preprocessed very differently may be refused and need
  the lab's unmodified SeSAMe/minfi β values.
- Reviewer independence is a process rule (a separate subagent with the reviewer contract); the harness can check the
  review's form and binding, not who wrote it.
- Reference change values always use CVA = 0.5 × CVI (EFLM planning value); there is no input yet for a lab's own CVA.
- WGBS/EM-seq coverage files are registered but no clock reads them yet; methylation clocks need array β values.
- nf-core containers run under amd64 emulation on Apple Silicon; the `arm` profile is not used.
- Methods with undeclared inputs in longevity-skills (`manual`) are never auto-run.

## Layout

```
skills/longevity-analyst/
  SKILL.md                 router
  workflows/01..06-*.md    stage procedures
  references/              subagent contracts (system analyst, reviewer), PDF transcription, plan schema
  scripts/la.py            harness CLI (typed commands, observe/validate)
  scripts/lalib/           intake, preflight, pipelines, methods, native, integrate, evidence, twin, report
  data/                    pinned pipelines + benchmarks, platform gate, licence policy, systems map, GMHI lists, menu
tests/                     pytest (gates each carry a negative test)
examples/case-A, case-B    two demo members built from public data (see examples/README.md)
```

## Tests

```bash
LONGEVITY_SKILLS_HOME=... python -m pytest tests -q
```

## Licence

MIT for this repository. The GMHI species lists come from jaeyunsung/GMHI_2020, and their citation is kept in
`data/gmhi.json`. Methods from longevity-skills keep their own licences.
