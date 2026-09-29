# Changelog

All notable changes to this skill are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · Versioning: [SemVer](https://semver.org/)

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

### Fixed

- Round-4 re-attack fixes (9 P1, 15 P2)
- β slope/intercept/deviation check; age gate covers all methods fed by the same β
- locale/TZ-independent pid start time; a blocked trace can never pass; live PMID check in prose
- literals cannot hold links or hide doses; APOE-declined spelling variants; stale deliverables deleted

## [0.4.0] - 2026-09-29

### Added

- Round-4: strict placeholder number binding, trace_id-bound review, reference-profile β check
- prose digits only inside {{r:}} / {{n:}} / {{pmid:}}; literals listed for the reviewer
- reviewer JSON carries trace_id; verdicts pass|block; robust severity parsing
- whole-state digest; report prints fixed Chinese status wording, never free-text reasons
- β values must correlate with the whole-blood reference profile; any implausible clock age stops the run
- verify refuses runs on excluded inputs; exclusion stops running pipelines; pid start-time identity

## [0.3.0] - 2026-09-29

### Added

- Round-3 hardening after re-attack (5 new P0, 10 new P1 closed)
- duplicate probes refused; detection-p shaped columns refused; clock plausibility gate
- gVCF reference calls need GQ>=20 and depth>=10; RefCall/LowQual not called
- excluded inputs void pipeline runs; init --force refuses live processes; lock waits
- state digest bound into trace/review/render; reviewer file freshness and reuse checks; stale deliver moved aside
- APOE disclosure consent; uncertain-provenance labelling; Weidner/Ying out of commercial list

## [0.2.0] - 2026-09-29

### Added

- Hardening after three adversarial reviews and two E2E rounds
- content-bound trace/review/render; enforced stage order and downstream invalidation
- sample-column/value-scale/assembly/primary judgments; derived single-sample files; identity exclusion
- commercial licence list narrowed and marked pending legal review; pyaging blocked
- APOE multi-record + gVCF; GMHI MetaPhlAn2-only; epiage coverage per clock; sarek DeepVariant gVCF
- trace normalisation (NFKC, entities, Chinese numerals), HTML escaping, live PMID check, workspace lock

