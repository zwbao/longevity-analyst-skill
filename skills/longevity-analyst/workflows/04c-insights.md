# 04c · Insights: genes vs labs, population position, causal projections, question board

## Use When
The organ table is done. This stage turns the member's data plus live public databases into insights a single-omics
report cannot give: which abnormal labs have a genetic basis, where the member stands among peers, whether lowering a
risk factor would plausibly help him, and a board of questions investigated one by one.

## Inputs Needed
- A workspace past `organs`.
- `genetic_disclosure=yes` and a usable VCF for the genomics step; otherwise skip it with the reason.
- Network access (GWAS Catalog, myvariant.info, Ensembl, EpiGraphDB, PubMed). A source that cannot be reached is
  reported as not retrieved; nothing is filled from memory.

## Workflow
1. Genes vs labs: `$LA insights genomics <ws>`. For every confirmed lab with a trait mapping it retrieves the genome-wide
   significant GWAS loci, reads the member's genotype there, places a risk-allele count in the East Asian distribution,
   and scans the member's own variants in the monogenic genes against ClinVar. Output: `work/insights/genotype_phenotype.json`.
   No VCF or genetic results declined: `$LA insights skip-genomics <ws> --reason "<why>"`.
   Read the result before using it:
   - `score_not_computed` says why no percentile exists. If most loci are `not_called` because the VCF lists variant
     sites only, decide whether an absent site can be read as reference for this delivery (the lab's note says every
     non-reference call is listed, whole-genome depth, joint calling). If yes, re-run with
     `--absent-as-ref "<what says so>"`; the report states that judgment. If you cannot tell, leave the score out.
   - `monogenic_scan.not_scanned` lists genes not fully checked: a source failed, or the VCF does not cover the gene
     region (a variant-only VCF without the absent-as-reference judgment). Carried variants in those genes were still
     matched, so a finding there is real, but "none found" may not be said for them.
   - `plp_not_unanimous` are ClinVar records that mention pathogenic without a unanimous, reviewed classification
     (conflicting, "risk factor", no criteria): mention them as "needs genetic counselling", never as a diagnosis,
     never as nothing.
   - `carrier_only: true` is a heterozygous finding in a recessive gene: the member is a carrier; it does not explain
     their lab value. `possible_compound_het` needs phasing before anyone calls it a cause.
   - `low_quality_not_counted` hits need orthogonal confirmation (Sanger) before they are mentioned as findings.
2. Position: `$LA insights position <ws>` → `work/insights/positions.json` (lab percentiles in the NHANES same-sex,
   same-age-band population; GMHI against healthy, East Asian healthy and non-healthy metagenomes; the several ages
   side by side).
3. Wearables, when a daily export exists: map its columns and run
   `$LA insights wearable <ws> --file F00x --map '{"date": "日期", "steps": "步数", "rhr": "静息心率", "hrv": "HRV_RMSSD_ms", "sleep_h": "睡眠时长_小时", "spo2_min": "血氧_最低"}'`.
   Map only columns whose meaning you can tell from the header and the lab's note.
4. Causal projections, for risk factors the member can change: retrieve published MR estimates
   `$LA insights mr <ws> --exposure "LDL cholesterol" --outcome "Coronary heart disease"`, choose one record (prefer IVW
   or weighted median, high `moescore`), then project it onto the member:
   `$LA insights project <ws> --mr-ref "<ref>" --analyte ldl --target <value in the table's unit> --baseline china-par-ascvd-risk.risk_10y_pct`.
   The target must lie inside the population's 1st–99th percentile. The baseline must be an absolute risk: a method
   readout named `*risk*_pct` / `*mortality*_pct` (percent) or an organ-table `organ.<organ>.risk.<n>`, whose outcome
   is as close as possible to the MR outcome; CHD in the MR study against a 10-year ASCVD baseline is a mismatch the
   report states. The exposure must be the same quantity as `--analyte`; a name outside the known list needs
   `--exposure-match "<why it is the same quantity>"`. MR estimates reflect lifelong differences, so the projected
   benefit of lowering the value in adulthood is an upper bound.
5. Question board. Read `genotype_phenotype.json`, `positions.json`, projections, readouts, labs, the questionnaire and
   the wearable summary. Write 5–10 questions specific to this member, the kind a physician-scientist would ask after
   seeing everything at once (for example: an abnormal lab with a matching gene finding; ages that disagree; a symptom
   in the questionnaire that the wearable data speaks to; a family history that a variant could explain). Save
   `{"questions": [{"id": "Q1", "title_zh": "…", "hypothesis_zh": "…", "basis": ["<member data ids>"], "why_zh": "…"}]}`
   and run `$LA board questions <ws> --file <json>`.
6. Dispatch one researcher subagent per question, in parallel, with [`../references/researcher.md`](../references/researcher.md),
   the question, the insight files and the output path `work/insights/board/Q<n>.json`, plus: "You do not know the
   member or the company's preferences." Register each: `$LA board finding <ws> --id Q<n>`; rejections name what to
   fix. A question that cannot be studied with the data at hand: `$LA board skip <ws> --id Q<n> --reason "<why>"`.
7. The stage is done when position ran, genomics ran or was skipped, and every question has a finding or a skip.
   If an upstream change voids this stage (for example an organ estimate registered again), `work/insights` is moved
   to `work/insights.previous`. Run genomics / position / wearable / projections again (they are quick), then copy
   `board/questions.json` and the researchers' `board/Q<n>.json` back and register them again: they are checked again
   against the new data. Never rewrite a researcher's finding yourself; a finding that no longer passes goes back to
   a researcher.
   The plan (workflow 05) may then target insight readouts and cite `{"type": "proj", "ref": "<projection ref>"}` or
   `{"type": "board", "ref": "Q<n>"}`.

## Output Format
1. Genes vs labs: per abnormal lab, the genetic percentile with coverage and any ClinVar P/LP finding.
2. Position: the percentiles that matter and the age comparison, with which reference they come from.
3. Projections: exposure → outcome, current → target, baseline → projected risk.
4. The board: question, verdict, confidence, next step.

## Guardrails
- A genetic percentile is a tendency, not a diagnosis; say that a common-variant score explains part of a trait at
  most. A ClinVar P/LP finding is reported as "needs clinical confirmation and genetic counselling", never as a
  diagnosis, and only when the member agreed to genetic results.
- NHANES is a US population: for a Chinese member it gives direction, not a Chinese percentile. Say so.
- An MR projection assumes the population causal effect applies to the member; it is not a prediction for him.
- Researchers cite only records retrieved in this workspace; their verdicts are `supported`, `not_supported` or
  `insufficient`, with confidence `low` or `moderate`.
