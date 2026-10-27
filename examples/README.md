# Demo members

Each demo member is a synthetic person spliced together from public data. The files for one demo member come from
different people, so do not read any relationship between them as a finding.

| Case | File | Source | Real or synthetic |
|---|---|---|---|
| A (female, 67) | `methylation_450k_betas.csv.gz` | GEO GSM989827 (Hannum et al. 2013, GSE40279), whole blood, 67 y F | real |
| A | `wgs_NA18525.chr19_APOE.vcf.gz` | 1000 Genomes 30× phased panel (20220422), NA18525 (CHB), chr19:44.89–44.92 Mb | real, one region only |
| A | `gut_metaphlan_profile.tsv` | GMHI_2020 `4347_final_relative_abundances.txt`, sample Healthy_641 (SAMEA104142355, China, 67 y M) | real |
| A | `体检化验_2026-09.csv` | written for the demo | synthetic |
| B (male, 73) | `甲基化_beta值.csv` | GEO GSM989833 (Hannum), whole blood, 73 y M | real |
| B | `HG00096.apoe.vcf.gz` | 1000 Genomes, HG00096 (GBR), same region | real, one region only |
| B | `stool_profile.txt` | GMHI_2020, Healthy_650 (SAMEA104142342, China, 70 y M) | real |
| B | `plasma_proteome_DIA_matrix.tsv` | 30 real UniProt IDs, random log2 intensities | synthetic (gate test) |
| B | `2026年度体检报告_检验部分.pdf` | generated with reportlab | synthetic |
| B | `unknown_run_R{1,2}.fastq.gz` | 200 random reads | synthetic (modality judgment + pipeline wiring test) |

## case-C: one member, four visits (PLAN appendix C)

A synthetic Chinese man, 49 at the first visit, tested on 2023-09-04, 2024-09-09, 2025-09-15 and 2026-09-15. Each
visit is a folder `examples/case-C/<yyyy-mm>/` with `raw/` (what the member hands in), `plan.json` (the plan the
agent registers at that visit) and, from the second visit on, `exposures.json` (what the member did since, as the
coach's member file reports it). `tests/test_cli_chain.py::test_full_chain_case_c` runs all four visits.

| File | Source | Real or synthetic |
|---|---|---|
| `2023-09/raw/HG00472.caseC_regions.vcf.gz` | 1000 Genomes 30× phased panel (20220422, `1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`, chrX `.v2`), sample HG00472 (CHS, male), read with pysam over the remote tabix index on 2026-10-24; small regions only: ±1 kb around the 40 urate and 39 ALT lead loci the genes-vs-labs step scores, APOE (chr19:44,905,000–44,910,000), SLCO1B1 (chr12:21,177,615–21,179,615) and ±200 bp around each pharmacogenomic and nutrigenetic site; 4,475 records, 0\|0 rows kept | real, regions only |
| `*/raw/体检化验_<yyyy-mm>.csv` | written for the demo: fasting glucose 4.6 / 5.0 / 5.3 / 5.6 mmol/L, LDL-C 3.1 / 3.3 / 3.0 / 3.2, urate in range at visits 1–2 and above from visit 3, ALT high at visit 2 only, Lp(a) at visit 1 only, homocysteine at visits 1 and 4, clinic blood pressure at visits 1–3 | synthetic |
| `2026-09/raw/家庭血压_2026-07至09.csv` | written for the demo: 14 morning readings from 2026-07-06 (mean 138/86 mmHg), salt reduction adopted 2026-07-20, 14 readings from 2026-08-31 after six weeks (mean 131/82) | synthetic |
| `*/plan.json`, `*/exposures.json` | written for the demo (each plan item cites a row of the library's effects table; check-ins weekly) | synthetic |

Why HG00472: among the 208 CHB and CHS samples of the panel, 7 have a urate genetic percentile of at least 86 (band
"elevated", away from the 80 edge) and an ALT percentile between 30 and 70 ("typical"); HG00472 is the one man among
them who carries SLCO1B1 rs4149056 C (decreased function), so the CPIC item has something to say. Its urate
percentile is 99, ALT 40, APOE ε2/ε3. One sample covers both designs, so there is no case-D. The selection script's
inputs and the per-sample scores are kept outside the repository (sources/integrate/case-c-vcf).

## case-C-swap: the third visit brings someone else's DNA

`examples/case-C-swap/2025-09/raw/` is case-C's 2025-09 folder with the VCF replaced by another real sample:

| File | Source | Real or synthetic |
|---|---|---|
| `NA18639.caseC_regions.vcf.gz` | same panel and regions as above, sample NA18639 (CHB, male): differs from HG00472 at APOE, SLCO1B1 and many lead loci | real, regions only |
| `体检化验_2025-09.csv` | copy of case-C's | synthetic |

The other visits are case-C's own folders. `test_case_c_swap_identity_alert` expects an identity alert for that visit
and its values left out of every trend.
