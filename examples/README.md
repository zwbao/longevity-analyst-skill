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
