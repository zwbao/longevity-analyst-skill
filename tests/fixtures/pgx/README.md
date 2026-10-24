# pgx fixtures

| File | Source | Real or synthetic |
|---|---|---|
| `NA18526.chr12_SLCO1B1.vcf.gz` | 1000 Genomes 30× phased panel (20220422), `1kGP_high_coverage_Illumina.chr12.filtered.SNV_INDEL_SV_phased_panel.vcf.gz` on ftp.1000genomes.ebi.ac.uk, sample NA18526 (CHB), chr12:21177615–21179615 (GRCh38), read with pysam over the tabix index on 2026-10-24; 53 records. rs4149056 (chr12:21178615 T>C, SLCO1B1 c.521T>C) = 1\|0 | real, one region only |
| `NA18526_rs4149056_ensembl.23andme.txt` | one 23andMe-format row written from Ensembl REST `/variation/human/rs4149056?genotypes=1` (retrieved 2026-10-24): `1000GENOMES:phase_3:NA18526` = C\|T | real genotype, file written for the test |

`tests/test_pgx_adapter.py` runs `cpic-pharmacogenomics` on the adapter's output from the VCF and on the Ensembl-derived
file and asserts the skill's outputs are identical (SLCO1B1 Decreased Function).
