# Reference data fitted or assembled by this library

These are not published methods. Each file says where its numbers come from; the engine labels every result that
uses them as this library's own method.

| File | Used by | Source | Licence |
| --- | --- | --- | --- |
| `nhanes_age_reference.json` | `estimate-age` | Fitted by `tools/reference/fit_nhanes_age.py` on NHANES 2011-2012, 2013-2014, 2015-2016 (adults, not pregnant); held out 2017-2018 | NHANES: US public domain |
| `olink_age_direction_ckb.csv.gz` | `compare-age` (proteomics) | Iona et al., Eur J Epidemiol 2025, doi 10.1007/s10654-025-01311-z, eTable 5 (sex All, main model, exposure Age); values copied, rounded to 5 significant digits | CC BY 4.0 |
| `dnam_age_direction_gs.csv.gz` | `compare-age` (methylation) | McCartney et al., Genome Med 2020, doi 10.1186/s13073-019-0693-z; Edinburgh DataShare doi 10.7488/ds/2709, Supplementary Table 1 (autosomes); kept CpGs whose discovery and replication t have the same sign and discovery P < 1e-5; z = (t1 + t2) / sqrt(2) | CC BY 4.0 |
| `expression_age_gene_sets.json` | `age-evidence --modality expression` | Curated literature priors (SenMayo subset, inflammaging/complement, mito OXPHOS decline, ECM remodeling); **not** fitted on LongevityBench or GTEx eval | assembled here |

None of the methylation/proteomics/expression reference sources is a LongevityBench training source for the vote tables. Checked against the benchmark's own source labels on
2026-11-16: its clinical items are NHANES 1999-2000 participants; its methylation items come from 16 named GEO
series; its proteomic items come from small public Olink Explore cohorts (not UK Biobank). Generation Scotland and
China Kadoorie Biobank individual data are not on GEO and are not among those cohorts.
