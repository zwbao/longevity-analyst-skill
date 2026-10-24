# case-C public-data fixture

`pubdata/` holds the answers of the public databases that `la.py insights genomics` asks for case-C, in the format of
pubdata's on-machine cache (`~/.cache/longevity-analyst/pubdata`, file name = sha1 of `source|query`). The test points
`pubdata.SHARED` here and turns the network off, so the chain runs offline and gives the same result every time.
Each file keeps the URL and the time it was retrieved; the harness copies both into the workspace's evidence ledger
exactly as it would from a live run.

| Source | What | Retrieved |
|---|---|---|
| GWAS Catalog REST v2 | associations, p ≤ 5e-8, for the six case-C analytes (EFO_0004761 urate, OBA_2050062 ALT, OBA_VT0000188 glucose, EFO_0004611 LDL, EFO_0006925 Lp(a), EFO_0004578 homocysteine) | 2026-10-24 |
| myvariant.info | position, alleles and gnomAD East Asian frequency of the lead variants (hg38) | 2026-10-24 |
| Ensembl REST | spans of the monogenic genes in `data/trait_map.json` | 2026-10-24 |
| NCBI ClinVar (E-utilities) | pathogenic / likely pathogenic / conflicting records per monogenic gene | 2026-10-24 |

The ClinVar answers were shortened after retrieval: only records within 2 kb of a record of the two case-C VCFs
(`examples/case-C/2023-09/raw/HG00472…`, `examples/case-C-swap/2025-09/raw/NA18639…`) are kept (each file says how
many were retrieved in `filtered_for_fixture`). The scan only looks up the member's own variants (and indels within
100 bp), so for these two VCFs the result is the same as with the full answers; for any other VCF it is not.

To record again (network): `LA_CASE_C_RECORD=<empty dir> pytest tests/test_cli_chain.py -k case_c`, then copy the
directory here and shorten the ClinVar files the same way.
