#!/usr/bin/env Rscript
# IDAT -> beta values with SeSAMe (prep "QCDPB": QC mask, channel inference,
# dye bias, pOOBAH detection masking, noob background). Writes a long CSV
# (cpg,beta) for the first array and one column per array otherwise.
# Usage: Rscript idat_to_beta.R <out.csv> <idat_prefix> [<idat_prefix> ...]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: idat_to_beta.R <out.csv> <idat_prefix>...")
suppressMessages(library(sesame))
out <- args[1]
prefixes <- args[-1]
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
betas <- sapply(prefixes, function(p) getBetas(prepSesame(readIDATpair(p), "QCDPB")))
if (is.null(dim(betas))) betas <- matrix(betas, ncol = 1, dimnames = list(names(betas), basename(prefixes)))
colnames(betas) <- basename(prefixes)
df <- data.frame(cpg = rownames(betas), betas, check.names = FALSE)
if (ncol(betas) == 1) colnames(df) <- c("cpg", "beta")
write.csv(df, out, row.names = FALSE, na = "")
cat(sprintf("wrote %d probes x %d arrays to %s\n", nrow(betas), ncol(betas), out))
