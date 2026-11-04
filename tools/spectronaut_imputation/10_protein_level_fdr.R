# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(data.table))
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; G <- X$global; fill <- is.na(N$q); flag <- G$imp; qg <- G$qval
pg <- N$info$PG.ProteinGroups
idf <- fill & qg <= 0.01           # identified but blanked (13,573)
dt <- data.table(pg = rep(pg, 9), run = rep(1:9, each = length(pg)), fill = as.vector(fill), idf = as.vector(idf), flag = as.vector(flag), q = as.vector(qg))
pr <- dt[, .(n = .N, nfill = sum(fill), min_q = min(q, na.rm = TRUE)), by = .(pg, run)]
pr[, all_filled := nfill == n]
dt <- merge(dt, pr, by = c("pg", "run"))
cat(sprintf("Identified-but-blanked cells (%d): their protein group has EVERY precursor blanked in that run: %.1f%%\n", sum(dt$idf), 100 * dt[idf == TRUE, mean(all_filled)]))
cat(sprintf("  ... protein-group size (precursors) median %.0f; for comparison all cells median %.0f\n", dt[idf == TRUE, median(n)], dt[, median(n)]))
cat(sprintf("Filled cells with q > 0.01 (%d): protein group fully blanked in that run: %.1f%%\n", dt[fill & !idf & !is.na(q), .N], 100 * dt[fill & !idf & !is.na(q), mean(all_filled)]))
# protein-run groups fully blanked although best precursor q <= 0.01
fb <- pr[all_filled == TRUE]
cat(sprintf("\nProtein-group x run combinations fully blanked: %d; of these, best precursor q <= 0.01: %d (%.0f%%)\n", nrow(fb), sum(fb$min_q <= 0.01), 100 * mean(fb$min_q <= 0.01)))
cat("  best-precursor q for fully blanked groups with an identified precursor: quantiles", signif(quantile(fb[min_q <= 0.01, min_q], c(.05, .5, .95)), 3), "\n")
# not-fully-blanked: identified-but-blanked inside groups that still have measured precursors
cat(sprintf("Identified-but-blanked in groups that still have measured precursors in that run: %d\n", dt[idf == TRUE & all_filled == FALSE, .N]))
