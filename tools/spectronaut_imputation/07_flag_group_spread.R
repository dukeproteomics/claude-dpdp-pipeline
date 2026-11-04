# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(data.table))
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; G <- X$global
fill <- is.na(N$q); flag <- G$imp; qg <- G$qval
fm <- flag & !fill   # flagged but measured
cat("cells: filled", sum(fill), "| flagged", sum(flag), "| flagged & measured", sum(fm), "| flagged & filled", sum(flag & fill), "| filled & not flagged", sum(fill & !flag), "\n")
cat("no-impute export: quantity missing but q/RT present:", sum(fill & !is.na(N$qval)), " (q values there:", paste(round(quantile(N$qval[fill & !is.na(N$qval)], c(0, .5, 1)), 4), collapse = "/"), ")\n")
# H1/H2: group propagation within a run
test_group <- function(key, label) {
  dt <- data.table(grp = rep(key, 9), run = rep(1:9, each = length(key)), fill = as.vector(fill), flag = as.vector(flag))
  dt[, gfill := any(fill), by = .(grp, run)]
  cat(sprintf("%-22s flagged & measured cells whose %s has a filled precursor in that run: %.1f%% | unflagged measured cells with that: %.1f%%\n",
      label, label, 100 * dt[flag & !fill, mean(gfill)], 100 * dt[!flag & !fill, mean(gfill)]))
  dt
}
d1 <- test_group(N$info$PG.ProteinGroups, "protein group")
d2 <- test_group(N$info$PEP.StrippedSequence, "peptide")
d3 <- test_group(N$info$EG.ModifiedSequence, "modified peptide")
# H3: identification quality of flagged-but-measured cells in that run
cat("\nrun-level q-value (global export) of measured cells: flagged median", signif(median(qg[fm]), 3), "| unflagged median", signif(median(qg[!flag & !fill]), 3), "\n")
cat("  share with q > 0.001: flagged", round(100 * mean(qg[fm] > 0.001), 1), "% | unflagged", round(100 * mean(qg[!flag & !fill] > 0.001), 1), "%\n")
gq <- N$info$EG.GlobalPrecursorQvalue
