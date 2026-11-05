# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input/11488")           # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(data.table))
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; B <- X$background; G <- X$global; fill <- is.na(N$q); q <- G$qval
nob <- fill & is.na(B$q)
cat("Gaps background could not fill:", sum(nob), "\n")
cat("  their q-value (global export): NA", sum(nob & is.na(q)), "| <= 0.01", sum(nob & q <= 0.01, na.rm = TRUE), "| > 0.01", sum(nob & q > 0.01, na.rm = TRUE), "\n")
cat("  RT reported in background export for them:", sum(nob & !is.na(B$rt)), "\n")
pg <- N$info$PG.ProteinGroups
dt <- data.table(pg = rep(pg, 9), run = rep(1:9, each = length(pg)), fill = as.vector(fill), nob = as.vector(nob), idf = as.vector(fill & q <= 0.01))
dt[, all_filled := all(fill), by = .(pg, run)]
cat(sprintf("  share in protein-group x runs that are fully blanked: %.1f%% (all gaps: %.1f%%)\n", 100 * dt[nob == TRUE, mean(all_filled)], 100 * dt[fill == TRUE, mean(all_filled)]))
cat(sprintf("  share that are identified-but-blanked (q <= 0.01): %.1f%%\n", 100 * dt[nob == TRUE, mean(idf, na.rm = TRUE)]))
# flagged & filled vs background value presence
cat("\nGaps background DID fill:", sum(fill & !is.na(B$q)), "| of the identified-but-blanked (13,573):", sum(fill & q <= 0.01 & !is.na(B$q), na.rm = TRUE), "filled by background\n")
