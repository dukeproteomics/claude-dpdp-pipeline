# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input/11488")           # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports.rds"))
N <- X$none; miss <- is.na(N$q); lm <- log2(N$q)
ends <- function(M, j) { iv <- log2(M$q[miss[, j], j]); iv <- iv[!is.na(iv)]; c(min(iv), max(iv)) }
cat("RUN-WISE: percentile of the run's measured log2 values at the imputation range ends (and at the range ends on the linear scale)\n")
for (j in 1:9) { e <- ends(X$runwise, j); mv <- lm[!miss[, j], j]
  cat(sprintf("  run %d: range %.3f-%.3f  -> measured percentile %.3f%% - %.3f%%  | width %.3f | mean of measured %.3f, sd %.3f, (lo-mean)/sd %.3f, (hi-mean)/sd %.3f\n", j, e[1], e[2],
      100 * mean(mv <= e[1]), 100 * mean(mv <= e[2]), e[2] - e[1], mean(mv), sd(mv), (e[1] - mean(mv)) / sd(mv), (e[2] - mean(mv)) / sd(mv))) }
eg <- c(min(sapply(1:9, function(j) ends(X$global, j)[1])), max(sapply(1:9, function(j) ends(X$global, j)[2]))); am <- lm[!miss]
cat(sprintf("GLOBAL: range %.3f-%.3f -> experiment measured percentile %.3f%% - %.3f%% | mean %.3f sd %.3f | (lo-mean)/sd %.3f (hi-mean)/sd %.3f\n",
    eg[1], eg[2], 100 * mean(am <= eg[1]), 100 * mean(am <= eg[2]), mean(am), sd(am), (eg[1] - mean(am)) / sd(am), (eg[2] - mean(am)) / sd(am)))
# is it percentiles of each precursor's... try percentiles of the run's measured values of precursors that are complete? Try the median-of-runs distribution
pmed <- apply(lm, 1, median, na.rm = TRUE); cat(sprintf("GLOBAL vs precursor-median distribution: percentile %.3f%% - %.3f%%\n", 100 * mean(pmed <= eg[1], na.rm=TRUE), 100 * mean(pmed <= eg[2], na.rm=TRUE)))
