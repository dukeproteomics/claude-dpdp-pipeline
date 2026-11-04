# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports.rds"))
N <- X$none; miss <- is.na(N$q); lm <- log2(N$q); grp <- list(A = 1:3, B = 4:6, SPQC = 7:9)
lb <- log2(X$background$q); lg <- log2(X$global$q)
pm <- rowMeans(lm, na.rm = TRUE); idx <- which(miss, arr.ind = TRUE)
okb <- !is.na(lb[idx]) & is.finite(pm[idx[, 1]])
cat("Background value vs the precursor's mean in its measured runs: cor =", round(cor(lb[idx][okb], pm[idx[, 1]][okb]), 3),
    "| median (background - precursor mean) =", round(median(lb[idx][okb] - pm[idx[, 1]][okb]), 2), "log2\n")
# group context of each missing value
gm <- sapply(grp, function(cc) rowMeans(lm[, cc], na.rm = TRUE)); gn <- sapply(grp, function(cc) rowSums(!miss[, cc]))
gid <- c(rep(1, 3), rep(2, 3), rep(3, 3))[idx[, 2]]
others <- gn[cbind(idx[, 1], gid)]; gmean <- gm[cbind(idx[, 1], gid)]
mar <- others == 2   # the other two replicates of the group are measured
cat(sprintf("\nMissing values whose other 2 group replicates are measured: %d of %d (%.0f%%)\n", sum(mar), length(mar), 100 * mean(mar)))
for (k in c("background", "global")) { v <- if (k == "background") lb else lg
  d <- v[idx][mar] - gmean[mar]
  cat(sprintf("  %-10s fill minus the group's measured mean (log2): median %.2f, 25%%-75%% [%.2f, %.2f]; %.0f%% of fills more than 2 log2 (4x) below\n", k, median(d, na.rm = TRUE),
      quantile(d, .25, na.rm = TRUE), quantile(d, .75, na.rm = TRUE), 100 * mean(d < -2, na.rm = TRUE))) }
cat("  (the group's measured mean for these: median log2", round(median(gmean[mar]), 2), ")\n")
# SPQC technical replicate CV for precursors with a fill in SPQC, by method
sp <- 7:9; hasfill <- rowSums(miss[, sp]) == 1
cvf <- function(M) { x <- M[hasfill, sp]; ok <- rowSums(is.na(x)) == 0; apply(x[ok, ], 1, function(r) sd(r) / mean(r) * 100) }
cat(sprintf("\nSPQC precursors missing in exactly 1 of 3 injections (%d): median CV after filling: background %.0f%%, global %.0f%%, runwise %.0f%%  (fully measured SPQC precursors: %.0f%%)\n",
    sum(hasfill), median(cvf(X$background$q)), median(cvf(X$global$q)), median(cvf(X$runwise$q)),
    median(apply(N$q[rowSums(miss[, sp]) == 0, sp], 1, function(r) sd(r) / mean(r) * 100))))
