# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; G <- X$global; B <- X$background; R <- X$runwise
fill <- is.na(N$q); flag <- G$imp; meas <- !fill
both <- meas & !is.na(N$qval)
for (k in c("background", "global", "runwise")) { M <- X[[k]]
  same <- abs(M$qval[both] - N$qval[both]) < 1e-12
  cat(sprintf("%-10s q-value equal to the no-impute export where both are measured: %.2f%% | RT equal: %.2f%%\n", k, 100 * mean(same), 100 * mean(abs(M$rt[both] - N$rt[both]) < 1e-9)))
  cat(sprintf("           measured cells with q > 0.01 in this export: %d (flagged among them: %d)\n", sum(meas & M$qval > 0.01, na.rm = TRUE), sum(meas & M$qval > 0.01 & M$imp, na.rm = TRUE)))
}
r <- G$qval[both] / N$qval[both]; cat("global/none q ratio quantiles (1,5,50,95,99%):", signif(quantile(r, c(.01, .05, .5, .95, .99), na.rm = TRUE), 3), "\n")
fm <- flag & meas
cat("flagged-measured: q in no-impute export median", signif(median(N$qval[fm], na.rm = TRUE), 3), " vs in global export", signif(median(G$qval[fm], na.rm = TRUE), 3), "\n")
cat("Global vs runwise flags identical:", identical(G$imp, R$imp), " | background flags all FALSE:", !any(B$imp), "\n")
# filled cells: q-values in the imputing exports
cat("filled cells: q in global export: quantiles", signif(quantile(G$qval[fill], c(0, .05, .5, .95, 1), na.rm = TRUE), 3), "\n")
cat("filled cells with q <= 0.01 in global export:", sum(fill & G$qval <= 0.01, na.rm = TRUE), "\n")
