# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports.rds"))
N <- X$none; miss <- is.na(N$q); lm <- log2(N$q)
qs <- c(0, .01, .02, .05, .1, .25, .5)
for (k in c("global", "runwise")) { M <- X[[k]]; li <- log2(M$q)
  cat("\n=====", k, "=====\n")
  allm <- lm[!miss]
  cat("imputed log2 overall: min/quantiles(5,25,50,75,95)/max:", round(c(min(li[miss], na.rm=TRUE), quantile(li[miss], c(.05,.25,.5,.75,.95), na.rm = TRUE), max(li[miss], na.rm=TRUE)), 2), "\n")
  cat("measured log2 percentiles 0,1,2,5,10,25,50 (experiment):", round(quantile(allm, qs), 2), "\n")
  cat("per run: n imputed | imputed min / median / max | measured 1st / 5th / 10th pct | exact copies of a measured value in same run / any run\n")
  for (j in seq_len(ncol(lm))) { iv <- M$q[miss[, j], j]; iv <- iv[!is.na(iv)]; mv <- N$q[!miss[, j], j]
    cat(sprintf("  run %d: %5d | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %4.1f%% / %4.1f%%\n", j, length(iv), min(log2(iv)), median(log2(iv)), max(log2(iv)),
      quantile(log2(mv), .01), quantile(log2(mv), .05), quantile(log2(mv), .10), 100 * mean(iv %in% mv), 100 * mean(iv %in% N$q[!miss]))) }
  # does the imputed value depend on the precursor? correlation with the precursor's mean measured intensity, and with background
  pm <- rowMeans(lm, na.rm = TRUE); idx <- which(miss, arr.ind = TRUE); ok <- is.finite(pm[idx[, 1]]) & !is.na(li[idx])
  cat("cor(imputed, precursor mean of measured runs):", round(cor(li[idx][ok], pm[idx[, 1]][ok]), 3), "\n")
  B <- log2(X$background$q); okb <- ok & !is.na(B[idx]); cat("cor(imputed, background value at same position):", round(cor(li[idx][okb], B[idx][okb]), 3), "\n")
  d <- li[idx][ok]; h <- hist(d, breaks = 30, plot = FALSE); cat("imputed log2 histogram (counts):", h$counts, "\n  breaks:", round(range(h$breaks), 1), "\n")
}
