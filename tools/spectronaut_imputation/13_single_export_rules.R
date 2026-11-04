# Part of tools/spectronaut_imputation (see README.md). Can ONE export tell measured from filled? Scored against the no-imputation export.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; truth <- is.na(N$q)                      # TRUE = Spectronaut filled this value
score <- function(pred, lab) { tp <- sum(pred & truth); fp <- sum(pred & !truth); fn <- sum(!pred & truth)
  cat(sprintf("  %-58s filled found %6d of %6d (%.2f%%) | measured wrongly called filled %6d (%.3f%% of measured)\n", lab, tp, sum(truth), 100 * tp / sum(truth), fp, 100 * fp / sum(!truth))) }
for (k in c("global", "runwise")) { M <- X[[k]]; q <- M$qval; l <- log2(M$q); f <- M$imp
  cat("\n==", k, "export ==\n")
  score(f, "EG.IsImputed alone")
  score(f & (is.na(q) | q > 0.01), "IsImputed & q > 0.01")
  sure <- f & (is.na(q) | q > 0.01)   # certainly fills: use them to learn the fill range
  if (k == "global") { r <- range(l[sure], na.rm = TRUE); inr <- l >= r[1] - 1e-9 & l <= r[2] + 1e-9 }
  else { inr <- matrix(FALSE, nrow(l), ncol(l)); for (j in 1:ncol(l)) { r <- range(l[sure[, j], j], na.rm = TRUE); inr[, j] <- l[, j] >= r[1] - 1e-9 & l[, j] <= r[2] + 1e-9 } }
  score(f & (is.na(q) | q > 0.01 | inr), "IsImputed & (q > 0.01 or value inside the learned fill range)")
  amb <- f & !is.na(q) & q <= 0.01 & inr & !truth
  cat(sprintf("  ambiguous: flagged measured values that happen to sit inside the fill range: %d\n", sum(amb, na.rm = TRUE)))
}
cat("\n== background export ==\n")
B <- X$background; q <- B$qval
score(is.na(q) | q > 0.01, "q > 0.01 (or no q)")
score(!is.na(B$q) & (is.na(q) | q > 0.01), "q > 0.01 and a background value is present")
# the identified-but-blanked background values: do they behave like measured values?
l <- log2(B$q); clean <- !truth & !X$global$imp
ref <- apply(ifelse(clean, l, NA), 1, median, na.rm = TRUE); dev <- l - ref
ibb <- truth & !is.na(q) & q <= 0.01 & !is.na(B$q)
cat(sprintf("  identified-but-blanked background values (%d): vs precursor's own level median %+.2f, IQR [%+.2f, %+.2f]  (clean measured: [%+.2f, %+.2f]; filled with q > 0.01: median %+.2f)\n",
    sum(ibb), median(dev[ibb], na.rm = TRUE), quantile(dev[ibb], .25, na.rm = TRUE), quantile(dev[ibb], .75, na.rm = TRUE),
    quantile(dev[clean], .25, na.rm = TRUE), quantile(dev[clean], .75, na.rm = TRUE), median(dev[truth & q > 0.01], na.rm = TRUE)))
