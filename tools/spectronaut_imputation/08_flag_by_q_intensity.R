# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; G <- X$global; fill <- is.na(N$q); flag <- G$imp; qg <- G$qval; meas <- !fill
nmiss <- matrix(rowSums(fill), nrow(fill), 9)
rate <- function(cond, by, brks, lab) { b <- cut(by[cond], brks, include.lowest = TRUE); t <- tapply(flag[cond], b, mean); n <- table(b)
  cat(lab, "\n"); print(round(data.frame(n = as.integer(n), flagged_pct = 100 * as.numeric(t), row.names = names(t)), 1)) }
rate(meas, qg, c(0, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 0.01), "Measured values: % flagged by run q-value")
rate(meas, nmiss, c(-0.5, 0.5, 1.5, 2.5, 4.5, 8.5), "Measured values: % flagged by number of runs the precursor is missing in (no-impute)")
rate(meas, log2(N$q), c(-Inf, 4, 6, 8, 10, 12, Inf), "Measured values: % flagged by log2 intensity")
cat("\nPer precursor: are the flags all-or-nothing across measured runs?\n")
fr <- rowSums(flag & meas); mr <- rowSums(meas); some <- fr > 0
print(table(all_measured_runs_flagged = (fr == mr)[some]))
# within precursor: flagged runs vs unflagged runs of the same precursor - q-value comparison
idx <- which(rowSums(flag & meas) > 0 & rowSums(!flag & meas) > 0)
d <- sapply(idx, function(i) { a <- qg[i, flag[i, ] & meas[i, ]]; b <- qg[i, !flag[i, ] & meas[i, ]]; c(max(b), min(a), median(a) > median(b)) })
cat("Precursors with both flagged and unflagged measured runs:", length(idx), "| flagged runs have the higher median q:", round(100 * mean(d[3, ]), 1), "%",
    "| every flagged run has higher q than every unflagged run:", round(100 * mean(d[2, ] > d[1, ]), 1), "%\n")
