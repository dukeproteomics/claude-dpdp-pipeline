# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input")              # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports2.rds"))
N <- X$none; G <- X$global; fill <- is.na(N$q); flag <- G$imp; meas <- !fill; l <- log2(N$q)
clean <- meas & !flag
ref <- apply(ifelse(clean, l, NA), 1, median, na.rm = TRUE)          # precursor's median over its clean (unflagged measured) runs
dev <- l - ref
ok <- is.finite(dev)
cat("value minus the precursor's median of its unflagged measured runs (log2):\n")
for (lab in c("flagged & measured", "unflagged & measured")) { sel <- ok & meas & (if (lab == "flagged & measured") flag else !flag)
  cat(sprintf("  %-22s n=%7d  median %+.2f  25%%-75%% [%+.2f, %+.2f]  share > 1 log2 below: %.1f%%\n", lab, sum(sel), median(dev[sel]), quantile(dev[sel], .25), quantile(dev[sel], .75), 100 * mean(dev[sel] < -1))) }
# flag rate by deviation bin
b <- cut(dev[ok & meas], c(-Inf, -3, -2, -1, -0.5, 0, 0.5, 1, Inf)); t <- tapply(flag[ok & meas], b, mean)
cat("\n% flagged by deviation from the precursor's own median:\n"); print(round(100 * t, 1))
# and with intensity held fixed (log2 8-10)
sel <- ok & meas & l > 8 & l <= 10; b2 <- cut(dev[sel], c(-Inf, -2, -1, -0.5, 0, 0.5, Inf)); cat("\n... only values at log2 8-10:\n"); print(round(100 * tapply(flag[sel], b2, mean), 1))
