# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input/11488")           # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
X <- readRDS(file.path(WORK, "exports.rds"))
N <- X$none; B <- X$background; G <- X$global; R <- X$runwise
miss <- is.na(N$q); qv <- B$qval
cat("No-impute missing:", sum(miss), " of which q-value NA:", sum(miss & is.na(qv)), " q>0.01:", sum(miss & qv > 0.01, na.rm = TRUE), " q<=0.01:", sum(miss & qv <= 0.01, na.rm = TRUE), "\n")
cat("Present in no-impute:", sum(!miss), " of which q<=0.01:", sum(!miss & qv <= 0.01, na.rm = TRUE), " q>0.01:", sum(!miss & qv > 0.01, na.rm = TRUE), " q NA:", sum(!miss & is.na(qv)), "\n\n")
for (k in c("global", "runwise")) { M <- X[[k]]
  cat(k, ": imputed flag TRUE where no-impute is missing:", sum(M$imp & miss), " | where no-impute has a value:", sum(M$imp & !miss), "\n")
  same <- !miss & !M$imp; cat("   not imputed and present in no-impute: values identical?", isTRUE(all.equal(M$q[same], N$q[same])), "\n")
  r <- M$q[M$imp & !miss] / N$q[M$imp & !miss]; cat("   flagged-but-measured: imputed/measured ratio quantiles:", round(quantile(r, c(.05, .5, .95), na.rm = TRUE), 3), "\n")
  cat("   flagged-but-measured: measured log2 quantiles:", round(quantile(log2(N$q[M$imp & !miss]), c(.05, .5, .95), na.rm = TRUE), 2), "| all measured:", round(quantile(log2(N$q[!miss]), c(.05, .5, .95), na.rm = TRUE), 2), "\n")
}
cat("\nBackground: values where no-impute missing:", sum(!is.na(B$q) & miss), " still NA:", sum(is.na(B$q)), "\n")
same <- !miss & !is.na(B$q); cat("   background equals no-impute where both present?", isTRUE(all.equal(B$q[same], N$q[same])), "\n")
cat("   background-filled log2 quantiles:", round(quantile(log2(B$q[miss & !is.na(B$q)]), c(.05, .25, .5, .75, .95)), 2), "\n")
cat("   measured log2 quantiles:          ", round(quantile(log2(N$q[!miss]), c(.05, .25, .5, .75, .95)), 2), "\n")
# per precursor: how many runs flagged in global
pr <- rowSums(G$imp); cat("\nglobal: precursors with 0..9 runs flagged:", tabulate(pr + 1, 10), "\n")
pm <- rowSums(miss);  cat("no-impute: precursors with 0..9 runs missing:", tabulate(pm + 1, 10), "\n")
cat("cross-tab (rows: runs missing in no-impute 0..9, cols: runs flagged in global) for rows 0-3:\n"); print(table(missing = pmin(pm, 9), flagged = pr)[1:4, ])
