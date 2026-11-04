# Part of tools/spectronaut_imputation (see README.md). What changes downstream if "measured" = EG.IsImputed FALSE (global export)
# instead of "has a value in the no-imputation export"? One worked example plus totals for 11488.
suppressPackageStartupMessages(library(data.table))
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")
X <- readRDS(file.path(WORK, "exports2.rds")); N <- X$none; G <- X$global
lab <- c("A_1","A_2","A_3","B_1","B_2","B_3","SPQC_1","SPQC_2","SPQC_3"); grp <- list(A = 1:3, B = 4:6, SPQC = 7:9)
m_empty <- !is.na(N$q)          # measured by the empty-value rule
m_flag  <- !G$imp               # measured by EG.IsImputed (global export)
# --- 1. a precursor that is measured in all 9 runs, but flagged in 2 of the 3 B runs ---
cand <- which(rowSums(m_empty) == 9 & rowSums(!m_flag[, 4:6]) == 2 & rowSums(!m_flag[, c(1:3, 7:9)]) == 0)
i <- cand[order(-rowMeans(log2(N$q[cand, ])))][1]
cat("EXAMPLE PRECURSOR:", N$info$EG.PrecursorId[i], "| protein", N$info$PG.ProteinGroups[i], "\n")
print(data.frame(sample = lab, value_no_impute = round(N$q[i, ]), value_global = round(G$q[i, ]), q_value = signif(G$qval[i, ], 2),
                 IsImputed = G$imp[i, ], measured_empty_rule = m_empty[i, ], measured_flag_rule = m_flag[i, ]), row.names = FALSE)
cat("-> same numbers in both exports; by the flag rule B has 1 of 3 measured instead of 3 of 3.\n")
# --- 2. precursor filters (as in Step 5): at least 2 measured overall, and >= 50% measured in the best group ---
pass <- function(m) rowSums(m) >= 2 & apply(sapply(grp, function(cc) rowSums(m[, cc]) / length(cc)), 1, max) >= 0.5
pe <- pass(m_empty); pf <- pass(m_flag)
cat(sprintf("\nPRECURSOR FILTERS (min 2 measured, 50%% of a group): pass with empty rule %d, with flag rule %d -> %d more precursors removed\n", sum(pe), sum(pf), sum(pe & !pf)))
# --- 3. protein level: missing factor for B vs A (app: share of samples with >= 1 measured precursor, best of the two groups) ---
pg <- N$info$PG.ProteinGroups
prot <- function(m) { d <- data.table(pg = pg, m); d[, lapply(.SD, function(x) as.integer(any(x))), by = pg] }
pe_ <- prot(m_empty); pf_ <- prot(m_flag)
mf <- function(p) pmax(rowMeans(as.matrix(p[, 2:4])), rowMeans(as.matrix(p[, 5:7])))   # A = V1-V3, B = V4-V6
mfe <- mf(pe_); mff <- mf(pf_)
cat(sprintf("PROTEINS removed by missing factor < 0.6 (B vs A): empty rule %d, flag rule %d\n", sum(mfe < 0.6), sum(mff < 0.6)))
one <- function(p) { a <- rowSums(as.matrix(p[, 2:4])); b <- rowSums(as.matrix(p[, 5:7])); (a == 0 & b > 0) | (b == 0 & a > 0) }
oe <- one(pe_); of <- one(pf_)
cat(sprintf("PROTEINS measured in one group only (B vs A): empty rule %d, flag rule %d (%d new)\n", sum(oe), sum(of), sum(of & !oe)))
# --- 4. an example protein whose status changes ---
k <- which(mfe >= 0.6 & mff < 0.6)
if (length(k)) { j <- k[1]; p <- pe_$pg[j]; rows <- which(pg == p)
  cat("\nEXAMPLE PROTEIN:", p, "-", N$info$PEP.StrippedSequence[rows[1]], "...", length(rows), "precursors\n")
  cat("  samples with >= 1 measured precursor (A_1..B_3): empty rule", paste(unlist(pe_[j, 2:7]), collapse = " "), " -> mf", mfe[j], "(kept)\n")
  cat("                                                   flag rule ", paste(unlist(pf_[j, 2:7]), collapse = " "), " -> mf", mff[j], "(removed)\n")
  print(data.frame(sample = lab[1:6], t(sapply(1:6, function(s) c(values = sum(m_empty[rows, s]), flagged_imputed = sum(G$imp[rows, s] & m_empty[rows, s])))), row.names = NULL)) }
