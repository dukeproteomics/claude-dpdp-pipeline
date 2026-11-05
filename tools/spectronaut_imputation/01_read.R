# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input/11488")           # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
# Read the four Spectronaut exports into aligned matrices (quantity, qvalue, imputed flag) and save as RDS.
suppressPackageStartupMessages(library(data.table))
f <- c(none = "20261104_140253_11488_Enkbayar_103126_Report.tsv", background = "20261104_140341_11488_Enkbayar_103126_background_Report.tsv",
       global = "20261104_135116_11488_Enkbayar_103126_impute_global_Report.tsv", runwise = "20261104_135110_11488_Enkbayar_103126_impute_runwise_Report.tsv")
out <- list()
for (k in names(f)) {
  d <- fread(file.path(INPUT, f[[k]]), sep = "\t", na.strings = c("", "NA", "NaN", "Filtered"), showProgress = FALSE)
  nm <- names(d)
  qc <- grep("EG.TotalQuantity", nm, value = TRUE); qv <- grep("EG.Qvalue", nm, value = TRUE); im <- grep("EG.IsImputed", nm, value = TRUE)
  run <- sub("^\\[(\\d+)\\].*", "\\1", qc)
  raw_q <- d[, ..qc]
  out[[k]] <- list(id = d$EG.PrecursorId, acc = d$PG.ProteinAccessions, genes = d$PG.Genes, desc = d$PG.ProteinDescriptions,
                   q = as.matrix(raw_q[, lapply(.SD, function(x) suppressWarnings(as.numeric(x)))]),
                   qval = as.matrix(d[, ..qv][, lapply(.SD, function(x) suppressWarnings(as.numeric(x)))]),
                   imp = as.matrix(d[, ..im][, lapply(.SD, function(x) toupper(as.character(x)) == "TRUE")]),
                   runs = run, qcols = qc,
                   raw_tokens = head(sort(table(unlist(lapply(raw_q, function(x) x[is.na(suppressWarnings(as.numeric(x))) & !is.na(x)]))), decreasing = TRUE), 5))
  cat(sprintf("%-10s rows %6d  unique precursors %6d  runs %d  NA quantities %6d  zeros %d  imputed TRUE %6d  q>0.01 %6d\n", k, nrow(d), uniqueN(d$EG.PrecursorId), length(qc),
              sum(is.na(out[[k]]$q)), sum(out[[k]]$q == 0, na.rm = TRUE), sum(out[[k]]$imp, na.rm = TRUE), sum(out[[k]]$qval > 0.01, na.rm = TRUE)))
}
saveRDS(out, file.path(WORK, "exports.rds"))
cat("column order same in all:", length(unique(lapply(out, function(x) x$runs))) == 1, "\n")
cat("same precursor sets:", all(sapply(out, function(x) setequal(x$id, out$none$id))), "\n")
