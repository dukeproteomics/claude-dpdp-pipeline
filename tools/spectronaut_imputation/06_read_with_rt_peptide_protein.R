# Part of tools/spectronaut_imputation (see README.md). Analysis of the 11488 Spectronaut exports, 2026-11-04.
WORK <- Sys.getenv("SPN_WORK", "/home/user/Claude_Shiny/data/spectronaut_compare")   # intermediate .rds files (git-ignored data/)
INPUT <- Sys.getenv("SPN_INPUT", "/home/user/Claude_Shiny/data/input/11488")           # the four Spectronaut exports
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(data.table))
f <- c(none = "20261104_140253_11488_Enkbayar_103126_Report.tsv", background = "20261104_140341_11488_Enkbayar_103126_background_Report.tsv",
       global = "20261104_135116_11488_Enkbayar_103126_impute_global_Report.tsv", runwise = "20261104_135110_11488_Enkbayar_103126_impute_runwise_Report.tsv")
num <- function(d, cols) as.matrix(d[, ..cols][, lapply(.SD, function(x) suppressWarnings(as.numeric(x)))])
out <- list()
for (k in names(f)) {
  d <- fread(file.path(INPUT, f[[k]]), sep = "\t", na.strings = c("", "NA", "NaN", "Filtered"), showProgress = FALSE)
  nm <- names(d)
  if (k == "none") cat("all columns:\n", paste(nm, collapse = " | "), "\n\n")
  out[[k]] <- list(info = d[, c("PG.ProteinGroups", "PG.Qvalue", "PEP.StrippedSequence", "EG.ModifiedSequence", "EG.PrecursorId"), with = FALSE],
                   q = num(d, grep("EG.TotalQuantity", nm, value = TRUE)), qval = num(d, grep("EG.Qvalue", nm, value = TRUE)),
                   rt = num(d, grep("EG.ApexRT", nm, value = TRUE)),
                   imp = as.matrix(d[, grep("EG.IsImputed", nm, value = TRUE), with = FALSE][, lapply(.SD, function(x) toupper(as.character(x)) == "TRUE")]))
}
saveRDS(out, file.path(WORK, "exports2.rds"))
for (k in names(out)) cat(k, ": q-value matrix identical to none?", isTRUE(all.equal(out[[k]]$qval, out$none$qval)), " | RT identical?", isTRUE(all.equal(out[[k]]$rt, out$none$rt)),
                          " | RT NA:", sum(is.na(out[[k]]$rt)), " q NA:", sum(is.na(out[[k]]$qval)), "\n")
