# Build (or rebuild) a step's review page: <run_dir>/review/<step>.html, and the run's review/index.html.
# Usage: Rscript pipeline/review_page.R --db=<db> --step=04_noise
# Steps with a review page: 01_load, 03_parameters, 04_noise, 05_filter, 06_normalize, 07_impute, 08_rollup, 09_qc, 10_stats

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "common.R"))
source(file.path(here, "review.R"))
source(file.path(here, "contaminants.R"))
args <- parse_args()
need(args, c("db", "step"))
db_path <- normalizePath(args$db, mustWork = TRUE)
f <- file.path(here, "reviews", paste0(args$step, ".R"))
if (!file.exists(f)) stop("No review page for step ", args$step, " yet. Available: ",
                          paste(sub("\\.R$", "", list.files(file.path(here, "reviews"), pattern = "\\.R$")), collapse = ", "))
use_app("Shiny_File.R")
source(file.path(here, "reviews", "01_load.R"))   # shared accept helpers
source(f)
special <- c("01_load" = "load", "03_parameters" = "parameters")
fun <- get(paste0("review_", if (args$step %in% names(special)) special[[args$step]] else sub("^[0-9]+_", "", args$step)))
out <- fun(db_path)
cat("Review page:", out, "\n")
