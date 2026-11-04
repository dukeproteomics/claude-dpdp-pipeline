# Record the analyst's review of a review-only step and rebuild its page.
# Usage: Rscript pipeline/accept.R --db=<db> --step=01_load|03_parameters --decided_by=user [--analyst=<name>] [--note="<why>"]
# (QC and Stats have their own: 09_qc.R --action=accept, 10_stats.R --action=accept.)

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "common.R"))
source(file.path(here, "review.R"))
args <- parse_args()
need(args, c("db", "step"))
if (!args$step %in% c("01_load", "03_parameters")) stop("--step must be 01_load or 03_parameters (use 09_qc.R / 10_stats.R --action=accept for those).")
if (!identical(args$decided_by, "user")) stop("accept needs --decided_by=user (after the analyst has reviewed the page).")
db_path <- normalizePath(args$db, mustWork = TRUE)
use_app("Shiny_File.R")
source(file.path(here, "reviews", "01_load.R")); source(file.path(here, "reviews", "03_parameters.R"))
log_decision(dirname(db_path), args$step, list(accepted = TRUE), note = args$note %||% "", analyst = args$analyst)
page <- if (args$step == "01_load") review_load(db_path) else review_parameters(db_path)
cat(args$step, "accepted\nReview page:", page, "\n")
