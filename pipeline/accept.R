# Record the analyst's review of a review-only step and rebuild its page.
# Usage: Rscript pipeline/accept.R --db=<db> --step=01_load|03_parameters --decided_by=user [--from_page] [--analyst=<name>] [--note="<why>"]
# 03_parameters: accepts the settings of the last run of the step. --from_page takes the values the analyst sent with "Accept these
# settings"; option values that differ from the last run are refused (re-run the step with them first).
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
run_dir <- dirname(db_path)
acc <- list(accepted = TRUE)
if (args$step == "03_parameters") {
  source(file.path(here, "options.R")); source(file.path(here, "options", "03_parameters.R"))
  summ <- fromJSON(file.path(run_dir, "steps", "03_parameters", "summary.json"), simplifyVector = FALSE)
  ls <- file.path(run_dir, "steps", "01_load", "summary.json")
  opts <- parameters_options(get_params(db_path), read_table("design", db_path),
                             if (file.exists(ls)) fromJSON(ls, simplifyVector = FALSE)$detected else NULL)
  spec <- opts_spec(opts); all <- opts_byname(opts)
  args <- page_args(args, run_dir, "03_parameters", opts, "apply")
  ran <- resolve_settings(spec, summ$settings[intersect(names(spec), names(summ$settings))])
  given <- resolve_settings(spec, modifyList(summ$settings, args[intersect(names(spec), names(args))]))
  diff <- Filter(function(k) k %in% names(args) && opt_in_use(all[[k]], given, all) && !isTRUE(all.equal(given[[k]], ran[[k]])), names(spec))
  if (length(diff)) stop("These settings differ from the last run of Parameters: ", paste(sprintf("%s (run %s, now %s)", diff,
    vapply(diff, function(k) as.character(ran[[k]]), ""), vapply(diff, function(k) as.character(given[[k]]), "")), collapse = "; "),
    ". Re-run Parameters with them first.", call. = FALSE)
  acc$settings <- ran
}
log_decision(run_dir, args$step, acc, note = args$note %||% "", analyst = args$analyst)
page <- if (args$step == "01_load") review_load(db_path) else review_parameters(db_path)
cat(args$step, "accepted\nReview page:", page, "\n")
