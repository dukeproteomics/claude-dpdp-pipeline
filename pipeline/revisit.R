# Go back to a finished step to change it (user, 2026-11-13). The run is changed in place, with a backup.
#
#   Rscript pipeline/revisit.R --db=<db> --step=<NN_step> [--mode=auto|step] [--from_page] [--analyst=<name>]
#     Backs up the run (<run>/Backup/revisit_<time>_<step>/: the .db, steps/, review/, QC/, decisions.jsonl), writes
#     <run>/revisit.json, and re-runs the step's preview with the settings the analyst applied, so its page is editable again.
#     mode auto (default): after the changed step is applied, Claude runs pipeline/replay.R, which redoes the later steps with the
#     earlier settings. mode step: the analyst goes through the later steps again (each previewed with the earlier settings).
#     --from_page: the analyst pressed "Change this step" on the done page (choice "revisit", field mode).
#   Rscript pipeline/revisit.R --db=<db> --cancel [--from_page]
#     "Keep what I had": restores the run from the backup (undoes the change and anything redone since) and drops the backup.
#     Works while the change is open, and after the automatic redo until the next change (the backup is kept until then).
# One revisit at a time. Load can't be revisited: start a new run. A new revisit drops the previous (finished) revisit's backup.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
root <- normalizePath(file.path(here, ".."))
source(file.path(here, "common.R"))
source(file.path(here, "review.R"))
source(file.path(here, "redo.R"))
args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app("Shiny_File.R")
rv <- revisit_info(run_dir)
rv_open <- !is.null(rv) && identical(rv$status, "open")
from_page <- isTRUE(args$from_page) || identical(args$from_page, "TRUE")
copy_dir <- function(from, to) { if (dir.exists(from)) { dir.create(to, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(from, full.names = TRUE, all.files = TRUE, no.. = TRUE), to, recursive = TRUE, copy.date = TRUE)
  if (!all(ok)) stop("Copy failed: ", from) } }

if (isTRUE(args$cancel) || identical(args$cancel, "TRUE")) {
  if (is.null(rv) || !rv$status %in% c("open", "done") || isTRUE(rv$backup_removed))
    stop("Nothing to undo: no revisit with a backup (", file.path(run_dir, "revisit.json"), ").", call. = FALSE)
  if (identical(rv$status, "done") && length(Filter(function(r) identical(r$step, "10_stats") && isTRUE(r$settings$accepted) &&
      (r$time %||% "") > (rv$replay$time %||% rv$finished %||% ""), read_decisions(run_dir))))
    cat("Note: the redone Stats results were already accepted; undoing restores the run as it was before the change.\n")
  if (from_page) { pend <- file.path(run_dir, "review", "pending_decision.json")
    pst <- if (file.exists(pend)) fromJSON(pend)$step else rv$step
    take_page_decision(run_dir, pst, "cancel_revisit") }
  b <- rv$backup
  if (!dir.exists(b)) stop("The backup is missing: ", b, call. = FALSE)
  cat("Restoring the run from", b, "\n")
  file.copy(file.path(b, basename(db_path)), db_path, overwrite = TRUE, copy.date = TRUE)
  for (d in c("steps", "review", "QC")) { unlink(file.path(run_dir, d), recursive = TRUE); copy_dir(file.path(b, d), file.path(run_dir, d)) }
  file.copy(file.path(b, "decisions.jsonl"), file.path(run_dir, "decisions.jsonl"), overwrite = TRUE)
  rv$status <- "cancelled"; rv$cancelled <- format(Sys.time(), "%Y-%m-%d %H:%M:%S"); rv$backup_removed <- TRUE
  write_revisit(run_dir, rv)
  log_revisit_event(run_dir, list(revisit_cancelled = TRUE, step = rv$step), note = "Keep what I had: run restored from the backup",
                    analyst = args$analyst)
  unlink(b, recursive = TRUE)
  write_review_index(run_dir)
  cat("Restored: the run is as it was before", rv$step, "was opened for changes. Backup removed.\n")
  quit(save = "no")
}

need(args, "step")
step <- args$step
if (step == "01_load") stop("Load can't be changed in place: start a new run with the corrected inputs.", call. = FALSE)
if (!step %in% STEP_ORDER) stop("--step must be one of: ", paste(STEP_ORDER[-1], collapse = ", "), call. = FALSE)
if (rv_open && !identical(rv$step, step))
  stop(sprintf("%s is already being changed. Finish it (apply it, then the later steps) or undo it (--cancel) first.", rv$step), call. = FALSE)
dec <- read_decisions(run_dir)
if (!step_was_applied(dec, step)) stop(step, " was not applied yet: change it on its page as usual.", call. = FALSE)
mode <- args$mode %||% "auto"
analyst <- args$analyst
if (from_page) {
  pd <- take_page_decision(run_dir, step, "revisit")
  m <- pd$fields$mode %||% ""
  mode <- if (grepl("go through", m)) "step" else "auto"
  analyst <- analyst %||% (if (nzchar(pd$analyst %||% "")) pd$analyst)
}
if (!mode %in% c("auto", "step")) stop("--mode must be auto or step")

if (!rv_open) {
  # The step the analyst had reached: the last step with a summary or a decision.
  reached <- STEP_ORDER[max(which(vapply(STEP_ORDER, function(s) nzchar(step_time(run_dir, s, dec)), TRUE)))]
  # Keep one backup per run: drop the previous (finished) revisit's backup.
  if (!is.null(rv) && !isTRUE(rv$backup_removed) && dir.exists(rv$backup %||% "")) { unlink(rv$backup, recursive = TRUE); cat("Removed the previous revisit's backup:", rv$backup, "\n") }
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  b <- file.path(run_dir, "Backup", paste0("revisit_", ts, "_", step))
  dir.create(b, recursive = TRUE)
  cat("Backing up the run to", b, "...\n")
  if (!file.copy(db_path, file.path(b, basename(db_path)), copy.date = TRUE)) stop("Backup of the .db failed (disk space?)")
  for (d in c("steps", "review", "QC")) copy_dir(file.path(run_dir, d), file.path(b, d))
  file.copy(file.path(run_dir, "decisions.jsonl"), file.path(b, "decisions.jsonl"))
  rv <- list(step = step, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), backup = normalizePath(b), mode = mode,
             analyst = analyst %||% Sys.getenv("USER"), status = "open", reached = reached)
  write_revisit(run_dir, rv)
  log_revisit_event(run_dir, list(revisit = TRUE, step = step, mode = mode, reached = reached),
                    note = paste("opened", step, "for changes; backup", b), analyst = analyst)
  cat(sprintf("Backup: %s (%.0f MB)\n", b, sum(file.info(list.files(b, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE) / 1e6))
} else { rv$mode <- mode; write_revisit(run_dir, rv) }

# Re-run the step's preview with the applied settings (from the backup), so its page shows them and is editable.
b_db <- file.path(rv$backup, basename(db_path))
cmd <- redo_cmds(step, db_path, b_db, read_decisions(rv$backup), rv$backup)[[1]]
status <- run_redo_cmd(cmd, root)
if (status != 0) stop("The preview of ", step, " failed (see above). The run is unchanged otherwise; --cancel restores it.", call. = FALSE)
write_review_index(run_dir)
cat(sprintf("\n%s is open for changes (later steps: %s). Its page shows the settings applied before.\n", step,
            if (mode == "auto") "redone automatically after you apply it" else "you go through them again"))
cat("Review page:", file.path(run_dir, "review", paste0(step, ".html")), "\n")
