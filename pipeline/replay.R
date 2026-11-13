# Redo the later steps with the analyst's earlier settings after a changed step was applied (user, 2026-11-13).
#
#   Rscript pipeline/replay.R --db=<db> --decided_by=user [--analyst=<name>] [--step=<NN_step> --from_page]
#
# Needs an open revisit (pipeline/revisit.R) whose step was applied after it was opened. Redoes, in order, every out-of-date step
# after it up to the step the analyst had reached: preview + apply with the settings from the backup (--settings_from), for
# steps the analyst had applied; preview only for the step they were on. Stats: preview, apply, graphs and output, then it stops
# for the analyst to accept the (changed) results. Stops at the first step that fails (e.g. an earlier setting that no longer fits
# the changed data); that step's page then waits for the analyst.
# The analyst approved this once: mode auto ("redo them with my earlier settings", chosen with Change this step), or
# "Accept the rest" on a later step's page (--step=<that step> --from_page).
# Prints, and saves in revisit.json, the Stats calls before and after.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
root <- normalizePath(file.path(here, ".."))
source(file.path(here, "common.R"))
source(file.path(here, "review.R"))
source(file.path(here, "redo.R"))
args <- parse_args()
need(args, "db")
if (!identical(args$decided_by, "user")) stop("replay needs --decided_by=user (the analyst chose to redo the later steps).")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app("Shiny_File.R")
rv <- revisit_info(run_dir)
if (is.null(rv) || !identical(rv$status, "open")) stop("No open revisit: use pipeline/revisit.R to change a step first.", call. = FALSE)
analyst <- args$analyst %||% rv$analyst
if (isTRUE(args$from_page) || identical(args$from_page, "TRUE")) {
  need(args, "step")
  pd <- take_page_decision(run_dir, args$step, "accept_rest")
  analyst <- args$analyst %||% (if (nzchar(pd$analyst %||% "")) pd$analyst) %||% analyst
  log_revisit_event(run_dir, list(accept_rest = TRUE, step = args$step), note = "Accept the rest: redo the remaining steps with the earlier settings", analyst = analyst)
} else if (!identical(rv$mode, "auto") && is.null(args$step))
  stop("This revisit is in step mode: the analyst goes through the later steps, or presses Accept the rest on one of them.", call. = FALSE)

dec <- read_decisions(run_dir)
s_dec <- Filter(function(r) identical(r$step, rv$step) && (r$time %||% "") >= rv$time, dec)
if (!step_was_applied(s_dec, rv$step))
  stop(rv$step, " has not been applied since it was opened for changes: apply it first.", call. = FALSE)

b_db <- file.path(rv$backup, basename(db_path)); b_dec <- read_decisions(rv$backup)
i0 <- match(rv$step, STEP_ORDER); i1 <- match(rv$reached, STEP_ORDER)
todo <- if (i1 > i0) STEP_ORDER[(i0 + 1):i1] else character(0)
# Redo every later step that hasn't been applied (accepted) since the changed step was applied.
apply_time <- function(s) { r <- Filter(function(x) identical(x$step, s) && (!s %in% c("03_parameters", "09_qc") || isTRUE(x$settings$accepted)) &&
                                          !(s == "10_stats" && isTRUE(x$settings$accepted)), dec)
  max(c("", vapply(r, function(x) x$time %||% "", ""))) }
t_s <- apply_time(rv$step)
todo <- todo[vapply(todo, apply_time, "") < t_s]
note <- sprintf("redone with the earlier settings after %s was changed (revisit %s)", rv$step, rv$time)
extra <- c("--decided_by=user", paste0("--analyst=", analyst %||% Sys.getenv("USER")), paste0("--note=", note))
done <- character(0); stopped <- NULL
for (s in todo) {
  applied <- step_was_applied(b_dec, s)
  cat(sprintf("\n==== %s: %s ====\n", s, if (applied) "redo with the earlier settings" else "preview (the step you were on)"))
  for (cmd in redo_cmds(s, db_path, b_db, b_dec, rv$backup, applied)) {
    if (run_redo_cmd(cmd, root, extra) != 0) { stopped <- list(step = s, script = cmd$script); break }
  }
  if (!is.null(stopped)) break
  done <- c(done, s)
  if (!applied) { stopped <- list(step = s, reason = "this is the step you were on: decide it on its page"); break }
}

# Stats before / after (calls per comparison)
cmp <- function(f) { if (!file.exists(f)) return(NULL); s <- fromJSON(f, simplifyVector = FALSE)
  lapply(s$comparisons %||% list(), function(c) list(comparison = paste0(c$N, "_v_", c$D), features = c$features, up = c$up, down = c$down)) }
ba <- list(before = cmp(file.path(rv$backup, "steps", "10_stats", "summary.json")),
           after = if ("10_stats" %in% done) cmp(file.path(run_dir, "steps", "10_stats", "summary.json")))
rv$replay <- list(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), redone = as.list(done), stopped = stopped, stats = ba)
if (is.null(stopped) || identical(stopped$reason, "this is the step you were on: decide it on its page")) rv$status <- "done"
write_revisit(run_dir, rv)
# Stats page again, now with the before/after (and Keep what I had next to it)
if ("10_stats" %in% done) system2(file.path(root, "pipeline", "run_step.sh"), c(file.path(root, "pipeline", "review_page.R"), shQuote(paste0("--db=", db_path)), "--step=10_stats"), stdout = FALSE)
write_review_index(run_dir)

cat("\n==== Replay ====\nRedone with the earlier settings:", if (length(done)) paste(done, collapse = ", ") else "none", "\n")
if (!is.null(stopped$script)) {
  cat(sprintf("STOPPED at %s (%s failed, see above). The earlier settings may not fit the changed data: decide this step on its page. The revisit stays open (--cancel restores the run).\n", stopped$step, stopped$script))
} else if (!is.null(stopped)) cat(sprintf("Now at %s: %s.\n", stopped$step, stopped$reason))
if (length(ba$before) && length(ba$after)) {
  cat("\nStats calls (before -> after):\n")
  for (k in seq_along(ba$after)) { a <- ba$after[[k]]; b <- Filter(function(x) identical(x$comparison, a$comparison), ba$before)
    b <- if (length(b)) b[[1]] else list(features = NA, up = NA, down = NA)
    cat(sprintf("  %s: tested %s -> %s, up %s -> %s, down %s -> %s\n", a$comparison, b$features, a$features, b$up, a$up, b$down, a$down)) }
  cat("Stats results were redone, not accepted: the analyst reviews them on the Stats page.\n")
}
