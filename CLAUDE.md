# Claude Proteomics Pipeline (from shiny_DPDP)

Read `PLAN.md` first: it has the plan, the status of each step, and the test dataset.

- `shiny_DPDP/` is a **git submodule**: the fork `dukeproteomics/shiny_DPDP` (of `gmw12/shiny_DPDP`), pinned to the app commit the
  pipeline is built against. Treat it as read-only reference code; the pipeline calls its functions through
  `pipeline/common.R::use_app()`. Clone with `git clone --recurse-submodules`. Moving to a newer app commit is a deliberate commit
  (`cd shiny_DPDP && git fetch && git checkout <sha>`, then re-check results). App fixes, if the user wants any, go in the fork.
- `pipeline/steps/NN_<step>.R`: one CLI script per nav-menu step. Each writes
  `data/<prefix>/steps/NN_<step>/summary.json` and appends to `data/<prefix>/decisions.jsonl`.
- `pipeline/compare.R`: compares tables/params against a Shiny reference `.db`.
- `.claude/skills/dpdp-<step>/SKILL.md`: how to run and read each step, including the user's rules. `dpdp-run`: start to end.
- Every step builds an HTML review page in `data/<prefix>/review/` (`pipeline/reviews/NN_<step>.R`) with a decision form. Serve the pages
  with `pipeline/review_server.py` (localhost:8765); the analyst clicks **Send to Claude**, and `pipeline/wait_decision.py` (run in the
  background) returns the choice. Or ask in the chat. Then apply with `--decided_by=user`. One decision per answer. See skill `dpdp-run`.
- `tools/build_plugin.sh` builds the shareable plugin in `dist/dpdp-marketplace/`. Edit sources here, not in `dist/`.
- Run R with plain `Rscript`, never `--vanilla` (packages are in `/home/R_Packages`).
- Work one step at a time. Show results, recommend options, and let the user decide before
  running the next step.
- Data stays under `data/`. Never commit it or push it to GitHub.
- Logs go in the run's own `Error/` folder, never loose in `data/`. The app's logs already do (`run_bg` -> `Error/error_<name>.txt`).
  To keep a step's stderr, run it through `pipeline/run_step.sh <step.R> <args>` (shows stderr and appends it to
  `<run>/Error/<step>_stderr.log`); don't redirect with `2> data/...`.

## Git
- Repo: https://github.com/dukeproteomics/claude-dpdp-pipeline (public). Commit as the repo-local identity `dukeproteomics`
  (already set in `.git/config`), never the global identity, which belongs to the app's GitHub account.
- Pushing uses `gh` (logged in as dukeproteomics) through a repo-local credential helper in `.git/config`. The global helper
  (`store`, gmw12) is for the user's other repos: don't change it. If a push gets 403 for gmw12, the local helper is missing.
- `.gitignore` keeps out `data/`, the old app copy `shiny_DPDP-main/`, `dist/`, and data file types.
  Check `git status` before committing; never force-add those.
- Commit after each finished piece of work and push to `origin main` when the user asks to back up.
