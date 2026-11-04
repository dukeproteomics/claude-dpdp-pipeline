# Claude Proteomics Pipeline (from shiny_DPDP)

Read `PLAN.md` first: it has the plan, the status of each step, and the test dataset.

- `shiny_DPDP/` is a clone of the user's Shiny app. Treat it as read-only reference code; the
  pipeline calls its functions through `pipeline/common.R::use_app()`.
- `pipeline/steps/NN_<step>.R`: one CLI script per nav-menu step. Each writes
  `data/<prefix>/steps/NN_<step>/summary.json` and appends to `data/<prefix>/decisions.jsonl`.
- `pipeline/compare.R`: compares tables/params against a Shiny reference `.db`.
- `.claude/skills/dpdp-<step>/SKILL.md`: how to run and read each step, including the user's rules. `dpdp-run`: start to end.
- Every step builds an HTML review page in `data/<prefix>/review/` (`pipeline/reviews/NN_<step>.R`). Give the analyst the path, ask
  the decision with AskUserQuestion, then apply with `--decided_by=user --analyst=<name>`. One decision per answer.
- `tools/build_plugin.sh` builds the shareable plugin in `dist/dpdp-marketplace/`. Edit sources here, not in `dist/`.
- Run R with plain `Rscript`, never `--vanilla` (packages are in `/home/R_Packages`).
- Work one step at a time. Show results, recommend options, and let the user decide before
  running the next step.
- Data stays under `data/`. Never commit it or push it to GitHub.
