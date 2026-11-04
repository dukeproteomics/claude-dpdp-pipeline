# dpdp: Claude-driven DPDP proteomics pipeline

The shiny_DPDP app's nav-menu steps as Claude Code skills and R scripts. Claude runs each step, builds an HTML **review page**
(all graphs, interactive explorers, the options in plain words), asks the analyst for the decision, and records it in
`decisions.jsonl`. The run is stored in the app's own SQLite `.db` format, so a run can be opened in the Shiny app.

## Requirements
- R with the shiny_DPDP packages. Run plain `Rscript`, never `--vanilla`.
- A clone of the shiny_DPDP app. Set `DPDP_APP` to its path (default `/home/user/Claude_Shiny/shiny_DPDP`).
- `DPDP_DATA`: where run folders are written (default `/home/user/Claude_Shiny/data`). Data never leaves this folder.
- python3 (for the decision hook) and a web browser for the review pages.

R packages: `renv/renv.lock` (363 packages, R 4.6.1, Bioconductor 3.23); see `renv/README.md` to set up a new machine.

## Install
```bash
claude plugin marketplace add <path or git URL of the dpdp-marketplace folder>
claude plugin install dpdp@dpdp-marketplace
```
Then, in Claude Code: "process project <name> from <report.tsv> and <samples.xlsx>" (skill `dpdp-run`).

## What's included
- `skills/`: `dpdp-run` (start to end) and one skill per step: load, parameters, noise, filter, normalize, impute, rollup, qc, stats.
- `pipeline/`: `common.R`, `steps/NN_<step>.R`, `review.R` + `reviews/NN_<step>.R` (review pages), `review_page.R`, `accept.R`,
  `compare.R` (compare a run with a Shiny reference `.db`), `noise_keep.R`.
- `hooks/require_answer.py`: blocks a command that applies an analyst decision (`--decided_by=user`) unless the analyst has
  answered since the last one. One decision per answer.

## Not built yet
TIC (needs raw files), Pathway, Phos, Report; other input types (protein, peptide, phospho, Proteome Discoverer). Cluster runs.
