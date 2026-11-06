# dpdp: Claude-driven DPDP proteomics pipeline

The shiny_DPDP app's nav-menu steps as Claude Code skills and R scripts. Claude runs each step, builds an HTML **review page**
(all graphs, interactive explorers, the options in plain words), asks the analyst for the decision, and records it in
`decisions.jsonl`. The run is stored in the app's own SQLite `.db` format, so a run can be opened in the Shiny app.

## Requirements
- R with the shiny_DPDP packages. Run plain `Rscript`, never `--vanilla`.
- A clone of the shiny_DPDP app (`git clone https://github.com/dukeproteomics/shiny_DPDP`).
- Two environment variables, set before starting Claude Code (e.g. in `~/.bashrc`):
  - `DPDP_APP`: path to the shiny_DPDP clone.
  - `DPDP_DATA`: where run folders are written. Data never leaves this folder.

  The defaults (`/home/user/Claude_Shiny/shiny_DPDP` and `/home/user/Claude_Shiny/data`) only fit the machine the plugin was built on.
- python3 (for the decision hook) and a web browser for the review pages.

R packages: `renv/renv.lock` (363 packages, R 4.6.1, Bioconductor 3.23); see `renv/README.md` to set up a new machine.

## Install
```bash
claude plugin marketplace add <path of the dpdp-marketplace folder>
claude plugin install dpdp@dpdp-marketplace
```
Update to a new version: `claude plugin marketplace update dpdp-marketplace`, then `claude plugin update dpdp@dpdp-marketplace`,
then restart Claude Code.

## Use
Start Claude Code in any folder and ask, e.g. "process project 11540 from <report.tsv> and <samples.xlsx>" (skill `dpdp-run`).
Claude runs one step at a time: Load, Parameters, Noise, Filter, Normalize, Impute, Rollup, QC, Stats. After each step it opens the
review page (served on localhost:8765). Choose there and click **Send to Claude**, or answer in the chat. Nothing is applied
without your answer.

## Input and output
- Input: Spectronaut precursor report (`.tsv`) and the sample description Excel file.
- Output level: **Protein**, or **Peptide** for phospho / PTM sets (precursors roll up to peptides, Stats tests peptides).

## What's included
- `skills/`: `dpdp-run` (start to end) and one skill per step: load, parameters, noise, filter, normalize, impute, rollup, qc, stats.
- `pipeline/`: `common.R`, `steps/NN_<step>.R`, `review.R` + `reviews/NN_<step>.R` (review pages), `review_page.R`, `accept.R`,
  `compare.R` (compare a run with a Shiny reference `.db`), `noise_keep.R`.
- `hooks/require_answer.py`: blocks a command that applies an analyst decision (`--decided_by=user`) unless the analyst has
  answered since the last one. One decision per answer.

## Not built yet
Steps: TIC (needs raw files), Pathway, Phos, Report. Inputs: Proteome Discoverer, FragPipe, protein- or peptide-level exports.
TMT and multi-batch sets. Cluster runs.
