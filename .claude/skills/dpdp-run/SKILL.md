---
name: dpdp-run
description: Start-to-end guide for the DPDP proteomics pipeline (Claude version of the shiny_DPDP app). Use when an analyst wants to process a new project ("process project 11510", "run the proteomics pipeline", "start from my Spectronaut report and sample sheet"), resume a run, or asks what the steps are. Walks the nav-menu steps in order (Load, Parameters, Noise, Filter, Normalize, Impute, Rollup, QC, Stats), building a review page at each step and asking the analyst for every decision.
---

# DPDP pipeline: start to end

The analyst makes every decision. At each step: run the preview, give the analyst the **review page** path (offer `xdg-open <path>`),
summarize what it shows in neutral terms, ask with **AskUserQuestion** (no option marked as recommended unless they ask for an opinion),
then apply with `--decided_by=user --analyst=<name>` and the reason in `--note`. **One decision per answer.** The plugin's hook blocks a
second `--decided_by=user` command if the analyst has not answered in between.

## Setup (once per machine)
- R with the app's packages. On the lab Linux PC they are in `/home/R_Packages`. **Run plain `Rscript`, never `--vanilla`.**
- `DPDP_APP`: path to a clone of the shiny_DPDP app (default `/home/user/Claude_Shiny/shiny_DPDP`). The pipeline calls its functions.
- `DPDP_DATA`: folder for run outputs (default `/home/user/Claude_Shiny/data`). **Data stays local: never commit or upload it.**
- Each run gets `<DPDP_DATA>/<prefix>/` with the `.db` (the same format as the Shiny app, so a run can be opened in Shiny), `QC/` plots,
  `steps/NN_*/summary.json`, `decisions.jsonl` (who decided what and when), and `review/` (one HTML page per step, plus `index.html`).

## Inputs the analyst provides
1. The search report (Spectronaut precursor report `.tsv`, so far).
2. The sample sheet (`.xlsx`) with groups and the QC spike level column.
3. A run prefix (e.g. `11488_110326`).
Ask for anything missing. Ask which spike-ins the project uses if they aren't the lab defaults (casein P02662/P02663 by level, ADH P00330/P00331 constant).

## Steps (skill to follow for each)
| # | Step | Skill | Decision |
|---|---|---|---|
| 1 | Load | dpdp-load | design matches data (accept) |
| 2 | TIC | - | skipped (needs raw files) |
| 3 | Parameters | dpdp-parameters | output level, peptide selection (accept) |
| 4 | Noise | dpdp-noise | none / fixed / dynamic, noise_keep |
| 5 | Filter | dpdp-filter | filter settings |
| 6 | Normalize | dpdp-normalize | exclude rule, method |
| 7 | Impute | dpdp-impute | duke / bottomx, cutoffs |
| 8 | Rollup | dpdp-rollup | sum / iq_maxlfq |
| 9 | QC | dpdp-qc | accept, flag samples |
| 10 | Stats | dpdp-stats | settings + comparisons, then accept results |
| 11-13 | Pathway, Phos, Report | - | not built yet |

To reuse a finished run's settings: most steps take `--settings_from=<that run's .db>`. The analyst still approves each step.

## Resuming
Read `<run>/decisions.jsonl` and `<run>/review/index.html`. Continue from the first step without a recorded decision.
The status badge on each page says whether it is waiting or done.

## Lab rules (from the lab lead)
- Casein: the sample sheet is correct; a mismatch is a sample problem. ADH: judge without normalization.
- Normalization exclude: only added material (digestion trypsin, human contaminant keratins, bovine alpha/beta casein).
- noise_keep: noise values are removed for the filter decision, then put back for the precursors that pass the filter.
- Report known app problems neutrally (each step's skill lists them). Don't change the app.
