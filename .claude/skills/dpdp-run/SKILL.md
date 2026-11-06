---
name: dpdp-run
description: Start-to-end guide for the DPDP proteomics pipeline (Claude version of the shiny_DPDP app). Use when an analyst wants to process a new project ("process project 11510", "run the proteomics pipeline", "start from my Spectronaut report and sample sheet"), resume a run, or asks what the steps are. Walks the nav-menu steps in order (Load, Parameters, Noise, Filter, Normalize, Impute, Rollup, QC, Stats), building a review page at each step and asking the analyst for every decision.
---

# DPDP pipeline: start to end

The analyst makes every decision. At each step: run the preview, give the analyst the **review page** path (offer `xdg-open <path>`),
summarize what it shows in neutral terms, let the analyst decide on the page (below) or ask in the chat with **AskUserQuestion**
(no option marked as recommended unless they ask for an opinion), then apply with `--decided_by=user --analyst=<name>` and the reason in `--note`. **One decision per answer.** The plugin's hook blocks a
second `--decided_by=user` command if the analyst has not answered in between.

## In-page decisions (default)
The analyst can decide **on the review page** (radio buttons, fields, a note, and **Send to Claude**) or in the chat. Both work.
1. Once per session, start the server in the background: `python3 pipeline/review_server.py` (http://127.0.0.1:8765, localhost only).
2. Open each page through it: `xdg-open http://127.0.0.1:8765/<prefix>/review/<step>.html`. Opened as a file, the buttons are disabled.
3. Start the wait in the background: `python3 pipeline/wait_decision.py --run=<run_dir> --step=<step>`. Then tell the analyst
   "choose on the page, or reply here" and **end your turn**. Don't call AskUserQuestion at the same time: it blocks the session.
4. When the wait finishes, read its JSON (`choice`, `fields`, `note`, `analyst`) and run the matching command below with
   `--decided_by=user`. The analyst and note come from the page automatically. `log_decision()` archives the page decision to
   `review/decisions_done/` and logs `source: "review page"`; the page then reloads itself.
5. If the analyst replies in the chat instead, stop the background wait and use their answer.
6. A choice that needs discussion ("change": free text) goes back to the chat. Show what you would run, and confirm first.

| Step | choice -> command |
|---|---|
| 01_load | `accept` -> `accept.R --step=01_load`; `fix` -> stop, discuss `fields.problem` |
| 03_parameters | `accept` -> `accept.R --step=03_parameters`; `change` -> re-run 03 with the change, then ask again |
| 04_noise | `none` / `dynamic` / `fixed` -> `04_noise.R --action=apply --noise_type=<choice>` [`--noise_baseline_value`] [`--noise_keep=TRUE` if `fields.noise_keep` is true] |
| 05_filter | `apply` -> `05_filter.R --action=apply` with the previewed settings; `change` -> re-preview with the change |
| 06_normalize | `default_exclude` -> default rule (PTM runs: PTM peptides only); `all_peptides` (PTM runs) -> default exclude + `--norm_include=not_used`; `current_exclude` -> current `norm_exclude_grep`; `custom_exclude` -> preview with `fields.norm_exclude_grep` first, then apply; `preset_adh` / `preset_carboxylases` -> `--norm_preset=adh` / `carboxylases`; `custom_include` -> `--norm_include=accession --norm_include_grep=<fields.norm_include_grep>` (presets and include: re-preview and re-run 06b first). Method = `fields.norm_type` |
| 07_impute | `duke` / `bottomx` / `perseus` / `minprob` / `qrilc` / `limpa` -> `07_impute.R --action=apply --impute_type=<choice>` plus that choice's fields (`--missing_cutoff --bottom_x`, `--perseus_shift --perseus_width`, `--minprob_q`); re-preview first if they changed |
| 08_rollup | `iq_maxlfq_scaled` -> `--rollup_method=iq_maxlfq --maxlfq_scale=TRUE`; `iq_maxlfq` -> `... --maxlfq_scale=FALSE`; `sum` -> `--rollup_method=sum`; `limpa` -> `--rollup_method=limpa --limpa_background=<field>` (the only choice after limpa in Impute) |
| 09_qc | `accept_flag` -> `09_qc.R --action=accept --flag_samples=<fields.flag_samples>`; `accept` -> `--action=accept`; `stop` -> stop |
| 10_stats (settings) | `run` -> `10_stats.R --action=apply`, then 10b and 10c; `run_adjp` -> the same with `--checkbox_filter_adjpval=TRUE`; `change` -> discuss |
| 10_stats (results) | `accept` -> `10_stats.R --action=accept`; `change` -> discuss, then re-run; `plots` -> `10c --proteins=<fields.proteins>`, then ask again |

## Logs
The app's logs go to `<run>/Error/error_<name>.txt` (`run_bg`). To keep a step's own stderr, run it as
`pipeline/run_step.sh pipeline/steps/NN_<step>.R <args>`: stderr still shows, and is appended to `<run>/Error/<step>_stderr.log`.
Never redirect logs into `data/` itself.

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
| 6 | Normalize | dpdp-normalize | exclude rule, method (after the 06b method comparison) |
| 7 | Impute | dpdp-impute | duke / bottomx / perseus / minprob / qrilc, settings |
| 8 | Rollup | dpdp-rollup | sum / iq_maxlfq (Peptide output: sum); limpa after limpa |
| 9 | QC | dpdp-qc | accept, flag samples |
| 10 | Stats | dpdp-stats | settings + comparisons, then accept results |
| 11-13 | Pathway, Phos, Report | - | not built yet |

**Phospho / PTM runs (Peptide output):** the result is per peptide. Precursors roll up to peptides (charge states), Stats tests
peptides, and the protein level isn't used except to name a peptide's protein and for the QC spike checks. Say "peptide" to the analyst.

To reuse a finished run's settings: most steps take `--settings_from=<that run's .db>`. The analyst still approves each step.

## Resuming
Read `<run>/decisions.jsonl` and `<run>/review/index.html`. Continue from the first step without a recorded decision.
The status badge on each page says whether it is waiting or done.

## Lab rules (from the lab lead)
- Casein: the sample sheet is correct; a mismatch is a sample problem. ADH: judge without normalization.
- Normalization exclude: only added material (digestion trypsin, human contaminant keratins, bovine alpha/beta casein).
- noise_keep: noise values are removed for the filter decision, then put back for the precursors that pass the filter.
- Report known app problems neutrally (each step's skill lists them). Don't change the app.
