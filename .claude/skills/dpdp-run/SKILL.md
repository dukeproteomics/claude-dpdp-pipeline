---
name: dpdp-run
description: Start-to-end guide for the DPDP proteomics pipeline (Claude version of the shiny_DPDP app). Use when an analyst wants to process a new project ("process project 11510", "run the proteomics pipeline", "start from my Spectronaut report and sample sheet"), resume a run, or asks what the steps are. Walks the nav-menu steps in order (Load, Parameters, Noise, Filter, Normalize, Impute, Rollup, QC, Stats), building a review page at each step and asking the analyst for every decision.
---

# DPDP pipeline: start to end

The analyst makes every decision. At each step: run the preview, give the analyst the **review page** path (offer `xdg-open <path>`),
summarize what it shows in neutral terms, let the analyst decide on the page (below) or ask in the chat with **AskUserQuestion**
(no option marked as recommended unless they ask for an opinion), then apply with `--decided_by=user --analyst=<name>` and the reason in `--note`. **One decision per answer.** The plugin's hook blocks a
second `--decided_by=user` command if the analyst has not answered in between.

**Go straight on (user, 2026-11-13):** once a step's decision is applied, start the next step right away: run its preview (or
calculation) and build its page, without asking whether to continue. Tell the analyst in one or two lines what was applied and
that the next page is ready. The analyst goes back to a finished step by clicking it (web app step bar or `review/index.html`).
Stop instead, and say why, only when the applied step reports a stop condition or warning the analyst must see first (e.g.
`sample_id_alert`, every precursor removed, a failed command), or when the analyst asked to pause.

## In-page decisions (default)
The analyst can decide **on the review page** (radio buttons, fields, a note, and **Send to Claude**) or in the chat. Both work.
**Web app sessions (`DPDP_WEB=1` in the environment):** skip steps 1 and 2. The DPDP web app serves the pages and opens each new
one in its left panel by itself; say "the <step> page is ready on the left". Runs go in `$DPDP_DATA`, uploads are in
`$DPDP_DATA/input/<name>/`. Steps 3 to 6 are the same.

1. Once per session, start the server in the background: `python3 pipeline/review_server.py` (http://127.0.0.1:8765, localhost only).
2. Open each page through it: `xdg-open http://127.0.0.1:8765/<prefix>/review/<step>.html`. Opened as a file, the buttons are disabled.
3. Start the wait in the background: `python3 pipeline/wait_decision.py --run=<run_dir> --step=<step>`. Then tell the analyst
   "choose on the page, or reply here" and **end your turn**. Don't call AskUserQuestion at the same time: it blocks the session.
4. When the wait finishes, read its JSON (`choice`, `fields`, `note`, `analyst`) and run the matching command below with
   `--decided_by=user`. The analyst and note come from the page automatically. `log_decision()` archives the page decision to
   `review/decisions_done/` and logs `source: "review page"`; the page then reloads itself.
5. If the analyst replies in the chat instead, stop the background wait and use their answer.
   Requests that no step command records ("other", "plots", "fix", "stop") are set aside by `wait_decision.py` once read (noted in
   `decisions.jsonl` as step `request`), so the next wait doesn't return the same request again.
   **After a restart** (or when resuming), first look for a waiting page decision in the run (`review/pending_decision.json`): a
   decision sent while no wait was running is still there; read it and act on it. (The web app also says so when Claude starts.)
6. A choice that needs discussion ("change": free text) goes back to the chat. Show what you would run, and confirm first.

| Step | choice -> command |
|---|---|
| 01_load | `accept` -> `accept.R --step=01_load`; `fix` -> stop, discuss `fields.problem` |
| any done step | `revisit` -> `revisit.R --step=<step> --from_page`; `cancel_revisit` -> `revisit.R --cancel --from_page`; `accept_rest` -> `replay.R --decided_by=user --step=<step> --from_page` (see Going back) |
| 03_parameters | Options panel: `preview` ("Run with these settings") -> `03_parameters.R --from_page` (re-runs the step; then wait again); `apply` ("Accept these settings") -> `accept.R --step=03_parameters --from_page --decided_by=user`; `other` -> discuss `fields.request` |
| 04_noise | Options panel: `preview` -> `04_noise.R --action=preview --from_page` (then wait again); `apply` -> `04_noise.R --action=apply --from_page --decided_by=user`; `other` -> discuss `fields.request` |
| 05_filter | Options panel: `preview` -> `05_filter.R --action=preview --from_page` (then wait again); `apply` -> `05_filter.R --action=apply --from_page --decided_by=user`; `other` -> discuss `fields.request` |
| 06_normalize | `default_exclude` -> default rule (PTM runs: PTM peptides only); `all_peptides` (PTM runs) -> default exclude + `--norm_include=not_used`; `current_exclude` -> current `norm_exclude_grep`; `custom_exclude` -> preview with `fields.norm_exclude_grep` first, then apply; `preset_adh` / `preset_carboxylases` -> `--norm_preset=adh` / `carboxylases`; `custom_include` -> `--norm_include=accession --norm_include_grep=<fields.norm_include_grep>` (presets and include: re-preview and re-run 06b first). Method = `fields.norm_type` |
| 07_impute | `duke` / `bottomx` / `perseus` / `minprob` / `qrilc` / `limpa` -> `07_impute.R --action=apply --impute_type=<choice>` plus that choice's fields (`--missing_cutoff --bottom_x`, `--perseus_shift --perseus_width`, `--minprob_q`); re-preview first if they changed |
| 08_rollup | `iq_maxlfq_scaled` -> `--rollup_method=iq_maxlfq --maxlfq_scale=TRUE`; `iq_maxlfq` -> `... --maxlfq_scale=FALSE`; `sum` -> `--rollup_method=sum`; `limpa` -> `--rollup_method=limpa --limpa_background=<field>` (the only choice after limpa in Impute) |
| 09_qc | `accept_flag` -> `09_qc.R --action=accept --flag_samples=<fields.flag_samples>`; `accept` -> `--action=accept`; `stop` -> stop |
| 10_stats | Options panel (also on the results page): `preview` -> `10_stats.R --action=preview --from_page` (then wait again); `apply` -> `10_stats.R --action=apply --from_page --decided_by=user`, then 10b and 10c; `other` -> discuss `fields.request` |
| 10_stats (results) | `accept` -> `10_stats.R --action=accept`; `plots` -> `10c --proteins=<fields.proteins>`, then ask again |

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

## Going back to an earlier step (user, 2026-11-13)
Each step works on the result of the one before, so changing a step makes every later step out of date. The run is changed in
place, with a backup. Pages of out-of-date steps show a red **Out of date** banner and the web app's step bar marks them.
1. **Change this step** (done page, choice `revisit`, field `mode`): `Rscript pipeline/revisit.R --db=<db> --step=<step> --from_page`.
   It backs up the run (`<run>/Backup/revisit_<time>_<step>/`), writes `<run>/revisit.json` (mode `auto` = "redo them with my earlier
   settings", the default; `step` = "let me go through them again") and previews the step with the settings applied before. Then
   wait on that page as usual. From the chat: `--mode=auto|step` (ask which). Load can't be changed (new run). One change at a time.
2. The analyst changes and applies the step as usual.
3. Then, without asking again:
   - mode `auto`: `Rscript pipeline/replay.R --db=<db> --decided_by=user` (the hook allows it: the analyst chose this with Change
     this step). It previews and applies each later step with the earlier settings (`--settings_from=<backup .db>`), previews the step
     the analyst had been on, and redoes the Stats graphs/output, then stops for the analyst to accept the (changed) Stats results.
     Tell the analyst what was redone and give the Stats calls before -> after it prints (also on the Stats page).
   - mode `step`: go on as usual, but preview each next step with `--settings_from=<backup .db>` (`revisit.json` `backup`), so its
     page starts from the earlier settings; the analyst applies each one.
4. **Accept the rest** (choice `accept_rest` on a later step's page, step mode): `Rscript pipeline/replay.R --db=<db> --decided_by=user --step=<that step> --from_page`.
5. **Keep what I had** (choice `cancel_revisit`, on the changed step and later pages, and on the Stats results page after an
   automatic redo): `Rscript pipeline/revisit.R --db=<db> --cancel --from_page`. Restores the run exactly as it was and drops the backup.
6. **If the redo stops** ("STOPPED at <step>"): an earlier setting no longer fits (e.g. Rollup limpa after Impute was changed to
   perseus). Run that step's preview without `--settings_from`, explain what didn't fit, and wait for the analyst. After they apply
   it, run `replay.R --db=<db> --decided_by=user` again: it continues with the steps not yet redone.
The backup is kept until the next change (one per run), so the before/after and Keep what I had stay available.

## Resuming
Read `<run>/decisions.jsonl` and `<run>/review/index.html`. Continue from the first step without a recorded decision (run its
preview if its page is missing or out of date, then wait for the analyst).
The status badge on each page says whether it is waiting or done.

## Lab rules (from the lab lead)
- Casein: the sample sheet is correct; a mismatch is a sample problem. ADH: judge without normalization.
- Normalization exclude: only added material (digestion trypsin, human contaminant keratins, bovine alpha/beta casein).
- noise_keep: noise values are removed for the filter decision, then put back for the precursors that pass the filter.
- Report known app problems neutrally (each step's skill lists them). Don't change the app.
