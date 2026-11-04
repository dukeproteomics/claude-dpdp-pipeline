---
name: dpdp-filter
description: Step 5 of the DPDP proteomics pipeline (same as the Shiny app's Filter tab). Use after dpdp-noise. Previews the misaligned intensity cutoff on the histogram, then applies the precursor filters (minimum measured, % per group, CV, misaligned, precursor quality, SPQC ratio, duplicates) and shows the filtered bar and box plots. The user approves the settings and reviews every plot.
---

# DPDP Step 5 — Filter

**Show the user every graph this step makes** and let them decide. Present settings and results;
only apply settings the user has approved.

## Settings (app widget → CLI `--name=value`)
| Setting | App default | 11488 reference |
|---|---|---|
| `filter_min_measured_all` (min values across all samples) | 2 | 2 |
| `filter_x_percent` / `filter_x_percent_value` (measured in ≥ X% of a group) | TRUE / 50 | TRUE / 50 |
| `filter_cv` / `filter_cv_group` / `filter_cv_value` | FALSE / SPQC / 99 | off |
| `filter_ptm` (keep only `ptm_grep` sequences) | FALSE | off |
| `checkbox_misaligned` / `misaligned_cutoff` / `misaligned_target` (group, dataset) | FALSE / 50 / dataset | **TRUE / 51 / dataset** |
| `intensity_cutoff_sd` (misaligned cutoff = mean + this × SD) | 0.5 | 0.5 |
| `precursor_quality` / `_sd` / `_intensity` / `_min` | FALSE / 50 / 500 / 300 | **TRUE** / 50 / 500 / 300 |
| `precursor_spqc_ratio` / `_accuracy` / `_intensity` | FALSE / 50 / 500 | **TRUE** / 50 / 500 |

To reuse a previous run's settings: `--settings_from=<that run's .db>`. Explicit `--name=value` args override it.

## Run
```bash
cd /home/user/Claude_Shiny
# 1) preview: saves settings, redraws histogram with the misaligned intensity cutoff
Rscript pipeline/steps/05_filter.R --db=<db> --action=preview [--settings_from=<db>] [--name=value ...]
# 2) after the analyst approves:
Rscript pipeline/steps/05_filter.R --db=<db> --action=apply  [--settings_from=<db>] [--name=value ...] \
  --decided_by=user --analyst=<name> [--note="<reason>"]
```
Apply refuses to run without `--decided_by=user`. Both actions (re)build **`data/<prefix>/review/05_filter.html`**
(rebuild any time: `Rscript pipeline/review_page.R --db=<db> --step=05_filter`).

## Review page and decision
Give the analyst the page path after preview (offer `xdg-open <path>`). The page has:
- **Settings table:** each filter, on/off, its values, and what it does.
- **Graphs:** the noise-filtered histogram (misaligned intensity cutoff). After apply, also the start vs filtered bar and box plots.
- **What-if explorers:** one curve per filter (minimum measured, X% per group, misaligned with an intensity-cutoff list, CV with a group
  list, SPQC ratio with an intensity list). Each shows the precursors that filter would remove **on its own** at any setting; the current
  setting is the orange dashed line. Precursor quality has no explorer.
- **After apply:** counts before → after, a removed-per-stage bar chart, missing % per sample, any `noise_keep` values put back, warnings.

Then ask with **AskUserQuestion** (header "Filter"): "Apply these filter settings to <prefix>?". Options:
- "Apply as shown"
- "Change settings": ask which filter and value, re-run preview with `--name=value`, show the page again, ask again.
Don't mark an option as recommended. Ask for an optional one-line reason for `--note`.

## Show the user
From `steps/05_filter/summary.json`, and **open every file in `plots`**:
1. **Histogram** (`Precursor_NoiseFiltered_Histogram.png`): green line = misaligned intensity cutoff; `total_na_by_group` and `total_misaligned_by_group`.
2. **Filter cascade:** proteins, peptides and precursors removed at each stage, then the final counts.
   Use `counts_final`, not `app_reported_counts`. The app stores its counts before the SPQC ratio and duplicate stages.
3. **Bar plot** (`Precursor_Filter_barplot.png`) and **box plot** (`Precursor_Filter_boxplot.png`): compare with the Step 3 start plots.
4. **Missing values per sample after filtering** (`missing_after_filter_pct`): feeds the Impute step.
Then ask the user to accept, or change settings and re-run.

## Known app behavior (reproduced as is, so results match Shiny)
- The "reapply minimum" stage doesn't reapply the % per group filter (result of `filter_min_group()` isn't assigned). It always removes 0 for that part.
- `meta_*_filter` params are taken after precursor quality, before SPQC ratio and duplicate removal.
- With `checkbox_misaligned=TRUE` and `misaligned_target=dataset`, **zero misaligned rows would delete all rows** in the app.
  The pipeline stops before this happens.
- The SPQC ratio filter has the same flaw (`df[-bad_final,]` with no failing rows deletes every row). The pipeline stops if the
  filter leaves 0 precursors. The SPQC ratio explorer shows whether the setting would remove anything.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db --tables=precursor_filter,missing_values,missing_values_plots \
  --params=intensity_cutoff,total_na,total_misaligned,meta_precursor_filter,meta_peptide_filter,meta_protein_filter
```

## Next
Step 6 — Normalize.

## Background exports
Background values count as not detected for every filter rule; after the filter they are put back for the precursors that passed
(not where the misaligned filter blanked a group, not below the noise baseline unless noise_keep). Report the summary's `background` counts.
