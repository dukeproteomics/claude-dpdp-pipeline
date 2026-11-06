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
| `checkbox_misaligned` / `misaligned_cutoff` / `misaligned_target` (group, dataset) | FALSE / 50 / dataset (**pipeline default: TRUE / 50 / dataset**) | **TRUE / 51 / dataset** |
| `intensity_cutoff_sd` (misaligned cutoff = mean + this × SD) | 0.5 | 0.5 |
| `precursor_quality` / `_sd` / `_intensity` / `_min` | FALSE / 50 / 500 / 300 | **TRUE** / 50 / 500 / 300 |
| `precursor_quality_method` (pipeline: `app` or `profile`) | app | - |
| `precursor_spqc_ratio` / `_accuracy` / `_intensity` | FALSE / 50 / 500 | **TRUE** / 50 / 500 |

**Misaligned filter (user decision 2026-11-06):** on by default with `misaligned_target=dataset`. It flags a precursor when a group is
more than 50% missing (not 100%) and the values it has are above the intensity cutoff, then removes that precursor. On the triple-spike it
flagged 0.19% of precursors, mostly wrong lone values, and gave slightly more true and fewer false calls (PLAN.md). `target=group` blanks the
flagged group instead; imputation then treats that group as low/absent and it added false calls, so **limpa is refused after target=group**
(Step 7 stops) and the other imputations warn. Recommend dataset; only use group if the analyst asks. Groups of 2 are never flagged
(1 of 2 missing is exactly 50%); when nothing is flagged, the pipeline skips the filter for that run and says so.

**Precursor quality filter (tested 2026-11-06, PLAN.md):** the app's method (`precursor_quality_method=app`) takes each precursor's
share of itself + the weaker precursors across all samples; missing values shrink that sum, so it removes good precursors of proteins
that change (triple-spike at 50%: 12-15% of E. coli / yeast precursors, more false calls). `profile` (pipeline/filter_quality.R) scores
each precursor by its residual to the protein's median-polish profile (CV %, same cutoff meaning) and removed 1-3% of each species with a
small gain. If the analyst turns the filter on, recommend `profile`; the app method only to match a Shiny run. `profile` runs after the
app's SPQC ratio and duplicate stages.
**Peptide output (phospho):** the quality filter groups by **peptide**, not protein: phosphopeptides of one protein change
independently, so comparing them would remove real site changes (on 11517, grouping by protein would have removed 14,644 precursors;
by peptide it tests the 3,475 precursors of peptides with 3+ charge states and removes 56). The app's method does nothing for
peptides, so the step refuses `precursor_quality_method=app` there.

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
  The pipeline turns the misaligned filter off for that run instead (NOTE in the output and on the page).
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
