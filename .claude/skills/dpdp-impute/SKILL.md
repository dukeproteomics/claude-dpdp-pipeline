---
name: dpdp-impute
description: Step 7 of the DPDP proteomics pipeline (same as the Shiny app's Impute tab). Use after dpdp-normalize. Builds per-group intensity bins and a seeded random table, shows missing value plots, imputes every normalized table (duke or bottom-X), applies TMM after imputation for sltmm, and shows the imputed bar plots. The user approves the settings and reviews every plot.
---

# DPDP Step 7 — Impute

**Show the user every graph** and the fill plan; the user decides.

## Settings
| Setting | App default | 11488 reference | Meaning |
|---|---|---|---|
| `impute_type` | duke | duke | `duke` or `bottomx` (only these two are implemented in the app's `impute_apply_bg2()`, although the dropdown lists more) |
| `bottom_x` | 2 | 2 | bottom-X percentile range used for bottom-X fills |
| `missing_cutoff` | 50 | 50 | % measured in a group needed for a duke fill |

**Duke rule (per group):** a precursor missing at most `Count × (100 − missing_cutoff)/100` values in the group
gets duke fills (group mean + random × bin SD). Precursors with more missing values get bottom-X fills (random value in the lowest `bottom_x`%).
With 3 samples per group and cutoff 50, a precursor missing 1 of 3 gets duke, and one missing 2 or 3 of 3 gets bottom-X.
The random numbers are seeded (`set.seed(123)`), so results are reproducible.

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/07_impute.R --db=<db> --action=preview [--settings_from=<db>] [--impute_type=duke --bottom_x=2 --missing_cutoff=50]
Rscript pipeline/steps/07_impute.R --db=<db> --action=apply   [--settings_from=<db>] ... --decided_by=user --analyst=<name> [--note="<reason>"]
```
Apply imputes every table in `norm_type` (e.g. `precursor_norm_sltmm` → `precursor_impute_sltmm`, plus the impute-only table).
For `sltmm`, **TMM is applied here**, after imputation, with factors from **measured values only** (user rule 2026-11-04; the app uses
the imputed table). The summary has both sets of factors (`tmm.factors_measured`, `tmm.factors_app_from_imputed`); report the largest
difference. **Background exports:** background values are kept and only empty values are imputed; `precursor_missing` /
`protein_missing` count background as not detected, so the stats missing factor uses identified values.
Duke runs in parallel (`detectCores() - 2` workers). Keep this in mind on shared or cluster nodes.

## Review page and decision
Both actions (re)build **`data/<prefix>/review/07_impute.html`** (rebuild: `Rscript pipeline/review_page.R --db=<db> --step=07_impute`).
`apply` refuses to run without `--decided_by=user`. Give the analyst the page path after preview. The page has:
- **The options in plain words:** duke uses the group mean ± the SD at that intensity, for groups that are mostly measured; bottom-X fills
  everything else at the low end. Fills are seeded, so re-runs give the same values.
- **Fill plan explorer:** duke vs bottom-X counts per group as `missing_cutoff` (slider) and `impute_type` change.
  With 3 samples per group, any cutoff from 34 to 66 gives the same split (duke only for precursors missing 1 of 3).
- **Bottom-X range explorer:** the fill band on the intensity histogram as `bottom_x` changes.
- **Duke SD by intensity**, per group.
- **The app's missing-value plots.**
- **After apply:** a histogram of measured vs duke vs bottom-X values (bottom-X should sit at the low end), filled % per sample,
  the TMM factors applied here for sltmm, and the imputed bar plots.

Ask with **AskUserQuestion** (header "Impute"), no option marked as recommended: "How should <prefix> be imputed?". Options:
"duke, cutoff <mc>%, bottom-X <bx>% (as shown)" / "bottomx only, <bx>%" / "Change values" (ask for missing_cutoff and/or bottom_x,
re-run preview, show the page again).

## Show the user
From `steps/07_impute/summary.json`, and **open every file in `plots`**:
1. **Missing values bar plot** (`missing_bar_plot.png`) and **percent plot** (`missing_percent_plot.png`): per sample.
2. **Fill plan** per group: total missing, how many get duke vs bottom-X fills, and rows by missing count.
   A group with mostly bottom-X fills has many precursors missing in most of its samples. Mention it, because it affects the stats.
3. **Imputed bar plots** (`Precursor_impute_<norm>_barplot.png`; added by the pipeline, not shown in the app):
   For sltmm, totals don't have to be equal after TMM, because TMM aligns the trimmed middle of the distribution.
   If one group stands out, check whether it comes from TMM factors or from imputed intensity before drawing conclusions.
Then ask the user to accept, or change settings and re-run.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db \
  --tables=impute_bin_<G>...,random,precursor_impute_<norm>,precursor_impute_impute,precursor_missing,protein_missing \
  --params=impute_type,bottom_x,missing_cutoff,meta_impute_na
```

## Next
Step 8 — Rollup. For Peptide output, the app runs rollup and QC automatically after imputation; the pipeline doesn't do that yet.
