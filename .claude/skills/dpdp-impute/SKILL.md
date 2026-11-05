---
name: dpdp-impute
description: Step 7 of the DPDP proteomics pipeline (same as the Shiny app's Impute tab). Use after dpdp-normalize. Builds per-group intensity bins and a seeded random table, shows missing value plots, imputes every normalized table (app: duke or bottom-X; pipeline: perseus, minprob, qrilc), applies TMM after imputation for sltmm, and shows the imputed bar plots. The user approves the settings and reviews every plot.
---

# DPDP Step 7 — Impute

**Show the user every graph** and the fill plan; the user decides.

## Settings
| Setting | Pipeline default | 11488 reference | Meaning |
|---|---|---|---|
| `impute_type` | **limpa** (raw exports), **perseus** (background exports) | duke | app: `duke`, `bottomx` (the only two in `impute_apply_bg2()`); pipeline: `perseus`, `minprob`, `qrilc` |
| `perseus_shift`, `perseus_width` | 1.8, 0.3 | - | Perseus: per sample N(mean - shift x SD, (width x SD)^2) of the measured log2 values |
| `minprob_q` | 0.01 | - | MinProb: quantile of each sample's measured values the fills centre on |
| `bottom_x` | 2 | 2 | bottom-X percentile range used for bottom-X fills |
| `missing_cutoff` | 50 | 50 | % measured in a group needed for a duke fill |

**Duke rule (per group):** a precursor missing at most `Count × (100 − missing_cutoff)/100` values in the group
gets duke fills (group mean + random × bin SD). Precursors with more missing values get bottom-X fills (random value in the lowest `bottom_x`%).
With 3 samples per group and cutoff 50, a precursor missing 1 of 3 gets duke, and one missing 2 or 3 of 3 gets bottom-X.
The random numbers are seeded (`set.seed(123)`), so results are reproducible.

**Pipeline methods** (`pipeline/impute_methods.R`): on log2, per sample, the fill distribution from **measured values only**
(background blanked); only empty values get fills, background values are never overwritten. All three treat every empty value as low.
**Evidence** (`tools/impute_benchmark`, triple-spike, 10 comparisons; false call rate mean / worst): perseus 3.7 / 9.0%,
minprob 4.0 / 10.5%, qrilc 3.7 / 9.8%, duke + bottom-X 8.3 / 22.5% (most true calls, but group-mean fills inflate false calls),
background export + perseus 3.5 / 6.8%. Show these on the page; the analyst chooses. **`limpa`**: no imputation; the gaps stay for DPC-Quant (Rollup) and dpcDE (Stats). Its weights carry the benefit
(4.2 / 8.3% without them vs 3.1 / 5.1%). The page shows the detection curve (limpa::dpc, measured values) for every run.

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
