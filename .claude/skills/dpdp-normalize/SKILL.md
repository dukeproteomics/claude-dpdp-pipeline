---
name: dpdp-normalize
description: Step 6 of the DPDP proteomics pipeline (same as the Shiny app's Normalize tab). Use after dpdp-filter. Builds the subset used to calculate normalization factors (exclude/include by accession, description, name, genes or sequence), compares every normalization method on the dataset (06b: CVs, intensity and RT bias, fold-change spread), applies the chosen method (sltmm, sl, tmm, median, medianratio, rlr, loess, rt, rt_loess, vsn, quantile, none), and shows the NormData and normalized bar plots. The user approves the settings and reviews every plot.
---

# DPDP Step 6 — Normalize

**Show the user every graph** and the list of excluded proteins; the user decides.

## Settings
| Setting | Choices | App default | 11488 reference |
|---|---|---|---|
| `norm_exclude` + `norm_exclude_grep` | not_used, accession, description, name, genes, sequence + regex | pipeline: description + user-rule regex (below); app: not_used | 11488 Shiny run: **description**, `trypsin\|keratin\|casein` (old, too broad) |
| `norm_include` + `norm_include_grep` | same | not_used; **PTM runs: sequence + `ptm_grep`** (e.g. Phospho) | not_used |
| `norm_type` (method) | sltmm, sl, tmm, median, medianratio, rlr, loess, rt, rt_loess, vsn, quantile, none | — | **sltmm** |
| `norm_preset` | adh (yeast ADH1 spike P00330), carboxylases (endogenous biotinylated carboxylases) | — | — |

- Exclude/include only changes which precursors are used to **calculate** the factors (`precursor_normdata`).
  All filtered precursors are normalized and kept.
- **Factors from measured values only.** Background exports: the subset is built from `precursor_filter_measured` (background
  blanked). Every method is fitted on the subset and applied to all values, background included, so all methods work on background exports. Spike-in designs where half the proteome
  changes (e.g. triple-spike): normalize on the constant proteome with `--norm_include=description --norm_include_grep="org=Homo sapiens$"`.
- **PTM runs (Parameters `ptm` = TRUE; user rule 2026-11-05):** the samples are enriched, so the factors come from the PTM precursors only
  (`--norm_include=sequence --norm_include_grep=<ptm_grep>`, the default), plus the usual exclude. Non-PTM precursors stay in the data,
  are scaled with the same factors, and are left out of Stats. The TMM part of sltmm uses the same subset. Page choice `all_peptides`
  -> `--norm_include=not_used`.
- `sltmm`: SL is applied at this step; **TMM is applied after imputation** (`Shiny_Impute.R`). So `precursor_norm_sltmm` holds SL-only values.

## Two separate choices (user, 2026-11-05: methods for every kind of dataset)
1. **What is constant** = the subset (exclude/include, preset). 2. **Which bias is corrected** = the method. Methods are in
`pipeline/normalize.R` (`NORM_METHODS` has each one's description and when to use it); sl and sltmm still run through the app.

| Method | Corrects | When |
|---|---|---|
| sltmm, sl | one factor per sample (total) | most proteins unchanged (lab default) |
| tmm | one factor, trimmed log ratios | part of the proteome changes, maybe one-sided. **Fixed**: the app divided by the TMM factor only, so loading was never corrected |
| median / medianratio | one factor | median: the common standard; medianratio: robust when samples differ in missingness |
| rlr | intensity bias, a line | small subsets (one spike-in protein) where loess would overfit. Rebuilt from the app's broken `lr` |
| loess | intensity bias, a curve | saturation, different loads or columns. Now fitted on the subset (app: whole table) |
| rt, rt_loess | drift along the gradient (+ intensity) | long runs, spray drift; needs `EG.ApexRT` in the export (10658 and 11488 have it) |
| vsn | variance vs intensity (glog) | cohorts; **shrinks fold changes** (spread vs SL 0.36 on the triple-spike); fitted on the subset |
| quantile | whole distribution | large cohorts of similar samples; keeps the sample names now (the app's version lost them) |
| none | nothing | IP / AP-MS, proximity labeling, secretome: total signal is the biology |

Subsets for designs where everything changes: `--norm_preset=adh` (constant yeast ADH spike), `--norm_preset=carboxylases`
(PC, PCCA/B, MCCC1/2, ACACA/B: endogenously biotinylated, for streptavidin pulldowns), or the bait: `--norm_include=accession --norm_include_grep=<acc>`.
App names removed: `ai`, `ti` (SL with another target), `mti` (no function), `mi` -> `median`, `lr` -> `rlr`; `directlfq` is a protein rollup method (Step 8, not built).

## Method comparison (06b) - run after the preview, before asking
```bash
Rscript pipeline/steps/06b_normalize_compare.R --db=<db> [--methods=sl,loess,...]
```
Fits every available method on the current subset and scores it on measured values: median CV within groups and in the pool, ADH CV
(only when yeast is a spike; information only: a constant spike varies after a correct loading correction, so it is lowest with none), sample spread by intensity and by RT (SD between samples within deciles; 0 = samples agree), fold-change
spread vs SL (1 = same size of differences; below 1 = shrinks them, which also lowers CVs), median log2 FC per comparison, and a
known-truth error when the group names are species mixes (E<n>H<n>Y<n>). Page: `review/06b_normalize_compare.html`. No single score: show
the analyst the table and the trend charts, point out where methods disagree, and re-run it when the subset changes.
- `--settings_from=<db>` copies the exclude/include settings (an old Shiny run brings back `trypsin|keratin|casein`; the script warns) and the first `norm_type` from a previous run.

## Run
```bash
cd /home/user/Claude_Shiny
# 1) preview: builds the normalization subset + NormData bar plot
Rscript pipeline/steps/06_normalize.R --db=<db> --action=preview [--settings_from=<db>] [--norm_exclude=description --norm_exclude_grep='<regex>']
# 2) after the analyst approves:
Rscript pipeline/steps/06_normalize.R --db=<db> --action=apply --norm_type=<method> [--settings_from=<db>] --decided_by=user --analyst=<name> [--note="<reason>"]
```
The first apply also creates `precursor_norm_impute` (impute-only, no normalization), as the app does.

## Review page and decisions
Both actions (re)build **`data/<prefix>/review/06_normalize.html`** (rebuild: `Rscript pipeline/review_page.R --db=<db> --step=06_normalize`).
`apply` refuses to run without `--decided_by=user`; add `--analyst=<name>` and optionally `--note="<reason>"`.
The page shows:
- **The exclude rule in use** and how many proteins and precursors it leaves out. If it isn't the pipeline default rule, a warning
  gives what the default would leave out.
- **Scaling factors per sample** for SL, SL + TMM (estimate), TMM and median, from the normalization subset (a list switches to the
  default-rule subset). The SL + TMM estimate was within about 1% of the real applied factors for 11488.
- **Graphs:** the filtered and NormData bar plots, the normalized bar plot after apply, and log2 box plots per sample, before vs after.
- **The full list of left-out proteins,** each marked current only / default only / both.

In the chat, point out samples whose factors are far from 1 (11488: B_1 about 0.82 to 0.84, so it was loaded heavier), and whether the
methods agree. Then ask with **AskUserQuestion** (two questions, no option marked as recommended):
1. "Which proteins should be left out of the factor calculation?" (header "Exclude"). Options: "Pipeline default (added material only)" /
   "Keep current" / Other (they type a regex).
2. "Which normalization method?" (header "Method"). Options: "sltmm" / "sl" / the method(s) that 06b points to / Other (any method in the table).
If the exclude rule changes, re-run `--action=preview --norm_exclude=description --norm_exclude_grep='<regex>'` and show the page again
before applying.
## Show the user
From `steps/06_normalize/summary.json`, and **open every file in `plots`**:
1. **Normalization subset:** precursors used vs excluded, and `top_excluded_descriptions`.
   Point out any matches the regex probably didn't intend.
   **Exclude rules (user, 2026-11-04):** exclude only what was *added* to the sample, never the sample's own proteins:
   - **trypsin** added for digestion (pig, cow, or another source), not endogenous trypsins or trypsin inhibitors
   - **keratins from human contamination** only (`org=Homo sapiens`), not the sample organism's keratins
   - **bovine alpha-S1, alpha-S2 and beta casein** (the spike), not casein kinases
   Default description regex (grep is case-insensitive):
   `trypsin, org=(sus scrofa|bos taurus)|keratin, type [^,]*, org=homo sapiens|(alpha-s1|alpha-s2|beta)-casein, org=bos taurus`
   For 11488 (mouse) this excludes 18 proteins / 415 precursors, against 45 / 1,014 with `trypsin|keratin|casein`.
   - If the digestion trypsin came from another species, add it to `(sus scrofa|bos taurus)`; check `top_excluded_descriptions`.
   - **Human samples:** `org=homo sapiens` would also exclude the sample's own keratins. Ask the user before using it.
   With `keratin`, intracellular cytokeratins (K8/K18/K7/K19) match along with skin keratins.
2. **NormData bar plot** (`Precursor_NormData_barplot.png`): totals of the subset before normalization.
3. **SL factors** per sample (`sl_factors`): how much each sample is scaled.
4. **Normalized bar plot** (`Precursor_<method>_barplot.png`): totals should be about equal. Report `total_cv_pct_before` and `_after`.
Then ask the user to accept, or change settings and re-run.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db --tables=precursor_normdata,precursor_norm_<method>,precursor_norm_impute \
  --params=norm_exclude,norm_exclude_grep,norm_include,norm_type
```

## Next
Step 7 — Impute.
