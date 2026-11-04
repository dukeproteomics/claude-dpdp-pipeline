---
name: dpdp-normalize
description: Step 6 of the DPDP proteomics pipeline (same as the Shiny app's Normalize tab). Use after dpdp-filter. Builds the subset used to calculate normalization factors (exclude/include by accession, description, name, genes or sequence), applies the chosen normalization method (sl, tmm, sltmm, quantile, loess, vsn, lr, mi, ai, directlfq), and shows the NormData and normalized bar plots. The user approves the settings and reviews every plot.
---

# DPDP Step 6 — Normalize

**Show the user every graph** and the list of excluded proteins; the user decides.

## Settings
| Setting | Choices | App default | 11488 reference |
|---|---|---|---|
| `norm_exclude` + `norm_exclude_grep` | not_used, accession, description, name, genes, sequence + regex | pipeline: description + user-rule regex (below); app: not_used | 11488 Shiny run: **description**, `trypsin\|keratin\|casein` (old, too broad) |
| `norm_include` + `norm_include_grep` | same | not_used | not_used |
| `norm_type` (method) | sl, tmm, sltmm, directlfq, quantile, lr, loess, vsn, mi, ai | — | **sltmm** |

- Exclude/include only changes which precursors are used to **calculate** the factors (`precursor_normdata`).
  All filtered precursors are normalized and kept.
- **Factors from measured values only.** Background exports: the subset is built from `precursor_filter_measured` (background
  blanked), and only `sl`, `tmm`, `sltmm` are allowed (the others transform the whole table). Spike-in designs where half the proteome
  changes (e.g. triple-spike): normalize on the constant proteome with `--norm_include=description --norm_include_grep="org=Homo sapiens$"`.
- `sltmm`: SL is applied at this step; **TMM is applied after imputation** (`Shiny_Impute.R`). So `precursor_norm_sltmm` holds SL-only values.
- `mti` is in the app's dropdown, but `mti_normalize()` doesn't exist. The pipeline refuses it.
- `directlfq` needs a Python path that isn't set up yet.
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
2. "Which normalization method?" (header "Method"). Options: "sltmm" / "sl" / "tmm" / Other (quantile, loess, vsn, lr).
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
