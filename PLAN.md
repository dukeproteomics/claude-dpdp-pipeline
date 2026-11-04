# Claude Proteomics Pipeline — Plan

Goal: move the shiny_DPDP workflow (https://github.com/gmw12/shiny_DPDP) to a Claude-driven
pipeline, one nav-menu step at a time. The deterministic R code stays in R; Claude runs
steps, reads the plots and summaries, suggests options, and the user decides.

## Layout
- `shiny_DPDP/` — clone of the Shiny app (source of the R logic; don't modify in place)
- `data/input/` — starting files (data report + sample description). Kept out of git.
- `PLAN.md` — this file; update the status as each step is finished.
- `pipeline/` — `common.R`, `compare.R`, `steps/NN_<step>.R`
- `.claude/skills/dpdp-<step>/SKILL.md` — one skill per step
- `CLAUDE.md` — short orientation for new sessions

## How each step runs
1. Run: the step's CLI calls the existing `*_bg` R function with the chosen options.
2. Show: it saves the same plots as the app (PNG) plus `summary.json` with the key numbers.
3. Recommend: Claude reads the plots and summary and suggests options for the next step.
4. Decide: the user approves or changes them; every choice is logged.
State stays in the app's SQLite `.db` format, so a Claude run can be opened in Shiny.

## How each step is built
1. Read that step's Shiny code (options on screen, compute function, plots).
2. User explains how they choose options from the plots; this goes into the skill.
3. Build the CLI wrapper and `SKILL.md`.
4. Run on the test dataset and check it matches the Shiny output.
5. Move to the next step.

## Environment notes
- R 4.6.1. Packages live in the shared library `/home/R_Packages` (set in `~/.Renviron`
  and `/etc/R/Renviron.site`). Use plain `Rscript`, never `--vanilla`.
- Not installed: `rmotifx`, `PTMphinder` (Phos motif step only); `pca3d` is loaded from
  `shiny_DPDP/Packages/`.
- Duke cluster: later, after the local pipeline works (renv or an Apptainer image built
  from `gregwaitt/dpmsr_proteomics`).

## Test dataset (project 11488)
- Data: `data/input/20261031_081024_11488_Enkbayar_103126_Report.tsv`
  (Spectronaut precursor report, 131,975 precursors, 9 samples, EG.TotalQuantity)
- Design: `data/input/11488_SampleList_110226.xlsx`
  (groups A ×3, B ×3, SPQC ×3; has a QC spike level column)
- Parameters findings: total intensity even (0.91–1.20× median, B_1 highest). Histogram is a single peak
  with no noise shoulder. Group A has more missing values (A total_na 32,753 vs about 19.5–20k)
  with normal intensity, so fewer IDs rather than low load.
- Path: precursor → protein
- Load check (2026-11-03): detected SP / precursor / short format, no PTM → Protein output.
  All 9 design IDs match the data columns in Raw_Order. 7,245 protein groups (7,218 mouse;
  a few human, bovine, yeast, pig: contaminants or spike-in?). 99,941 peptides; no duplicate
  precursor IDs. Missing values: A_2 and A_3 about 12k vs about 6.5–9k for the other samples.
  Look at this in TIC and QC.
- Note: the app's ID check strips `_N` suffixes, so SPQC_1/2/3 all match `154038` and can't
  catch an SPQC order swap.
- Shiny reference run: `data/reference_shiny/11488_110226.db` (unzipped from `data/input/11488_110226.zip`).
  Full run through Stats (B_v_A). Its `params` table records the user's choices at each step.
- Claude run: `data/11488_110326/` (prefix `11488_110326`)
- User confirmed: yeast = spike-in; SPQC = pooled QC (order within SPQC doesn't matter).
- Spike-ins: bovine casein P02662/P02663 at design `QC Spike Level` (varies); yeast ADH P00330/P00331 at the same level in all samples.
- QC finding: casein is high only in A_2 (about 20x the others); B_1/B_3 (level 2) are at background; SPQC is about A_2/5.
  ADH1 CV 4.3% before normalization, 10.2% after SL+TMM (B_1 lowest).
- Stats B_v_A: 221 Up/Down (133/88) by raw p (checkbox_filter_adjpval off); 0 pass FDR (min adj p 0.104).
  57 of the 221 have |FC| >= 10, mostly detected in B and imputed in A. PCA: PC1 separates A from B; SPQC tight in the
  centre; A_2 sits next to A_3 (not a proteome outlier); B samples spread out (B_3 near the centre/SPQC).

## Steps (precursor → protein path first)
| # | Step | Status |
|---|------|--------|
| 1 | Load | done 2026-11-03: `pipeline/steps/01_load.R`, skill `dpdp-load`; matches Shiny ref |
| 2 | TIC | skipped for now (needs raw files; the reference run didn't use it) |
| 3 | Parameters | done 2026-11-03: `03_parameters.R`, skill `dpdp-parameters`; matches Shiny ref (except random colors) |
| 4 | Noise | done 2026-11-03: `04_noise.R`, skill `dpdp-noise`. **Human decision step.** User chose none; matches Shiny ref |
| 5 | Filter | done 2026-11-03: `05_filter.R`, skill `dpdp-filter`; ref settings via `--settings_from`; matches Shiny ref (127,756 precursors / 7,155 proteins) |
| 6 | Normalize | done 2026-11-03: `06_normalize.R`, skill `dpdp-normalize`; sltmm (SL now, TMM after impute); matches Shiny ref |
| 7 | Impute | done 2026-11-03: `07_impute.R`, skill `dpdp-impute`; duke, bottom_x 2; TMM applied here for sltmm; matches Shiny ref |
| 8 | Rollup | done 2026-11-03: `08_rollup.R`, skill `dpdp-rollup`; iq_maxlfq + scale; also runs QC CV stats like the app; matches Shiny ref |
| 9 | QC | done 2026-11-03: `09_qc.R`, skill `dpdp-qc`; casein spike by level + ADH per sample. casein deviates from the design (sample problem, per user) |
| 10 | Stats (Setup, Comparisons, Graphs, Data, Protein Plots) | Setup/Comparisons/Graphs done 2026-11-03: `10_stats.R`, `10b_stats_graphs.R`, skill `dpdp-stats`; matches Shiny ref. Excel + Protein Plots done: `10c_stats_output.R` (own protein plot prep, app version broken for group 'A') |
| 11 | Pathway (Wiki, GO Profile/Analysis/Volcano, STRING, STRING Enrich) | skipped for now (user choice 2026-11-03) |
| 12 | Phos (Setup, Motif), PTM datasets only | not applicable to 11488 (no PTM) |
| 13 | Report | skipped for now (user choice 2026-11-04); requirements kept below for later |

Later: other data types (protein, peptide, phospho), then the Duke cluster.

## Reference run choices (from its params table), for later steps
noise_type=none · filter: min_measured_all=2, x_percent 50, misaligned 51 (dataset), precursor_quality on,
spqc_ratio on · norm_exclude description `trypsin|keratin|casein` · norm_type=sltmm · impute=duke ·
rollup=iq_maxlfq (maxlfq_scale on) · stats: B_v_A, fdr adj p<0.05, FC 1.5, missing_factor 0.6, peptide_refilter on

## Packaging for other analysts (agreed 2026-11-04)
- Claude Code plugin (skills + step scripts + review template), shared from a private repo; app commit and R packages pinned
  (renv); the same setup later becomes the Apptainer image for the cluster.
- Each step builds a self-contained review page `data/<prefix>/review/<step>.html` (graphs embedded, interactive parts inline),
  plus `review/index.html` with all pages and the decision log. Claude asks each decision as a multiple-choice question
  (AskUserQuestion); the decision is logged to `decisions.jsonl` with `--analyst` and `--note`.
- **Done 2026-11-04:** review pages for every built step (01, 03-10): `pipeline/review.R` (shared page code), `pipeline/reviews/NN_<step>.R`,
  `pipeline/review_page.R` (rebuild any page), `pipeline/accept.R` (accept 01/03). Decision steps refuse `apply` without
  `--decided_by=user` and record `--analyst`/`--note`; QC and Stats have `--action=accept`. `08_rollup.R` and `10_stats.R` now
  need `--action`. Skill `dpdp-run` covers start to end.
- **Plugin:** `tools/build_plugin.sh` builds `dist/dpdp-marketplace/` (marketplace + plugin `dpdp`: 10 skills, pipeline, hook) from
  `pipeline/`, `.claude/skills/dpdp-*` and `tools/plugin_src/`. Validated with `claude plugin validate`, and installed + run in a scratch
  HOME. Hook `require_answer.py` blocks a second `--decided_by=user` command without an analyst answer in between.
  The hook is only active where the plugin is installed (not in this project's own settings).
- **Packaging is paused (user, 2026-11-04):** the plugin, hook and lockfile exist, but are parked until the pipeline covers more
  sample/project types and the user has tested it. Don't extend them unprompted.
- Was open for packaging: author/owner name and a LICENSE (needed only for the public plugin directory), where to host the
  marketplace (a private Duke git repo), then the Apptainer image for the cluster.
- **renv lockfile done 2026-11-04:** `tools/renv/renv.lock` (363 packages, R 4.6.1, Bioc 3.23), `tools/renv/snapshot.R` to regenerate,
  `tools/renv/README.md` (restore steps, Ubuntu system libraries). Copied into the plugin as `renv/`. Partial restore tested.
  Not locked: rmotifx, PTMphinder (Phos only, not installed).

## App issues found (not fixed; the app is reference only)
- `table_exists()` uses unqualified `dbConnect` (works only because the app loads RSQLite).
- `create_dir()` deletes an existing folder.
- `parameter_create_plots()`: the wait loop never ends if the bar plot PNG isn't written; `sesion` typo.
- `update_widgets_parameters()`: `params$primary_goup` typo.
- `noise_remove_bg()`: `noise_keep` is **inverted** (confirmed with the user 2026-11-04). Intended: rows with ≥ x values above
  noise keep their measured low values; rows below that are cleaned (all-noise precursors are dropped later). The code cleans
  rows with `noise_above >= noise_keep_value`. **Fixed in the pipeline 2026-11-04:** step 4 removes all noise, then step 5 puts the
  measured values back for precursors that pass the filter (`pipeline/noise_keep.R`). `noise_keep_value` is unused. The app itself is unchanged.
- `filter_data_bg()`: reapply stage doesn't assign `filter_min_group()`; `meta_*_filter` taken before SPQC ratio and dedupe
  (UI shows 128,647 precursors, real count 127,756); misaligned target=dataset with 0 misaligned rows deletes all rows;
  `precursor_spqc_quality()` also deletes all rows when no precursor fails (`df[-bad_final,]`). Pipeline stops on 0 rows.
- Normalize: `mti` is in the dropdown but `mti_normalize()` doesn't exist; `norm_apply_widget_save` saves `norm_ptp` (typo).
  `norm_apply_bg()` calls `ai_normalize()` / `mi_normalize()` with 2 arguments (they need 3), so ai and mi fail; the pipeline refuses them.
  `vsn_normalize()` has `data_out < data.frame(data_out)` (typo for `<-`, harmless).
- Stats Protein Plots: `create_stats_oneprotein_plots_bg` finds sample columns with grep(group name), so group 'A' matches
  'Accession' and it crashes; it also has a leftover `save(..., file="z1b24")` debug line. Excel Table1 Sample Info is always empty.
- Impute: dropdown lists floor/min/avg/knn/lls/... but `impute_apply_bg2()` only implements duke and bottomx.
  Duke SD bins are built from `precursor_filter` (before normalization) but looked up with the normalized group mean; small mismatch.
- Rollup: only sum and iq_maxlfq work. `rollup_median`/`rollup_mean` group without `Name` and then call sum/median on the character
  `Name` column (error); `rollup_mean` is called with 1 argument (needs 2) and uses an undefined `df`; `rollup_median_polish` calls
  `replace_na` without `tidyr::`; `rollup_topn` is passed `params` as `topn_count` and logs the Precursors column. Pipeline refuses them.
- Start histogram: "Misaligned Intensity Cutoff" label is cut off at the right edge.

## Working rules from the user
- Every step shows the user ALL graphs it makes; the human makes the choice (first said for Noise).
- First pass: match the Shiny run's settings and show the plots. Add detailed decision rules later.

## Report step needs (from Shiny_Report*.R)
- Template: `/mnt/DocShare/Proteomics_Report_Templates/base_templates/proteomics_report_protein_SP_v1.docx` (+ .xlsx key map)
- Raw file path (TIC step) gives the project number (`/RawData/<num>/`) and instrument method via `rawrr`
- `get_google_data()`: authenticated Google Sheet (resuspend volume, trypsin, sample amount, order ID)
- `get_pamlims_data()`: PAM LIMS login with stored credentials (customer, PI, project name)
- Report tab inputs: staff roles, report type, searched databases, Duke Box folder
- Don't call Google/LIMS without the user's explicit OK.


## Open items (answered by the user 2026-11-04)
- Casein spike: the description sheet is correct. Mismatches (A_2 high, B_1/B_3 at background) are a sample problem.
- ADH: added by volume before injection, loaded on the EvoSep tip with the sample. Judge it without normalization
  (CV 4.3% un-normalized is the injection check; 10.2% after SL+TMM is expected).
- Norm exclude: only added material: digestion trypsin (pig/cow/other), human contaminant keratins, bovine alpha/beta casein.
  New default regex in the dpdp-normalize skill: 18 proteins / 415 precursors for 11488, against 45 / 1,014 in the Shiny run.
  The current 11488 runs (Shiny and Claude) used the old regex.
- noise_keep: fixed in the pipeline (see App issues). All later test runs use the fix.
- Credentials file: the user says it's fine for now.
- Contaminants: **flag, don't remove** (user 2026-11-04). `pipeline/contaminants.R`; shown on the stats page, in the summary, and as a
  'Flagged proteins' Excel sheet. 11488: 7 of the 221 calls are human keratins, all higher in A (A has more keratin contamination).
  The 11488 Excel hasn't been regenerated yet (that needs 10c, which also remakes the protein plots).
