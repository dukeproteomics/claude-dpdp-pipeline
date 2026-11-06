# Claude Proteomics Pipeline — Plan

Goal: move the shiny_DPDP workflow (https://github.com/gmw12/shiny_DPDP) to a Claude-driven
pipeline, one nav-menu step at a time. The deterministic R code stays in R; Claude runs
steps, reads the plots and summaries, suggests options, and the user decides.

## Layout
- `shiny_DPDP/` — clone of the Shiny app (source of the R logic; don't modify in place)
- `data/input/<project>/` — starting files (data report + sample description), one folder per project. Kept out of git.
  - `10658_triple_spike/` — three-species spike-in (known ratios); designs `_101926` (with SPQC), `_noSPQC`, `_noSPQC_mix`; raw + background exports, old `_MOD` files.
  - `11485_phos/` — phospho, groups 1050 / DMSO / SPQC x3; `_RePrep_Report.tsv` (raw) + `_ImputeBackground_Report.tsv`. Not run yet.
  - `11488/` — main test set (below); four Spectronaut exports (none / background / global / run-wise) + Shiny reference `.zip`.
  - `11517_phos/` — phospho, WT/KO x 30/37 C + SPQC; `_Report.tsv` (raw) + `_preimpute.tsv` (background).
  - `11531/` — groups CreNeg/CrePos x Biotin/NoBiotin + SPQC (CrePos_Biotin x6, others x3); raw + `_Impute_Report.tsv`. Not run yet.
  - `11538/` — groups KD / NTC / SPQC x3; raw + `_preimpute.tsv` (background). Not run yet.
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
- Data: `data/input/11488/20261104_140253_11488_Enkbayar_103126_Report.tsv` (same as the original `20261031_081024_..._Report.tsv`, now removed)
  (Spectronaut precursor report, 131,975 precursors, 9 samples, EG.TotalQuantity)
- Design: `data/input/11488/11488_SampleList_110226.xlsx`
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
- Shiny reference run: `data/reference_shiny/11488_110226.db` (unzipped from `data/input/11488/11488_110226.zip`).
  Full run through Stats (B_v_A). Its `params` table records the user's choices at each step.
- Claude run: `data/11488_110326/` (prefix `11488_110326`; folder deleted 2026-11-05 at the user's request; re-run from the inputs if needed)
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

## Spectronaut imputation path (discussion 2026-11-04, next major piece of work)
- The lab finds Spectronaut's imputation gives better data quality than the app's duke/bottom-X. Plan: use Spectronaut's filled
  values (not replicate them) and track which values were measured vs filled.
- Spectronaut options: **background signal** (signal extracted at the expected RT/m/z in runs where the precursor wasn't identified;
  not imputation, so `EG.IsImputed` is False for these, confirmed by the user), **global** and **run-wise** imputation (random draws from
  low-abundance signals of the experiment / the run; these do set `EG.IsImputed`, per the user).
- Tracking: global/run-wise -> `EG.IsImputed`; background -> per-run `EG.Qvalue` (identified = q <= 0.01). Check the q-value rule once
  against a paired export of the same analysis with no imputation and q-value filtering (filled = present in A, empty in B).
- Order: factors from identified values only, applied to all values, so "impute first, normalize second" is safe. Check whether
  Spectronaut cross-run normalization was on (11488 looks un-normalized: medians 9.53-9.91, totals 0.91-1.20x) to avoid normalizing twice.
- Pipeline changes: Load reads quantity + IsImputed/Qvalue per run and stores a measured/filled mask; Noise mostly redundant
  (q-value is the noise criterion); Filter and Stats (missing factor, Detected_Imputed, one-group-only) use the mask; Impute skipped.
- Then compare on 11488: app imputation vs Spectronaut background / global / run-wise, scored by SPQC CV, ADH CV, casein level
  agreement, one-group-only calls and the p-value histogram.
- **Exports received 2026-11-04** (data/input/11488, no Spectronaut normalization): `..._Report.tsv` (no imputation; identical to the original
  11488 report), `_background_`, `_impute_global_`, `_impute_runwise_`. Same 131,975 precursors and run order in all four.
- Findings:
  - `EG.IsImputed` (revisited after the user's pushback): TRUE for all 72,232 filled values and also for 79,547 measured values whose
    quantity is unchanged. Those flagged-but-measured values are **not random**: flag rate rises with run q-value (2% -> 21%) and falls with
    intensity (40% below log2 4, 0.1% above 12); they sit at their precursor's own level (median -0.06 log2) but are noisier (IQR +/-0.5 vs
    +/-0.2). The flag is set run by run, not per peptide or protein. Best explanation: Spectronaut treats them as not fully meeting its
    quantification/FDR criteria (weak or partial peaks) and flags them, but keeps the measured precursor quantity. Untested guesses: a
    fragment-level fill, or a protein-level use. The q-values and RTs are identical in all four exports.
  - Why identified values (q <= 0.01) are blanked in the no-imputation export (13,573): ~15% are protein-group x run combinations where the
    whole protein group is blanked (1,499 of 2,256 such combinations still had an identified precursor, so likely a run-wise protein FDR
    cutoff). The other ~11,500 need another rule (candidates: precursor PEP cutoff, minimum fragments for quantification).
  - Background can't fill 3,956 gaps although RT and q exist for 3,943 of them: no usable signal at the expected position.
  - **Mask for the pipeline: filled = empty in the no-imputation export (exact). Keep `EG.IsImputed` as a second flag for "weak measured value".**
  - **User requirement (2026-11-04): one Spectronaut analysis, one TSV.** Single-export rules (script 13): `q > 0.01 (or no q)` finds 81% of
    fills with zero measured values misclassified; the other 13,573 fills (identified, q <= 0.01, but blanked by Spectronaut) can't be told
    apart from flagged measured values in a global/run-wise export (2% false calls if the fill range is used). In a background export those
    13,573 hold real extracted signal (median 0.58 log2 below the precursor's level, between measured and q > 0.01 fills).
  - To settle the rest: Spectronaut analysis settings (identification cutoffs incl. PEP and run-wise protein q; quantification filtering,
    minimum fragments, interference correction), per-run EG.PEP, per-run PG.Quantity in the none vs global exports, fragment-level columns
    for a subset.
  - Global: uniform on log2 between the ~1.5th and ~9th percentile of all measured values (4.43-6.58), the same in every run.
    Run-wise: the same rule per run (1.48-1.54% / 8.98-9.19%). Fills are unrelated to the precursor (cor 0.00). **Replicable.**
  - Background: fills 68,289 of 72,245 gaps (3,956 stay empty); median log2 6.3 (measured 9.6); **correlates with the precursor's level in
    its measured runs (r = 0.69)**, about 0.84 log2 below it. It needs the raw chromatograms, so **not replicable**: import it.
  - 46% of gaps (33,022) have both other group replicates measured (missing at random). Global puts these a median 1.95 log2 below the
    group mean (49% more than 4x below); background a median 0.87 below (24%). SPQC precursors missing 1 of 3 injections: median CV after
    filling 52% background vs 71% global/run-wise (fully measured: 15%).
  - The app's bottom-X 2% fills lower (min to 2nd percentile, log2 0-4.7) than Spectronaut (1.5-9%); duke fills missing-at-random gaps
    near the group mean, which Spectronaut's global/run-wise don't do.

## Decision: how to tell measured from filled in a Spectronaut export (user, 2026-11-04)
One Spectronaut analysis, one TSV, with per-run `EG.TotalQuantity`, `EG.Qvalue` and `EG.IsImputed`.
- **Background signal export:** filled = `EG.Qvalue > 0.01` (or no q-value). Values with q <= 0.01 count as measured (this includes the
  identified-but-blanked ones, which hold real extracted signal). Gaps with no value at all are still missing and the pipeline fills them.
- ~~Global or run-wise imputation export: filled = EG.IsImputed~~ **Dropped (user, 2026-11-04):** the flag has false positives (it marks
  79,547 measured, unchanged values in 11488, some strong and well identified), and rule-bending workarounds aren't acceptable.
  **Accepted inputs: background-signal exports, or raw exports (no background, no imputation). For raw input the pipeline imputes, with
  several improved, best-practice methods to choose from (to be designed and tested).**
- No-imputation exports keep working as now (empty = missing; the app's imputation step fills them).
- **Built in Load (2026-11-04):** `pipeline/spectronaut.R::sp_detect_fill()`; Load refuses imputed exports (IsImputed True), reads the type
  from run-wise q-values, needs `--sp_fill=raw` when there are none (background without q refused), writes `precursor_fill` + param
  `sp_fill_mode`. 11488: raw 0 filled / 72,245 empty; background 55,705 filled / 3,956 empty; global and run-wise refused (151,779 flagged).
  Triple-spike `_MOD` exports have no q-values (raw loads with `--sp_fill=raw`; the background one can't be loaded).
- **Steps 3-7 use the mask (built 2026-11-04, `pipeline/background.R`; user rule: background = signal, not identification):**
  Step 3 moves background values out of `precursor_start` (kept in `precursor_start_bg`), so Parameters, Noise and Filter see measured
  values only; Step 5 puts them back for precursors that pass (not where the misaligned filter blanked, not below the noise baseline
  unless noise_keep); Step 6 builds the factor subset from `precursor_filter_measured` and allows only sl/tmm/sltmm on background runs;
  Step 7 imputes only empty values and writes `precursor_missing`/`protein_missing` with background = not detected (stats missing factor).
- **sltmm, all runs:** the app computes TMM after imputation on the imputed table; the pipeline now computes it from measured values
  (Step 7, `tmm_measured_factors`). 11488 raw: factors differ by up to 1.14%; triple-spike background: 3.83%.
- Tested end to end on the new triple-spike background export (10658, 15 runs, no SPQC): all steps ran, measured and background values
  only rescaled, only empty values imputed. Normalize on human only for this design (`--norm_include=description
  --norm_include_grep="org=Homo sapiens$"`): SL on all proteins pushes human E45/E5 to +0.67 log2 because total signal isn't constant.
- **No SPQC:** Stats defaults to no pool (`comp_spqc=0`) when there's no SPQC group; Noise/QC/07b already handle it. A one-sample group
  breaks the app's filter (Shiny_Filter.R:274), so Load now warns STOP. New sheets: `10658_SampleList_noSPQC.xlsx` (true groups, E30
  alone) and `_noSPQC_mix.xlsx` (E20 x2 + E30 x1 as group Mix_E20_E30).

## Imputation benchmark (2026-11-04)
- 11488 (step 07b, no known truth) and the triple-spike known-ratio data (`tools/impute_benchmark/triple_spike.R`, data
  `data/input/10658_triple_spike`, 10658, design after Guzman et al. 2024; the "SPQC" columns in the `_MOD` files are the E20H50Y30 x2 + E30H50Y20 x1 runs, renamed because the app needs an SPQC; dropped).
- Triple-spike, 6 comparisons, adj p <= 0.05 and FC >= 1.5: limpa observed FDR 3.2% mean / 4.8% worst, background 3.5 / 6.2,
  Perseus 3.9 / 9.4, duke+bottom-X 8.6 / 22.7 (most true calls), kNN and group-blind mixed compress ratios (dropped).
- Rerun on the new exports (2026-11-04, `data/10658_triplespike_bench2`, 5 groups, E30 dropped as a one-sample group, 10 comparisons):
  same ranking. Observed FDR mean / worst: limpa 3.1 / 5.1, background 3.5 / 6.8, Perseus 3.7 / 9.0, duke+bottom-X 8.3 / 22.5.
- **Open question (identified-but-blanked):** 84,753 values are empty in the raw export but have q <= 0.01 in the background export
  (Load counts them as measured). Their accuracy is between measured and background (spread vs true level: measured MAD 0.43-0.65 log2,
  these 1.18-1.41, background 1.55-1.91; S/N medians 18.7 / 6.4 / 2.1). Nothing in one background export separates them from
  measured values (IsImputed all False, PG.Qvalue global). Ask the user whether a Spectronaut setting/column explains the blanking.
- **Impute step options built (2026-11-04):** perseus, minprob, qrilc next to the app's duke/bottomx (measured-only fill distribution,
  background never overwritten). limpa not yet: DPC-Quant values + plain limma lose most of its advantage (triple-spike 4.2 / 8.3% vs
  3.1 / 5.1% with dpcDE), so it needs Rollup (dpcQuant) and Stats (dpcDE p-values) changes.
- **limpa built into Impute / Rollup / Stats (2026-11-04, `pipeline/limpa.R`).** Tested end to end on 11488 (raw) and triple-spike
  (raw and background). Pipeline triple-spike E45/E5, E10/E5 observed FDR: raw 4.0 / 5.6%; background export with background left
  missing 4.9 / 5.6% (default `limpa_background=missing`), with background used 6.9 / 9.2%. Remaining gap raw vs background is most
  likely the identified-but-blanked values (open question above). Detected_Imputed: the app adds it in Stats, in the peptide
  refilter path (precursor_refilter_rollup, refilter is on by default), so it is in every normal data return. With refilter off, Stats
  reads it from protein_<norm>_final, which the Rollup step doesn't write, and fails (tested on 11488 with duke + iq_maxlfq). limpa
  runs need refilter off, so the limpa rollup writes the column itself.
- **Peptide refilter removed (user decision 2026-11-04):** Stats never re-rolls a comparison's samples; what the refilter would remove is
  flagged (`pipeline/refilter_flags.R`). Found and handled two app bugs in the no-refilter path (Detected_Imputed missing from
  protein_<norm>_final; CV columns read as samples). The earlier limpa test tables had wrong sample/CV columns (FC and p were right).
- **Phospho / peptide output (built 2026-11-04, user away; test runs zz_phos_bg / zz_phos_raw, decisions logged as code tests).** Data:
  `data/input/11517_phos`, 11517 (Cryptococcus, WT/KO x 30/37 C + SPQC, 148,530 precursors, 93,815 phospho). `_Report.tsv` = raw,
  `_preimpute.tsv` = background (214,754 fills). Load/Parameters (PTM localization, ~7 min), Filter, Normalize, Impute work unchanged;
  Rollup now allows Peptide output (app sum per modified peptide, or limpa per peptide), QC runs on summed proteins, Stats / refilter flags /
  limpa key on the peptide (Accession, Sequence, PTM_Loc, Protein_PTM_Loc, PeptidePosition, ...), 10c plots top peptides.
  Results: app t-test finds 0 adjusted calls (WT_37 vs WT_30 has 17,441 raw p < 0.05, 3x chance); limpa finds 3,410 (background export) /
  3,372 (raw); null splits 0.09-0.34% p < 0.01 (expect 1%). Raw, limpa: WT 37 vs 30 3,372; KO 37 vs 30 2; KO vs WT 37 C 40; 30 C 0.
  Not done: 10b graphs on the raw run (stopped: system low on memory); no phospho benchmark with known truth.
- **PTM runs: phosphopeptides only for normalization and output (user rule 2026-11-05: enriched samples, non-phospho peptides are not
  supposed to be there).** Parameters `ptm` = TRUE makes the defaults: Normalize `norm_include=sequence` + `ptm_grep` (factors from PTM
  precursors; TMM part of sltmm too; non-PTM rows stay in the data), Stats `checkbox_report_ptm` + `peptide_report_grep=ptm_grep`.
  Differs from the app: the app adjusts p over all peptides and then hides the non-PTM rows; the pipeline drops them **before** testing
  (app path and limpa dpcDE), so adj p and limma's prior cover PTM peptides only. Filter's `filter_ptm` stays off (rows are kept).
  zz_phos_raw rerun (raw, limpa, sltmm): 92,775 of 146,982 precursors used for factors; 75,813 phosphopeptides tested (42,079 not).
  Adjusted calls, all peptides -> phospho only: KO vs WT 30 C 0 -> 0; WT 37 vs 30 3,372 -> 5,340 (2,669 up / 2,671 down);
  KO vs WT 37 C 40 -> 358; KO 37 vs 30 2 -> 31. Null splits (30 C, mixed genotype): 0.40-0.49% p < 0.01, 0 calls.
  ADH (non-phospho standard) CV on sltmm rose 14-15% -> 17-18%, expected when factors come from phosphopeptides.
  Bovine caseins (phospho spike, levels 1/1.5/2 by design) are significant in WT 37 vs 30 and KO vs WT 37: flagged, not removed.
- **Live run 11517_110526 (phospho, raw export; 2026-11-05, complete, accepted).** Decisions: noise none, filter as shown
  (146,982 precursors), normalize phospho-only default + sltmm, impute limpa, rollup limpa, QC accepted (no flags), stats run with raw-p
  calls. Raw-p / adjusted calls: KO_30 v WT_30 5,360 / 0; WT_37 v WT_30 19,514 / 5,340; KO_37 v WT_37 11,740 / 358; KO_37 v KO_30
  12,643 / 31. Raw-p calls include many chance hits (~3,700 expected at p < 0.05 among ~74k tests). Caseins significant everywhere (flagged).
  Same session: all review pages got interactive charts (app PNGs collapsed); 10b/10c write chart data (heatmap = 1,000 most variable rows).
- **Impute defaults (user decision 2026-11-04): limpa for raw exports, perseus for background exports.** The form preselects the run's
  default. Rollup and Stats go by the imputation that was applied (Step 7 summary), not the param a preview writes.
- Proposed Impute options (decided above): raw default limpa; background default Spectronaut background + Perseus for
  empty gaps; Perseus / MinProb baselines; duke+bottom-X for comparison with old runs.

## Normalization rework (2026-11-05)
- **User rule:** methods for every kind of dataset (nothing changes ... everything changes); remove a method only if it has no value on
  any dataset type, never because of one benchmark. Two separate choices: what is constant (subset / preset) and which bias is corrected (method).
- `pipeline/normalize.R`: every method is fitted on the subset (measured values) and applied to all values, so all methods work on
  background exports. sl / sltmm still run through the app. Fixed: tmm (app never corrected loading), loess / vsn / quantile (ignored
  the subset; quantile lost the sample names; vsn had a /10 hack), mi -> median, lr -> rlr. New: medianratio, rt and rt_loess
  (need `EG.ApexRT`; 10658 and 11488 exports have it), none, presets `--norm_preset=adh|carboxylases`. Removed: ai, ti (SL duplicates),
  mti (never existed). directlfq: belongs in Rollup, not built (Python package not installed).
- Step 6b `06b_normalize_compare.R` + page: CVs (groups, pool), sample spread by intensity and by RT, fold-change spread vs SL,
  median FC per comparison, ADH CV (information only), known-truth error for species-mix designs.
- Tested on six datasets (code-test runs `data/zz_norm_*`, Steps 1-5 defaults; deleted 2026-11-06):
  - TMM / medianratio / sltmm agree everywhere; plain sl and median leave more intensity bias (11488: 0.08 / 0.045 vs 0.02).
  - loess flattens the intensity bias in every dataset (0.004-0.027) with FC spread 0.98-1.01 and CVs equal to sltmm.
  - rt / rt_loess remove a gradient drift (11488 RT spread 0.059 -> 0.006) and lower CVs (11488 groups 19.6 -> 19.0%, SPQC 14.3 -> 13.5%);
    FC spread 0.93-0.97.
  - vsn has the lowest CVs everywhere but shrinks fold changes (spread 0.36-0.80; triple-spike truth error 0.60 vs 0.08).
  - quantile: no advantage on these datasets (no large cohort available).
  - ADH CV is lowest with none in 5 of 6 (constant spike, loads differ), so it can't rank methods.
  - 11531 (design unknown to the user; mouse, albumin / fibrinogen top, totals vary 7x, within-group CV 57-66% for every method): looks
    like a pulldown from plasma/serum where everything changes. No method fixes the replicate spread; ADH-only factors made it worse (71%).
- **Through Stats (2026-11-06; limpa impute + rollup + dpcDE, test runs zz_norm_11488 and zz_norm_10658, since deleted; 10658 = raw triple-spike export,
  human-only subset, 6 comparisons of the pure groups):** the 06b gains barely change the calls.
  - 11488 B_v_A: sltmm / loess / rt_loess 224 / 225 / 224 raw-p calls, 1 adjusted call each, same p-value distribution.
  - 10658, adj p <= 0.05 and FC >= 1.5, mean over 6 comparisons: true calls 4,527 / 4,526 / 4,518; false calls 151 / 171 / 139;
    observed FDR 4.17 / 4.53 / 4.01% (worst 9.1 / 9.3 / 9.7%, E10 vs E5); human log2 FC MAD 0.17 / 0.17 / 0.16.
  - Reading: the biases these methods remove (0.04-0.06 log2) are far below the FC cutoff (0.585), so they rarely move a call.
    rt_loess is a small, consistent gain when RT is exported; loess is not better than sltmm in Stats on these data.
- Exports without `EG.ApexRT`: 11485, 11517, 11531, 11538 (11531 and 11538 also lack run-wise `EG.Qvalue`). Asked the user to add both to the report schema.

## Misaligned filter test (2026-11-06, triple-spike raw export)
- Runs zz_misal_off / _dataset / _group (deleted after the test): same settings (Noise none, Filter defaults, sltmm human subset, limpa Impute + Rollup, dpcDE,
  6 comparisons); only the misaligned filter differs (cutoff 50%, intensity mean + 0.5 SD = 1,674). Off reproduces the earlier baseline exactly.
- Flags 664 group-rows in 565 precursors (0.19%); dataset target removes them (3 one-precursor human proteins lost, one is trypsin);
  group target blanks the group. Read-only check on 10658_110426: flagged human lone values sit +0.87 log2 (median) above the row's
  level in fully measured groups; 52% are off by more than 2x vs 4% for normal human precursors.
- Mean over 6 comparisons (off / dataset / group): true calls 4,527 / 4,537 / 4,534; false 151.3 / 149.8 / 153.7; observed FDR
  4.17 / 4.11 / 4.26% (worst 9.09 / 8.86 / 9.41%); human MAD 0.17 for all. Protein-level, summed over comparisons vs off:
  dataset +67 / -7 true, +13 / -16 false; group +44 / -5 true, +16 / -8 false.
- Reading: dataset target is a small, consistent gain (more true calls, slightly fewer false); group target adds false calls (blanking the
  group makes limpa treat it as low). App detail: the step-5 re-apply of the X%-per-group rule is not assigned (`filter_min_group(df, ...)`
  without `df <-`, Shiny_Filter.R), so groups the misaligned filter blanks are not re-checked.
- Decision (user, 2026-11-06): misaligned filter on by default, target=dataset. Step 7 refuses limpa after target=group (other imputations
  warn). With nothing flagged (e.g. groups of 2), the filter is skipped for that run instead of the app dropping every row.

## Live run 10658_110426 (triple-spike background export, complete 2026-11-04)
- Load accept; Parameters accept; Noise none; Filter as shown (294,715 precursors / 15,114 proteins; 311,693 background values back);
  Normalize SL, human only; Impute perseus (81,509 empty values); Rollup iq_maxlfq scaled; QC accept, no spikes (the 11488 ADH/casein
  defaults don't apply), no flags; Stats 6 comparisons, calls by adjusted p; results accepted.
- Scored against the true ratios: observed FDR 0.9-2.9%, 22,624 true calls; E. coli FC compressed when E5 is the low group (+2.59 vs
  +3.17): app FC is the ratio of linear means, inflated by noisy background values. Open question for the user: app t-test / ratio-of-means
  vs limma or limpa as the Stats default.

## Live run 11488_110426b (paused 2026-11-04 at Step 6 Normalize; deleted 2026-11-05)
- Paused before the Spectronaut/limpa/phospho changes; run folder deleted at the user's request.

## In-page decisions (option B, built 2026-11-04)
- `pipeline/review_server.py` (localhost only, per-start token, writes only `review/pending_decision.json`), a form on every page
  (`pipeline/review.R`: `decision_spec()` + `DECISION_JS`), `pipeline/wait_decision.py`, `log_decision()` archives the page decision and logs
  `source`. The plugin hook accepts a pending page decision for the same step. Tested end to end with curl + headless Chrome; not yet
  clicked by the user.

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
