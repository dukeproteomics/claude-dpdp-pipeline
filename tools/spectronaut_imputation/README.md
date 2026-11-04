# Spectronaut imputation analysis (11488, 2026-11-04)

Scripts behind the findings in `PLAN.md` ("Spectronaut imputation path"). They compare four Spectronaut exports of the same analysis
(no Spectronaut normalization): no imputation (`..._Report.tsv`), background signal, global imputation and run-wise imputation.

Run from the project root with plain `Rscript` (needs data.table). Inputs: `$SPN_INPUT` (default `data/input`). Intermediate files:
`$SPN_WORK` (default `data/spectronaut_compare/`, git-ignored). The export file names are hard-coded for 11488 in `01` and `06`.

| Script | What it checks | Result for 11488 |
|---|---|---|
| 01_read.R | reads the four exports into aligned matrices -> exports.rds | same 131,975 precursors and run order |
| 02_flags_vs_noimpute.R | EG.IsImputed vs values missing in the no-imputation export | flag TRUE for 72,232 filled + 79,547 unchanged measured values |
| 03_how_global_runwise_fill.R | distribution of global / run-wise fills | uniform on log2, unrelated to the precursor (r = 0.00) |
| 04_fill_range_percentiles.R | fill range vs percentiles of measured values | ~1.5th-9th percentile; experiment-wide (global) or per run (run-wise) |
| 05_background_vs_global_quality.R | background vs global for missing-at-random gaps; SPQC CV | background tracks the precursor (r = 0.69); SPQC CV 52% vs 71% |
| 06_read_with_rt_peptide_protein.R | re-reads with RT, q-value, peptide, protein -> exports2.rds | q/RT blanked where filtered in the no-imputation export |
| 07_flag_group_spread.R | does the flag spread within a peptide or protein? | no |
| 08_flag_by_q_intensity.R | flag rate by run q-value, intensity, missingness | rises with q (2% -> 21%), falls with intensity (40% -> 0.1%) |
| 09_qvalue_rt_identical.R | q-values / RTs across exports | identical in all four |
| 10_protein_level_fdr.R | identified (q <= 0.01) but blanked values | ~15% from whole protein groups blanked in a run (protein-level FDR) |
| 11_flag_vs_own_level.R | flagged values vs their precursor's own level | centred (-0.06 log2) but noisier (IQR +/-0.5 vs +/-0.2) |
| 12_background_gaps.R | gaps background could not fill | 3,956; RT/q known for 3,943 (no usable signal) |
| 13_single_export_rules.R | can ONE export tell measured from filled? (scored against the no-imputation export) | q > 0.01 finds 81% of fills with no false calls; the 13,573 identified-but-blanked fills can't be separated exactly |
