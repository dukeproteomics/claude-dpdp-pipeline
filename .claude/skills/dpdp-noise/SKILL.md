---
name: dpdp-noise
description: Step 4 of the DPDP proteomics pipeline (same as the Shiny app's Noise tab). HUMAN DECISION STEP. Use after dpdp-parameters. Calculates the noise inflection point and plot, shows them to the user, waits for the user to choose None/Fixed/Dynamic, then applies it and makes the noise-filtered histogram.
---

# DPDP Step 4 — Noise (the user decides)

**The user makes this decision from the graph.** Claude runs the calculation, shows the plot and
numbers, describes what is visible, and waits. Don't recommend a noise type unless the user asks
for an opinion. Never run `--action=apply` until the user has stated their choice.

## 1. Calculate inflection
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/04_noise.R --db=data/<prefix>/<prefix>.db --action=inflection [--noise_baseline_value=1]
```
Then read `steps/04_noise/summary.json` and **open `QC/Inflection_Point.png`**.

### What the plot is
Every data point in the dataset (all samples, every 100th value plotted), sorted from highest to lowest log2 intensity.
- **Vertical line ("Dynamic"):** the inflection point the app found where the tail curve gets steep.
  Values below it count as dynamic noise.
- **Horizontal line ("Custom"):** log2 of the fixed baseline (`noise_baseline_value`). Baseline 1 → log2 0, which removes nothing.

### Show the user: the review page
The inflection run also builds **`data/<prefix>/review/04_noise.html`** (one self-contained file; `review/index.html` lists all steps
and decisions). Give the analyst the full path and offer to open it (`xdg-open <path>` on a desktop session).
It holds the options with exact counts, the app's Inflection and Start histogram plots, and a **cutoff explorer**:
a slider over log2 intensity that shows the values below the cutoff per sample, how many precursors have a value below it,
how many have every value below it, and how many would be left with too few values for the filter's minimum-measured rule.
Rebuild it at any time with `Rscript pipeline/review_page.R --db=<db> --step=04_noise`.

In the chat, also give the key numbers: inflection point (raw and log2), total data points, and the dynamic and fixed counts and %.
Describe what is visible neutrally (where the curve bends, whether the histogram has a low-end shoulder or second peak).

### Ask for the decision (multiple choice)
Use the **AskUserQuestion** tool so the analyst picks from buttons (they can always type "Other"). Don't mark an option as recommended.
1. Question "Which noise setting for <prefix>?", header "Noise", options:
   - `none`: keep every value (0 removed)
   - `dynamic`: below the inflection <value> (<count>, <pct>)
   - `fixed`: below a value the analyst gives (now <baseline>, <count>). If they pick this, ask for the value, re-run
     `--action=inflection --noise_baseline_value=<n>`, show the new count and page, and confirm.
2. If fixed or dynamic: a second question "Use noise_keep?" (header "noise_keep"). Options: on (put measured values back for
   precursors that pass the filter) or off (removed values are imputed later).
3. Ask for a one-line reason if they want one (it goes in `--note`).

## 2. Apply the analyst's choice
```bash
Rscript pipeline/steps/04_noise.R --db=<db> --action=apply --noise_type=<none|fixed|dynamic> \
  [--noise_baseline_value=<n>] [--noise_keep=TRUE|FALSE] --decided_by=user --analyst=<name> [--note="<reason>"]
```
This writes `precursor_noise` and `QC/Precursor_NoiseFiltered_Histogram.png`, and logs the choice to `decisions.jsonl` with the
analyst's name (default `$USER`). It also rebuilds the review page, which now shows the decision and the values removed per sample.
Show the analyst the page path again, plus the counts of values removed.

## noise_keep (user's rule, 2026-11-04; pipeline fix)
1. Step 4 removes every value below the noise baseline, so the filter treats noise as missing.
2. Step 5 filters as usual. A precursor that is mostly noise fails the filter and is removed from the analysis.
3. For every precursor that **passes** the filter, Step 5 puts the measured values removed in step 1 back
   (`pipeline/noise_keep.R`). This way they aren't imputed: removing a measured value and imputing noise in its place makes no sense.
   Groups blanked by the misaligned filter (target = group) stay NA.

The app's own `noise_keep` branch in `noise_remove_bg()` cleans the wrong rows (`noise_above >= noise_keep_value`), so the
pipeline doesn't use it, and `noise_keep_value` is not used. Report `noise_keep` from `steps/05_filter/summary.json`:
values put back, per sample, and how many precursors were affected.
Test (11488, dynamic 14.93): 10,034 noise values removed; 920 precursors with noise values were dropped by the filter;
7,454 values put back in 4,739 kept precursors.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db --tables=precursor_noise \
  --params=noise_type,noise_baseline_value,noise_inflection,dynamic_noise_count,custom_noise_count,noise_total,noise_keep,noise_keep_value
```

## Next
Step 5 — Filter.
