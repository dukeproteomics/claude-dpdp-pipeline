#!/usr/bin/env python3
"""Kinase Library (Johnson et al., Nature 2023: 303 Ser/Thr kinases; Yaron-Barir et al., Nature 2024: Tyr kinases) for Step 12 Phospho.

Input  (--sites, TSV): comp, site, window (15-mer, phospho-residue in the middle, lower case; '_' padding), label (1 Up, -1 Down,
       0 tested but not called). Output (--out, JSON): per comparison and kinase type (ser_thr for S/T sites, tyrosine for Y sites)
       the enrichment of each kinase's predicted substrates among the Up and among the Down sites against the unchanged sites
       (Fisher's exact test, BH; the package's differential-phosphorylation design), and for every site its top kinases by
       percentile. A site counts as a kinase's substrate when that kinase ranks in its top --rank kinases (percentile rank;
       the package's recommended method, default 15).
Runs locally (matrices ship with the package); no web service.
"""
import argparse, json, sys, warnings, io, contextlib
warnings.filterwarnings("ignore")
import pandas as pd
import kinase_library as kl

ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--rank", type=int, default=15); ap.add_argument("--top", type=int, default=3)
a = ap.parse_args()
sites = pd.read_csv(a.sites, sep="\t", dtype={"comp": str, "site": str, "window": str})
quiet = lambda: contextlib.redirect_stdout(io.StringIO())

def kin_types(df):
    res = df["window"].str[7].str.upper()
    return [("ser_thr", df[res.isin(["S", "T"])]), ("tyrosine", df[res == "Y"])]

def family(k):
    try:
        return kl.get_kinase_family(k)
    except Exception:
        return ""

def table(e, direction):
    sfx = "_" + direction
    cols = ["fg_counts", "fg_total", "bg_counts", "bg_total", "fg_percent", "log2_freq_factor", "fisher_pval", "fisher_adj_pval"]
    t = e[[c + sfx for c in cols]].copy(); t.columns = cols
    t = t.reset_index().rename(columns={"index": "kinase"})
    if "kinase" not in t.columns: t = t.rename(columns={t.columns[0]: "kinase"})
    t["family"] = [family(k) for k in t["kinase"]]
    t = t.sort_values(["fisher_adj_pval", "log2_freq_factor"], ascending=[True, False])
    return t.round({"fg_percent": 2, "log2_freq_factor": 3}).to_dict(orient="records")

out = {"version": getattr(kl, "__version__", "?"), "method": "percentile_rank", "rank": a.rank, "comparisons": {}, "top_kinases": {}}
# every site's top kinases (unique windows, all comparisons)
for kt, d in kin_types(sites.drop_duplicates("window")):
    if not len(d): continue
    with quiet():
        pc = kl.PhosphoProteomics(d[["window"]].rename(columns={"window": "seq"}).reset_index(drop=True), seq_col="seq").percentile(kin_type=kt)
    kcols = [c for c in pc.columns if c not in ("seq", "phos_res", "Sequence")]
    for i, row in pc.iterrows():
        top = row[kcols].astype(float).nlargest(a.top)
        out["top_kinases"][row["seq"].upper()] = [[k, round(float(v), 1)] for k, v in top.items()]
for comp, d0 in sites.groupby("comp"):
    r = {}
    for kt, d in kin_types(d0):
        n_up, n_dn = int((d["label"] > 0).sum()), int((d["label"] < 0).sum())
        entry = {"up_sites": n_up, "down_sites": n_dn, "unchanged_sites": int((d["label"] == 0).sum()), "up": [], "down": []}
        if (n_up >= 3 or n_dn >= 3) and (d["label"] == 0).sum() >= 10:
            try:
                with quiet():
                    dp = kl.DiffPhosData(d.rename(columns={"window": "seq"}), lfc_col="label", lfc_thresh=0.5, seq_col="seq", suppress_warnings=True)
                    res = dp.kinase_enrichment(kin_type=kt, kl_method="percentile_rank", kl_thresh=a.rank)
                e = res.combined_enrichment_results
                if n_up >= 3: entry["up"] = table(e, "upreg")
                if n_dn >= 3: entry["down"] = table(e, "downreg")
            except Exception as ex:
                entry["error"] = str(ex)[:300]
        else:
            entry["note"] = "too few sites (3 Up or Down and 10 unchanged needed)"
        r[kt] = entry
    out["comparisons"][comp] = r
json.dump(out, open(a.out, "w"))
print("kinase library: %d sites, %d comparisons" % (len(sites), len(out["comparisons"])))
