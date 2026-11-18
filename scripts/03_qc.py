"""Step 3 - Cell-level quality control for gut mucosal biopsies.

Choices (biological rationale):
  * min genes/cell = 200           : floor for an intact cell; intentionally NOT an adaptive (MAD) lower bound on UMI/genes, because
                                     resting lymphocytes and neutrophils are naturally RNA-poor and the project's target is immune cells.
  * max mitochondrial fraction 25% : higher than the usual 10-20% because absorptive epithelium (enterocytes/colonocytes) is mitochondria-rich;
                                     damaged cells are identified by high-MT *together with* low gene counts (see density figure).
  * doublets = Scrublet per sample : tissue contains very different cell sizes, so an upper UMI cut-off would remove large epithelial /
                                     plasma cells; Scrublet is run on each sample separately (doublets only form within a 10x lane).
  * ambient RNA: deposited matrices contain called cells only, so SoupX/CellBender cannot be applied; flagged as a limitation.
Output: $PROC/02_qc.h5ad (kept cells, raw counts), out/qc_summary_per_sample.csv, out/qc_thresholds.json, out/fig03*.png
"""
import argparse, json, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from common import *
from qc_figures import make_figures
style()
ap = argparse.ArgumentParser()
ap.add_argument("--min-genes", type=int, default=200)
ap.add_argument("--max-pct-mt", type=float, default=25.0)
ap.add_argument("--doublet-rate", type=float, default=0.06)
ap.add_argument("--plots-only", action="store_true", help="re-draw fig03a-c from the saved per-cell QC table and per-sample summary (no recomputation)")
ap.add_argument("--summary", default=f"{OUT}/qc_summary_per_sample.csv")
a = ap.parse_args()
if a.plots_only:
    make_figures(pd.read_csv(f"{PROC}/qc_flags_all_cells.csv.gz", index_col=0), pd.read_csv(a.summary), a.min_genes, a.max_pct_mt, None); sys.exit(0)
adata = sc.read_h5ad(f"{PROC}/01_raw.h5ad")
print("loaded", adata, flush=True)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]\d")
adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDEGMQZ]\d?$")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
adata.obs["log10_umi"] = np.log10(adata.obs["total_counts"])
print("MT genes:", int(adata.var.mt.sum()), "| ribosomal:", int(adata.var.ribo.sum()), "| haemoglobin:", int(adata.var.hb.sum()), flush=True)

# ---- flags
o = adata.obs
o["fail_low_genes"] = o["n_genes_by_counts"] < a.min_genes
o["fail_high_mt"] = o["pct_counts_mt"] > a.max_pct_mt
o["fail_basic"] = o["fail_low_genes"] | o["fail_high_mt"]
o["mt_20_25"] = (o["pct_counts_mt"] > 20) & (o["pct_counts_mt"] <= a.max_pct_mt)        # sensitivity: cells that a 20% rule would remove

# ---- Scrublet per sample on cells that pass the basic filters.
# Scrublet's automatic threshold (valley of the simulated-doublet score histogram) is unreliable in tissue dominated by one cell type: in
# test runs it exceeded 0.6 in large epithelial samples and called ~0 doublets (0.1% of cells, far below the 3-12% expected for these loadings).
# The threshold actually applied is therefore the automatic one capped at 0.25 (Scrublet's usual manual default) and floored at 0.20. Doublets
# distort immune-state structure more than they cost in cells, so the conservative (inclusive) side is chosen. Both thresholds and the per-sample
# score histograms are reported (fig03d); clusters with mixed lineage signatures are screened again in step 5a.
THR_LO, THR_HI = 0.20, 0.25
o["doublet_score"] = np.nan; o["predicted_doublet"] = False
thr_rows, diag = [], {}
for s_id in o["sample_id"].cat.categories:
    idx = np.where((o["sample_id"] == s_id).values & ~o["fail_basic"].values)[0]
    sub = adata[idx].copy(); sub.X = sp.csr_matrix(sub.X)
    auto, how = np.nan, "auto"
    try:
        sc.pp.scrublet(sub, expected_doublet_rate=a.doublet_rate, n_prin_comps=30, random_state=0, verbose=False)
        auto = float(sub.uns["scrublet"].get("threshold", np.nan))
    except Exception as e:
        how = f"auto failed ({type(e).__name__}); fixed {THR_HI}"
        sc.pp.scrublet(sub, expected_doublet_rate=a.doublet_rate, n_prin_comps=30, threshold=THR_HI, random_state=0, verbose=False)
    score = sub.obs["doublet_score"].values.astype(float)
    used = THR_HI if np.isnan(auto) else float(np.clip(auto, THR_LO, THR_HI))
    pred = score >= used
    o.iloc[idx, o.columns.get_loc("doublet_score")] = score
    o.iloc[idx, o.columns.get_loc("predicted_doublet")] = pred
    diag[s_id] = (score, np.asarray(sub.uns["scrublet"].get("doublet_scores_sim", [])), auto, used)
    thr_rows.append(dict(sample_id=s_id, n_scored=len(idx), threshold_auto=auto, threshold_used=used, method=how,
                         doublet_rate=float(pred.mean()), doublet_rate_if_auto=float((score >= auto).mean()) if not np.isnan(auto) else np.nan))
    print(f"scrublet {s_id}: {len(idx)} cells, auto thr {auto:.3f} -> used {used:.3f}, predicted doublets {pred.mean():.3%} (auto would give {thr_rows[-1]['doublet_rate_if_auto']:.3%}) [{how}]", flush=True)
thr = pd.DataFrame(thr_rows)
o["keep"] = ~o["fail_basic"] & ~o["predicted_doublet"]

# ---- per-sample summary
g = o.groupby("sample_id", observed=True)
k = o[o["keep"]].groupby("sample_id", observed=True)
summ = pd.DataFrame({
    "status": g["status"].first(), "site": g["site"].first(), "n_raw": g.size(),
    "n_low_genes": g["fail_low_genes"].sum(), "n_high_mt": g["fail_high_mt"].sum(),
    "n_fail_basic": g["fail_basic"].sum(), "n_doublet": g["predicted_doublet"].sum(), "n_kept": g["keep"].sum(),
    "n_mt_20_25_kept": (o["mt_20_25"] & o["keep"]).groupby(o["sample_id"], observed=True).sum(),
    "median_umi_kept": k["total_counts"].median().round(0),
    "median_genes_kept": k["n_genes_by_counts"].median().round(0),
    "median_pct_mt_kept": k["pct_counts_mt"].median().round(2),
}).reset_index()
summ["pct_kept"] = (100 * summ.n_kept / summ.n_raw).round(1)
summ = summ.merge(thr[["sample_id", "threshold_auto", "threshold_used", "doublet_rate", "doublet_rate_if_auto", "method"]], on="sample_id")
summ["flag"] = np.where(summ.n_kept < 1000, "FEW_CELLS", np.where(summ.pct_kept < 70, "HIGH_LOSS", np.where(summ.doublet_rate > 0.15, "HIGH_DOUBLETS", "")))
summ.to_csv(f"{OUT}/qc_summary_per_sample.csv", index=False)
from scipy.stats import mannwhitneyu
_rows = []
for lab, d_ in [("all samples", summ), ("terminal ileum", summ[summ.site == "TI"]), ("ascending colon", summ[summ.site == "AC"])]:
    a_, b_ = d_[d_.status == "non_inflamed"].pct_kept, d_[d_.status == "inflamed"].pct_kept
    _rows.append(dict(subset=lab, n_non_inflamed=len(a_), n_inflamed=len(b_), median_pct_kept_non_inflamed=a_.median(), median_pct_kept_inflamed=b_.median(),
                      mannwhitney_p=mannwhitneyu(a_, b_).pvalue if len(a_) > 0 and len(b_) > 0 else np.nan))
pd.DataFrame(_rows).round(3).to_csv(f"{OUT}/qc_loss_by_status.csv", index=False)   # does QC remove a different share of cells in inflamed vs non-inflamed biopsies?
json.dump(dict(min_genes=a.min_genes, max_pct_mt=a.max_pct_mt, doublet_rate=a.doublet_rate, n_raw=int(len(o)), n_kept=int(o["keep"].sum()),
               n_low_genes=int(o.fail_low_genes.sum()), n_high_mt=int(o.fail_high_mt.sum()), n_doublet=int(o.predicted_doublet.sum()),
               n_samples_flagged=int((summ.flag != "").sum())), open(f"{OUT}/qc_thresholds.json", "w"), indent=1)
print(summ[["sample_id", "n_raw", "n_kept", "pct_kept", "doublet_rate", "doublet_rate_if_auto", "flag"]].to_string(), flush=True)
print(f"TOTAL raw {len(o)} -> kept {int(o.keep.sum())} ({100*o.keep.mean():.1f}%) | low genes {int(o.fail_low_genes.sum())} | high MT {int(o.fail_high_mt.sum())} | doublets {int(o.predicted_doublet.sum())}", flush=True)
o.drop(columns=[c for c in ["mt_20_25"] if c in o]).to_csv(f"{PROC}/qc_flags_all_cells.csv.gz")

# ---- figures (code in qc_figures.py so they can be regenerated from the saved per-cell table)
make_figures(o, summ, a.min_genes, a.max_pct_mt, diag)

adata = adata[o["keep"].values].copy()
adata.obs = adata.obs.drop(columns=["fail_low_genes", "fail_high_mt", "fail_basic", "predicted_doublet", "keep"])
adata.write_h5ad(f"{PROC}/02_qc.h5ad", compression="lzf")
print("wrote", f"{PROC}/02_qc.h5ad", adata.shape, flush=True)
