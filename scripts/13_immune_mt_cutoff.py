"""13_immune_mt_cutoff.py - is a stricter mitochondrial cut-off for IMMUNE cells biologically sound, and what would it cost?

Input : 05_immune.h5ad (immune cells that already passed QC: >=200 genes, <=25% MT, Scrublet; layers['counts'] = raw UMI counts)
Part A (descriptive, 25% baseline set): MT% distribution per immune lineage; cell quality by MT bin (UMI, genes, own-lineage marker score,
        PTPRC detection, epithelial-marker share of counts as an ambient-RNA / doublet indicator); retention by inflammation status.
Part B (sensitivity): re-run the two immune analyses of slides 6-7 at MT caps 25, 20, 15, 10 %:
        (1) immune-lineage shares of immune cells per biopsy (median by status, Mann-Whitney, site-adjusted OLS),
        (2) pseudobulk PCA (same recipe as 05c) for immune cells, T/NK, myeloid, plasma: association of PC1-3 with status (eta2, site-adjusted p),
            and agreement of the PC scores with the 25 % version.
Biopsy is the unit; everything exploratory. Outputs (CSV) in $GSE_OUT. No data are modified.
"""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.formula.api as smf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import PROC, OUT, HVG_EXCLUDE_REGEX, SEX_GENES, LINEAGE_MARKERS, eta2

CAPS = [25, 20, 15, 10]; MIN_CELLS = 20
IMM = ["T_NK", "plasma", "myeloid", "mast", "B"]
path = os.environ.get("GSE_IMMUNE", f"{PROC}/05_immune.h5ad")
a = ad.read_h5ad(path); o = a.obs.copy()
if "counts" in a.layers: counts = a.layers["counts"]
else:
    b = ad.read_h5ad(f"{PROC}/03_norm.h5ad", backed="r"); counts = b[a.obs_names, :].to_memory().layers["counts"]
counts = sp.csr_matrix(counts); genes = np.asarray(a.var_names)
o["sample_id"] = o["sample_id"].astype(str); o["lin"] = o["coarse_lineage"].astype(str)
assert set(o["lin"].unique()) <= set(IMM), set(o["lin"].unique())
print("immune cells:", a.shape, "| biopsies:", o.sample_id.nunique(), "| max MT%:", round(float(o.pct_counts_mt.max()), 2), flush=True)

# ------------------------------------------------ Part A1: MT distribution per lineage
def q(x, p): return float(np.percentile(x, p))
rows = []
for lin in ["all immune"] + IMM:
    d = o if lin == "all immune" else o[o.lin == lin]; m = d.pct_counts_mt.values
    rows.append(dict(lineage=lin, n_cells=len(d), mean_mt=m.mean(), median_mt=q(m, 50), p75=q(m, 75), p90=q(m, 90), p95=q(m, 95),
                     pct_gt10=100 * (m > 10).mean(), pct_gt15=100 * (m > 15).mean(), pct_gt20=100 * (m > 20).mean()))
pd.DataFrame(rows).round(3).to_csv(f"{OUT}/immune_mt_distribution_by_lineage.csv", index=False)
edges = np.arange(0, 26, 1.0); H = []
for lin in IMM:
    for stt in ["non_inflamed", "inflamed"]:
        d = o[(o.lin == lin) & (o.status == stt)]; h, _ = np.histogram(np.clip(d.pct_counts_mt, 0, 24.999), bins=edges)
        for lo, n in zip(edges[:-1], h): H.append(dict(lineage=lin, status=stt, mt_bin_lo=lo, n=int(n), frac=n / max(len(d), 1)))
pd.DataFrame(H).to_csv(f"{OUT}/immune_mt_histogram.csv", index=False)

# ------------------------------------------------ Part A2: cell quality by MT bin
epi_idx = np.where(np.isin(genes, LINEAGE_MARKERS["epithelial"]))[0]; ptprc = np.where(genes == "PTPRC")[0]
tot = np.asarray(counts.sum(1)).ravel()
o["epi_share"] = np.asarray(counts[:, epi_idx].sum(1)).ravel() / np.maximum(tot, 1)
o["ptprc_det"] = (np.asarray(counts[:, ptprc].sum(1)).ravel() > 0) if len(ptprc) else np.nan
o["own_score"] = np.nan
for l in IMM:
    if f"score_{l}" in o.columns: o.loc[o.lin == l, "own_score"] = o.loc[o.lin == l, f"score_{l}"].astype(float)
o["mt_bin"] = pd.cut(o.pct_counts_mt, [-0.01, 5, 10, 15, 20, 25.01], labels=["0-5", "5-10", "10-15", "15-20", "20-25"])
g = o.groupby(["lin", "mt_bin"], observed=True)
bq = g.agg(n=("total_counts", "size"), median_umi=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"), median_ribo_pct=("pct_counts_ribo", "median"),
           mean_epi_marker_share=("epi_share", "mean"), ptprc_detected=("ptprc_det", "mean"), own_marker_score=("own_score", "mean")).reset_index()
bq.round(4).to_csv(f"{OUT}/immune_quality_by_mt_bin.csv", index=False)
allq = o.groupby("mt_bin", observed=True).agg(n=("total_counts", "size"), median_umi=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"),
        mean_epi_marker_share=("epi_share", "mean"), ptprc_detected=("ptprc_det", "mean")).reset_index()
allq.round(4).to_csv(f"{OUT}/immune_quality_by_mt_bin_all.csv", index=False)

# ------------------------------------------------ Part A3: retention by status
bio = o.groupby("sample_id").agg(status=("status", "first"), site=("site", "first"), patient_id=("patient_id", "first"), biologic=("biologic", "first")).reset_index()
rows, perb = [], []
for cap in CAPS:
    k = o.pct_counts_mt <= cap
    n25 = o.groupby("sample_id").size(); nk = o[k].groupby("sample_id").size().reindex(n25.index).fillna(0)
    d = bio.set_index("sample_id").assign(n_base=n25, n_kept=nk); d["pct_kept"] = 100 * d.n_kept / d.n_base; d["cap"] = cap
    for lin in IMM:
        nb = o[o.lin == lin].groupby("sample_id").size().reindex(n25.index).fillna(0); nl = o[k & (o.lin == lin)].groupby("sample_id").size().reindex(n25.index).fillna(0)
        d[f"pct_kept_{lin}"] = 100 * nl / nb.replace(0, np.nan); d[f"n_kept_{lin}"] = nl
    perb.append(d.reset_index())
    a_, b_ = d[d.status == "non_inflamed"].pct_kept, d[d.status == "inflamed"].pct_kept
    r = dict(cap=cap, cells_kept=int(k.sum()), pct_cells_kept=100 * k.mean(), median_pct_kept_non=a_.median(), median_pct_kept_inf=b_.median(),
             mwu_p_status=mannwhitneyu(a_, b_).pvalue, min_cells_per_biopsy=int(d.n_kept.min()), n_biopsies_lt100=int((d.n_kept < 100).sum()), n_biopsies_lt200=int((d.n_kept < 200).sum()))
    for lin in IMM:
        x, y = d[d.status == "non_inflamed"][f"pct_kept_{lin}"].dropna(), d[d.status == "inflamed"][f"pct_kept_{lin}"].dropna()
        r[f"{lin}_median_pct_kept_non"], r[f"{lin}_median_pct_kept_inf"] = x.median(), y.median(); r[f"{lin}_mwu_p"] = mannwhitneyu(x, y).pvalue if len(x) > 2 and len(y) > 2 else np.nan
    rows.append(r)
pd.concat(perb).round(3).to_csv(f"{OUT}/immune_retention_per_biopsy.csv", index=False)
pd.DataFrame(rows).round(4).to_csv(f"{OUT}/immune_retention_by_cap.csv", index=False)
# fraction of cells >15 % by lineage and status (pooled)
f15 = o.assign(gt15=o.pct_counts_mt > 15).groupby(["lin", "status"]).gt15.mean().mul(100).round(2).reset_index(); f15.to_csv(f"{OUT}/immune_pct_gt15_by_lineage_status.csv", index=False)
print(pd.DataFrame(rows)[["cap", "cells_kept", "pct_cells_kept", "median_pct_kept_non", "median_pct_kept_inf", "mwu_p_status", "min_cells_per_biopsy", "n_biopsies_lt100"]].round(3).to_string(), flush=True)

# ------------------------------------------------ Part B1: composition at each cap
rows = []
for cap in CAPS:
    k = (o.pct_counts_mt <= cap).values
    ct = pd.crosstab(o.sample_id[k], o.lin[k]).reindex(columns=IMM, fill_value=0); ct = ct[ct.sum(1) >= MIN_CELLS]
    sh = 100 * ct.div(ct.sum(1), axis=0); m = bio.set_index("sample_id").loc[sh.index]
    for lin in IMM:
        d = pd.DataFrame({"y": sh[lin], "status": m.status, "site": m.site}); x, y = d[d.status == "non_inflamed"].y, d[d.status == "inflamed"].y
        f = smf.ols("y ~ C(status, Treatment('non_inflamed')) + C(site)", data=d).fit()
        rows.append(dict(cap=cap, lineage=lin, n_biopsies=len(d), median_non=x.median(), median_inf=y.median(), mwu_p=mannwhitneyu(x, y).pvalue,
                         siteadj_diff_pp=f.params.filter(like="status").iloc[0], siteadj_p=f.pvalues.filter(like="status").iloc[0]))
cs = pd.DataFrame(rows); cs.round(5).to_csv(f"{OUT}/immune_composition_by_cap.csv", index=False)
print(cs.pivot(index="lineage", columns="cap", values="mwu_p").round(4).to_string(), flush=True)

# ------------------------------------------------ Part B2: pseudobulk PCA at each cap (recipe of 05c)
excl = pd.Series(genes).str.match(HVG_EXCLUDE_REGEX).values | np.isin(genes, SEX_GENES)
sid_all = o.sample_id.values
def pseudobulk(mask):
    sid = sid_all[mask]; ids = pd.Index(sorted(set(sid))); n_per = pd.Series(sid).value_counts().reindex(ids); ids = ids[n_per.values >= MIN_CELLS]
    keep = np.isin(sid, ids); codes = pd.Categorical(sid[keep], categories=ids).codes
    M = sp.csr_matrix((np.ones(keep.sum()), (codes, np.arange(keep.sum()))), shape=(len(ids), keep.sum()))
    return pd.DataFrame((M @ counts[np.where(mask)[0][keep]]).toarray(), index=ids, columns=genes), n_per.reindex(ids)
def pca_samples(pb, n_top=2000):
    cpm = pb.div(pb.sum(1), axis=0) * 1e6; ok = ((cpm >= 1).mean(0) >= 0.25).values & ~excl
    L = np.log2(cpm.loc[:, ok] + 1); top = L.var(0).sort_values(ascending=False).index[:n_top]
    X = L[top].values - L[top].values.mean(0); U, S, Vt = np.linalg.svd(X, full_matrices=False); ev = (S ** 2) / (S ** 2).sum()
    return pd.DataFrame(U[:, :3] * S[:3], index=L.index, columns=["PC1", "PC2", "PC3"]), ev[:3]
def ols_p(y, X, j=1):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ beta; dof = len(y) - X.shape[1]; cov = (r @ r / dof) * np.linalg.inv(X.T @ X)
    return float(beta[j]), float(2 * stats.t.sf(abs(beta[j] / np.sqrt(cov[j, j])), dof))
SETS = {"immune cells": np.ones(len(o), bool), "T/NK cells": (o.lin == "T_NK").values, "myeloid cells": (o.lin == "myeloid").values, "plasma cells": (o.lin == "plasma").values}
BASE_PC = {"immune cells": "PC2", "T/NK cells": "PC2", "myeloid cells": "PC2", "plasma cells": "PC1"}      # PCs used on slide 6 / in the axis analysis (25 % run)
meta = bio.set_index("sample_id"); res, rows = {}, []
for cap in CAPS:
    for name, sm in SETS.items():
        pb, npc = pseudobulk(sm & (o.pct_counts_mt <= cap).values)
        sc_, ev = pca_samples(pb); res[(cap, name)] = sc_; m = meta.loc[sc_.index]
        stv = (m.status == "inflamed").astype(float).values; siv = (m.site == "AC").astype(float).values
        base = res.get((25, name)); common_ids = sc_.index.intersection(base.index)
        for i in range(3):
            pc = f"PC{i+1}"; y = sc_[pc].values; yz = (y - y.mean()) / y.std()
            b, p = ols_p(yz, np.c_[np.ones(len(yz)), stv, siv])
            rho = spearmanr(sc_.loc[common_ids, pc], base.loc[common_ids, BASE_PC[name]]).statistic
            rows.append(dict(cap=cap, set=name, n_biopsies=len(sc_), median_cells_per_biopsy=int(npc.median()), PC=pc, var_pct=round(ev[i] * 100, 1), eta2_status=eta2(y, m.status.values),
                             eta2_site=eta2(y, m.site.values), status_gap_SD_siteadj=b, p_siteadj=p, rho_with_baseline_PC=rho, baseline_PC=BASE_PC[name]))
pb_ = pd.DataFrame(rows); pb_["q_BH_within_cap_set"] = pb_.groupby(["cap", "set"]).p_siteadj.transform(lambda p: stats.false_discovery_control(p.values, method="bh"))
pb_.round(6).to_csv(f"{OUT}/immune_pseudobulk_pca_by_cap.csv", index=False)
best = pb_.loc[pb_.groupby(["cap", "set"]).p_siteadj.idxmin()][["cap", "set", "PC", "var_pct", "eta2_status", "p_siteadj", "rho_with_baseline_PC", "n_biopsies"]]
print(best.round(6).to_string(index=False), flush=True)
print("DONE", flush=True)
