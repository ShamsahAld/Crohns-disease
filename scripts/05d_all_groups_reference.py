"""Step 5d - Reference embedding of ALL tissue states: non-IBD control, non-inflamed CD, adjacent-to-inflamed CD, inflamed CD.

Question: do the four groups separate, and in what order? Two views, because cells from one biopsy are not independent replicates:
  A. cell level   - PCA/UMAP of every QC-passed cell, coloured by group (cell identity is expected to dominate; group effects are subtle);
  B. sample level - pseudobulk PCA of the biopsies (all cells and per lineage), coloured by group, marker = site;
  C. held-out check - an 'inflammation axis' is defined using ONLY inflamed vs non-inflamed biopsies (a PC of their pseudobulk profiles, chosen for its
     association with inflammation after adjusting for site); control and adjacent biopsies play no part in defining it and are projected on to it.
     Inflamed vs non-inflamed separation on this axis is therefore true by construction; where control and adjacent fall is the evidence.
Run with GSE_STATUS_ORDER=control,non_inflamed,adjacent,inflamed and GSE_PROC pointing at the all-group processed folder.
Outputs: out/fig05d_*.png and out/*all_groups*.csv
"""
import argparse, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
from common import *
style()
ap = argparse.ArgumentParser()
ap.add_argument("--n-pcs", type=int, default=30); ap.add_argument("--n-neighbors", type=int, default=30)
ap.add_argument("--leiden-res", type=float, default=0.5); ap.add_argument("--min-dist", type=float, default=0.3)
ap.add_argument("--skip-cells", action="store_true", help="skip the cell-level embedding (part A); reuse the saved lineage labels")
a = ap.parse_args()
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
MIN_CELLS = int(os.environ.get("MIN_CELLS", 20)); MIN_SAMPLES = int(os.environ.get("MIN_SAMPLES", 8)); MIN_TRAIN = int(os.environ.get("MIN_TRAIN", 8))
G = [g for g in STATUS_ORDER]; assert len(G) == 4, f"expected four groups in GSE_STATUS_ORDER, got {G}"
adata = sc.read_h5ad(f"{PROC}/03_norm.h5ad"); print("loaded", adata, flush=True)

# ================================================================== A. cell level
if a.skip_cells:
    _lab = pd.read_csv(f"{PROC}/cell_lineage_all_groups.csv.gz", index_col=0).reindex(adata.obs_names)
    for _c in _lab.columns: adata.obs[_c] = _lab[_c].astype(str).values
else:
    h = adata[:, adata.var["highly_variable"].values].copy(); sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=50, svd_solver="arpack", random_state=0); vr = np.asarray(h.uns["pca"]["variance_ratio"])
    sc.pp.neighbors(h, n_neighbors=a.n_neighbors, n_pcs=a.n_pcs, random_state=0); sc.tl.umap(h, min_dist=a.min_dist, random_state=0)
    sc.tl.leiden(h, resolution=a.leiden_res, flavor="igraph", n_iterations=2, directed=False, random_state=0, key_added="leiden")
    adata.obsm["X_umap"] = h.obsm["X_umap"]; adata.obsm["X_pca"] = h.obsm["X_pca"]; adata.obs["leiden"] = h.obs["leiden"].astype(str).values
    tab = assign_lineage(adata); tab.round(3).to_csv(f"{OUT}/leiden_lineage_scores_all_groups.csv")
    adata.obs[["coarse_lineage", "compartment", "leiden"]].to_csv(f"{PROC}/cell_lineage_all_groups.csv.gz")      # reused by --skip-cells
    o = adata.obs; xy = adata.obsm["X_umap"]
    FACT = ["coarse_lineage", "sample_id", "site", "status", "platform_id", "biologic"]
    et = eta2_table(h.obsm["X_pca"], o, FACT, n_pcs=20); et.to_csv(f"{OUT}/pca_eta2_all_groups_cells.csv")
    print(f"cells {adata.n_obs}; PC1 {vr[0]*100:.1f}%; clusters {o.leiden.nunique()}; lineage counts:", o.coarse_lineage.value_counts().to_dict(), flush=True)
    print("max eta2 over PC1-20:", et.max().round(2).to_dict(), flush=True)
    pd.crosstab(o["coarse_lineage"], o["status"], normalize="columns").mul(100).round(1).to_csv(f"{OUT}/lineage_share_by_group_cells.csv")

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9)); axes = axes.ravel()
    embed_cat(axes[0], xy, o["coarse_lineage"], LINEAGE_COL, s=0.2, label_groups=True, title="Provisional cell lineage")
    embed_cat(axes[1], xy, o["status"], STATUS_COL, s=0.2, legend=True, legend_ncol=2, title="Tissue state (shuffled)")
    for ax_, g in zip(axes[2:], G):
        embed_highlight(ax_, xy, (o["status"] == g).values, STATUS_COL[g], f"{DISPLAY[g][0].upper() + DISPLAY[g][1:]} cells only", s=0.2)
    for ax_, l in zip(axes, "abcdef"): panel_letter(ax_, l)
    fig.subplots_adjust(wspace=0.1, hspace=0.38); savefig(fig, "fig05d_all_groups_umap.png")
    LABEL = {"coarse_lineage": "cell lineage", "sample_id": "sample", "site": "tissue site", "status": "tissue state", "platform_id": "platform", "biologic": "biologic therapy"}
    fig, ax = plt.subplots(figsize=(7.2, 2.6)); heat_with_values(ax, et.rename(columns=LABEL).iloc[:15].T, vmax=1.0, title="Share of each PC's variance explained by each factor (eta², all four tissue states)")
    savefig(fig, "fig05d_all_groups_pca_eta2.png")


# ================================================================== B. sample level (pseudobulk)
o = adata.obs.copy(); counts = adata.layers["counts"].tocsr()
meta_s = o[["sample_id", "patient_id", "site", "status", "platform_id", "biologic"]].astype({"sample_id": str}).drop_duplicates("sample_id").set_index("sample_id")
excl = adata.var_names.str.match(HVG_EXCLUDE_REGEX) | adata.var_names.isin(SEX_GENES)
def pseudobulk(mask):
    sid = o["sample_id"].astype(str).values[mask]; n_per = pd.Series(sid).value_counts()
    ids = pd.Index(sorted(n_per.index[n_per >= MIN_CELLS])); keep = np.isin(sid, ids); codes = pd.Categorical(sid[keep], categories=ids).codes
    M = sp.csr_matrix((np.ones(keep.sum()), (codes, np.arange(keep.sum()))), shape=(len(ids), keep.sum()))
    return pd.DataFrame((M @ counts[np.where(mask)[0][keep]]).toarray(), index=ids, columns=adata.var_names)
def prep(pb):
    cpm = pb.div(pb.sum(1), axis=0) * 1e6; ok = ((cpm >= 1).mean(0) >= 0.25).values & ~excl
    return np.log2(cpm.loc[:, ok] + 1)
def fit_pca(L, k=4, n_top=2000):
    top = L.var(0).sort_values(ascending=False).index[:n_top]; X = L[top].values; mu = X.mean(0); U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    k = min(k, len(S)); return dict(top=top, mu=mu, Vt=Vt[:k], sc=U[:, :k] * S[:k], ev=(S ** 2 / (S ** 2).sum())[:k], idx=L.index)
def project(L, f): return (L[f["top"]].values - f["mu"]) @ f["Vt"].T
def ols(y, X):
    """OLS; returns coefficient and two-sided p per column (NaN if rank-deficient)."""
    if np.linalg.matrix_rank(X) < X.shape[1] or len(y) - X.shape[1] < 1: return np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan)
    b = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ b; dof = len(y) - X.shape[1]; cov = (r @ r / dof) * np.linalg.inv(X.T @ X)
    return b, 2 * stats.t.sf(np.abs(b / np.sqrt(np.diag(cov))), dof)
def design(ms, ref):                   # intercept + group dummies (all but ref) + site
    cols = [g for g in G if g != ref]; X = np.c_[np.ones(len(ms)), *[(ms.status == g).astype(float).values for g in cols], (ms.site == "AC").astype(float).values]
    return X, cols

SETS = {"all cells": np.ones(len(o), bool), "immune cells": (o["compartment"] == "immune").values, "T/NK cells": (o["coarse_lineage"] == "T_NK").values,
        "myeloid cells": (o["coarse_lineage"] == "myeloid").values, "plasma cells": (o["coarse_lineage"] == "plasma").values, "epithelial cells": (o["coarse_lineage"] == "epithelial").values}
res, rows_eta, rows_grp, axis_rows, axis_tests = {}, [], [], [], []
for name, m in SETS.items():
    L = prep(pseudobulk(m))
    if len(L) < MIN_SAMPLES: print("skip", name, len(L)); continue
    ms = meta_s.loc[L.index]; f = fit_pca(L); sc_ = pd.DataFrame(f["sc"], index=L.index, columns=[f"PC{i+1}" for i in range(f["sc"].shape[1])]).join(ms); res[name] = (sc_, f["ev"])
    print(f"{name}: {len(L)} biopsies; PC1-3 variance {np.round(f['ev'][:3]*100,1)}%", flush=True)
    for i in range(min(3, f["sc"].shape[1])):
        y = f["sc"][:, i]; y = (y - y.mean()) / y.std()
        for fac in ["status", "site", "biologic", "platform_id"]: rows_eta.append(dict(set=name, n_samples=len(L), PC=f"PC{i+1}", var_pct=round(f["ev"][i] * 100, 1), factor=fac, eta2=eta2(y, ms[fac].values)))
        X, cols = design(ms, G[0]); b, p = ols(y, X)               # each group vs control, adjusted for site
        for j, g in enumerate(cols): rows_grp.append(dict(set=name, PC=f"PC{i+1}", group=g, vs=G[0], diff_SD=b[1 + j], p_adj_site=p[1 + j]))
    # ---- C. held-out inflammation axis
    tr = ms.status.isin(["inflamed", "non_inflamed"]).values
    if tr.sum() >= MIN_TRAIN and ms.status[tr].nunique() == 2:
        ft = fit_pca(L[tr], k=3); st = (ms.status[tr] == "inflamed").astype(float).values; si = (ms.site[tr] == "AC").astype(float).values
        best = None
        for i in range(ft["sc"].shape[1]):
            y = ft["sc"][:, i]; b, p = ols((y - y.mean()) / y.std(), np.c_[np.ones(len(y)), st, si])
            if np.isfinite(p[1]) and (best is None or p[1] < best[1]): best = (i, p[1], np.sign(b[1]))
        if best is None: continue
        i, p_sel, sgn = best; sc_all = project(L, ft)[:, i] * sgn; sd = (ft["sc"][:, i] * sgn).std(); mu0 = (ft["sc"][:, i] * sgn)[st == 0].mean()
        z = (sc_all - mu0) / sd                                   # 0 = mean of non-inflamed training biopsies; units = SD of training scores
        d_ = ms.assign(score=z, set=name, axis_PC=f"PC{i+1}", axis_var_pct=round(ft["ev"][i] * 100, 1)); axis_rows.append(d_.reset_index().rename(columns={"index": "sample_id"}))
        for ref in ["non_inflamed", "inflamed"]:
            X, cols = design(ms, ref); b, p = ols(z, X)
            for j, g in enumerate(cols): axis_tests.append(dict(set=name, axis_PC=f"PC{i+1}", reference=ref, group=g, diff_SD=b[1 + j], p_adj_site=p[1 + j],
                                                                 circular=bool({g, ref} == {"inflamed", "non_inflamed"})))
pd.DataFrame(rows_eta).to_csv(f"{OUT}/pseudobulk_all_groups_eta2.csv", index=False); pd.DataFrame(rows_grp).to_csv(f"{OUT}/pseudobulk_all_groups_vs_control.csv", index=False)
if not axis_rows: raise SystemExit(f"no held-out axis could be fitted (needs >= {MIN_TRAIN} inflamed + non-inflamed biopsies per set)")
ax_df = pd.concat(axis_rows); ax_df.to_csv(f"{OUT}/inflammation_axis_scores.csv", index=False)
at = pd.DataFrame(axis_tests); at["q_BH_noncircular"] = np.nan; ok = at.circular.eq(False) & at.p_adj_site.notna()
if ok.any(): at.loc[ok, "q_BH_noncircular"] = stats.false_discovery_control(at.loc[ok, "p_adj_site"].values, method="bh")
at.to_csv(f"{OUT}/inflammation_axis_group_tests.csv", index=False)
med = ax_df.groupby(["set", "status", "site"]).score.agg(["median", "count"]).round(2); med.to_csv(f"{OUT}/inflammation_axis_medians.csv"); print(at.round(4).to_string(), flush=True)

# ---- figure: unsupervised sample-level PCA, all four groups
names = list(res); MK = {"TI": "o", "AC": "^"}
fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.4)); axes = axes.ravel()
for ax_, nm in zip(axes, names):
    sc_, ev = res[nm]
    for g in G:
        for site, mk in MK.items():
            d_ = sc_[(sc_.status == g) & (sc_.site == site)]; ax_.scatter(d_.PC1, d_.PC2, s=20, c=STATUS_COL[g], marker=mk, edgecolors="white", linewidths=0.4, zorder=3)
    ax_.axhline(0, color="0.85", lw=0.5, zorder=0); ax_.axvline(0, color="0.85", lw=0.5, zorder=0); ax_.margins(0.08)
    ax_.set_xlabel(f"PC1 ({ev[0]*100:.0f}%)"); ax_.set_ylabel(f"PC2 ({ev[1]*100:.0f}%)"); ax_.set_title(f"{nm[0].upper() + nm[1:]} (n={len(sc_)})", loc="left")
for ax_ in axes[len(names):]: ax_.axis("off")
hs = [mpl.lines.Line2D([], [], marker="o", ls="", color=STATUS_COL[g], label=DISPLAY[g], ms=4.5) for g in G] + [mpl.lines.Line2D([], [], marker=m, ls="", color="0.4", label=l, ms=4.5) for l, m in [("terminal ileum", "o"), ("ascending colon", "^")]]
fig.legend(handles=hs, loc="lower center", ncol=6, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, -0.01))
for ax_, l in zip(axes[:len(names)], "abcdef"): panel_letter(ax_, l)
fig.subplots_adjust(hspace=0.5, wspace=0.4, bottom=0.12); savefig(fig, "fig05d_pseudobulk_pca_all_groups.png")

# ---- figure: held-out inflammation axis, by group, split by site
sets_ax = [s_ for s_ in ["all cells", "immune cells", "T/NK cells", "plasma cells", "myeloid cells", "epithelial cells"] if s_ in set(ax_df["set"])]
fig, axes = plt.subplots(len(sets_ax), 2, figsize=(7.2, 1.55 * len(sets_ax) + 0.6), sharex=True, sharey="row", squeeze=False)   # sharey="row": same y-scale for ileum and colon within a row
rng = np.random.default_rng(0)
for r_, nm in enumerate(sets_ax):
    d_ = ax_df[ax_df["set"] == nm]
    for c_, site in enumerate(["TI", "AC"]):
        ax_ = axes[r_, c_]; dd = d_[d_.site == site]
        for k, g in enumerate(G):
            v = dd[dd.status == g].score.values; held = g in ("control", "adjacent")
            ax_.scatter(k + rng.uniform(-0.18, 0.18, len(v)), v, s=16, marker="o", facecolors="none" if held else STATUS_COL[g], edgecolors=STATUS_COL[g], linewidths=0.9 if held else 0.3, zorder=3)
            if len(v): ax_.hlines(np.median(v), k - 0.3, k + 0.3, color="k", lw=1.2, zorder=4)
        ax_.axhline(0, color="0.85", lw=0.5, zorder=0); ax_.margins(x=0.08)
        if c_ == 0: ax_.set_ylabel(f"{nm}\n(SD units)", fontsize=6.5)
        if r_ == 0: ax_.set_title({"TI": "Terminal ileum", "AC": "Ascending colon"}[site], loc="left")
WRAP = {"control": "non-IBD\ncontrol", "non_inflamed": "non-\ninflamed", "adjacent": "adjacent\nto inflamed", "inflamed": "inflamed"}
for ax_ in axes[-1]: ax_.set_xticks(range(4)); ax_.set_xticklabels([WRAP[g] for g in G], fontsize=5.5)
fig.suptitle("Position of each biopsy on the inflamed-vs-non-inflamed axis (open circles: groups that did not define the axis)", x=0.01, ha="left", fontsize=7.5, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965)); savefig(fig, "fig05d_inflammation_axis_by_group.png")
