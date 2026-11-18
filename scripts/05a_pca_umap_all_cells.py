"""Step 5a - PCA, neighbours, UMAP and Leiden on ALL QC-passed cells (whole-biopsy view).

Biological expectation being tested (and how it is judged):
  1. PCs/UMAP of a whole intestinal biopsy should be dominated by CELL IDENTITY: epithelial vs immune vs stromal first, then immune lineages.
     -> coloured by marker-based coarse lineage; PC variance explained (eta^2) by lineage should be far larger than by any other factor.
  2. Technical covariates (library size, MT%, doublet score, platform) should NOT define clusters.
  3. Site (terminal ileum vs ascending colon) mainly separates EPITHELIAL cells (region-specific enterocytes/colonocytes); immune cells should mix.
  4. Individual samples/patients should mix within a lineage; residual donor structure is expected and is what MrVI (step 6) models later.
  5. Inflamed vs non-inflamed is NOT expected to be a dominant axis across all cells; it is a subtle shift inside particular cell types
     (5b/5c look at this within immune cells and at sample level).
Outputs: $PROC/04_embedded.h5ad ; out/pca_eta2_all_cells.csv, out/leiden_lineage_scores.csv, out/lineage_summary.csv, out/knn_mixing_by_lineage.csv ; out/fig05a_*.png
"""
import argparse, json, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, seaborn as sns
from common import *
style()
ap = argparse.ArgumentParser()
ap.add_argument("--n-pcs", type=int, default=30); ap.add_argument("--n-neighbors", type=int, default=30)
ap.add_argument("--leiden-res", type=float, default=0.5); ap.add_argument("--min-dist", type=float, default=0.3)
a = ap.parse_args()
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
adata = sc.read_h5ad(f"{PROC}/03_norm.h5ad")
print("loaded", adata, flush=True)

# ---------------- PCA on scaled HVGs
h = adata[:, adata.var["highly_variable"].values].copy()
sc.pp.scale(h, max_value=10)                      # z-score each gene, clip at 10 so a few extreme cells do not define a PC
sc.tl.pca(h, n_comps=50, svd_solver="arpack", random_state=0)
vr = np.asarray(h.uns["pca"]["variance_ratio"])
def knee(v):                                       # elbow = point of the scree curve farthest from the chord first->last PC
    x = np.arange(len(v)) / (len(v) - 1); y = (v - v.min()) / (v.max() - v.min())
    return int(np.argmax(np.abs((y[-1] - y[0]) * x - (x[-1] - x[0]) * y + x[-1] * y[0] - y[-1] * x[0]))) + 1
k = knee(vr); print(f"PCA: PC1 {vr[0]*100:.1f}%, first 10 PCs {vr[:10].sum()*100:.1f}%, first {a.n_pcs} PCs {vr[:a.n_pcs].sum()*100:.1f}% ; scree elbow at PC{k}; using {a.n_pcs} PCs", flush=True)
sc.pp.neighbors(h, n_neighbors=a.n_neighbors, n_pcs=a.n_pcs, random_state=0)
sc.tl.umap(h, min_dist=a.min_dist, random_state=0)
sc.tl.leiden(h, resolution=a.leiden_res, flavor="igraph", n_iterations=2, directed=False, random_state=0, key_added="leiden")
adata.obsm["X_pca"] = h.obsm["X_pca"]; adata.obsm["X_umap"] = h.obsm["X_umap"]
adata.obs["leiden"] = h.obs["leiden"].astype(str).values
adata.uns["pca_variance_ratio"] = vr
print("leiden clusters:", adata.obs["leiden"].nunique(), flush=True)

# ---------------- provisional coarse lineage from marker scores (cluster level)
for name, genes in LINEAGE_MARKERS.items():
    gl = [g for g in genes if g in adata.var_names]
    print(f"score {name}: {len(gl)}/{len(genes)} markers present", flush=True)
    sc.tl.score_genes(adata, gl, score_name=f"score_{name}", random_state=0)
S = adata.obs[[f"score_{n}" for n in LINEAGE_MARKERS]].copy(); S.columns = list(LINEAGE_MARKERS)
cm = S.groupby(adata.obs["leiden"].values, observed=True).mean()
z = (cm - cm.mean()) / cm.std(ddof=0).replace(0, 1)
assign = {}
for kk in cm.index:
    best = z.loc[kk].idxmax(); rawbest = cm.loc[kk, best]; other = cm.loc[kk].drop(best)
    assign[kk] = "unassigned" if rawbest < 0.1 else ("mixed" if (other.max() >= 0.8 * rawbest and other.max() >= 0.5) else best)
cmo = cm.copy(); cmo["assigned"] = pd.Series(assign); cmo["n_cells"] = adata.obs["leiden"].value_counts().reindex(cmo.index).values
cmo.round(3).to_csv(f"{OUT}/leiden_lineage_scores.csv")
adata.obs["coarse_lineage"] = adata.obs["leiden"].map(assign).astype(str)
adata.obs["compartment"] = adata.obs["coarse_lineage"].map(COMPARTMENT).astype(str)
def gene(g): return np.asarray(adata[:, g].X.todense()).ravel() if g in adata.var_names else None
ptprc, epcam = gene("PTPRC"), gene("EPCAM")
ls = adata.obs.assign(ptprc_pos=ptprc > 0, epcam_pos=epcam > 0).groupby("coarse_lineage").agg(
    n_cells=("sample_id", "size"), pct_PTPRC_pos=("ptprc_pos", lambda x: 100 * x.mean()), pct_EPCAM_pos=("epcam_pos", lambda x: 100 * x.mean()),
    median_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median"), median_doublet_score=("doublet_score", "median")).round(2)
ls["pct_of_cells"] = (100 * ls.n_cells / ls.n_cells.sum()).round(1); ls.to_csv(f"{OUT}/lineage_summary.csv"); print(ls.to_string(), flush=True)

# ---------------- how much of each PC is explained by which factor
FACT = ["coarse_lineage", "sample_id", "patient_id", "site", "status", "platform_id", "biologic"]
et = eta2_table(h.obsm["X_pca"], adata.obs, FACT, n_pcs=20); et.to_csv(f"{OUT}/pca_eta2_all_cells.csv")
print("eta^2 (PC1-5):\n", et.head(5).round(3).to_string(), flush=True)

# ---------------- neighbourhood mixing: observed / expected fraction of kNN neighbours sharing a label, per lineage.
# 1.0 = neighbours are drawn at random; larger = cells sit next to cells with the same label. Status and site are properties of whole BIOPSIES, so
# part of their ratio is simply donor/sample structure. The null therefore re-assigns status (or site) at random ACROSS SAMPLES (keeping every
# sample intact and the group sizes fixed, 200 permutations): a ratio inside the null range means no structure beyond sample identity.
dist = h.obsp["distances"].tocsr(); n = h.n_obs; kk_ = np.diff(dist.indptr)
rows = []; rng = np.random.default_rng(0); N_PERM = 200
if kk_.min() == kk_.max():
    nn = dist.indices.reshape(n, -1)
    smp = pd.Categorical(adata.obs["sample_id"].astype(str)).codes
    sinfo = adata.obs[["sample_id", "status", "site"]].astype(str).drop_duplicates("sample_id").sort_values("sample_id")
    samp_lab = {"status": pd.Categorical(sinfo["status"]).codes, "site": pd.Categorical(sinfo["site"]).codes}
    for lin in [l for l in pd.unique(adata.obs["coarse_lineage"]) if (adata.obs["coarse_lineage"] == l).sum() >= 500]:
        m = (adata.obs["coarse_lineage"] == lin).values
        r = {"lineage": lin, "n_cells": int(m.sum())}
        for f in ["sample_id", "patient_id", "platform_id"]:
            lab = pd.Categorical(adata.obs[f].astype(str)).codes; p = np.bincount(lab[m]) / m.sum()
            r[f] = round(float((lab[nn[m]] == lab[m][:, None]).mean() / (p ** 2).sum()), 2)
        cs, ns = smp[m], smp[nn[m]]
        def ratio(lab_s):
            p = np.bincount(lab_s[cs]) / len(cs)
            return float((lab_s[ns] == lab_s[cs][:, None]).mean() / (p ** 2).sum())
        for f in ["status", "site"]:
            ob = ratio(samp_lab[f]); nl = np.array([ratio(rng.permutation(samp_lab[f])) for _ in range(N_PERM)])
            r[f] = round(ob, 2); r[f + "_null_mean"] = round(float(nl.mean()), 2); r[f + "_null_p95"] = round(float(np.percentile(nl, 95)), 2)
            r[f + "_perm_p"] = round((1 + int((nl >= ob).sum())) / (1 + N_PERM), 3)
        rows.append(r)
pd.DataFrame(rows).to_csv(f"{OUT}/knn_mixing_by_lineage.csv", index=False); print(pd.DataFrame(rows).to_string(), flush=True)

adata.write_h5ad(f"{PROC}/04_embedded.h5ad", compression="lzf"); print("wrote 04_embedded.h5ad", flush=True)

# ================================================================= figures
xy = adata.obsm["X_umap"]; pc = adata.obsm["X_pca"]; o = adata.obs
LABEL = {"coarse_lineage": "cell lineage", "sample_id": "sample", "patient_id": "patient", "site": "tissue site", "status": "inflamed vs non-inflamed",
         "platform_id": "sequencing platform", "biologic": "biologic therapy"}
# --- PCA overview
fig = plt.figure(figsize=(7.2, 6.0)); gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15], width_ratios=[1, 1, 1, 0.42], hspace=0.45, wspace=0.35)
ax_a = fig.add_subplot(gs[0, 0]); ax_a.plot(np.arange(1, 51), vr * 100, "o-", ms=2, lw=0.6, color="#333333")
ax_a.axvline(a.n_pcs, color="#D55E00", lw=0.7, ls="--"); ax_a.text(a.n_pcs + 1, vr[0] * 100 * 0.8, f"{a.n_pcs} PCs used", fontsize=6, color="#D55E00")
ax_a.set_xlabel("Principal component"); ax_a.set_ylabel("Variance explained (%)"); ax_a.set_yscale("log"); ax_a.set_title("Variance explained per PC", loc="left")
ax_b = fig.add_subplot(gs[0, 1]); embed_cat(ax_b, pc[:, [0, 1]], o["coarse_lineage"], LINEAGE_COL, s=0.25, xlab=f"PC1 ({vr[0]*100:.1f}%)", ylab=f"PC2 ({vr[1]*100:.1f}%)", title="PC1 vs PC2")
ax_c = fig.add_subplot(gs[0, 2]); embed_cat(ax_c, pc[:, [2, 3]], o["coarse_lineage"], LINEAGE_COL, s=0.25, xlab=f"PC3 ({vr[2]*100:.1f}%)", ylab=f"PC4 ({vr[3]*100:.1f}%)", title="PC3 vs PC4")
ax_l = fig.add_subplot(gs[0, 3]); ax_l.axis("off")
present = set(o["coarse_lineage"].astype(str))
ax_l.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", ms=4, color=LINEAGE_COL[g], label=g) for g in LINEAGE_COL if g in present],
            loc="center left", fontsize=6, title="Lineage\n(provisional)", title_fontsize=6, handletextpad=0.2, borderaxespad=0)
ax_d = fig.add_subplot(gs[1, :]); et_plot = et.rename(columns=LABEL)
heat_with_values(ax_d, et_plot.iloc[:15].T, vmax=1.0, title="Share of each PC's variance explained by each factor (eta², 1 = fully explained)")
for ax_, l in zip([ax_a, ax_b, ax_c, ax_d], "abcd"): panel_letter(ax_, l)
savefig(fig, "fig05a_all_pca_overview.png")

# --- lineage + markers
fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.0)); axes = axes.ravel()
embed_cat(axes[0], xy, o["coarse_lineage"], LINEAGE_COL, s=0.25, label_groups=True, title="Provisional cell lineage")
for ax_, g in zip(axes[1:], ["PTPRC", "EPCAM", "CD3E", "MS4A1", "JCHAIN", "LYZ", "COL1A1", "PECAM1"]):
    v = gene(g)
    if v is not None: embed_cont(ax_, xy, v, title=g, italic=True, s=0.25)
for ax_, l in zip(axes, "abcdefghi"): panel_letter(ax_, l)
fig.subplots_adjust(wspace=0.3, hspace=0.3); savefig(fig, "fig05a_all_umap_lineage_markers.png")

# --- technical covariates
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9)); axes = axes.ravel()
embed_cont(axes[0], xy, np.log10(o["total_counts"]), title="UMI per cell (log10)", cmap="viridis", s=0.25)
embed_cont(axes[1], xy, np.log10(o["n_genes_by_counts"]), title="Genes per cell (log10)", cmap="viridis", s=0.25)
embed_cont(axes[2], xy, o["pct_counts_mt"].values, title="Mitochondrial reads (%)", cmap="viridis", s=0.25, vmax=25)
embed_cont(axes[3], xy, o["doublet_score"].values, title="Doublet score (Scrublet)", cmap="viridis", s=0.25)
embed_cat(axes[4], xy, o["platform_id"], PLATFORM_COL, s=0.25, legend=True, title="Sequencing platform", legend_ncol=2)
lv = sorted(o["leiden"].unique(), key=int); pal = dict(zip(lv, sns.color_palette("husl", len(lv)).as_hex()))
embed_cat(axes[5], xy, o["leiden"], pal, s=0.25, label_groups=True, title=f"Leiden clusters (n={len(lv)})")
for ax_, l in zip(axes, "abcdef"): panel_letter(ax_, l)
fig.subplots_adjust(wspace=0.3, hspace=0.38); savefig(fig, "fig05a_all_umap_technical.png")

# --- condition, site, sample
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9)); axes = axes.ravel()
embed_cat(axes[0], xy, o["status"], STATUS_COL, s=0.25, legend=True, title="Inflamed vs non-inflamed")
embed_cat(axes[1], xy, o["site"], SITE_COL, s=0.25, legend=True, title="Tissue site")
sp_ = dict(zip(sorted(o["sample_id"].unique()), sns.color_palette("husl", o["sample_id"].nunique()).as_hex()))
embed_cat(axes[2], xy, o["sample_id"], sp_, s=0.25, title=f"Sample ({o['sample_id'].nunique()} biopsies)")
embed_cat(axes[3], xy, o["biologic"], BIOLOGIC_COL, s=0.25, legend=True, title="Biologic therapy")
embed_highlight(axes[4], xy, (o["status"] == "non_inflamed").values, STATUS_COL["non_inflamed"], "Non-inflamed cells only", s=0.25)
embed_highlight(axes[5], xy, (o["status"] == "inflamed").values, STATUS_COL["inflamed"], "Inflamed cells only", s=0.25)
for ax_, l in zip(axes, "abcdef"): panel_letter(ax_, l)
fig.subplots_adjust(wspace=0.1, hspace=0.38); savefig(fig, "fig05a_all_umap_condition_site_sample.png")
