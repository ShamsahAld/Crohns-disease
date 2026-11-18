"""Step 5b - PCA/UMAP restricted to IMMUNE cells (the cells the research question is about).

In a whole biopsy the epithelial/immune split uses up the leading PCs (5a). Re-doing HVG selection, scaling and PCA inside the immune
compartment lets the PCs describe immune-cell state: T/NK vs B vs plasma vs myeloid vs mast first, then states within them, and only then
(if at all) differences between inflamed and non-inflamed tissue. Cells in clusters flagged 'mixed'/'unassigned' are excluded (possible doublets).
Outputs: $PROC/05_immune.h5ad ; out/immune_pca_eta2.csv, out/immune_eta2_within_lineage.csv, out/immune_composition_by_sample.csv ; out/fig05b_*.png
"""
import argparse, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, seaborn as sns
from common import *
style()
ap = argparse.ArgumentParser()
ap.add_argument("--n-pcs", type=int, default=30); ap.add_argument("--n-neighbors", type=int, default=30)
ap.add_argument("--leiden-res", type=float, default=0.6); ap.add_argument("--min-dist", type=float, default=0.3)
a = ap.parse_args()
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
adata = sc.read_h5ad(f"{PROC}/04_embedded.h5ad")
n_all = adata.n_obs
adata = adata[adata.obs["compartment"] == "immune"].copy()
print(f"immune cells: {adata.n_obs} of {n_all} ({100*adata.n_obs/n_all:.1f}%)", flush=True)
for c_ in ["status", "site", "sample_id", "patient_id", "platform_id", "biologic"]: adata.obs[c_] = adata.obs[c_].astype("category")
adata.obs["immune_lineage"] = adata.obs["coarse_lineage"].astype(str)
adata.obsm.pop("X_pca", None); adata.obsm.pop("X_umap", None)

# ---- re-select HVGs within immune cells (same exclusions as step 4)
excl = adata.var_names.str.match(HVG_EXCLUDE_REGEX) | adata.var_names.isin(SEX_GENES)
tmp = adata[:, ~excl].copy(); tmp.X = tmp.layers["counts"].copy()
sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=3000, batch_key="sample_id", layer="counts", span=0.3)
adata.var["highly_variable_immune"] = False; adata.var.loc[tmp.var_names[tmp.var["highly_variable"].values], "highly_variable_immune"] = True
h = adata[:, adata.var["highly_variable_immune"].values].copy(); sc.pp.scale(h, max_value=10)
sc.tl.pca(h, n_comps=50, svd_solver="arpack", random_state=0); vr = np.asarray(h.uns["pca"]["variance_ratio"])
print(f"immune PCA: PC1 {vr[0]*100:.1f}%, first {a.n_pcs} PCs {vr[:a.n_pcs].sum()*100:.1f}%", flush=True)
sc.pp.neighbors(h, n_neighbors=a.n_neighbors, n_pcs=a.n_pcs, random_state=0)
sc.tl.umap(h, min_dist=a.min_dist, random_state=0)
sc.tl.leiden(h, resolution=a.leiden_res, flavor="igraph", n_iterations=2, directed=False, random_state=0, key_added="leiden_immune")
adata.obsm["X_pca"] = h.obsm["X_pca"]; adata.obsm["X_umap"] = h.obsm["X_umap"]; adata.obs["leiden_immune"] = h.obs["leiden_immune"].astype(str).values
o = adata.obs; xy = adata.obsm["X_umap"]; pc = adata.obsm["X_pca"]

# ---- variance explained: overall, and WITHIN each immune lineage
FACT = ["immune_lineage", "sample_id", "patient_id", "site", "status", "platform_id", "biologic"]
et = eta2_table(pc, o, FACT, n_pcs=15); et.to_csv(f"{OUT}/immune_pca_eta2.csv"); print(et.head(6).round(3).to_string(), flush=True)
rows = []
for lin in IMMUNE:
    m = (o["immune_lineage"] == lin).values
    if m.sum() < 1000: continue
    for f in ["sample_id", "patient_id", "site", "status", "platform_id", "biologic"]:
        for i in range(5):
            rows.append(dict(lineage=lin, n_cells=int(m.sum()), factor=f, PC=f"PC{i+1}", eta2=eta2(pc[m, i], o.loc[m, f].values)))
ew = pd.DataFrame(rows); ew.to_csv(f"{OUT}/immune_eta2_within_lineage.csv", index=False)
comp = pd.crosstab(o["sample_id"], o["immune_lineage"], normalize="index").round(4); comp.to_csv(f"{OUT}/immune_composition_by_sample.csv")
adata.write_h5ad(f"{PROC}/05_immune.h5ad", compression="lzf"); print("wrote 05_immune.h5ad", flush=True)

# ================================================================= figures
def gene(g): return np.asarray(adata[:, g].X.todense()).ravel() if g in adata.var_names else None
LABEL = {"immune_lineage": "immune lineage", "sample_id": "sample", "patient_id": "patient", "site": "tissue site", "status": "inflamed vs non-inflamed",
         "platform_id": "sequencing platform", "biologic": "biologic therapy"}
lvl = sorted(o["leiden_immune"].unique(), key=int); pal = dict(zip(lvl, sns.color_palette("husl", len(lvl)).as_hex()))
fig, axes = plt.subplots(3, 4, figsize=(7.2, 5.6)); axes = axes.ravel()
embed_cat(axes[0], xy, o["immune_lineage"], {k: v for k, v in LINEAGE_COL.items() if k in IMMUNE}, s=0.4, label_groups=True, title="Immune lineage")
embed_cat(axes[1], xy, o["leiden_immune"], pal, s=0.4, label_groups=True, title=f"Leiden clusters (n={len(lvl)})")
for ax_, g in zip(axes[2:], ["CD3E", "CD4", "CD8A", "FOXP3", "NKG7", "MS4A1", "JCHAIN", "CD14", "C1QA", "TPSAB1"]):
    v = gene(g)
    if v is not None: embed_cont(ax_, xy, v, title=g, italic=True, s=0.4)
for ax_, l in zip(axes, "abcdefghijkl"): panel_letter(ax_, l)
fig.subplots_adjust(wspace=0.3, hspace=0.3); savefig(fig, "fig05b_immune_umap_lineage_markers.png")

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9)); axes = axes.ravel()
embed_cat(axes[0], xy, o["status"], STATUS_COL, s=0.4, legend=True, title="Inflamed vs non-inflamed")
embed_cat(axes[1], xy, o["site"], SITE_COL, s=0.4, legend=True, title="Tissue site")
sp_ = dict(zip(sorted(o["sample_id"].astype(str).unique()), sns.color_palette("husl", o["sample_id"].nunique()).as_hex()))
embed_cat(axes[2], xy, o["sample_id"], sp_, s=0.4, title=f"Sample ({o['sample_id'].nunique()} biopsies)")
embed_cat(axes[3], xy, o["biologic"], BIOLOGIC_COL, s=0.4, legend=True, title="Biologic therapy")
embed_highlight(axes[4], xy, (o["status"] == "non_inflamed").values, STATUS_COL["non_inflamed"], "Non-inflamed immune cells only", s=0.4)
embed_highlight(axes[5], xy, (o["status"] == "inflamed").values, STATUS_COL["inflamed"], "Inflamed immune cells only", s=0.4)
for ax_, l in zip(axes, "abcdef"): panel_letter(ax_, l)
fig.subplots_adjust(wspace=0.1, hspace=0.38); savefig(fig, "fig05b_immune_umap_condition_site_sample.png")

fig = plt.figure(figsize=(7.2, 3.6)); gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.55)
ax = fig.add_subplot(gs[0]); heat_with_values(ax, et.rename(columns=LABEL).iloc[:12].T, vmax=1.0, title="All immune cells")
ax = fig.add_subplot(gs[1]); lins = [l for l in IMMUNE if (ew.lineage == l).any()]
facs = [("status", "inflamed vs\nnon-inflamed"), ("site", "tissue\nsite"), ("sample_id", "sample")]
mat = pd.DataFrame({f"{f}_{p}": [ew[(ew.lineage == l) & (ew.factor == f) & (ew.PC == p)].eta2.iloc[0] for l in lins] for f, _ in facs for p in ["PC1", "PC2", "PC3"]},
                   index=[f"{l} (n={int(ew[ew.lineage == l].n_cells.iloc[0]):,})" for l in lins])
mat.columns = ["PC1", "PC2", "PC3"] * len(facs)
heat_grouped(ax, mat, [(lab, 3) for _, lab in facs], vmax=0.5, title="Within one lineage")
for ax_, l in zip(fig.axes[:2], "ab"): panel_letter(ax_, l, dx=-0.12)
fig.suptitle("Variance explained (eta²) by each factor", x=0.01, ha="left", fontsize=8, y=1.03)
savefig(fig, "fig05b_immune_pca_eta2.png")
