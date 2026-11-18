"""Step 4 - Normalisation and feature selection.

  * raw UMI counts are kept in layers['counts'] (needed later by MrVI/scVI and for pseudobulk)
  * X = log1p( counts / library size * 1e4 )  (library-size normalisation; makes each cell's total equal before log)
  * highly variable genes: 3,000, Seurat-v3 variance-stabilised ranking on raw counts, selected PER SAMPLE and then combined
    (batch_key = sample_id): genes that are variable in many biopsies are kept, genes variable in a single biopsy are not.
    Excluded from HVG selection because they encode technical, clone- or sex-specific signal rather than cell state:
    mitochondrial and ribosomal-protein genes, haemoglobin, immunoglobulin / TCR V-D-J segment genes (patient-specific clones), sex-linked genes.
  * no regression of total counts / MT fraction (would remove real biology: epithelial cells legitimately differ in both).
Output: $PROC/03_norm.h5ad, out/hvg_list.csv, out/fig04_normalisation.png
"""
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from common import *
style()
adata = sc.read_h5ad(f"{PROC}/02_qc.h5ad")
n0 = adata.n_vars
sc.pp.filter_genes(adata, min_cells=10)
print(f"genes: {n0} -> {adata.n_vars} after min_cells=10", flush=True)
adata.layers["counts"] = adata.X.copy()
sf = np.asarray(adata.layers["counts"].sum(1)).ravel()
sc.pp.normalize_total(adata, target_sum=1e4)
rs = np.asarray(adata.X[:2000].sum(1)).ravel()
print("check after normalize_total: per-cell sums mean %.1f sd %.4f (expected 10000)" % (rs.mean(), rs.std()), flush=True)
sc.pp.log1p(adata)
rs = np.asarray(np.expm1(adata.X[:2000]).sum(1)).ravel()
print("check after log1p: expm1 sums mean %.1f sd %.4f (expected 10000)" % (rs.mean(), rs.std()), flush=True)
adata.obs["library_size"] = sf

# ---- HVGs
excl = adata.var_names.str.match(HVG_EXCLUDE_REGEX) | adata.var_names.isin(SEX_GENES)
print("genes excluded from HVG selection:", int(excl.sum()), flush=True)
tmp = adata[:, ~excl].copy(); tmp.X = tmp.layers["counts"].copy()
sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=3000, batch_key="sample_id", layer="counts", span=0.3)
adata.var["highly_variable"] = False; adata.var["highly_variable_nbatches"] = 0; adata.var["hvg_rank"] = np.nan
sel = tmp.var_names[tmp.var["highly_variable"].values]
adata.var.loc[sel, "highly_variable"] = True
adata.var.loc[tmp.var_names, "highly_variable_nbatches"] = tmp.var["highly_variable_nbatches"].values
adata.var.loc[tmp.var_names, "hvg_rank"] = tmp.var["highly_variable_rank"].values
print("HVGs selected:", int(adata.var.highly_variable.sum()), flush=True)
adata.var.loc[adata.var.highly_variable].sort_values("hvg_rank")[["gene_ids", "highly_variable_nbatches", "hvg_rank"]].to_csv(f"{OUT}/hvg_list.csv")
adata.uns["normalization"] = dict(method="normalize_total(target_sum=1e4)+log1p", hvg="seurat_v3, n=3000, batch_key=sample_id, excluded MT/ribo/Hb/IG-TR V(D)J/sex genes", counts_layer="counts")
adata.write_h5ad(f"{PROC}/03_norm.h5ad", compression="lzf")

# ---- figure
o = adata.obs
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
axes[0].hist(np.log10(o["library_size"]), bins=80, color="#777777", linewidth=0)
axes[0].set_xlabel("UMI per cell (log10), before normalisation"); axes[0].set_ylabel("Cells"); axes[0].set_title("Library size per cell", loc="left")
for st in STATUS_ORDER:
    axes[1].hist(np.log10(o.loc[o.status == st, "library_size"]), bins=60, histtype="step", color=STATUS_COL[st], lw=1, density=True, label=DISPLAY[st])
axes[1].set_xlabel("UMI per cell (log10)"); axes[1].set_ylabel("Density"); axes[1].legend(fontsize=6); axes[1].set_title("By tissue state", loc="left")
nb = adata.var.loc[adata.var.highly_variable, "highly_variable_nbatches"]
axes[2].hist(nb, bins=np.arange(0.5, adata.obs["sample_id"].nunique() + 1.5, 1), color="#555555", linewidth=0)
axes[2].set_xlabel("Number of samples in which a HVG is variable"); axes[2].set_ylabel("Highly variable genes"); axes[2].set_title("HVG recurrence", loc="left")
for ax_, l in zip(axes, "abc"): panel_letter(ax_, l)
fig.tight_layout(); savefig(fig, "fig04_normalisation.png")
