"""17b_mrvi_mixing_diagnostics.py - the diagnostics part of step 6d that crashed in 17 (Arrow-backed string arrays cannot be 2-D indexed in pandas 3).
Reads the MrVI latent table written by 17 (mrvi_immune_cells_latent.csv.gz) + the immune input counts; computes kNN mixing of biopsies and lineage purity in
MrVI u-space vs plain PCA; writes knn_mixing_mrvi_vs_pca.csv and 06b_immune_mrvi_latent.h5ad (u + UMAP from step 4 on the 3,000-gene input)."""
import os, numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
PROC = os.path.expanduser("~/.claude-science-scratch/gse266616/processed"); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); lat = pd.read_csv("mrvi_immune_cells_latent.csv.gz").set_index("barcode")
assert a.n_obs == len(lat) and (a.obs_names == lat.index).all(), "latent table order differs from the input"
u = lat[[c for c in lat.columns if c.startswith("u") and c[1:].isdigit()]].values; sid = np.asarray(lat["sample_id"].astype(str), dtype=object); lin = np.asarray(lat["immune_lineage"].astype(str), dtype=object)
a.obsm["X_mrvi_u"] = u; a.obsm["X_umap"] = lat[["umap1", "umap2"]].values
for c_ in ["status", "site", "biologic", "immune_lineage", "sample_id"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
q = ad.AnnData(sp.csr_matrix(a.X), obs=a.obs[["sample_id"]].copy()); sc.pp.normalize_total(q, target_sum=1e4); sc.pp.log1p(q); sc.pp.scale(q, max_value=10); sc.tl.pca(q, n_comps=30, svd_solver="arpack", random_state=0)
a.obsm["X_pca_hvg"] = q.obsm["X_pca"]
def knn_idx(X, k=30): return NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]
rows = []
for nm, X in [("MrVI u (sample-unaware)", u), ("PCA of 3,000 HVGs", a.obsm["X_pca_hvg"])]:
    nn = knn_idx(X); same_b = (sid[nn] == sid[:, None]).mean(1); same_l = (lin[nn] == lin[:, None]).mean(1)
    for L in ["all"] + sorted(set(lin)):
        m = np.ones(len(lin), bool) if L == "all" else lin == L
        cnt = pd.Series(sid[m]).value_counts() / m.sum(); exp = pd.Series(sid[m]).map(cnt).values.astype(float)
        rows.append(dict(space=nm, lineage=L, n_cells=int(m.sum()), same_biopsy_neighbour_frac=same_b[m].mean(), expected_if_mixed=exp.mean(), ratio_observed_to_expected=same_b[m].mean() / exp.mean(), same_lineage_neighbour_frac=same_l[m].mean()))
mix = pd.DataFrame(rows); mix.round(4).to_csv(f"{OUT}/knn_mixing_mrvi_vs_pca.csv", index=False); print(mix.round(3).to_string(index=False), flush=True)
a.write_h5ad(f"{PROC}/06b_immune_mrvi_latent.h5ad", compression="lzf"); print("DIAGNOSTICS DONE", flush=True)
