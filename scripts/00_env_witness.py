"""Kernel-dispatch witness: a seeded end-to-end mini pipeline that exercises every function the real pipeline uses."""
import numpy as np, scanpy as sc, anndata as ad, scipy.sparse as sp, igraph, leidenalg, skmisc, scrublet, umap
print("versions: scanpy", sc.__version__, "anndata", ad.__version__, "igraph", igraph.__version__,
      "leidenalg", leidenalg.__version__, "scrublet", getattr(scrublet, "__version__", "?"), "umap", umap.__version__)
rng = np.random.default_rng(0)
n, g = 1500, 2000
mu = rng.gamma(0.3, 2.0, size=g)
grp = rng.integers(0, 3, n)
X = rng.poisson(mu[None, :] * (1 + 0.8 * (rng.random((3, g)) < 0.1)[grp]) * 2).astype(np.float32)
a = ad.AnnData(sp.csr_matrix(X))
a.var_names = [f"G{i}" for i in range(g)]
a.obs["sample_id"] = rng.choice(["s1", "s2", "s3"], n)
a.layers["counts"] = a.X.copy()
sc.pp.scrublet(a, batch_key="sample_id")
sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=500, batch_key="sample_id", layer="counts")
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.scale(a, max_value=10)
sc.tl.pca(a, n_comps=20, mask_var="highly_variable")
sc.pp.neighbors(a, n_neighbors=15, n_pcs=15)
sc.tl.umap(a)
sc.tl.leiden(a, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
print("WITNESS OK pca", a.obsm["X_pca"].shape, "umap", a.obsm["X_umap"].shape, "leiden", a.obs["leiden"].nunique(),
      "doublet_score_nan", int(np.isnan(a.obs["doublet_score"]).sum()))
