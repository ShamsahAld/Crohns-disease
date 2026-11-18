"""24c_mrvi_de_aggregate.py - plan step 8b: MrVI per-cell differential expression (status, site, biologic as sample covariates), aggregated per lineage and per state.

MrVI's per-cell 'pvalue'/'padj' test a shift in the 30-dimensional latent space, not individual genes (95 % of cells are 'significant' for status in the probe), so they are not used for gene calls.
Gene-level outputs used: per-cell log fold change 'lfc' (covariate effect on expression of each gene, others held fixed), its posterior sd 'lfc_std' and 'pde' = posterior probability |lfc| > delta (0.3).
Per unit (lineage or state; non-excluded cells) and covariate we store, for each of the 3,000 model genes: mean lfc, sd of lfc across cells, mean pde, fraction of cells with lfc > 0, n cells.
Cells are processed in chunks of 4,000 (each cell's DE is independent of the others, so chunking is exact).  Model: step-4 MrVI, 61,538 immune cells, biopsy = sample.
Output: out/mrvi_de_unit_gene_stats.npz.
"""
import os, time, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, torch, scvi
from scvi.external import MRVI
PROC = os.path.expanduser("~/.claude-science-scratch/gse266616/processed"); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True); torch.set_num_threads(8); scvi.settings.seed = 0
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); a.layers["counts"] = a.X.copy()
for c_ in ["sample_id", "status", "site", "biologic", "immune_lineage"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"]); a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(["none", "anti_TNF", "anti_IL12_IL23"])
MRVI.setup_anndata(a, layer="counts", sample_key="sample_id"); model = MRVI.load(f"{PROC}/mrvi_immune_model", adata=a, accelerator="cpu")
lab = pd.read_csv("cell_state_labels.csv.gz").set_index("barcode").loc[a.obs_names]; keep = np.where((lab.confidence != "exclude").values)[0]; print("cells analysed:", len(keep), flush=True)
unit_names = ["lineage:" + L for L in ["T_NK", "plasma", "B", "myeloid", "mast"]] + ["state:" + s for s in sorted(lab.state[lab.confidence != "exclude"].unique())]
memb = {u: (lab.immune_lineage.values == u.split(":", 1)[1]) if u.startswith("lineage:") else (lab.state.values == u.split(":", 1)[1]) for u in unit_names}
U = len(unit_names); nG = a.n_vars; covs = None; S1 = S2 = SP = SPOS = None; N = np.zeros(U); CH = 4000; t0 = time.time()
for ci, start in enumerate(range(0, len(keep), CH)):
    idx = keep[start:start + CH]; sub = a[idx].copy()
    de = model.differential_expression(adata=sub, sample_cov_keys=["status", "site", "biologic"], store_lfc=True, use_vmap=False, mc_samples=50, delta=0.3)
    if covs is None: covs = [str(x) for x in de.coords["covariate"].values]; genes = np.asarray(de.coords["gene"].values, dtype=str); S1 = np.zeros((U, len(covs), nG)); S2 = np.zeros_like(S1); SP = np.zeros_like(S1); SPOS = np.zeros_like(S1)
    lfc = de["lfc"].values.astype(np.float64); pde = de["pde"].values.astype(np.float64); M = np.stack([memb[u][idx] for u in unit_names]).astype(np.float64)    # U x cells
    for k in range(len(covs)):
        S1[:, k] += M @ lfc[:, k, :]; S2[:, k] += M @ (lfc[:, k, :] ** 2); SP[:, k] += M @ pde[:, k, :]; SPOS[:, k] += M @ (lfc[:, k, :] > 0)
    N += M.sum(1); del de, lfc, pde; print(f"chunk {ci + 1}/{int(np.ceil(len(keep) / CH))} done, {time.time() - t0:.0f}s", flush=True)
n3 = N[:, None, None]; mean = S1 / n3; sd = np.sqrt(np.maximum(S2 / n3 - mean ** 2, 0)); np.savez_compressed(f"{OUT}/mrvi_de_unit_gene_stats.npz", units=np.asarray(unit_names), covariates=np.asarray(covs), genes=genes, n_cells=N, mean_lfc=mean.astype(np.float32), sd_lfc=sd.astype(np.float32), mean_pde=(SP / n3).astype(np.float32), frac_pos=(SPOS / n3).astype(np.float32))
print("saved; units:", U, "covariates:", covs, "total seconds:", round(time.time() - t0), flush=True)
