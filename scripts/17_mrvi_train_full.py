"""17_mrvi_train_full.py - roadmap step 6d: train MrVI on all immune cells (15 % MT cap) and write diagnostics.

Input : $GSE_PROC/06_immune_mrvi_input_hvg.h5ad  (61,538 cells x 3,000 vote-first HVGs, X = raw counts, obs: sample_id, status, site, biologic, immune_lineage ...)
Model : scvi.external.MRVI, sample_key = biopsy, no batch_key; batch_size 128, max_epochs 300 with early stopping on validation ELBO; seed 0; CPU (torch is the CPU wheel).
Saves : model (PROC/mrvi_immune_model and out/), latent u (sample-unaware) and z (sample-aware), UMAP of u, local sample distances per lineage,
        training history, neighbour-mixing table (MrVI u vs PCA), summary JSON.   Nothing outside scratch is written."""
import os, sys, time, json, traceback, numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, torch, scvi
from scvi.external import MRVI
from sklearn.neighbors import NearestNeighbors
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
MAX_EPOCHS = int(os.environ.get("MRVI_EPOCHS", 300)); BATCH = int(os.environ.get("MRVI_BATCH", 128)); SEED = 0
torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 8))); scvi.settings.seed = SEED; sc.settings.n_jobs = torch.get_num_threads()
S = dict(scvi=scvi.__version__, torch=torch.__version__, threads=torch.get_num_threads(), batch_size=BATCH, max_epochs=MAX_EPOCHS, seed=SEED)
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); print("input:", a.shape, flush=True)
a.layers["counts"] = a.X.copy()
for c_ in ["sample_id", "status", "site", "biologic", "immune_lineage"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"])
a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(["none", "anti_TNF", "anti_IL12_IL23"]); a.obs["all_immune"] = "all_immune"
MRVI.setup_anndata(a, layer="counts", sample_key="sample_id"); model = MRVI(a)
t0 = time.time(); model.train(max_epochs=MAX_EPOCHS, batch_size=BATCH, early_stopping=True, accelerator="cpu", check_val_every_n_epoch=1); S["train_seconds"] = time.time() - t0
hist = pd.concat([pd.Series(v.values.ravel(), name=k) for k, v in model.history.items()], axis=1); hist.to_csv(f"{OUT}/mrvi_full_history.csv", index=False)
ev = hist["elbo_validation"].dropna().values; S.update(epochs_run=int(len(ev)), best_epoch=int(np.argmin(ev)) + 1, val_elbo_first=float(ev[0]), val_elbo_last=float(ev[-1]), val_elbo_min=float(ev.min()),
        val_change_last10_pct=float(100 * (ev[-1] - ev[-11]) / abs(ev[-11])) if len(ev) > 11 else None, stopped_early=bool(len(ev) < MAX_EPOCHS))
print("TRAIN", json.dumps(S), flush=True)
model.save(f"{PROC}/mrvi_immune_model", overwrite=True, save_anndata=False); model.save(f"{OUT}/mrvi_model", overwrite=True, save_anndata=False)

# ---- latent spaces
u = model.get_latent_representation(give_z=False); z = model.get_latent_representation(give_z=True); a.obsm["X_mrvi_u"] = u; a.obsm["X_mrvi_z"] = z; print("u", u.shape, "z", z.shape, flush=True)
sc.pp.neighbors(a, use_rep="X_mrvi_u", n_neighbors=30, random_state=0); sc.tl.umap(a, min_dist=0.3, random_state=0)
o = a.obs; cells = pd.DataFrame(dict(barcode=o.index, sample_id=o.sample_id.values, immune_lineage=o.immune_lineage.values, status=o.status.values, site=o.site.values, biologic=o.biologic.values,
                                     umap1=a.obsm["X_umap"][:, 0], umap2=a.obsm["X_umap"][:, 1]))
for i in range(u.shape[1]): cells[f"u{i+1}"] = u[:, i]
cells.round(4).to_csv(f"{OUT}/mrvi_immune_cells_latent.csv.gz", index=False)

# ---- local sample distances per lineage (and over all immune cells)
try:
    for gb, name in [("immune_lineage", "by_lineage"), ("all_immune", "all_immune")]:
        d = model.get_local_sample_distances(batch_size=256, keep_cell=False, groupby=gb); print(name, "distance dataset:", {k: (tuple(v.shape), list(v.dims)) for k, v in d.data_vars.items()}, flush=True)
        for var in d.data_vars: d[var].to_dataframe(name="distance").reset_index().to_csv(f"{OUT}/mrvi_sample_distances_{name}_{var}.csv", index=False)
except Exception: traceback.print_exc()

# ---- neighbour mixing in MrVI u-space vs plain PCA of the same cells
def knn_idx(X, k=30): return NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]
q = ad.AnnData(sp.csr_matrix(a.layers["counts"]), obs=o[["sample_id", "immune_lineage"]].copy()); sc.pp.normalize_total(q, target_sum=1e4); sc.pp.log1p(q); sc.pp.scale(q, max_value=10); sc.tl.pca(q, n_comps=30, svd_solver="arpack", random_state=0)
sid, lin = np.asarray(o.sample_id.astype(str), dtype=object), np.asarray(o.immune_lineage.astype(str), dtype=object); rows = []      # object arrays: pandas 3 Arrow strings cannot be 2-D indexed
for nm, X in [("MrVI u (sample-unaware)", u), ("PCA of 3,000 HVGs", q.obsm["X_pca"])]:
    nn = knn_idx(X); same_b = (sid[nn] == sid[:, None]).mean(1); same_l = (lin[nn] == lin[:, None]).mean(1)
    for L in ["all"] + sorted(set(lin)):
        m = np.ones(len(lin), bool) if L == "all" else lin == L
        cnt = pd.Series(sid[m]).value_counts() / m.sum(); exp = pd.Series(sid[m]).map(cnt).values.astype(float)
        rows.append(dict(space=nm, lineage=L, n_cells=int(m.sum()), same_biopsy_neighbour_frac=same_b[m].mean(), expected_if_mixed=exp.mean(), ratio_observed_to_expected=same_b[m].mean() / exp.mean(), same_lineage_neighbour_frac=same_l[m].mean()))
mix = pd.DataFrame(rows); mix.round(4).to_csv(f"{OUT}/knn_mixing_mrvi_vs_pca.csv", index=False); print(mix.round(3).to_string(index=False), flush=True)
a.write_h5ad(f"{PROC}/06b_immune_mrvi_latent.h5ad", compression="lzf")
json.dump(S, open(f"{OUT}/mrvi_full_summary.json", "w"), indent=1, default=str); print("FULL TRAIN DONE", flush=True)
