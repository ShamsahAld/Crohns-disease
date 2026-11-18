"""28_mrvi_hvg_sensitivity.py - plan step 10a: does the MrVI abundance result depend on the number of highly variable genes?

For N in (2000, 5000) the HVG list is rebuilt with exactly the step-1 procedure (Seurat v3 'paper' vote-first rule on raw counts, batch_key = biopsy, same exclusions: MT/ribo/Hb/V(D)J/sex genes, all IG*/JCHAIN genes, epithelial-dominant genes);
the N = 3000 call is repeated first and must reproduce the existing list (assert).  Then MrVI is retrained as in step 4 (biopsy = sample, batch 128, early stopping, seed 0, CPU) and per-biopsy cell densities are exported for the differential-abundance test.
Input : $GSE_PROC/06_immune_mrvi_input.h5ad (61,538 immune cells, all genes, layers['counts']), common.py, epithelial_dominant_genes.csv, hvg_immune_mt15.csv.
Output: out/hvg<N>_da_persample_logprobs.npz, out/hvg<N>_train_summary.json, out/hvg<N>_gene_list.csv.
"""
import os, re, sys, json, time, numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, torch, scvi
from scvi.external import MRVI
sys.path.insert(0, "."); from common import HVG_EXCLUDE_REGEX, SEX_GENES
PROC = os.path.expanduser("~/.claude-science-scratch/gse266616/processed"); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True); torch.set_num_threads(8); scvi.settings.seed = 0
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input.h5ad"); print(a.shape, flush=True)
for c_ in ["sample_id", "status", "site", "biologic", "immune_lineage"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"]); a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(["none", "anti_TNF", "anti_IL12_IL23"])
epi = set(pd.read_csv("epithelial_dominant_genes.csv").iloc[:, 0].astype(str)); excl = a.var_names.str.match(HVG_EXCLUDE_REGEX) | a.var_names.isin(SEX_GENES) | a.var_names.str.match(r"^(IG[HKL]|JCHAIN$)") | a.var_names.isin(epi)
tmp = a[:, ~excl].copy(); tmp.X = tmp.layers["counts"].copy(); print("candidate genes:", tmp.n_vars, flush=True); ref = set(pd.read_csv("hvg_immune_mt15.csv").iloc[:, 0].astype(str))
def hvgs(n):
    t = tmp.copy(); sc.pp.highly_variable_genes(t, flavor="seurat_v3_paper", n_top_genes=n, batch_key="sample_id", layer="counts", span=0.3); return list(t.var_names[t.var["highly_variable"].values])
g3 = hvgs(3000); print("N=3000 reproduces existing list:", set(g3) == ref, "| overlap", len(set(g3) & ref), flush=True); assert set(g3) == ref
for N in (2000, 5000):
    t0 = time.time(); g = hvgs(N); pd.DataFrame(dict(gene=g)).to_csv(f"{OUT}/hvg{N}_gene_list.csv", index=False); print(f"N={N}: {len(g)} genes; overlap with 3000 list: {len(set(g) & ref)}", flush=True)
    b = a[:, g].copy(); b.layers["counts"] = a[:, g].layers["counts"].copy(); b.X = b.layers["counts"].copy(); MRVI.setup_anndata(b, layer="counts", sample_key="sample_id"); model = MRVI(b)
    model.train(max_epochs=300, batch_size=128, early_stopping=True, accelerator="cpu", check_val_every_n_epoch=1); ev = pd.concat([pd.Series(v.values.ravel(), name=k) for k, v in model.history.items()], axis=1)["elbo_validation"].dropna().values
    S = dict(n_hvg=N, epochs_run=int(len(ev)), best_epoch=int(np.argmin(ev)) + 1, val_elbo_first=float(ev[0]), val_elbo_last=float(ev[-1]), train_seconds=round(time.time() - t0)); print("TRAIN", json.dumps(S), flush=True)
    ds = model.differential_abundance(sample_cov_keys=None, batch_size=512); lp = ds["log_probs"]
    np.savez_compressed(f"{OUT}/hvg{N}_da_persample_logprobs.npz", log_probs=lp.values.astype(np.float32), barcodes=np.asarray(ds["cell_name"].values, dtype=str), samples=np.asarray(lp["sample"].values, dtype=str))
    u = model.get_latent_representation(give_z=False); np.save(f"{OUT}/hvg{N}_latent_u.npy", u.astype(np.float32)); S["seconds_total"] = round(time.time() - t0); json.dump(S, open(f"{OUT}/hvg{N}_train_summary.json", "w")); print("done", N, flush=True)
