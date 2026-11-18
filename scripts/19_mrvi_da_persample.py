"""19_mrvi_da_persample.py - plan step 7a: per-sample cell densities from the trained immune MrVI model, for differential abundance.

Input : $GSE_PROC/06_immune_mrvi_input_hvg.h5ad + trained model $GSE_PROC/mrvi_immune_model (step 4; nothing is retrained).
Method: model.differential_abundance(sample_cov_keys=[status, site, biologic], compute_log_enrichment=True). The dataset has, for every cell, log p(u | aggregated posterior of biopsy s)
        for each of the 33 biopsies (variable 'log_probs'), plus the library's own group-level log-probs/enrichments (own biopsy omitted from the aggregate posterior).
Output: out/da_persample_logprobs.npz (log_probs float32 cells x biopsies, barcodes, biopsy order), out/da_native_<var>.npy for every other variable, da_info.json, sample_info.csv.
        The local script recomputes the group-level scores from log_probs so that permutation nulls over biopsy labels are possible."""
import os, json, time, numpy as np, pandas as pd, anndata as ad, torch, scvi
from scvi.external import MRVI
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 8))); scvi.settings.seed = 0
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); a.layers["counts"] = a.X.copy(); print("input", a.shape, flush=True)
for c_ in ["sample_id", "status", "site", "biologic", "immune_lineage"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"])
a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(["none", "anti_TNF", "anti_IL12_IL23"])
MRVI.setup_anndata(a, layer="counts", sample_key="sample_id"); model = MRVI.load(f"{PROC}/mrvi_immune_model", adata=a, accelerator="cpu")
t0 = time.time(); ds = model.differential_abundance(sample_cov_keys=["status", "site", "biologic"], compute_log_enrichment=True, batch_size=512); sec = time.time() - t0
print("DA seconds", round(sec, 1)); print(ds, flush=True)
lp = ds["log_probs"]; np.savez_compressed(f"{OUT}/da_persample_logprobs.npz", log_probs=lp.values.astype(np.float32), barcodes=np.asarray(ds["cell_name"].values, dtype=str), samples=np.asarray(lp["sample"].values, dtype=str))
info = dict(seconds=sec, vars={}, coords={})
for k, v in ds.data_vars.items():
    if k == "log_probs": continue
    np.save(f"{OUT}/da_native_{k}.npy", np.asarray(v.values, dtype=np.float32)); info["vars"][k] = dict(dims=list(v.dims), shape=list(v.shape))
    for d in v.dims:
        if d not in ("cell_name",) and d in v.coords: info["coords"][f"{k}:{d}"] = [str(x) for x in v[d].values]
json.dump(info, open(f"{OUT}/da_info.json", "w"), indent=1); model.sample_info.to_csv(f"{OUT}/sample_info.csv"); print(json.dumps(info)[:1500])
