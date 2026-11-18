"""24a_mrvi_de_probe.py - read-only probe of MRVI.differential_expression on 1,000 mast cells: signature, output structure, timing, memory."""
import os, time, json, inspect, resource, numpy as np, pandas as pd, anndata as ad, torch, scvi
from scvi.external import MRVI
print(inspect.signature(MRVI.differential_expression)); print(inspect.getdoc(MRVI.differential_expression)[:5500], flush=True)
PROC = os.path.expanduser("~/.claude-science-scratch/gse266616/processed"); torch.set_num_threads(8)
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); a.layers["counts"] = a.X.copy()
for c_ in ["sample_id", "status", "site", "biologic", "immune_lineage"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"]); a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(["none", "anti_TNF", "anti_IL12_IL23"])
MRVI.setup_anndata(a, layer="counts", sample_key="sample_id"); model = MRVI.load(f"{PROC}/mrvi_immune_model", adata=a, accelerator="cpu")
sub = a[a.obs.immune_lineage == "mast"][:1000].copy(); t0 = time.time()
de = model.differential_expression(adata=sub, sample_cov_keys=["status", "site", "biologic"], store_lfc=True, use_vmap=False, mc_samples=50) if "mc_samples" in inspect.signature(MRVI.differential_expression).parameters else model.differential_expression(adata=sub, sample_cov_keys=["status", "site", "biologic"], store_lfc=True, use_vmap=False)
print("DE seconds for 1000 cells:", round(time.time() - t0, 1), "| peak RSS GB:", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2)); print(de, flush=True)
for k, v in de.data_vars.items(): print(k, v.dims, v.dtype, float(np.nanmin(v.values)), float(np.nanmax(v.values)))
print("coords:", {k: list(map(str, v.values[:6])) for k, v in de.coords.items() if k not in ("cell_name", "gene")})
cov = [c for c in de.coords["covariate"].values if "status" in str(c)][0]; pv = de["padj"].sel(covariate=cov).values; ef = de["effect_size"].sel(covariate=cov).values; print("status cov:", cov, "| cells x genes padj<0.05 share:", float((pv < 0.05).mean()), "| |effect_size| median:", float(np.nanmedian(np.abs(ef))))
