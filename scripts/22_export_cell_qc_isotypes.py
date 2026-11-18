"""22_export_cell_qc_isotypes.py - plan step 7d helper: per-cell QC values and immunoglobulin constant-region counts for the 61,538 immune cells (full-gene file, raw counts).
Output: out/immune_cell_qc_isotype_counts.csv.gz (barcode, QC columns present in obs, raw UMI of IG constant genes + JCHAIN, total heavy-chain UMI) and out/obs_columns.txt."""
import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input.h5ad", backed="r"); print(a.shape, list(a.obs.columns), flush=True); open(f"{OUT}/obs_columns.txt", "w").write(" ".join(a.obs.columns))
genes = ["IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGHD", "IGHE", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7", "JCHAIN"]; gi = [i for i, g in enumerate(a.var_names) if g in genes]; names = list(a.var_names[gi]); print("found", names, flush=True)
X = a.X[:, gi]; X = X.toarray() if sp.issparse(X) else np.asarray(X); layers = list(a.layers.keys()); print("layers", layers, "X dtype", X.dtype, "max", X.max(), flush=True)
df = pd.DataFrame(X, columns=names, index=a.obs_names); keep = [c for c in a.obs.columns if any(k in c.lower() for k in ["mt", "total", "n_genes", "umi", "doublet", "count"])]
df = pd.concat([a.obs[keep], df], axis=1); df.index.name = "barcode"; df.to_csv(f"{OUT}/immune_cell_qc_isotype_counts.csv.gz"); print(df.shape, df.columns.tolist(), flush=True)
