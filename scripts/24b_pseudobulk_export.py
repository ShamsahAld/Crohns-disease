"""24b_pseudobulk_export.py - plan step 8a: pseudobulk raw-count sums per biopsy for every lineage and for each state with enough cells (all 30,675 genes, counts layer).

Input : $GSE_PROC/06_immune_mrvi_input.h5ad (raw counts in layers['counts']), cell_state_labels.csv.gz (step 5; barcode -> state, confidence, lineage).
Units : 5 lineages (cells in non-excluded states) and every non-excluded state that has >= 10 biopsies with >= 20 cells and >= 800 cells overall.
Output: out/pseudobulk_counts.npz  (for each unit 'U::<name>': int32 matrix biopsies x genes; 'N::<name>': cells per biopsy), genes, samples, units table (unit, kind, lineage, state, confidence, n_cells, n_biopsies_ge20).
"""
import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input.h5ad"); X = a.layers["counts"].tocsr().astype(np.float64) if sp.issparse(a.layers["counts"]) else sp.csr_matrix(a.layers["counts"]); print(a.shape, "counts max", X.max(), "integer:", bool(np.all(X.data == np.round(X.data))), flush=True)
lab = pd.read_csv("cell_state_labels.csv.gz").set_index("barcode"); assert set(a.obs_names) == set(lab.index), "barcode mismatch"; lab = lab.loc[a.obs_names]
samples = sorted(lab.sample_id.unique()); sidx = pd.Index(samples).get_indexer(lab.sample_id.values); ok = (lab.confidence != "exclude").values; res = {}; rows = []
def pseudobulk(mask):
    m = np.where(mask)[0]; G = sp.csr_matrix((np.ones(len(m)), (sidx[m], np.arange(len(m)))), shape=(len(samples), len(m))); S = (G @ X[m]).toarray(); n = np.bincount(sidx[m], minlength=len(samples)); return S.astype(np.int32), n
for L in ["T_NK", "plasma", "B", "myeloid", "mast"]:
    S, n = pseudobulk(ok & (lab.immune_lineage.values == L)); res["U::lineage:" + L] = S; res["N::lineage:" + L] = n; rows.append(dict(unit="lineage:" + L, kind="lineage", lineage=L, state=L, confidence="", n_cells=int(n.sum()), n_biopsies_ge20=int((n >= 20).sum())))
for st, g in lab[ok].groupby("state"):
    n_ = np.bincount(sidx[(lab.state.values == st)], minlength=len(samples))
    if (n_ >= 20).sum() >= 10 and n_.sum() >= 800:
        S, n = pseudobulk(lab.state.values == st); res["U::state:" + st] = S; res["N::state:" + st] = n; rows.append(dict(unit="state:" + st, kind="state", lineage=g.immune_lineage.iloc[0], state=st, confidence=g.confidence.iloc[0], n_cells=int(n.sum()), n_biopsies_ge20=int((n >= 20).sum())))
units = pd.DataFrame(rows); units.to_csv(f"{OUT}/pseudobulk_units.csv", index=False); print(units.to_string(), flush=True)
np.savez_compressed(f"{OUT}/pseudobulk_counts.npz", genes=np.asarray(a.var_names, dtype=str), samples=np.asarray(samples, dtype=str), **res); print("saved", len(units), "units", flush=True)
