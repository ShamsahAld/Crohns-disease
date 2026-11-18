"""29_pca_robustness.py - plan step 10b: robustness of the biopsy-level pseudobulk PCA 'inflammation' signal.

Replicates 05c_pseudobulk_pca.py (counts of the QC-passed cells of a biopsy summed -> log2 CPM -> most variable genes -> PCA over biopsies), for two cell sets:
 'all cells'    = all QC-passed cells (25 % MT cap), 'immune cells' = immune compartment with pct_counts_mt <= 15 (the 61,538 cells used for MrVI);
with variants: drop the extreme biopsy HA98AC_inf (about 0 % epithelium) from the PCA; number of top variable genes 1,000 / 2,000 (original) / 5,000.
The PCA is computed over all biopsies with >= 20 cells in the set (as in 05c: groups inflamed, non-inflamed, adjacent, control), the tests use the 33 primary biopsies (inflamed vs non-inflamed).
Test per PC1-3: OLS of the standardised PC score on status + site (primary); + biologic; + age + ileocolonic location; + patient-clustered standard errors (cluster-robust OLS).  eta^2 of status within the primary set as well.
Output: out/pca_robustness.csv, out/pca_robustness_scores.csv (scores of the primary biopsies for every variant).
"""
import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, statsmodels.api as sm
sys.path.insert(0, "."); from common import HVG_EXCLUDE_REGEX, SEX_GENES
PROC = os.path.expanduser("~/.claude-science-scratch/gse266616/processed"); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True); MIN_CELLS = 20
adata = sc.read_h5ad(f"{PROC}/04_embedded.h5ad"); o = adata.obs.copy(); counts = adata.layers["counts"].tocsr(); print(adata.shape, o.status.astype(str).value_counts().to_dict(), flush=True)
sid_all = o["sample_id"].astype(str).values; excl = (adata.var_names.str.match(HVG_EXCLUDE_REGEX) | adata.var_names.isin(SEX_GENES)); genes = adata.var_names
meta = o[["sample_id", "patient_id", "site", "status", "biologic", "age", "cd_history"]].astype({"sample_id": str}).drop_duplicates("sample_id").set_index("sample_id"); meta["ileocolonic"] = meta.cd_history.astype(str).str.startswith("ileocolonic").astype(float)
prim = meta.index[meta.status.astype(str).isin(["inflamed", "non_inflamed"])]; print("primary biopsies:", len(prim), flush=True)
masks = {"all cells": np.ones(len(o), bool), "immune cells (MT<=15%)": ((o["compartment"].astype(str) == "immune") & (o["pct_counts_mt"] <= 15)).values}
def pseudobulk(mask):
    sid = sid_all[mask]; ids = pd.Index(sorted(set(sid))); n_per = pd.Series(sid).value_counts().reindex(ids); ids = ids[n_per.values >= MIN_CELLS]; keep = np.isin(sid, ids); codes = pd.Categorical(sid[keep], categories=ids).codes
    M = sp.csr_matrix((np.ones(keep.sum()), (codes, np.arange(keep.sum()))), shape=(len(ids), keep.sum())); return pd.DataFrame((M @ counts[np.where(mask)[0][keep]]).toarray(), index=ids, columns=genes)
def pca(pb, n_top):
    cpm = pb.div(pb.sum(1), axis=0) * 1e6; ok = ((cpm >= 1).mean(0) >= 0.25).values & ~excl; L = np.log2(cpm.loc[:, ok] + 1); top = L.var(0).sort_values(ascending=False).index[:n_top]; X = L[top].values - L[top].values.mean(0)
    U, S, Vt = np.linalg.svd(X, full_matrices=False); ev = (S ** 2) / (S ** 2).sum(); return pd.DataFrame(U[:, :3] * S[:3], index=L.index, columns=["PC1", "PC2", "PC3"]), ev[:3]
def fit(y, cols, groups=None):
    X = sm.add_constant(pd.DataFrame(cols)); m = sm.OLS(y, X).fit(cov_type="cluster", cov_kwds={"groups": groups}) if groups is not None else sm.OLS(y, X).fit(); return float(m.params["infl"]), float(m.pvalues["infl"])
rows, scs = [], []
for sname, mask in masks.items():
    pb0 = pseudobulk(mask)
    for drop in ("none", "HA98AC_inf"):
        pb = pb0.drop(index=[drop], errors="ignore") if drop != "none" else pb0
        for n_top in (1000, 2000, 5000):
            sc_, ev = pca(pb, n_top); pr = [s for s in prim if s in sc_.index]; d = meta.loc[pr]; P = sc_.loc[pr]
            infl = (d.status.astype(str) == "inflamed").astype(float).values; AC = (d.site.astype(str) == "AC").astype(float).values; bio = pd.get_dummies(d.biologic.astype(str), drop_first=True).astype(float).values; age = ((d.age - d.age.mean()) / d.age.std()).values; ic = d.ileocolonic.values; pid = d.patient_id.astype(str).values
            for pcn in ("PC1", "PC2", "PC3"):
                y = P[pcn].values; y = (y - y.mean()) / y.std(); base = dict(infl=infl, AC=AC); b1, p1 = fit(y, base); b2, p2 = fit(y, {**base, **{f"b{i}": bio[:, i] for i in range(bio.shape[1])}}); b3, p3 = fit(y, {**base, "age": age, "ileocolonic": ic}); b4, p4 = fit(y, base, groups=pid)
                eta = float(((y[infl == 1].mean() - y.mean()) ** 2 * (infl == 1).sum() + (y[infl == 0].mean() - y.mean()) ** 2 * (infl == 0).sum()) / ((y - y.mean()) ** 2).sum())
                rows.append(dict(set=sname, dropped=drop, n_top_genes=n_top, n_biopsies_in_pca=len(sc_), n_primary=len(pr), PC=pcn, var_pct=round(100 * ev[int(pcn[2]) - 1], 1), beta_site=b1, p_site=p1, beta_site_biologic=b2, p_site_biologic=p2, beta_site_age_loc=b3, p_site_age_loc=p3, beta_patient_clustered=b4, p_patient_clustered=p4, eta2_status=eta))
            scs.append(P.assign(set=sname, dropped=drop, n_top_genes=n_top).reset_index())
pd.DataFrame(rows).to_csv(f"{OUT}/pca_robustness.csv", index=False); pd.concat(scs).to_csv(f"{OUT}/pca_robustness_scores.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(), flush=True)
