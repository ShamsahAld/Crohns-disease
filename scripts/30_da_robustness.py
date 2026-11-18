"""30_da_robustness.py - plan step 10c: robustness of the differential-abundance headline results (states and plasma IgG share).

Variants for the per-biopsy proportion tests (same primary test as 20_da_proportions.py: site-stratified van Elteren, 10,000 within-site permutations, BH over 45 states):
  base | drop HA98AC_inf (extreme biopsy, ~0 % epithelium) | drop HA27 biopsies (the only patient with both states) | additionally per state (all 33 biopsies): OLS-logit with site + biologic + age + ileocolonic disease location,
  OLS-logit with patient-clustered standard errors (cluster-robust), and MrVI density permutation test without HA98AC_inf (21_mrvi_da_permutation.py logic, 32 biopsies).
"""
import importlib.util, os, sys, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
spec = importlib.util.spec_from_file_location("da20", "20_da_proportions.py"); da = importlib.util.module_from_spec(spec); spec.loader.exec_module(da)
cs, bm = da.load(); d_all = da.design(bm); bm["ileocolonic"] = bm.cd_history.astype(str).str.startswith("ileocolonic").astype(float); d_all["ileocolonic"] = bm.ileocolonic
conf_ok = cs.confidence != "exclude"; C_all = da.counts_table(cs, conf_ok).reindex(d_all.index).fillna(0); states = list(C_all.columns); rows = []
variants = {"base": [], "drop HA98AC_inf": ["HA98AC_inf"], "drop HA27 (both biopsies)": [s for s in d_all.index if s.startswith("HA27")]}
for vn, drop in variants.items():
    keep = [s for s in d_all.index if s not in drop]; C = C_all.loc[keep]; N = C.sum(1); r = da.run_level(C, N, d_all.loc[keep], states, vn); rows.append(r[["level", "state", "log2fc_site_std", "p_primary", "q_primary", "p_ols_perm", "p_adj_biologic_age"]].rename(columns={"level": "variant"}))
res = pd.concat(rows, ignore_index=True)
# covariate-adjusted and cluster-robust OLS on all 33 biopsies
N = C_all.sum(1).values; Y = np.log((C_all.values + 0.5) / (N[:, None] - C_all.values + 0.5)); d = d_all; pid = d.patient.values; ex = []
for k, s in enumerate(states):
    y = Y[:, k]; X1 = sm.add_constant(d[["AC", "infl", "anti_TNF", "anti_IL", "age", "ileocolonic"]]); m1 = sm.OLS(y, X1).fit(); X2 = sm.add_constant(d[["AC", "infl"]]); m2 = sm.OLS(y, X2).fit(cov_type="cluster", cov_kwds={"groups": pid})
    ex.append(dict(state=s, beta_adj_biologic_age_loc=m1.params["infl"], p_adj_biologic_age_loc=m1.pvalues["infl"], beta_patient_cluster=m2.params["infl"], p_patient_cluster=m2.pvalues["infl"]))
ex = pd.DataFrame(ex); ex["q_adj_biologic_age_loc"] = multipletests(ex.p_adj_biologic_age_loc, method="fdr_bh")[1]; ex["q_patient_cluster"] = multipletests(ex.p_patient_cluster, method="fdr_bh")[1]
res.to_csv("robust_da_proportion_variants.csv", index=False); ex.to_csv("robust_da_adjusted_models.csv", index=False)
# plasma IgG share: drop HA98AC_inf and covariate-adjusted
iso = pd.read_csv("da_plasma_isotype_per_biopsy.csv", index_col=0); out_iso = []
for vn, drop in [("base", []), ("drop HA98AC_inf", ["HA98AC_inf"])]:
    idx = [s for s in d_all.index if s not in drop]; v = iso.loc[idx, "IgG"].values[:, None]; _, p = da.van_elteren(v, d_all.loc[idx]); out_iso.append(dict(variant=vn, measure="IgG share (pseudobulk)", p_van_elteren=float(p[0])))
yv = np.log(np.clip(iso.loc[d_all.index, "IgG"].values, 0.005, 0.995) / (1 - np.clip(iso.loc[d_all.index, "IgG"].values, 0.005, 0.995))); m1 = sm.OLS(yv, sm.add_constant(d_all[["AC", "infl", "anti_TNF", "anti_IL", "age", "ileocolonic"]])).fit(); m2 = sm.OLS(yv, sm.add_constant(d_all[["AC", "infl"]])).fit(cov_type="cluster", cov_kwds={"groups": pid})
out_iso += [dict(variant="OLS logit + biologic + age + ileocolonic", measure="IgG share (pseudobulk)", p_van_elteren=float(m1.pvalues["infl"])), dict(variant="OLS logit, patient-clustered SE", measure="IgG share (pseudobulk)", p_van_elteren=float(m2.pvalues["infl"]))]
pd.DataFrame(out_iso).to_csv("robust_igg_variants.csv", index=False); print(pd.DataFrame(out_iso).round(4).to_string(index=False))
head = ["Mast, activated (IEG-high)", "CD8 tissue-resident memory", "Mast, core (KIT+ MS4A2+)", "pDC", "gd/innate-like IEL", "Mast low-information", "Plasmablast (B-lineage cluster)", "Plasma low-UMI (ribosomal-high)", "Neutrophil", "Inflammatory macrophage/monocyte"]
w = res[res.state.isin(head)].pivot(index="state", columns="variant", values="p_primary").loc[head]; q = res[res.state.isin(head)].pivot(index="state", columns="variant", values="q_primary").loc[head]
print("p (van Elteren) by variant\n", w.round(4).to_string()); print("q by variant\n", q.round(3).to_string()); print(ex.set_index("state").loc[head][["p_adj_biologic_age_loc", "q_adj_biologic_age_loc", "p_patient_cluster", "q_patient_cluster"]].round(4).to_string())
