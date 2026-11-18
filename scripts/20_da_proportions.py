"""20_da_proportions.py - plan step 7b: differential abundance of immune cell states, inflamed vs non-inflamed, from per-biopsy proportions.

Inputs (workspace): cell_state_labels.csv.gz (step 5), biopsy_metadata_mrvi.csv (step 1), qc_summary_per_sample.csv (step 3: n_kept = all cells per biopsy at the 25 % MT cap).
Unit = biopsy (33: 17 inflamed, 16 non-inflamed; 11+6 vs 8+8 TI+AC).  Primary denominator = immune cells in non-excluded states of the same biopsy.
Primary test  = site-stratified Wilcoxon (van Elteren) rank statistic on per-biopsy proportions, p from 10,000 permutations of status within site (TI / AC). Chosen as primary because many states are sparse (zero counts in several biopsies).\nConcordance   = OLS of logit((count+0.5)/(other+0.5)) on status + site with the same within-site permutation (p_ols_perm).
Effect size   = log2 ratio of site-standardised mean proportions (mean of the TI and AC group means, inflamed / non-inflamed).
Also reported: site-adjusted quasi-binomial GLM (variance-weighted, Pearson-scaled), unadjusted Mann-Whitney, patient-clustered Gaussian GEE on the logit scale (exchangeable), OLS + biologic + age, leave-one-biopsy-out range, TI-only and AC-only direction.
Multiple testing: Benjamini-Hochberg across the analysed states (exploratory throughout).
"""
import os, json, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf, warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore"); RNG = np.random.default_rng(20261109); B = 10000; OUT = os.environ.get("DA_OUT", "."); 

def load():
    cs = pd.read_csv("cell_state_labels.csv.gz"); bm = pd.read_csv("biopsy_metadata_mrvi.csv").set_index("sample_id"); qc = pd.read_csv("qc_summary_per_sample.csv").set_index("sample_id")
    bm["n_all_cells"] = qc.n_kept.reindex(bm.index); return cs, bm

def design(bm):
    d = pd.DataFrame(index=bm.index); d["infl"] = (bm.status == "inflamed").astype(float); d["AC"] = (bm.site == "AC").astype(float)
    d["anti_TNF"] = (bm.biologic == "anti_TNF").astype(float); d["anti_IL"] = (bm.biologic == "anti_IL12_IL23").astype(float); d["age"] = (bm.age - bm.age.mean()) / bm.age.std(); d["patient"] = bm.patient_id; return d

def ols_t(Y, X, j=-1):
    """vectorised OLS t statistic of column j of X for each column of Y"""
    n, p = X.shape; XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ X.T @ Y; res = Y - X @ beta; s2 = (res ** 2).sum(0) / (n - p); se = np.sqrt(s2 * XtXi[j, j]); return beta[j], beta[j] / se

def perm_p(Y, d, B=B):
    """site-stratified permutation p for the status coefficient, with all states at once (Frisch-Waugh)"""
    n = len(d); Z = np.column_stack([np.ones(n), d.AC.values]); Pz = Z @ np.linalg.pinv(Z); yr = Y - Pz @ Y; y2 = (yr ** 2).sum(0); st = d.infl.values.copy()
    def tstat(s):
        sr = s - Pz @ s; ss = sr @ sr; b = (sr @ yr) / ss; rss = y2 - (sr @ yr) ** 2 / ss; return b / np.sqrt(rss / (n - 3) / ss)
    tobs = tstat(st); strata = [np.where(d.AC.values == k)[0] for k in (0, 1)]; ge = np.zeros(Y.shape[1])
    for _ in range(B):
        s = st.copy()
        for ix in strata: s[ix] = RNG.permutation(st[ix])
        ge += np.abs(tstat(s)) >= np.abs(tobs) - 1e-12
    return tobs, (ge + 1) / (B + 1)

def van_elteren(P, d, B=B):
    """site-stratified Wilcoxon (van Elteren) statistic with a permutation p (status shuffled within TI / AC); ties get average ranks, zeros are handled naturally"""
    n = len(d); strata = [np.where(d.AC.values == k)[0] for k in (0, 1)]; R = np.zeros_like(P); w = np.zeros(n)
    for ix in strata: R[ix] = stats.rankdata(P[ix], axis=0); w[ix] = 1.0 / (len(ix) + 1)
    st = d.infl.values; tobs = (st * w) @ R; mean0 = sum(st[ix].sum() / 2.0 for ix in strata); ge = np.zeros(P.shape[1])
    for _ in range(B):
        s_ = st.copy()
        for ix in strata: s_[ix] = RNG.permutation(st[ix])
        ge += np.abs((s_ * w) @ R - mean0) >= np.abs(tobs - mean0) - 1e-9
    return tobs - mean0, (ge + 1) / (B + 1)

def std_log2fc(P, d):
    """log2 ratio of site-standardised mean proportion (equal weight TI and AC); P = biopsies x states"""
    m = {}
    for st_ in (1.0, 0.0):
        m[st_] = np.mean([P[(d.infl.values == st_) & (d.AC.values == k)].mean(0) for k in (0, 1)], axis=0)
    return np.log2((m[1.0] + 1e-4) / (m[0.0] + 1e-4)), m[1.0], m[0.0]

def run_level(C, N, d, names, label):
    """C: biopsies x states counts, N: denominators per biopsy (Series)"""
    Pm = C.values / N.values[:, None]; Y = np.log((C.values + 0.5) / (N.values[:, None] - C.values + 0.5)); X = np.column_stack([np.ones(len(d)), d.AC.values, d.infl.values])
    beta, t = ols_t(Y, X); _, pperm = perm_p(Y, d); ve_stat, pve = van_elteren(Pm, d); lfc, mi, mn = std_log2fc(Pm, d); rows = []
    for k, s in enumerate(names):
        c, n_ = C.iloc[:, k].values, N.values; r = dict(level=label, state=s, n_cells=int(c.sum()), n_biopsies_present=int((c > 0).sum()), mean_prop_inflamed=mi[k], mean_prop_non_inflamed=mn[k], log2fc_site_std=lfc[k], ols_logit_beta=beta[k], ols_t=t[k], p_ols_perm=pperm[k], ve_stat=ve_stat[k], p_primary=pve[k])
        r["p_ols"] = 2 * stats.t.sf(abs(t[k]), len(d) - 3)
        r["mwu_p"] = stats.mannwhitneyu(Pm[d.infl.values == 1, k], Pm[d.infl.values == 0, k]).pvalue
        dd = d.assign(c=c, o=n_ - c)
        try:
            Xg = sm.add_constant(dd[["AC", "infl"]]); g = sm.GLM(c / n_, Xg, family=sm.families.Binomial(), var_weights=n_).fit(scale="X2"); r.update(glm_logOR=g.params["infl"], glm_se=g.bse["infl"], glm_p=g.pvalues["infl"], glm_dispersion=g.scale)
        except Exception: r.update(glm_logOR=np.nan, glm_se=np.nan, glm_p=np.nan, glm_dispersion=np.nan)
        try:
            ge = sm.GEE(Y[:, k], Xg, groups=dd.patient.values, family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Exchangeable()).fit(); r.update(gee_logit_beta=ge.params["infl"], gee_p=ge.pvalues["infl"])
        except Exception: r.update(gee_logit_beta=np.nan, gee_p=np.nan)
        X2 = sm.add_constant(dd[["AC", "infl", "anti_TNF", "anti_IL", "age"]]); o2 = sm.OLS(Y[:, k], X2).fit(); r.update(p_adj_biologic_age=o2.pvalues["infl"], beta_adj_biologic_age=o2.params["infl"])
        for site, kk in (("TI", 0), ("AC", 1)):
            m = d.AC.values == kk; r[f"log2fc_{site}_only"] = np.log2((Pm[m & (d.infl.values == 1), k].mean() + 1e-4) / (Pm[m & (d.infl.values == 0), k].mean() + 1e-4))
        loo = []
        for i in range(len(d)):
            keep = np.arange(len(d)) != i; b_, t_ = ols_t(Y[keep][:, [k]], X[keep]); loo.append((b_[0], 2 * stats.t.sf(abs(t_[0]), keep.sum() - 3)))
        loo = np.array(loo); r.update(loo_p_max=loo[:, 1].max(), loo_sign_consistent=bool(np.all(np.sign(loo[:, 0]) == np.sign(beta[k]))))
        rows.append(r)
    out = pd.DataFrame(rows); out["q_primary"] = multipletests(out.p_primary, method="fdr_bh")[1]; out["q_ols_perm"] = multipletests(out.p_ols_perm, method="fdr_bh")[1]; out["q_glm"] = multipletests(out.glm_p.fillna(1), method="fdr_bh")[1]; out["q_mwu"] = multipletests(out.mwu_p, method="fdr_bh")[1]
    out["direction"] = np.where(out.log2fc_site_std > 0, "higher in inflamed", "lower in inflamed"); return out

def counts_table(cs, keep_mask, level="state"):
    sub = cs[keep_mask]; return pd.crosstab(sub.sample_id, sub[level])

def one_per_patient(C, N, d, names, draws=300, B_=500):
    """sensitivity to the 5 patients with two biopsies (HA27 has both states): draw one biopsy per patient, 300 times; report the fraction of draws with p<0.05 and the median p (van Elteren permutation, B=500)"""
    pid = d.patient.values; Pm = C.values / N.values[:, None]; ps = []
    groups = {p: np.where(pid == p)[0] for p in np.unique(pid)}
    for _ in range(draws):
        idx = np.sort([RNG.choice(v) for v in groups.values()]); _, pv = van_elteren(Pm[idx], d.iloc[idx], B=B_); ps.append(pv)
    ps = np.array(ps); return pd.DataFrame(dict(state=names, op_median_p=np.median(ps, 0), op_frac_p_lt_0_05=(ps < 0.05).mean(0), op_n_biopsies=len(groups)))

if __name__ == "__main__":
    cs, bm = load(); d = design(bm); allres = []
    conf_ok = cs.confidence != "exclude"
    variants = {"primary: non-excluded states": conf_ok, "high+medium states only": cs.confidence.isin(["high", "medium"])}
    for vn, mask in variants.items():
        C = counts_table(cs, mask).reindex(d.index).fillna(0); N = C.sum(1); r = run_level(C, N, d, list(C.columns), vn); allres.append(r); print(vn, C.shape, "min denominator", int(N.min()), flush=True)
    # all immune denominator including excluded states: states = non-excluded, denominator = all immune cells
    Call = counts_table(cs, conf_ok).reindex(d.index).fillna(0); Nall = cs.groupby("sample_id").size().reindex(d.index); allres.append(run_level(Call, Nall, d, list(Call.columns), "denominator = all immune cells (incl. excluded)"))
    # all-cell denominator (immune state count / all cells kept at the 25 % cap): a different question (immune share of the whole biopsy)
    allres.append(run_level(Call, bm.n_all_cells, d, list(Call.columns), "denominator = all cells in biopsy (25% MT cap)"))
    # stricter MT cap: keep only labelled cells with pct_counts_mt <= 10 % (the earlier cap analysis rejected 10 % for the main analysis because it removes 14.6 % of cells and 31 % of mast cells; used here only as a robustness check)
    qcc = pd.read_csv("immune_cell_qc_isotype_counts.csv.gz", index_col=0, usecols=["barcode", "pct_counts_mt"]).pct_counts_mt; mt_ok = cs.barcode.map(qcc).values <= 10
    C10 = counts_table(cs, conf_ok & mt_ok).reindex(d.index).fillna(0); allres.append(run_level(C10, C10.sum(1), d, list(C10.columns), "MT cap 10% (stricter)"))
    # states as a fraction of their OWN lineage, biopsies with >= 20 cells of that lineage only (removes lineage-level composition shifts, e.g. the CD8 Trm share falling because mast/myeloid cells rise)
    lin_of = cs[conf_ok].drop_duplicates("state").set_index("state").immune_lineage; parts = []
    for L_ in sorted(lin_of.unique()):
        sts = [x for x in Call.columns if lin_of[x] == L_]
        if len(sts) < 2: continue
        tot = Call[sts].sum(1); ok = tot >= 20; Cs = Call.loc[ok, sts]; r_ = run_level(Cs, Cs.sum(1), d.loc[ok], sts, "denominator = own lineage"); r_["n_biopsies_in_test"] = int(ok.sum()); parts.append(r_)
    w_ = pd.concat(parts, ignore_index=True); w_["q_primary"] = multipletests(w_.p_primary, method="fdr_bh")[1]; w_["q_ols_perm"] = multipletests(w_.p_ols_perm, method="fdr_bh")[1]; allres.append(w_)
    # lineage level
    mask = np.ones(len(cs), bool); Cl = counts_table(cs, mask, "immune_lineage").reindex(d.index).fillna(0); allres.append(run_level(Cl, Cl.sum(1), d, list(Cl.columns), "lineage within immune cells"))
    C0 = counts_table(cs, conf_ok).reindex(d.index).fillna(0); op = one_per_patient(C0, C0.sum(1), d, list(C0.columns)); op.to_csv(f"{OUT}/da_one_biopsy_per_patient.csv", index=False)
    res = pd.concat(allres, ignore_index=True); res = res.merge(op.assign(level="primary: non-excluded states"), on=["level", "state"], how="left"); res.to_csv(f"{OUT}/da_proportion_tests.csv", index=False)
    # per-biopsy counts and proportions for reference
    C = counts_table(cs, conf_ok).reindex(d.index).fillna(0).astype(int); C.to_csv(f"{OUT}/da_state_counts_per_biopsy.csv")
    print(res[res.level.str.startswith("primary")].sort_values("p_primary").head(14)[["state", "n_cells", "n_biopsies_present", "log2fc_site_std", "p_primary", "q_primary", "p_ols_perm", "glm_p", "gee_p", "mwu_p", "p_adj_biologic_age", "loo_p_max"]].round(4).to_string())
