"""23_plasma_isotypes.py - plan step 7d: plasma-cell immunoglobulin heavy-chain isotype use, inflamed vs non-inflamed.

Inputs: immune_cell_qc_isotype_counts.csv.gz (job 5a9f3975; Ig values were exported from the log-normalised X and are converted back to raw UMI as expm1(X) * total_counts / 1e4, then rounded; 98.2 % of entries were within 0.01 of an integer before rounding),
        cell_state_labels.csv.gz, biopsy_metadata_mrvi.csv.
Cells: plasma-lineage cells in non-excluded states.  Classes: IgA = IGHA1+IGHA2, IgG = IGHG1..4, IgM = IGHM (IgD and IgE are negligible and omitted).
(1) Biopsy pseudobulk: class share of summed heavy-chain UMI over plasma cells (sensitive to ambient Ig RNA, which the cell-level call below only partly removes).
(2) Cell-level call: a cell with >= 10 heavy-chain UMI whose top class holds >= 70 % of them; per biopsy the share of called cells in each class (biopsies with >= 30 called cells).
Tests as in 20_da_proportions.py: site-stratified van Elteren with 10,000 within-site permutations (primary), OLS-logit permutation, unadjusted Mann-Whitney.  Effect = site-standardised difference in mean share (percentage points).
"""
import importlib.util, numpy as np, pandas as pd
from scipy import stats
spec = importlib.util.spec_from_file_location("da20", "20_da_proportions.py"); da = importlib.util.module_from_spec(spec); spec.loader.exec_module(da)
q = pd.read_csv("immune_cell_qc_isotype_counts.csv.gz", index_col=0); cs = pd.read_csv("cell_state_labels.csv.gz").set_index("barcode"); bm = pd.read_csv("biopsy_metadata_mrvi.csv").set_index("sample_id"); d = da.design(bm)
ig = ["IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"]; raw = pd.DataFrame(np.round(np.expm1(q[ig].values.astype(np.float64)) * q.total_counts.values[:, None] / 1e4), index=q.index, columns=ig)
m = (cs.immune_lineage == "plasma") & (cs.confidence != "exclude"); P = raw.loc[cs.index[m]]; P["sample_id"] = cs.sample_id[m].values
P["IgA"] = P.IGHA1 + P.IGHA2; P["IgG"] = P.IGHG1 + P.IGHG2 + P.IGHG3 + P.IGHG4; P["IgM"] = P.IGHM; cls = ["IgA", "IgG", "IgM"]; tot = P[cls].sum(1)
print("plasma cells:", len(P), "| with >=10 heavy-chain UMI:", int((tot >= 10).sum()), "| called (top class >=70%):", int(((tot >= 10) & (P[cls].max(1) / tot.replace(0, np.nan) >= 0.7)).sum()))
rows = []; bio = P.groupby("sample_id")[cls + ["IGHA1", "IGHA2"]].sum().reindex(d.index); hc = bio[cls].sum(1)
share = bio[cls].div(hc, axis=0); share["IgA2_within_IgA"] = bio.IGHA2 / (bio.IGHA1 + bio.IGHA2)
call = pd.Series("none", index=P.index); ok = (tot >= 10) & (P[cls].max(1) / tot.replace(0, np.nan) >= 0.7); call[ok] = P.loc[ok, cls].idxmax(1)
ct = pd.crosstab(P.sample_id, call).reindex(d.index).fillna(0); ncalled = ct[cls].sum(1); cshare = ct[cls].div(ncalled, axis=0); cshare.columns = [c + "_cellcall" for c in cls]; cshare[ncalled < 30] = np.nan
out = pd.concat([share, cshare, ncalled.rename("n_called_cells"), hc.rename("heavy_chain_UMI"), bm[["status", "site", "patient_id"]]], axis=1); out.to_csv("da_plasma_isotype_per_biopsy.csv")
def test(vals, label, level):
    v = vals.copy(); keep = v.notna(); dd = d[keep.values]; y = v[keep].values[:, None]; Y = np.log((np.clip(y, 0.005, 0.995)) / (1 - np.clip(y, 0.005, 0.995)))
    tve, pve = da.van_elteren(y, dd); _, pol = da.perm_p(Y, dd); lfc, mi, mn = da.std_log2fc(y, dd)
    return dict(level=level, measure=label, n_biopsies=int(keep.sum()), mean_share_inflamed=mi[0], mean_share_non_inflamed=mn[0], site_std_diff_pp=100 * (mi[0] - mn[0]), p_primary=pve[0], p_ols_perm=pol[0], mwu_p=stats.mannwhitneyu(v[keep][d.infl[keep.values] == 1], v[keep][d.infl[keep.values] == 0]).pvalue)
rows = [test(out[c], c, "biopsy pseudobulk") for c in cls + ["IgA2_within_IgA"]] + [test(out[c + "_cellcall"], c + " (cell calls)", "cell-level call") for c in cls]
res = pd.DataFrame(rows); from statsmodels.stats.multitest import multipletests; res["q_primary"] = multipletests(res.p_primary, method="fdr_bh")[1]; res.to_csv("da_plasma_isotype_tests.csv", index=False); print(res.round(4).to_string(index=False))
