"""21_mrvi_da_permutation.py - plan step 7c: MrVI differential abundance (inflamed vs non-inflamed) with a biopsy-label permutation null.

Inputs: da_persample_logprobs.npz (job 42f405db: log p(u_cell | aggregated posterior of biopsy s) for 61,538 cells x 33 biopsies), cell_state_labels.csv.gz, biopsy_metadata_mrvi.csv.
Score of a cell = log mean_{inflamed biopsies, own biopsy omitted} p - log mean_{non-inflamed biopsies, own omitted} p; this reproduces MrVI's own status_log_enrichs exactly (max |diff| = 0).
The raw score is heavy-tailed (cells in latent regions with no support in one group reach |score| > 10), so the PRIMARY state statistic is bounded: net fraction = mean over the state's cells of
  sign(score) * 1[|score| > ln 2]  (= share of cells >=2-fold denser in inflamed minus share >=2-fold denser in non-inflamed), cell-weighted; the biopsy-balanced version averages the per-biopsy net fractions
  over biopsies with >=10 cells in the state. Mean score (cell-mean and biopsy-balanced) and the median score are reported as secondary summaries; the mean has a permutation p as well.
Null: status labels permuted among biopsies within site (TI, AC), 10,000 times, own-biopsy omission redone for every permutation. Because the site mix of the inflamed group is fixed by design, the null is not centred on zero for site-dependent states;
the site-adjusted score is therefore observed minus null mean, and p is two-sided around the null mean.  BH across analysed states.
"""
import os, json, numpy as np, pandas as pd, scipy.sparse as sps
from statsmodels.stats.multitest import multipletests
RNG = np.random.default_rng(20261110); B = int(os.environ.get("DA_B", 10000)); OUT = os.environ.get("DA_OUT", ".")
DROP = [x for x in os.environ.get("DA_DROP", "").split(",") if x]; TAG = os.environ.get("DA_TAG", "da_mrvi"); z = np.load(os.environ.get("DA_NPZ", "da_persample_logprobs.npz")); L = z["log_probs"].astype(np.float64); samples = z["samples"].astype(str); bc = z["barcodes"].astype(str)
cs = pd.read_csv("cell_state_labels.csv.gz").set_index("barcode").loc[bc].reset_index()
if DROP:   # leave the named biopsies out of the aggregated posteriors and of the tested cells
    kc = ~cs.sample_id.isin(DROP).values; kb = ~np.isin(samples, DROP); L = L[kc][:, kb]; samples = samples[kb]; bc = bc[kc]; cs = cs[kc].reset_index(drop=True)
bm = pd.read_csv("biopsy_metadata_mrvi.csv").set_index("sample_id").loc[samples]
own = pd.Index(samples).get_indexer(cs.sample_id.values); infl0 = (bm.status.values == "inflamed"); ac = (bm.site.values == "AC"); nS = len(samples)
E = np.exp(L - L.max(1, keepdims=True)); E[np.arange(len(own)), own] = 0.0; OWN = np.zeros((len(own), nS), bool); OWN[np.arange(len(own)), own] = True
def score(infl):
    ind = np.column_stack([infl, ~infl]).astype(float); S = E @ ind; n = ind.sum(0)[None, :] - OWN @ ind; return np.log(S[:, 0] / n[:, 0]) - np.log(S[:, 1] / n[:, 1])
def weights(labels, min_cells=10):
    """sparse S x cells weight matrices: cell-mean and biopsy-balanced (biopsies with >= min_cells cells in the state)"""
    names = sorted(pd.unique(labels)); Wc, Wb, keep = [], [], []
    for nm in names:
        m = np.where(labels == nm)[0]; wc = np.full(len(m), 1.0 / len(m)); b = own[m]; cnt = np.bincount(b, minlength=nS); ok = cnt >= min_cells
        if ok.sum() >= 3:
            w = np.where(ok[b], 1.0 / np.maximum(cnt[b], 1), 0.0); w = w / w.sum()
        else: w = np.full(len(m), np.nan)
        Wc.append(sps.csr_matrix((wc, (np.zeros(len(m), int), m)), shape=(1, len(own)))); Wb.append(sps.csr_matrix((np.nan_to_num(w), (np.zeros(len(m), int), m)), shape=(1, len(own)))); keep.append(int(ok.sum()))
    return names, sps.vstack(Wc).tocsr(), sps.vstack(Wb).tocsr(), np.array(keep)
def run(labels, level):
    names, Wc, Wb, nb = weights(labels); s0 = score(infl0); oc, ob = Wc @ s0, Wb @ s0; net = lambda s_: np.sign(s_) * (np.abs(s_) > np.log(2)); n0 = net(s0); fc_, fb_ = Wc @ n0, Wb @ n0
    med = np.array([np.median(s0[labels == nm]) for nm in names]); nc = np.zeros((B, len(names))); nbp = np.zeros((B, len(names))); ntc = np.zeros((B, len(names))); ntb = np.zeros((B, len(names))); strata = [np.where(ac == k)[0] for k in (False, True)]
    for i in range(B):
        inf = infl0.copy()
        for ix in strata: inf[ix] = RNG.permutation(infl0[ix])
        s = score(inf); nc[i] = Wc @ s; nbp[i] = Wb @ s; n_ = net(s); ntc[i] = Wc @ n_; ntb[i] = Wb @ n_
    out = pd.DataFrame(dict(level=level, state=names, n_cells=np.asarray(Wc.getnnz(1)), n_biopsies_ge10=nb, score_cell=oc, null_mean_cell=nc.mean(0), score_biopsy=ob, null_mean_biopsy=nbp.mean(0)))
    out["median_score"] = med; out["net_frac_cell"] = fc_; out["net_null_mean_cell"] = ntc.mean(0); out["net_frac_cell_site_adj"] = fc_ - ntc.mean(0); out["net_frac_biopsy"] = fb_; out["net_frac_biopsy_site_adj"] = fb_ - ntb.mean(0)
    out["p_primary"] = [(1 + np.sum(np.abs(ntc[:, k] - ntc[:, k].mean()) >= abs(fc_[k] - ntc[:, k].mean()) - 1e-12)) / (B + 1) for k in range(len(names))]
    out["p_net_biopsy"] = [(1 + np.sum(np.abs(ntb[:, k] - np.nanmean(ntb[:, k])) >= abs(fb_[k] - np.nanmean(ntb[:, k])) - 1e-12)) / (B + 1) if nb[k] >= 3 else np.nan for k in range(len(names))]
    out["score_cell_site_adj"] = out.score_cell - out.null_mean_cell; out["score_biopsy_site_adj"] = out.score_biopsy - out.null_mean_biopsy
    out["p_cell"] = [(1 + np.sum(np.abs(nc[:, k] - nc[:, k].mean()) >= abs(oc[k] - nc[:, k].mean()) - 1e-12)) / (B + 1) for k in range(len(names))]
    out["p_biopsy"] = [(1 + np.sum(np.abs(nbp[:, k] - np.nanmean(nbp[:, k])) >= abs(ob[k] - np.nanmean(nbp[:, k])) - 1e-12)) / (B + 1) if nb[k] >= 3 else np.nan for k in range(len(names))]
    out["q_primary"] = multipletests(out.p_primary, method="fdr_bh")[1]; ok_ = out.p_net_biopsy.notna(); out.loc[ok_, "q_net_biopsy"] = multipletests(out.loc[ok_, "p_net_biopsy"], method="fdr_bh")[1]
    out["q_cell"] = multipletests(out.p_cell, method="fdr_bh")[1]; ok = out.p_biopsy.notna(); out.loc[ok, "q_biopsy"] = multipletests(out.loc[ok, "p_biopsy"], method="fdr_bh")[1]; return out
if __name__ == "__main__":
    conf_ok = (cs.confidence != "exclude").values; res = [run(cs.state.values[conf_ok] if False else np.where(conf_ok, cs.state.values, "EXCLUDED"), "primary: non-excluded states")]
    res[0] = res[0][res[0].state != "EXCLUDED"]; res.append(run(cs.immune_lineage.values.astype(object), "lineage within immune cells"))
    r = pd.concat(res, ignore_index=True); r.to_csv(f"{OUT}/{TAG}_state_scores.csv", index=False)
    pd.DataFrame(dict(barcode=bc, sample_id=cs.sample_id.values, state=cs.state.values, mrvi_da_score_inflamed_vs_non=score(infl0).round(4))).to_csv(f"{OUT}/{TAG}_cell_scores.csv.gz", index=False)
    print(r[r.level.str.startswith("primary")].sort_values("p_primary").head(14)[["state", "n_cells", "median_score", "net_frac_cell", "net_frac_cell_site_adj", "p_primary", "q_primary", "net_frac_biopsy_site_adj", "p_net_biopsy", "q_net_biopsy", "p_cell"]].round(4).to_string()); print(r[r.level.str.startswith("lineage")][["state", "n_cells", "median_score", "net_frac_cell_site_adj", "p_primary", "q_primary", "p_cell"]].round(4).to_string())
