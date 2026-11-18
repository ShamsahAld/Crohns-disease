"""Step 5c - SAMPLE-LEVEL (pseudobulk) PCA: do whole biopsies group by inflammation, site or patient?

Counts of all QC-passed cells of a sample (optionally only of one lineage) are summed -> one profile per biopsy -> log2 CPM ->
top 2,000 most variable genes -> PCA of the samples. Because the design is between patients (only HA27 contributes both states, from different
sites), this is the appropriate unit for a first look at inflamed vs non-inflamed: cells from one biopsy are not independent replicates.
Outputs: out/pseudobulk_pca_scores.csv, out/pseudobulk_pca_eta2.csv, out/pseudobulk_composition.csv, out/fig05c_*.png
"""
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, seaborn as sns
from scipy.stats import spearmanr
from common import *
style()
adata = sc.read_h5ad(f"{PROC}/04_embedded.h5ad")
o = adata.obs.copy(); counts = adata.layers["counts"].tocsr()
samples = pd.DataFrame({"sample_id": o["sample_id"].astype(str).unique()}).merge(
    o[["sample_id", "patient_id", "site", "status", "platform_id", "biologic", "age"]].astype({"sample_id": str}).drop_duplicates("sample_id"), on="sample_id").set_index("sample_id")
excl = adata.var_names.str.match(HVG_EXCLUDE_REGEX) | adata.var_names.isin(SEX_GENES)
MIN_CELLS = int(os.environ.get('MIN_CELLS', 20)); MIN_SAMPLES = int(os.environ.get('MIN_SAMPLES', 8))

def pseudobulk(mask):
    sid = o["sample_id"].astype(str).values[mask]; ids = pd.Index(sorted(set(sid)))
    n_per = pd.Series(sid).value_counts().reindex(ids)
    ids = ids[n_per.values >= MIN_CELLS]
    keep = np.isin(sid, ids); codes = pd.Categorical(sid[keep], categories=ids).codes
    M = sp.csr_matrix((np.ones(keep.sum()), (codes, np.arange(keep.sum()))), shape=(len(ids), keep.sum()))
    X = M @ counts[np.where(mask)[0][keep]]
    return pd.DataFrame(X.toarray(), index=ids, columns=adata.var_names), n_per.reindex(ids)

def pca_samples(pb, n_top=2000):
    cpm = pb.div(pb.sum(1), axis=0) * 1e6
    ok = ((cpm >= 1).mean(0) >= 0.25).values & ~excl
    L = np.log2(cpm.loc[:, ok] + 1); top = L.var(0).sort_values(ascending=False).index[:n_top]
    X = L[top].values - L[top].values.mean(0)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    k = min(6, len(S)); ev = (S ** 2) / (S ** 2).sum()
    return pd.DataFrame(U[:, :k] * S[:k], index=L.index, columns=[f"PC{i+1}" for i in range(k)]), ev[:k], len(top)

sets = {"all cells": np.ones(len(o), bool), "immune cells": (o["compartment"] == "immune").values,
        "T/NK cells": (o["coarse_lineage"] == "T_NK").values, "myeloid cells": (o["coarse_lineage"] == "myeloid").values, "plasma cells": (o["coarse_lineage"] == "plasma").values}
res, rows, scores_all = {}, [], []
for name, m in sets.items():
    pb, npc = pseudobulk(m)
    if len(pb) < MIN_SAMPLES: print("skip", name, len(pb)); continue
    sc_, ev, ng = pca_samples(pb); sc_ = sc_.join(samples); sc_["n_cells"] = npc.values; res[name] = (sc_, ev)
    print(f"{name}: {len(pb)} samples (>= {MIN_CELLS} cells), {ng} genes, PC1-3 variance {np.round(ev[:3]*100,1)}%", flush=True)
    for i in range(3):
        for f in ["status", "site", "platform_id", "biologic"]:      # patient_id omitted: 28 patients in 33 samples, eta^2 ~ 1 by construction
            rows.append(dict(set=name, n_samples=len(pb), PC=f"PC{i+1}", var_pct=round(ev[i] * 100, 1), factor=f, eta2=eta2(sc_[f"PC{i+1}"].values, sc_[f].values)))
    scores_all.append(sc_.assign(set=name).reset_index())
et = pd.DataFrame(rows); et.to_csv(f"{OUT}/pseudobulk_pca_eta2.csv", index=False)
pd.concat(scores_all).to_csv(f"{OUT}/pseudobulk_pca_scores.csv", index=False)

# composition (fraction of QC-passed cells per coarse lineage) and its correlation with all-cell pseudobulk PCs
comp = pd.crosstab(o["sample_id"].astype(str), o["coarse_lineage"], normalize="index"); comp.to_csv(f"{OUT}/pseudobulk_composition.csv")
sc_all, ev_all = res["all cells"]
cor = pd.DataFrame({l: [spearmanr(sc_all[f"PC{i+1}"], comp.loc[sc_all.index, l]).statistic for i in range(3)] for l in comp.columns if comp[l].std() > 0}, index=["PC1", "PC2", "PC3"]).round(2)
cor.to_csv(f"{OUT}/pseudobulk_pc_vs_composition_spearman.csv"); print("Spearman(PC, lineage fraction), all-cell pseudobulk:\n", cor.to_string(), flush=True)

# ---- is inflamed vs non-inflamed associated with the sample-level PCs?  OLS of the standardised PC score on status, adjusting for site (and biologic).
from scipy import stats
def ols_coef_p(y, X, j):
    if np.linalg.matrix_rank(X) < X.shape[1] or len(y) - X.shape[1] < 1: return np.nan, np.nan      # status collinear with a covariate (cannot be adjusted)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]; resid = y - X @ beta; dof = len(y) - X.shape[1]
    cov = (resid @ resid / dof) * np.linalg.inv(X.T @ X); t = beta[j] / np.sqrt(cov[j, j]); return float(beta[j]), float(2 * stats.t.sf(abs(t), dof))
trows = []
for nm, sc_, ev in [(k, v[0], v[1]) for k, v in res.items()]:
    stv = (sc_["status"] == "inflamed").astype(float).values; sitev = (sc_["site"] == "AC").astype(float).values
    bio = pd.get_dummies(sc_["biologic"], drop_first=True).astype(float).values
    for pcn in ["PC1", "PC2", "PC3"]:
        y = sc_[pcn].values; y = (y - y.mean()) / y.std(); one = np.ones(len(y))
        b1, p1 = ols_coef_p(y, np.c_[one, stv], 1); b2, p2 = ols_coef_p(y, np.c_[one, stv, sitev], 1); b3, p3 = ols_coef_p(y, np.c_[one, stv, sitev, bio], 1)
        trows.append(dict(set=nm, PC=pcn, n_samples=len(y), inflamed_minus_non_inflamed_SD=round(b1, 2), p_unadjusted=p1, p_adj_site=p2, p_adj_site_biologic=p3))
tt = pd.DataFrame(trows); tt["q_BH_adj_site"] = np.nan; _ok = tt["p_adj_site"].notna()
if _ok.any(): tt.loc[_ok, "q_BH_adj_site"] = stats.false_discovery_control(tt.loc[_ok, "p_adj_site"].values, method="bh")
tt.to_csv(f"{OUT}/pseudobulk_status_tests.csv", index=False); print(tt.to_string(), flush=True)

# ---- descriptive composition by status (median share per biopsy; not a test - differential abundance is roadmap step 8)
st_of = samples["status"]; comp_all = comp.join(st_of).groupby("status").median().T * 100
imm = o[o["compartment"] == "immune"]; comp_imm = pd.crosstab(imm["sample_id"].astype(str), imm["coarse_lineage"].astype(str), normalize="index").join(st_of).groupby("status").median().T * 100
pd.concat({"share_of_all_cells_pct": comp_all, "share_of_immune_cells_pct": comp_imm}, axis=1).round(1).to_csv(f"{OUT}/composition_median_by_status.csv")

# ================================================================= figures
def scatter(ax, sc_, ev, title, i=0, j=1):
    for st, col in [(g, STATUS_COL[g]) for g in STATUS_ORDER]:
        for site, mk in {"TI": "o", "AC": "^"}.items():
            d = sc_[(sc_.status == st) & (sc_.site == site)]
            ax.scatter(d[f"PC{i+1}"], d[f"PC{j+1}"], s=22, c=col, marker=mk, edgecolors="white", linewidths=0.4, zorder=3)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0); ax.axvline(0, color="0.85", lw=0.5, zorder=0)
    ax.set_xlabel(f"PC{i+1} ({ev[i]*100:.0f}%)"); ax.set_ylabel(f"PC{j+1} ({ev[j]*100:.0f}%)"); ax.set_title(f"{title} (n={len(sc_)})", loc="left"); ax.margins(0.08)
names = list(res)
fig = plt.figure(figsize=(7.2, 5.0)); gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.42)
for k, nm in enumerate(names[:5]):
    ax = fig.add_subplot(gs[k // 3, k % 3]); scatter(ax, *res[nm], title=nm[0].upper() + nm[1:])
axl = fig.add_subplot(gs[1, 2]); axl.axis("off")
hs = [mpl.lines.Line2D([], [], marker="o", ls="", color=STATUS_COL[s_], label=s_.replace("_", "-"), ms=4) for s_ in STATUS_ORDER] + \
     [mpl.lines.Line2D([], [], marker=m, ls="", color="0.4", label=l, ms=4) for l, m in [("terminal ileum", "o"), ("ascending colon", "^")]]
axl.legend(handles=hs, fontsize=6.5, loc="center left", handletextpad=0.3, title="Each point = one biopsy", title_fontsize=6.5)
for ax_, l in zip(fig.axes[:5], "abcde"): panel_letter(ax_, l)
savefig(fig, "fig05c_pseudobulk_pca.png")

# variance explained by factor, per pseudobulk set, PC1-3
facs = [("status", "inflamed vs\nnon-inflamed"), ("site", "tissue\nsite"), ("biologic", "biologic\ntherapy"), ("platform_id", "platform")]
mat = pd.DataFrame({f"{f}_{p}": [et[(et["set"] == n) & (et.factor == f) & (et.PC == p)].eta2.iloc[0] for n in names] for f, _ in facs for p in ["PC1", "PC2", "PC3"]},
                   index=[f"{n[0].upper() + n[1:]} (n={int(et[et['set'] == n].n_samples.iloc[0])})" for n in names])
mat.columns = ["PC1", "PC2", "PC3"] * len(facs)
fig, ax = plt.subplots(figsize=(7.2, 2.6)); heat_grouped(ax, mat, [(lab, 3) for _, lab in facs], vmax=0.6, fs=6, title="Share of each sample-level PC explained by each factor (eta²)")
savefig(fig, "fig05c_pseudobulk_variance_explained.png")

order = samples.assign(status=pd.Categorical(samples.status, ["non_inflamed", "inflamed"])).sort_values(["status", "site"]).index.tolist()
fig, ax = plt.subplots(figsize=(7.2, 2.9)); bottom = np.zeros(len(order))
for lin in LINEAGE_COL:
    if lin in comp.columns:
        v = comp.loc[order, lin].values * 100; ax.bar(range(len(order)), v, bottom=bottom, color=LINEAGE_COL[lin], width=0.8, label=lin); bottom += v
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=5.5)
for t, sid in zip(ax.get_xticklabels(), order): t.set_color(SITE_COL[samples.loc[sid, "site"]])
for i, sid in enumerate(order): ax.add_patch(mpl.patches.Rectangle((i - 0.4, 101.5), 0.8, 2.5, color=STATUS_COL[samples.loc[sid, "status"]], clip_on=False))
ax.set_ylim(0, 100); ax.set_ylabel("Cells passing QC (%)"); ax.margins(x=0.01)
ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=6)
ax.set_title("Exploratory cell-lineage composition per biopsy (top bar: blue = non-inflamed, orange = inflamed; label colour = site)", loc="left", fontsize=7, pad=14)
savefig(fig, "fig05c_composition_per_sample.png")
