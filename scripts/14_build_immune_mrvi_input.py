"""14_build_immune_mrvi_input.py (v2) - roadmap step 6a: immune-only input for MrVI, with an ambient-epithelial-RNA safeguard.

Input : $GSE_PROC/04_embedded.h5ad (all 136,284 QC-passed cells, lineage labels, raw counts)  and  $GSE_PROC/05_immune.h5ad (65,191 immune cells)
        pseudobulk_pca_scores.csv (staged; immune-cell biopsy-level PCs from step 5c)
Part 0: data-defined EPITHELIAL-DOMINANT genes = genes with >=8x higher pseudobulk CPM in epithelial cells than in immune cells and CPM >= 5 in epithelium.
        Such genes are what ambient RNA from the surrounding tissue puts into immune droplets.
Part 1: keep immune cells with pct_counts_mt <= MT_CAP (15; decision 2026-11-09).
Part 2: ambient index per cell = share of counts from epithelial-dominant genes; per biopsy (median over immune cells, and per lineage). Tests: inflamed vs non-inflamed
        (Mann-Whitney, site-adjusted OLS); Spearman with the immune biopsy-PCA scores; does the PC2~status association survive adjusting for the ambient index?
Part 3: 3,000 HVGs among the immune cells (Seurat v3 'paper' rule = rank genes by the NUMBER OF BIOPSIES in which they are among each biopsy's top-3,000 most variable genes, ties by median rank;
        scanpy's default flavor 'seurat_v3' instead sorts by median rank first, which lets genes variable in only 1-2 biopsies in; raw counts, batch_key = sample_id), excluding MT-, ribosomal, haemoglobin, V(D)J segments, sex genes (as before)
        PLUS all immunoglobulin-locus genes (IGH*/IGK*/IGL*, JCHAIN) and the epithelial-dominant genes from Part 0.
Part 4: write 06_immune_mrvi_input.h5ad (all genes, log-normalised X, layers['counts']) and 06_immune_mrvi_input_hvg.h5ad (HVG subset, X = raw counts).
Only files under ~/.claude-science-scratch are written; the user's own folders are untouched.
"""
import os, sys, json, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.formula.api as smf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import PROC, OUT, HVG_EXCLUDE_REGEX, SEX_GENES
MT_CAP = float(os.environ.get("MT_CAP", 15)); N_HVG = 3000; RATIO, MIN_CPM = 8.0, 5.0
IG_ALL = r"^(IG[HKL]|JCHAIN$)"
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
IMM = ["T_NK", "plasma", "myeloid", "mast", "B"]

# ---------------- Part 0: epithelial-dominant genes
full = ad.read_h5ad(f"{PROC}/04_embedded.h5ad"); lin = full.obs["coarse_lineage"].astype(str).values
cnt = sp.csr_matrix(full.layers["counts"]); genes = np.asarray(full.var_names)
def cpm(mask):
    v = np.asarray(cnt[np.where(mask)[0]].sum(0)).ravel(); return v / v.sum() * 1e6
cpm_epi, cpm_imm = cpm(lin == "epithelial"), cpm(np.isin(lin, IMM))
ratio = (cpm_epi + 1e-3) / (cpm_imm + 1e-3)
epi_dom = (ratio >= RATIO) & (cpm_epi >= MIN_CPM)
pd.DataFrame(dict(gene=genes, cpm_epithelial=cpm_epi.round(2), cpm_immune=cpm_imm.round(2), ratio=ratio.round(2), epithelial_dominant=epi_dom)).loc[epi_dom].sort_values("cpm_epithelial", ascending=False).to_csv(f"{OUT}/epithelial_dominant_genes.csv", index=False)
print("epithelial-dominant genes:", int(epi_dom.sum()), "of", len(genes), "| cumulative share of epithelial CPM:", round(float(cpm_epi[epi_dom].sum() / 1e6), 3), flush=True)
del full, cnt

# ---------------- Part 1: immune cells at the MT cap
a = sc.read_h5ad(f"{PROC}/05_immune.h5ad"); n0 = a.n_obs
assert "counts" in a.layers and a.obs["compartment"].eq("immune").all()
a = a[a.obs["pct_counts_mt"] <= MT_CAP].copy(); print(f"immune cells: {n0} -> {a.n_obs} at MT <= {MT_CAP}%", flush=True)
if MT_CAP == 15: assert a.n_obs == 61538, a.n_obs
for c_ in ["status", "site", "sample_id", "patient_id", "platform_id", "biologic", "cd_history"]:
    if c_ in a.obs: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["immune_lineage"] = a.obs["coarse_lineage"].astype(str); assert a.obs["sample_id"].nunique() == 33
assert np.array_equal(np.asarray(a.var_names), genes)

# ---------------- Part 2: ambient index
Xc = sp.csr_matrix(a.layers["counts"]); tot = np.asarray(Xc.sum(1)).ravel()
a.obs["epi_dom_share"] = np.asarray(Xc[:, np.where(epi_dom)[0]].sum(1)).ravel() / np.maximum(tot, 1)
o = a.obs; bio = o.groupby("sample_id", observed=True).agg(status=("status", "first"), site=("site", "first"), patient_id=("patient_id", "first"), n_cells=("status", "size"),
        ambient_median=("epi_dom_share", "median"), ambient_mean=("epi_dom_share", "mean")).reset_index()
for c_ in ["sample_id", "status", "site", "patient_id"]: bio[c_] = bio[c_].astype(str)      # categoricals would make .map() return a non-numeric categorical
for l in IMM: bio[f"ambient_median_{l}"] = bio.sample_id.map(o[o.immune_lineage == l].groupby("sample_id", observed=True).epi_dom_share.median().rename(index=str)).astype(float)
bio.to_csv(f"{OUT}/ambient_index_per_biopsy.csv", index=False)
rows = []
for col in ["ambient_median"] + [f"ambient_median_{l}" for l in IMM]:
    d = bio[["status", "site", col]].dropna().rename(columns={col: "y"}); x, y = d[d.status == "non_inflamed"].y, d[d.status == "inflamed"].y
    f = smf.ols("y ~ C(status, Treatment('non_inflamed')) + C(site)", data=d).fit()
    rows.append(dict(measure=col, n_biopsies=len(d), median_non=x.median(), median_inf=y.median(), mwu_p=mannwhitneyu(x, y).pvalue, siteadj_diff=f.params.filter(like="status").iloc[0], siteadj_p=f.pvalues.filter(like="status").iloc[0]))
pd.DataFrame(rows).round(6).to_csv(f"{OUT}/ambient_index_status_tests.csv", index=False); print(pd.DataFrame(rows).round(5).to_string(index=False), flush=True)
sco = pd.read_csv("pseudobulk_pca_scores.csv"); sco = sco[sco["set"] == "immune cells"].set_index("index")
m = bio.set_index("sample_id").join(sco[["PC1", "PC2", "PC3"]], how="inner"); print("biopsies with immune PCs:", len(m), flush=True)
rr = []
for pc in ["PC1", "PC2", "PC3"]:
    r_, p_ = spearmanr(m.ambient_median, m[pc]); rr.append(dict(PC=pc, spearman_with_ambient=r_, p=p_))
    d = m.assign(y=(m[pc] - m[pc].mean()) / m[pc].std(), amb=(m.ambient_median - m.ambient_median.mean()) / m.ambient_median.std())
    f0 = smf.ols("y ~ C(status, Treatment('non_inflamed')) + C(site)", data=d).fit(); f1 = smf.ols("y ~ C(status, Treatment('non_inflamed')) + C(site) + amb", data=d).fit()
    rr[-1].update(status_coef_site_adj=f0.params.filter(like="status").iloc[0], status_p_site_adj=f0.pvalues.filter(like="status").iloc[0],
                  status_coef_plus_ambient=f1.params.filter(like="status").iloc[0], status_p_plus_ambient=f1.pvalues.filter(like="status").iloc[0], ambient_coef=f1.params["amb"], ambient_p=f1.pvalues["amb"])
pd.DataFrame(rr).round(6).to_csv(f"{OUT}/ambient_vs_immune_pcs.csv", index=False); print(pd.DataFrame(rr).round(5).to_string(index=False), flush=True)

# ---------------- Part 3: HVGs
excl = a.var_names.str.match(HVG_EXCLUDE_REGEX) | a.var_names.isin(SEX_GENES) | a.var_names.str.match(IG_ALL) | epi_dom
print("genes excluded from HVG candidates:", int(excl.sum()), "| Ig-locus genes present:", int(a.var_names.str.match(IG_ALL).sum()), flush=True)
tmp = a[:, ~excl].copy(); tmp.X = tmp.layers["counts"].copy()
sc.pp.highly_variable_genes(tmp, flavor="seurat_v3_paper", n_top_genes=N_HVG, batch_key="sample_id", layer="counts", span=0.3)
old = a.var["highly_variable_immune"].copy() if "highly_variable_immune" in a.var else pd.Series(False, index=a.var_names)
allc = a.var["highly_variable"].copy() if "highly_variable" in a.var else pd.Series(False, index=a.var_names)
a.var["epithelial_dominant"] = epi_dom; a.var["highly_variable_mrvi"] = False; a.var["hvg_votes"] = 0; a.var["hvg_rank"] = np.nan
a.var.loc[tmp.var_names[tmp.var["highly_variable"].values], "highly_variable_mrvi"] = True
a.var.loc[tmp.var_names, "hvg_votes"] = tmp.var["highly_variable_nbatches"].values; a.var.loc[tmp.var_names, "hvg_rank"] = tmp.var["highly_variable_rank"].values
sel = a.var["highly_variable_mrvi"].values; assert sel.sum() == N_HVG
cols = [c for c in ["gene_ids", "hvg_votes", "hvg_rank"] if c in a.var]; a.var.loc[sel].sort_values("hvg_rank")[cols].to_csv(f"{OUT}/hvg_immune_mt15.csv")
v = a.var.loc[sel, "hvg_votes"]
cand = a.var.loc[~excl, "hvg_votes"]
vd = pd.DataFrame(dict(min_votes=[1, 3, 5, 10, 17, 25, 33])); vd["n_candidate_genes_with_at_least"] = [int((cand >= k).sum()) for k in vd.min_votes]; vd["n_selected_with_at_least"] = [int((v >= k).sum()) for k in vd.min_votes]
vd.to_csv(f"{OUT}/hvg_vote_distribution.csv", index=False); print(vd.to_string(index=False), "| fewest votes among selected:", int(v.min()), flush=True)
ov = pd.DataFrame([
    dict(comparison="new immune HVGs vs previous immune HVGs (MT<=25, Ig constant genes allowed)", n_new=int(sel.sum()), n_other=int(old.sum()), n_shared=int((a.var["highly_variable_mrvi"] & old).sum())),
    dict(comparison="new immune HVGs vs all-cell HVGs (step 4)", n_new=int(sel.sum()), n_other=int(allc.sum()), n_shared=int((a.var["highly_variable_mrvi"] & allc).sum())),
    dict(comparison="previous immune HVGs that were epithelial-dominant", n_new=0, n_other=int(old.sum()), n_shared=int((old & epi_dom).sum())),
    dict(comparison="previous immune HVGs that were Ig-locus genes", n_new=0, n_other=int(old.sum()), n_shared=int((old & a.var_names.str.match(IG_ALL)).sum())),
    dict(comparison="new HVGs nominated by all 33 biopsies", n_new=int(sel.sum()), n_other=33, n_shared=int((v == 33).sum())),
    dict(comparison="new HVGs nominated by <=2 biopsies", n_new=int(sel.sum()), n_other=2, n_shared=int((v <= 2).sum())),
    dict(comparison="new HVGs nominated by <=5 biopsies", n_new=int(sel.sum()), n_other=5, n_shared=int((v <= 5).sum()))])
ov.to_csv(f"{OUT}/hvg_overlap.csv", index=False); print(ov.to_string(index=False), flush=True)

# ---------------- tables and objects
ct = pd.crosstab(o["sample_id"].astype(str), o["immune_lineage"])
meta_cols = [c for c in ["gsm", "patient_id", "site", "status", "age", "biologic", "cd_history", "platform_id"] if c in o]
bm = o.groupby("sample_id", observed=True)[meta_cols].first(); bm["n_immune_cells"] = o.groupby("sample_id", observed=True).size(); bm.index = bm.index.astype(str); bm["ambient_median"] = bio.set_index("sample_id").ambient_median
bm.to_csv(f"{OUT}/biopsy_metadata_mrvi.csv"); ct.join(bm[["status", "site", "patient_id", "biologic"]]).to_csv(f"{OUT}/immune_cells_per_biopsy_by_lineage_mt15.csv")
summ = dict(mt_cap=MT_CAP, n_cells_before=int(n0), n_cells_after=int(a.n_obs), n_genes=int(a.n_vars), n_hvg=int(sel.sum()), n_biopsies=int(o.sample_id.nunique()),
            n_epithelial_dominant_genes=int(epi_dom.sum()), min_cells_per_biopsy=int(bm.n_immune_cells.min()), median_cells_per_biopsy=float(bm.n_immune_cells.median()),
            lineage_counts={k: int(x) for k, x in o["immune_lineage"].value_counts().items()})
json.dump(summ, open(f"{OUT}/step6a_summary.json", "w"), indent=1); print(json.dumps(summ), flush=True)
a.write_h5ad(f"{PROC}/06_immune_mrvi_input.h5ad", compression="lzf")
h = a[:, sel].copy(); h.X = h.layers["counts"].copy(); h.write_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad", compression="lzf")
print("written:", sorted(f for f in os.listdir(PROC) if f.startswith("06_")), flush=True)
