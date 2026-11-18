"""Step 2 - Load the per-sample 10x matrices, harmonise them into one AnnData of RAW counts, and check the design.

Checks written to out/: metadata_checks.txt, sample_summary_raw.csv, fig02_cells_per_sample.png
Output AnnData: $PROC/01_raw.h5ad  (X = raw UMI counts, float32; var_names = unique gene symbols; var['gene_ids'])
"""
import argparse, hashlib, sys, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
from common import *
style()
ap = argparse.ArgumentParser()
ap.add_argument("--meta", default="sample_metadata_inf_vs_non.csv")
ap.add_argument("--smoke", action="store_true", help="2 inflamed + 2 non-inflamed samples only (pipeline test)")
a = ap.parse_args()
meta = pd.read_csv(a.meta)
if a.smoke: meta = meta.groupby("status").head(2).reset_index(drop=True)
log = open(f"{OUT}/metadata_checks.txt", "w")
def P(*x):
    s = " ".join(str(i) for i in x); print(s, flush=True); log.write(s + "\n")

# ------------------------------------------------------------------ metadata / design checks
P("== DESIGN CHECKS ==")
P(f"samples: {len(meta)} | patients: {meta.patient_id.nunique()} | GSM unique: {meta.gsm.is_unique} | sample titles unique: {meta.title.is_unique}")
P("missing values per column:", meta[["patient_id", "site", "status", "age", "biologic", "platform_id"]].isna().sum().to_dict())
P("\nstatus x site:\n" + pd.crosstab(meta.status, meta.site).to_string())
P("\nstatus x biologic therapy:\n" + pd.crosstab(meta.status, meta.biologic).to_string())
P("\nstatus x platform:\n" + pd.crosstab(meta.status, meta.platform_id).to_string())
P("\nstatus x CD history:\n" + pd.crosstab(meta.status, meta.cd_history).to_string())
per_pat = meta.groupby("patient_id").agg(n_samples=("gsm", "size"), statuses=("status", lambda x: ",".join(sorted(set(x)))), sites=("site", lambda x: ",".join(sorted(set(x)))))
P(f"\npatients with >1 sample: {(per_pat.n_samples > 1).sum()} ; patients contributing BOTH statuses: {(per_pat.statuses.str.contains(',')).sum()}")
P(per_pat[per_pat.n_samples > 1].to_string())
P("\nage by status (median, range):", meta.groupby("status").age.agg(["median", "min", "max"]).to_dict("index"))

# ------------------------------------------------------------------ load matrices
adatas, rows, feat_hash = [], [], {}
for r in meta.itertuples():
    d = f"{DATA}/mtx/{r.sample_folder}"
    x = sc.read_10x_mtx(d, var_names="gene_ids", gex_only=True, make_unique=False)       # keep unique Ensembl IDs for safe concat
    x.var = x.var.rename(columns={"gene_symbols": "gene_symbol"})
    bc = pd.read_csv(f"{d}/barcodes.tsv.gz", header=None)[0]
    cnt = x.X.data
    feat_hash[r.title] = hashlib.md5("|".join(x.var_names).encode()).hexdigest()[:10]
    tot = np.asarray(x.X.sum(1)).ravel(); ng = np.asarray((x.X > 0).sum(1)).ravel()
    rows.append(dict(sample_id=r.title, gsm=r.gsm, status=r.status, site=r.site, patient_id=r.patient_id, platform=r.platform_id,
                     n_barcodes=x.n_obs, n_features=x.n_vars, dup_barcodes=int(bc.duplicated().sum()),
                     integer_counts=bool(np.all(cnt[:200000] == np.round(cnt[:200000]))), min_umi=int(tot.min()), median_umi=int(np.median(tot)),
                     median_genes=int(np.median(ng)), min_genes=int(ng.min()), cells_zero_umi=int((tot == 0).sum()), feature_hash=feat_hash[r.title]))
    x.obs_names = [f"{r.title}_{b}" for b in x.obs_names]
    for col, val in dict(gsm=r.gsm, sample_id=r.title, patient_id=r.patient_id, site=r.site, status=r.status, age=r.age,
                         biologic=r.biologic, cd_history=r.cd_history, platform_id=r.platform_id).items():
        x.obs[col] = val
    adatas.append(x)
    print(f"loaded {r.title}: {x.n_obs} cells x {x.n_vars} features", flush=True)

summ = pd.DataFrame(rows); summ.to_csv(f"{OUT}/sample_summary_raw.csv", index=False)
P("\n== MATRIX CHECKS ==")
P(f"total cells: {summ.n_barcodes.sum()} ; per-sample cells min/median/max: {summ.n_barcodes.min()}/{int(summ.n_barcodes.median())}/{summ.n_barcodes.max()}")
P("integer counts in all samples:", bool(summ.integer_counts.all()), "| duplicate barcodes:", int(summ.dup_barcodes.sum()), "| zero-UMI cells:", int(summ.cells_zero_umi.sum()))
P("distinct feature-list fingerprints:", summ.feature_hash.nunique(), "| n_features per sample:", sorted(summ.n_features.unique()))
P("min UMI per cell across samples (deposited-cell floor):", int(summ.min_umi.min()), "| min genes per cell:", int(summ.min_genes.min()))

# ------------------------------------------------------------------ concat on shared genes
common_genes = set(adatas[0].var_names)
for x in adatas[1:]: common_genes &= set(x.var_names)
union = set().union(*[set(x.var_names) for x in adatas])
P(f"genes: union {len(union)} | in all samples {len(common_genes)} | dropped by intersection {len(union) - len(common_genes)}")
keep = [g for g in adatas[0].var_names if g in common_genes]
adatas = [x[:, keep].copy() for x in adatas]
adata = ad.concat(adatas, join="inner", merge="same", index_unique=None)
adata.var = adatas[0].var.loc[adata.var_names, ["gene_symbol"]].copy()
adata.var["gene_ids"] = adata.var_names
sym = pd.Index(adata.var["gene_symbol"].astype(str).values)
dups = int(sym.duplicated().sum()); P("duplicated gene symbols made unique:", dups)
adata.var_names = ad.utils.make_index_unique(sym)
adata.var_names.name = None
for c_ in ["status", "site", "platform_id", "biologic", "sample_id", "patient_id"]: adata.obs[c_] = pd.Categorical(adata.obs[c_])
adata.obs["status"] = adata.obs["status"].cat.reorder_categories([g for g in STATUS_ORDER if g in set(adata.obs["status"])])
adata.X = sp.csr_matrix(adata.X).astype(np.float32)
P(f"\nFINAL raw AnnData: {adata.n_obs} cells x {adata.n_vars} genes ; X dtype {adata.X.dtype}")
mt = adata.var_names.str.startswith("MT-"); P("MT- genes found:", int(mt.sum()), "| MALAT1 present:", "MALAT1" in adata.var_names, "| PTPRC/EPCAM present:", "PTPRC" in adata.var_names, "EPCAM" in adata.var_names)
adata.write_h5ad(f"{PROC}/01_raw.h5ad", compression="lzf")
P("wrote", f"{PROC}/01_raw.h5ad")

# ------------------------------------------------------------------ figure: cells per sample
s = summ.copy(); s["status"] = pd.Categorical(s.status, STATUS_ORDER); s = s.sort_values(["status", "site", "sample_id"]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.2, 2.8))
ax.bar(s.index, s.n_barcodes / 1000, color=[STATUS_COL[x] for x in s.status], width=0.8)
ax.set_xticks(s.index); ax.set_xticklabels(s.sample_id, rotation=90, fontsize=5.5)
for t, site in zip(ax.get_xticklabels(), s.site): t.set_color(SITE_COL[site])
ax.set_ylabel("Cells per sample (thousands)"); ax.margins(x=0.01)
hs = [mpl.patches.Patch(color=STATUS_COL[k], label=DISPLAY[k]) for k in STATUS_ORDER] + [mpl.lines.Line2D([], [], color=SITE_COL[k], lw=2, label={"TI": "terminal ileum (label colour)", "AC": "ascending colon (label colour)"}[k]) for k in SITE_COL]
ax.legend(handles=hs, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, fontsize=6, borderaxespad=0.1)
what = "biopsies" if "control" in STATUS_ORDER else "Crohn's disease biopsies"
ax.set_title(f"{int(s.n_barcodes.sum()):,} deposited cells in {len(s)} {what}", loc="left", pad=16)
savefig(fig, "fig02_cells_per_sample.png"); log.close()
