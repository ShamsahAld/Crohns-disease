"""18_annotate_immune_states.py - roadmap step 7 (plan step 5): immune cell-state clusters + marker evidence on the MrVI latent space.

Input : $GSE_PROC/06_immune_mrvi_input.h5ad (61,538 cells x 30,675 genes, X = log1p(CP10k), layers['counts']) and mrvi_immune_cells_latent.csv.gz (MrVI u + UMAP from step 4)
Method: per lineage (T/NK, plasma, myeloid, B, mast) build a 20-NN graph on the sample-unaware MrVI latent u and run Leiden (resolution per lineage below);
        score curated marker sets per cell (scanpy score_genes), average per cluster, z-score across clusters of the same lineage, and propose a label (argmax, with margin);
        Wilcoxon top markers per cluster; cross-lineage contamination scores; plasma isotype fractions from raw heavy-chain counts.
Outputs (CSV): cluster_summary, cluster_marker_scores, cluster_top_markers, marker_dot_data, plasma_isotype_by_cluster, cell_state_clusters, cluster_by_biopsy_counts.   Labels are PROPOSALS; final names are set after review."""
import os, sys, json, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
RES = {"T_NK": 0.7, "plasma": 0.4, "myeloid": 0.7, "B": 0.5, "mast": 0.3}
SETS = {
 "T_NK": {"CD4_T": ["CD4", "IL7R", "CD40LG", "LTB", "MAL", "TNFRSF25"], "naive_CM_T": ["CCR7", "SELL", "LEF1", "TCF7", "SATB1"], "Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF4", "TNFRSF18", "BATF"],
          "Tfh_like": ["CXCL13", "PDCD1", "CXCR5", "ICOS", "BCL6", "MAF", "TOX2"], "Th17": ["IL17A", "IL17F", "RORC", "CCR6", "IL23R", "CTSH", "KLRB1"], "CD8_T": ["CD8A", "CD8B", "GZMK", "CCL5", "CD3D"],
          "CD8_cytotoxic": ["GZMB", "PRF1", "GZMH", "GNLY", "NKG7", "CST7", "FGFBP2"], "Trm": ["CD69", "ITGAE", "CXCR6", "ITGA1", "ZNF683", "CD101"], "gdT": ["TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "KLRC1"],
          "NK": ["NCAM1", "KLRF1", "KLRC1", "KLRD1", "XCL1", "XCL2", "FCER1G", "TYROBP"], "ILC3": ["KIT", "IL1R1", "IL4I1", "NCR2", "LST1", "AREG", "TNFRSF25", "RORC"], "cycling": ["MKI67", "TOP2A", "STMN1", "TYMS", "HMGB2"]},
 "B": {"naive_B": ["IGHD", "TCL1A", "FCER2", "IL4R", "CD72", "SELL"], "memory_B": ["CD27", "TNFRSF13B", "AIM2", "COCH", "SSPN", "CLECL1", "CRIP1"], "GC_B": ["BCL6", "AICDA", "RGS13", "MEF2B", "LMO2", "SERPINA9", "GCSAM", "MYBL1", "S1PR2"],
       "activated_B": ["CD83", "EGR1", "NR4A1", "CD69", "FOS", "JUN"], "atypical_B": ["FCRL5", "ITGAX", "TBX21", "FCRL3", "CD86"], "plasmablast_like": ["MZB1", "XBP1", "PRDM1", "JCHAIN"], "cycling": ["MKI67", "TOP2A", "STMN1", "TYMS"]},
 "plasma": {"IgA": ["IGHA1", "IGHA2"], "IgG": ["IGHG1", "IGHG2", "IGHG3", "IGHG4"], "IgM": ["IGHM"], "cycling_plasmablast": ["MKI67", "TOP2A", "STMN1", "TYMS", "HMGB2"], "mature_plasma": ["TNFRSF17", "SDC1", "PRDM1", "XBP1"], "immature_plasmablast": ["CD27", "CD74", "HLA-DRA", "MS4A1", "CD19"]},
 "myeloid": {"resident_macrophage": ["C1QA", "C1QB", "C1QC", "FOLR2", "LYVE1", "MRC1", "CD163", "SELENOP", "F13A1", "STAB1"], "inflammatory_macrophage": ["IL1B", "TNF", "CXCL8", "IL6", "SOD2", "OSM", "IL23A", "CCL3", "CCL4", "NLRP3", "PTGS2", "G0S2"],
             "classical_monocyte": ["CD14", "VCAN", "FCN1", "S100A8", "S100A9", "S100A12", "LYZ"], "nonclassical_monocyte": ["FCGR3A", "LILRB2", "CDKN1C", "MTSS1"], "cDC1": ["CLEC9A", "XCR1", "BATF3", "CADM1", "IDO1"],
             "cDC2": ["CD1C", "FCER1A", "CLEC10A", "CD1E", "HLA-DQA1", "CLEC4A"], "mature_DC": ["LAMP3", "CCR7", "FSCN1", "CCL22", "CCL19", "BIRC3"], "pDC": ["LILRA4", "IL3RA", "TCF4", "IRF7", "CLEC4C"],
             "neutrophil": ["CSF3R", "FCGR3B", "CXCR2", "S100A8", "S100A9", "G0S2", "NAMPT"], "cycling": ["MKI67", "TOP2A", "STMN1", "TYMS"]},
 "mast": {"mast_core": ["TPSAB1", "TPSB2", "CPA3", "KIT", "HDC", "MS4A2"], "mast_activated": ["IL1RL1", "HPGDS", "LTC4S", "GATA2", "CD69", "FCER1A"], "mast_CMA1_CTSG": ["CMA1", "CTSG"]}}
a = sc.read_h5ad(f"{PROC}/06_immune_mrvi_input.h5ad"); lat = pd.read_csv("mrvi_immune_cells_latent.csv.gz").set_index("barcode")
assert a.n_obs == len(lat) and (a.obs_names == lat.index).all(); u = lat[[c for c in lat.columns if c.startswith("u") and c[1:].isdigit()]].values
a.obsm["X_mrvi_u"] = u; a.obs["immune_lineage"] = a.obs["immune_lineage"].astype(str); a.obs["sample_id"] = a.obs["sample_id"].astype(str); print("loaded", a.shape, flush=True)
Xc = a.layers["counts"].tocsr() if sp.issparse(a.layers["counts"]) else sp.csr_matrix(a.layers["counts"]); genes = np.asarray(a.var_names)
a.obs["cluster"] = ""; rows_sum, rows_score, rows_mark, dot, iso_rows = [], [], [], [], []
score_cols = [c for c in a.obs.columns if c.startswith("score_")]
for L, res in RES.items():
    sub = a[a.obs.immune_lineage == L].copy(); sc.pp.neighbors(sub, use_rep="X_mrvi_u", n_neighbors=20, random_state=0); sc.tl.leiden(sub, resolution=res, flavor="igraph", n_iterations=2, directed=False, random_state=0, key_added="cl")
    sub.obs["cl"] = L + "_" + sub.obs["cl"].astype(str); a.obs.loc[sub.obs_names, "cluster"] = sub.obs["cl"].values
    present = {k: [g for g in v if g in sub.var_names] for k, v in SETS[L].items()}; missing = {k: [g for g in v if g not in sub.var_names] for k, v in SETS[L].items()}
    print(L, "clusters:", sub.obs.cl.nunique(), "| missing markers:", {k: v for k, v in missing.items() if v}, flush=True)
    for k, gl in present.items():
        if len(gl) >= 1: sc.tl.score_genes(sub, gl, score_name=f"s_{k}", random_state=0)
    sc_cols = [f"s_{k}" for k in present if f"s_{k}" in sub.obs]
    cm = sub.obs.groupby("cl", observed=True)[sc_cols].mean(); z = (cm - cm.mean()) / cm.std(ddof=0).replace(0, 1); z.columns = [c[2:] for c in z.columns]
    cm.columns = [c[2:] for c in cm.columns]
    for cl in cm.index:
        zs = z.loc[cl].sort_values(ascending=False); rows_score.append(dict(cluster=cl, lineage=L, **{f"mean_{k}": cm.loc[cl, k] for k in cm.columns}, **{f"z_{k}": z.loc[cl, k] for k in z.columns}))
        m = sub.obs.cl == cl; o = sub.obs[m]; cross = a.obs.loc[o.index, score_cols].mean()
        rows_sum.append(dict(cluster=cl, lineage=L, n_cells=int(m.sum()), n_biopsies=int(o.sample_id.nunique()), frac_inflamed_cells=float((o.status == "inflamed").mean()), median_umi=float(o.total_counts.median()), median_pct_mt=float(o.pct_counts_mt.median()),
                             top1=zs.index[0], top1_z=float(zs.iloc[0]), top2=zs.index[1] if len(zs) > 1 else "", top2_z=float(zs.iloc[1]) if len(zs) > 1 else np.nan, **{f"xl_{c[6:]}": float(cross[c]) for c in score_cols}))
    # top markers (Wilcoxon) on genes detected in >= 50 cells of this lineage
    det = np.asarray((sub.X > 0).sum(0)).ravel() >= 50; s2 = sub[:, det].copy(); sc.tl.rank_genes_groups(s2, "cl", method="wilcoxon", n_genes=20, pts=True, use_raw=False)
    for cl in s2.obs.cl.unique():
        r = sc.get.rank_genes_groups_df(s2, group=cl); r["cluster"] = cl; r["lineage"] = L; rows_mark.append(r.head(20))
    # dot-plot data for the curated genes
    gl_all = sorted({g for v in present.values() for g in v}); sub_g = sub[:, gl_all]; Xd = sub_g.X.toarray() if sp.issparse(sub_g.X) else np.asarray(sub_g.X)
    for cl in cm.index:
        m = (sub.obs.cl == cl).values
        for j, g in enumerate(gl_all): dot.append(dict(cluster=cl, lineage=L, gene=g, mean_expr=float(Xd[m, j].mean()), frac_expr=float((Xd[m, j] > 0).mean())))
    if L == "plasma":
        ig = {"IgA": ["IGHA1", "IGHA2"], "IgG": ["IGHG1", "IGHG2", "IGHG3", "IGHG4"], "IgM": ["IGHM"], "IgD_E": ["IGHD", "IGHE"]}; idx = np.where(a.obs.immune_lineage.values == L)[0]
        tot = {k: np.asarray(Xc[idx][:, np.where(np.isin(genes, v))[0]].sum(1)).ravel() for k, v in ig.items()}; T = sum(tot.values()) + 1e-9
        d = pd.DataFrame({k: v / T for k, v in tot.items()}); d["cluster"] = sub.obs["cl"].values; d["heavy_chain_umi"] = T
        iso = d.groupby("cluster").agg(**{f"frac_{k}": (k, "mean") for k in ig}, median_heavy_umi=("heavy_chain_umi", "median"), n=("IgA", "size")).reset_index(); iso_rows.append(iso)
pd.DataFrame(rows_sum).round(4).to_csv(f"{OUT}/cluster_summary.csv", index=False); pd.DataFrame(rows_score).round(4).to_csv(f"{OUT}/cluster_marker_scores.csv", index=False)
pd.concat(rows_mark).round(4).to_csv(f"{OUT}/cluster_top_markers.csv", index=False); pd.DataFrame(dot).round(4).to_csv(f"{OUT}/marker_dot_data.csv", index=False)
if iso_rows: pd.concat(iso_rows).round(4).to_csv(f"{OUT}/plasma_isotype_by_cluster.csv", index=False)
a.obs[["sample_id", "immune_lineage", "cluster", "status", "site"]].to_csv(f"{OUT}/cell_state_clusters.csv")
pd.crosstab(a.obs["sample_id"], a.obs["cluster"]).to_csv(f"{OUT}/cluster_by_biopsy_counts.csv"); print(pd.DataFrame(rows_sum)[["cluster", "n_cells", "n_biopsies", "frac_inflamed_cells", "top1", "top1_z", "top2"]].round(2).to_string(index=False), flush=True)
print("ANNOTATION CLUSTERING DONE", flush=True)
