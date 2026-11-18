"""26_de_calls.py - plan step 8d: combine pseudobulk (edgeR / limma-voom) and MrVI differential expression into gene calls.

Inputs : pbde/pseudobulk_de_results.csv (25_pseudobulk_de.R), mrvi_de_unit_gene_stats.npz (24c), pseudobulk_counts.npz (24b), epithelial_dominant_genes.csv (step 1).
Call   : a gene is CALLED in a unit (lineage or state) when ALL of
  1. edgeR quasi-likelihood FDR < 0.05 (design ~ site + status; biopsies with >= 20 cells of the unit) and |log2FC| >= 0.25;
  2. the gene is one of the 3,000 genes modelled by MrVI, MrVI mean cell-level LFC (status_inflamed, adjusted for site and biologic) has the same sign, |mean LFC| >= 0.05 and >= 90 % of the unit's cells share that sign;
  3. not guarded: immunoglobulin genes, epithelial-dominant genes (>=8x CPM epithelial vs immune), mitochondrial / MTRNR2L, ribosomal, haemoglobin, sex-linked genes;
  4. not ambient-risk: expression CPM in the unit's lineage < 20 % of the highest CPM in another immune lineage (>= 5 CPM there).
ROBUST = called + logCPM >= 1 + limma-voom with duplicateCorrelation(patient) FDR < 0.05 and edgeR with biologic + age FDR < 0.05, both with the same sign.
IEG / stress genes (FOS, JUN, HSPA1A ...) and TCR genes are kept but labelled in 'guard_class' because dissociation stress can alter them.
Thresholds in 2. were set from the pass rates of pseudobulk-significant vs non-significant genes (MrVI LFCs are shrunk 2-4x relative to pseudobulk), not from MrVI p-values, which are not gene-level.
"""
import re, numpy as np, pandas as pd
de = pd.read_csv("pbde/pseudobulk_de_results.csv"); mz = np.load("mrvi_de_unit_gene_stats.npz"); pz = np.load("pseudobulk_counts.npz"); units = pd.read_csv("pseudobulk_units.csv"); epi = set(pd.read_csv("epithelial_dominant_genes.csv").iloc[:, 0].astype(str))
munits = list(mz["units"].astype(str)); mgenes = mz["genes"].astype(str); mcov = list(mz["covariates"].astype(str)); ks = mcov.index("status_inflamed")
IEG = set("FOS FOSB JUN JUNB JUND EGR1 EGR2 EGR3 ATF3 IER2 IER3 IER5 DUSP1 DUSP2 DUSP5 DUSP6 NR4A1 NR4A2 NR4A3 HSPA1A HSPA1B HSPA6 HSPB1 HSPH1 DNAJB1 ZFP36 ZFP36L1 ZFP36L2 KLF2 KLF4 KLF6 BTG1 BTG2 SOCS3 CEBPB CEBPD PPP1R15A GADD45B MAFF MCL1 CXCL1 CXCL2 CXCL3 CCL3 CCL4 CCL3L1 CCL4L2".split()); SEX = set("XIST TSIX UTY DDX3Y KDM5D EIF1AY RPS4Y1 ZFY USP9Y".split())
def guard(g):
    if g.startswith("MT-"): return "mito"
    if g.startswith("MTRNR2L"): return "mito-like"
    if re.match(r"^(RPS|RPL|MRPS|MRPL)\d", g) or g in ("RPLP0", "RPLP1", "RPLP2", "RPSA"): return "ribo"
    if re.match(r"^HB[ABDEGMQZ]", g): return "hemoglobin"
    if re.match(r"^IG[HKL]", g) or g == "JCHAIN": return "Ig"
    if g in SEX: return "sex"
    if g in epi: return "epithelial-dominant"
    if re.match(r"^TR[ABGD][VJC]", g): return "TCR"
    if g in IEG: return "IEG/stress"
    return ""
de["guard_class"] = de.gene.map(guard); de["guarded"] = de.guard_class.isin({"Ig", "epithelial-dominant", "mito", "mito-like", "ribo", "hemoglobin", "sex"})
gidx = {g: i for i, g in enumerate(mgenes)}; uidx = {u: i for i, u in enumerate(munits)}; ok = de.gene.isin(gidx).values; ui = de.unit.map(uidx).values; gi = de.gene.map(gidx).values
for nm, key in [("mrvi_mean_lfc", "mean_lfc"), ("mrvi_frac_cells_pos", "frac_pos"), ("mrvi_mean_pde", "mean_pde"), ("mrvi_sd_lfc", "sd_lfc")]:
    v = np.full(len(de), np.nan); v[ok] = mz[key][ui[ok].astype(int), ks, gi[ok].astype(int)]; de[nm] = v
de["in_mrvi_genes"] = ok; de["mrvi_sign_consistency"] = np.where(de.mrvi_mean_lfc > 0, de.mrvi_frac_cells_pos, 1 - de.mrvi_frac_cells_pos)
LN = ["T_NK", "plasma", "B", "myeloid", "mast"]; genes = pz["genes"].astype(str); cpm = {}
for L in LN: S = pz["U::lineage:" + L].sum(0).astype(float); cpm[L] = S / S.sum() * 1e6
cpm = pd.DataFrame(cpm, index=genes); lin = lambda u: u.split(":", 1)[1] if u.startswith("lineage:") else units.set_index("unit").loc[u, "lineage"]; amb = {}
for u in de.unit.unique():
    L = lin(u); mx = cpm[[x for x in LN if x != L]].max(1); amb[u] = ((cpm[L] < 0.2 * mx) & (mx >= 5)).to_dict()
de["ambient_risk"] = [bool(amb[u].get(g, False)) for u, g in zip(de.unit, de.gene)]
de["pb_sig"] = (de.FDR < 0.05) & (de.logFC.abs() >= 0.25); de["mrvi_support"] = de.in_mrvi_genes & (np.sign(de.mrvi_mean_lfc) == np.sign(de.logFC)) & (de.mrvi_mean_lfc.abs() >= 0.05) & (de.mrvi_sign_consistency >= 0.9)
de["called"] = de.pb_sig & de.mrvi_support & ~de.guarded & ~de.ambient_risk; de["low_expr"] = de.logCPM < 1
de["robust"] = de.called & ~de.low_expr & (de.FDR_patient < 0.05) & (de.FDR_adj < 0.05) & (np.sign(de.logFC_patient) == np.sign(de.logFC)) & (np.sign(de.logFC_adj) == np.sign(de.logFC)); de["direction"] = np.where(de.logFC > 0, "up in inflamed", "down in inflamed")
# broader list for pathway sensitivity: pseudobulk-only genes (not covered by the 3,000-gene MrVI model) that pass the guards and both sensitivity models
de["pb_robust_not_mrvi_covered"] = de.pb_sig & ~de.in_mrvi_genes & ~de.guarded & ~de.ambient_risk & ~de.low_expr & (de.FDR_patient < 0.05) & (de.FDR_adj < 0.05) & (np.sign(de.logFC_patient) == np.sign(de.logFC)) & (np.sign(de.logFC_adj) == np.sign(de.logFC))
de.to_csv("de_pseudobulk_mrvi_all.csv.gz", index=False); de[de.called].sort_values(["unit", "FDR"]).to_csv("de_called_genes.csv", index=False)
s = de.groupby("unit").agg(genes_tested=("gene", "size"), pb_FDR05=("FDR", lambda x: int((x < 0.05).sum())), pb_sig=("pb_sig", "sum"), called=("called", "sum"), robust=("robust", "sum"), pb_robust_not_mrvi_covered=("pb_robust_not_mrvi_covered", "sum")).reset_index()
c = de[de.called]; s["called_up"] = s.unit.map(c[c.logFC > 0].groupby("unit").size()).fillna(0).astype(int); s["called_down"] = s.unit.map(c[c.logFC < 0].groupby("unit").size()).fillna(0).astype(int)
s["pb_sig_not_in_mrvi_genes"] = s.unit.map(de[de.pb_sig & ~de.in_mrvi_genes].groupby("unit").size()).fillna(0).astype(int); s["pb_sig_guarded_or_ambient"] = s.unit.map(de[de.pb_sig & de.in_mrvi_genes & (de.guarded | de.ambient_risk)].groupby("unit").size()).fillna(0).astype(int)
info = pd.read_csv("pbde/pseudobulk_de_unit_info.csv")[["unit", "n_biopsies", "n_inflamed", "n_non"]]; s = s.merge(info, on="unit", how="left"); s.to_csv("de_unit_summary.csv", index=False); print(s.to_string(index=False)); print("called total", int(de.called.sum()), "robust total", int(de.robust.sum()))
