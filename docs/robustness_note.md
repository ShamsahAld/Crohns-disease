# Robustness checks (plan step 10) - complete, including MrVI retraining at 2,000 and 5,000 genes

**Scope.** Headline results from steps 7-9 and the biopsy-level PCA were re-tested under variants that could plausibly change them. Summary table: `robustness_summary.csv` (91 checks); figure: fig10. MrVI was retrained with 2,000 and 5,000 variable genes (same selection procedure; as a control, the procedure at 3,000 genes reproduces the existing gene list exactly). Training: 2,000 genes 190 epochs (best 175); 5,000 genes 182 epochs (best 167); base 3,000 genes 157 epochs (best 142). The first 5,000-gene attempt (132 epochs) was stopped before its densities were written (exit code 143; most likely the 3.5 h client-side run limit I set, which that run exceeded after an unexplained ~2 h pause in its log); it was rerun alone with a longer limit and a saved model, and the rerun is the one reported. Stochastic CPU training means the two attempts did not stop at the same epoch.

## Variants
- **Extreme biopsy HA98AC_inf removed** (about 0 % epithelium) and **HA27 removed** (the only patient with both states).
- **Gene count for the biopsy PCA:** top 1,000 / 2,000 (original) / 5,000 variable genes.
- **Covariates:** site + biologic therapy + age + ileocolonic disease location (25 of 33 biopsies are ileocolonic).
- **Patient structure:** patient-clustered standard errors; one biopsy per patient (300 random draws); limma-voom with duplicateCorrelation for DE.
- **Other:** stricter 10 % MT cap; denominator = all cells in the biopsy; MrVI density test as a second method, with and without HA98AC_inf.

## Results
**Biopsy-level PCA.** The inflamed vs non-inflamed separation survives every variant (best-PC p between 1e-6 and 1.1e-4 in all cells, and between 4e-6 and 1.1e-3 in immune cells). Details: in all-cell pseudobulk the axis stays on PC2 (eta2 0.38-0.41 with all biopsies, 0.56-0.60 without HA98AC_inf). In immune-cell pseudobulk the axis is PC2 (eta2 0.29-0.32) with all biopsies, but when HA98AC_inf is removed PC2 weakens (p 0.01-0.02) and the inflammation axis appears as PC3 with a stronger association (eta2 0.37-0.47, p < 1e-4). The axis is therefore real and not caused by one biopsy, but its PC number is not stable. Choosing the best of three PCs per variant slightly flatters p-values; all three PCs per variant are in `pca_robustness.csv`.

**Differential abundance (per-biopsy proportions).**
- All six headline states keep BH q < 0.05 when HA98AC_inf or HA27 is removed (q 0.004-0.019), with the stricter 10 % MT cap (q 0.004-0.039), and in 100 % of one-biopsy-per-patient draws.
- Patient-clustered SEs: q < 0.02 for mast activated, CD8 Trm, pDC, gamma-delta IEL, plasmablast; **mast core weakens (q 0.13)**.
- Adding biologic + age + disease location (6 covariates on 33 biopsies): CD8 Trm, pDC, gamma-delta IEL and plasmablast stay (q 0.011-0.027); **mast activated weakens (p 0.009, q 0.066) and mast core is lost (p 0.081)**. This model costs power; the direction does not change.
- **CD8 Trm falls only relative to other immune cells:** with all cells in the biopsy as the denominator it is lost (q 0.18).
- Plasma IgG share: p 0.005 without HA98AC_inf, 0.001 with covariates or patient clustering.
- Second method (MrVI density): mast activated, mast core and CD8 Trm stay at q 0.031 without HA98AC_inf; pDC (q 0.14) and gamma-delta IEL (q 0.10) are weaker in MrVI than in the proportions test; **plasmablast is not supported by MrVI** (q 0.84), so it rests on one method. Inflammatory macrophage/monocyte falls from p 0.023 to 0.055 without HA98AC_inf.
- **MrVI retrained with 2,000 and 5,000 genes (new):** net-fraction effects across the 45 states correlate with the 3,000-gene run (Spearman 0.92 at 2,000 genes, same sign for 84% of states; 0.89 at 5,000 genes, same sign for 91%). Headline states keep their size (net fraction, 3,000 / 2,000 / 5,000 genes: CD8 Trm -0.60 / -0.59 / -0.59; mast activated +0.72 / +0.72 / +0.67; mast core +0.75 / +0.64 / +0.66; gamma-delta IEL -0.63 / -0.68 / -0.72). Significance depends on the gene set: with 5,000 genes CD8 Trm (q 0.018), mast activated (0.028), mast core (0.028) and gamma-delta IEL (0.028) pass q < 0.05, as the first three do with 3,000 genes; with 2,000 genes the same states have p 0.002-0.007 but BH q 0.077, so none passes. pDC (p 0.11 at 2,000, 0.071 at 5,000) and plasmablast (net fraction changes sign; p 0.46 / 0.15) lose MrVI support in both retrainings.

**Differential expression (called genes).** Direction is preserved for 100 % of called genes in all variants. FDR < 0.05 is retained for 96 % (T/NK) and 90 % (plasma) of called genes with the full covariate model, 68 % when HA98AC_inf is dropped (it removes one of 17 inflamed biopsies, so power falls; 90-92 % remain at FDR < 0.1), and 56 % (T/NK) / 74 % (plasma) with patient-aware limma-voom. B-lineage calls are fragile (20 % under limma-voom; 3 robust genes in all).

## What this changes in the interpretation
- Most reliable: CD8 tissue-resident memory lower, activated mast cells higher, plasma IgG share higher, pDC higher, gamma-delta/innate-like IEL lower, T/NK proliferation and interferon/TNF programmes (pathway step), and the existence of an inflammation axis in biopsy-level expression.
- MrVI support for the abundance changes depends on the gene set: it passes q < 0.05 at 3,000 and 5,000 genes but only reaches q 0.077 at 2,000 genes, with direction and size unchanged. Treat MrVI as supportive, not confirmatory, and rely on the per-biopsy proportion tests. pDC and plasmablast have no MrVI support at any gene count.
- Use caution: mast core (not significant with covariates or patient clustering), plasmablast (single method), CD8 Trm as an absolute tissue change (only relative to immune cells), all B-lineage DE calls, all mast DE and pathway results.
- Thresholds in the summary figure (green = q or p < 0.05, DE >= 80 % retained) are conventions for display, not tests.
