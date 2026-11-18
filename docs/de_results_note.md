# Within-cell-type differential expression, inflamed vs non-inflamed (plan step 8)

**Design.** Two independent methods per cell unit (5 lineages and 15 states with at least 5 biopsies in each status group; low-confidence states are not analysed): (1) pseudobulk, raw UMI summed per biopsy over the unit's cells, biopsies with at least 20 cells of the unit, edgeR quasi-likelihood with TMM and design ~ site + status (all expressed genes, log2 fold change inflamed / non-inflamed, BH within unit), plus two sensitivity models: limma-voom with duplicateCorrelation(patient) as the patient-aware version, and edgeR with biologic therapy and age added; (2) MrVI per-cell differential expression on the 3,000 model genes (status effect with site and biologic fitted together), averaged over the unit's cells.
MrVI's own cell-level p-values test a shift in the 30-dimensional latent space (95 % of cells are 'significant'), so they are not used. The gene-level MrVI quantities (mean log fold change, share of cells with the same sign) have no calibrated error, so they act as an agreement filter, and pseudobulk provides the inference.

**A gene is CALLED in a unit when** edgeR FDR < 0.05 and |log2FC| >= 0.25; MrVI has the same sign, |mean LFC| >= 0.05 and at least 90 % of the unit's cells agree; and the gene is not an immunoglobulin gene, epithelial-dominant gene (at least 8x CPM epithelial vs immune), mitochondrial/MTRNR2L, ribosomal, haemoglobin or sex gene, nor ambient-risk (another immune lineage expresses it at least 5x higher, so ambient RNA from that lineage could explain it). **ROBUST** = called, logCPM >= 1, and FDR < 0.05 with the same sign in both sensitivity models. IEG/stress genes (5 among called genes) and TCR genes (2) are kept but labelled.
MrVI LFCs are shrunk 2 to 4-fold relative to pseudobulk log2FC, so the MrVI size threshold was chosen from pass rates among pseudobulk-significant versus non-significant genes (not from MrVI p-values).

## Agreement between methods (fig08b)
Spearman correlation of MrVI mean LFC with pseudobulk log2FC on the 3,000 shared genes: plasma 0.81, T/NK 0.61, myeloid 0.54, B 0.48, mast 0.18. Among pseudobulk-significant genes the MrVI sign matches in 99 % (plasma), 97 % (T/NK), 95 % (myeloid), 93 % (B) and only 67 % (mast, 21 genes). Mast calls are therefore low reliability.

## Genes called (fig08c; table de_unit_summary.csv)
| Unit | biopsies infl / non | pseudobulk sig. | called (up / down) | robust | pseudobulk-only robust, outside MrVI genes |
|---|---|---|---|---|---|
| Plasma | 17 / 16 | 971 | 167 (113 / 54) | 115 | 218 |
| T/NK | 17 / 16 | 788 | 156 (124 / 32) | 86 | 301 |
| B | 14 / 13 | 92 | 15 | 3 | 15 |
| Myeloid | 17 / 16 | 41 | 10 | 2 | 9 |
| Mast | 14 / 9 | 45 | 6 | 0 | 3 |
States with most calls: IgA plasma core 97 (81 robust), IgA plasma stress-response (JUN/FOS) 96 (65), Treg-like CD4 17 (7), IL7R+ memory CD4 13 (10), MAIT/Th17-like 11 (6), IgA plasma secretory 16 (6). CD8 tissue-resident memory has 7 (3), NK 0, resident macrophage 0. Units with fewer than 5 biopsies in a group (mast activated, memory B TNFRSF13B+, C1Q+ macrophage low-UMI) were not tested. Neutrophils, inflammatory macrophage/monocytes, pDC and core mast cells are too sparse per biopsy for pseudobulk.

## What the calls show (exploratory)
- Most called genes are modestly expressed lineage-specific genes and non-coding transcripts; canonical inflammatory effector genes are mostly absent from the called set. Among the strongest: T/NK CXCL13 up (log2FC +4.4, FDR 0.003; adjusted-model FDR 0.006 but patient-aware FDR 0.20, so called and not robust), GBP5 up (+0.7, interferon-induced; robust), EZH2 and BCL6 up (robust), CD96 down (-1.1, robust); TNFSF9 down (called, not robust). IFNG is not significant (FDR 0.27).
- Plasma cells: the largest set. Up: CCND2, TLE1, DCC, NEAT1, PIK3R3; down: CECR2, NRG3, NPAS3, UCHL1, TTLL7. Some of these may reflect plasma-cell maturation state or anatomical origin rather than inflammation per se.
- 20 genes are called in both T/NK and plasma (for example CCNE2, HIST1H2BD/BN, HIST1H3J, MIR4435-2HG, SLFN11, GADD45A, SLC7A1, SPCS1), which suggests a shared proliferation/inflammatory-environment signal across lineages.
- B-lineage plasma-cell genes (MZB1, PRDM1, FKBP11) are up in the B lineage by pseudobulk but not called (they are flagged ambient-risk because plasma cells express them far more); this is consistent with the plasmablast-like cells in that lineage that are more abundant in inflamed biopsies (see differential abundance), but this explanation has not been tested separately.

## Caveats
- Calls are restricted to the 3,000 MrVI genes, which excludes Ig and epithelial-dominant genes by design. This is conservative: 67 % of the pseudobulk-significant T/NK genes lie outside the model, hence the extra pseudobulk-only 'robust' list for the pathway step.
- Lineage-level pseudobulk mixes composition changes with within-state changes (T/NK and plasma composition differ by status); the MrVI counterfactual estimate is within-cell-state, which is why the intersection is stricter than pseudobulk alone.
- All FDRs are within a unit, not across units; 33 biopsies from 28 patients; five patients have two biopsies (duplicateCorrelation handles this; in the mast subset there are no repeats, correlation set to 0).
- Dissociation stress genes can change with tissue handling; the 'stress-response (JUN/FOS)' plasma state is itself defined by them.
- Thresholds (|log2FC| >= 0.25; MrVI |LFC| >= 0.05; 90 % agreement; 5x ambient rule) were set once from the data distributions described above; the full table lets any other threshold be applied.
