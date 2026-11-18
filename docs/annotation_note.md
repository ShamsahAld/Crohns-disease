# Immune cell-state annotation (plan step 5)

**Input:** 61,538 immune cells (15% MT cap) from 33 inflamed/non-inflamed biopsies, MrVI latent space (u, 10 dimensions).
**Method:** per-lineage Leiden clustering on the MrVI u-space (T/NK 13, plasma 7, myeloid 18, B 9, mast 5 = 52 clusters), then manual annotation of each cluster from marker-score tables and the dot plot (fig06f). Evidence per cluster is in `cluster_to_state_labels.csv`; cell-level labels are in `cell_state_labels.csv.gz`. Two neutrophil clusters (myeloid_14, 208 cells; myeloid_15, 853 cells) are merged, which gives 51 named states.

**Confidence classes (cells):** high 21,339; medium 30,284; low 9,235; exclude 680.
Confidence reflects the marker evidence for each cluster, not its size or biopsy count (high clusters occur in 12-33 biopsies, medium and low in 24-33, exclude in 9-29).
- high = lineage and subtype supported by several specific markers.
- medium = identity clear but the subtype call rests on fewer markers, or on a graded rather than a specific marker.
- low = technical or low-information cluster (low UMI, MT-high, ribosomal-high, stress/IEG); kept in the data but not interpreted biologically.
- exclude (680 cells) = epithelial-contaminated (myeloid clusters 12 and 17, 242 cells), low-quality plasma (plasma_4, 344 cells), or ambiguous (mast_1, myeloid_5, myeloid_16; 94 cells). Removed from abundance/DE and reported separately.

## Caveats (state these when presenting)
1. **Plasma isotype does not split clusters.** IgA dominates every plasma cluster (69-83% of heavy-chain counts), so plasma states are separated by quality/stress/secretory programme, not by isotype. Isotype fractions will be analysed separately at cell level (DA step).
2. **'Treg-like CD4' (cluster T_NK_3, 3,211 cells) is only weakly FOXP3+** (small dot in fig06f; mean 0.31 in the evidence column). It is CTLA4+ ICOS+ TIGIT+ with BATF, MAF and IL2RA, so it is best read as activated/regulatory-like CD4 until FOXP3 is checked at pseudobulk level.
3. **Many fine clusters are technical** (ribosomal, MT-high, IEG/stress, low-UMI). They are kept as states labelled low/exclude and abundance will be reported with and without them.
4. **Residual epithelial contamination:** 242 cells sit in two clusters labelled epithelial-contaminated (AGR2/PIGR; Paneth genes), both in the exclude class. Smaller contamination within other clusters is not removed (see the ambient-gene analysis in step 1).
5. **Labels are marker-based and not validated against the authors' 79 cell types**, which are not public (data 'available upon request'; email sent). If the authors send labels we compare and report agreement (adjusted Rand index / confusion matrix) rather than replace ours.
6. **Cell-level inflamed share is not enrichment.** Inflamed biopsies contribute more immune cells (61.1% of all immune cells), so a state with 61% inflamed cells is *not* enriched. Enrichment is tested at biopsy level in the next step.

## Early signals to test (not results)
| State | cells | biopsies | % cells from inflamed biopsies |
|---|---|---|---|
| Inflammatory macrophage/monocyte | 1,135 | 16 | 98.3 |
| Neutrophil | 1,061 | 18 | 93.6 |
| pDC | 63 | 17 | 95.2 |
| Mature DC (LAMP3+) | 71 | 17 | 90.1 |
| CD8 tissue-resident memory | 5,365 | 33 | 42.6 |
| gd/innate-like IEL | 1,402 | 30 | 35.5 |

Baseline for comparison is 61.1%. Biopsy-level tests (site-adjusted, BH-corrected) follow in 'Differential abundance'.
