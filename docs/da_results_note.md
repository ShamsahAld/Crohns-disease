# Differential abundance, inflamed vs non-inflamed (plan step 7)

**Scope.** 33 primary biopsies (17 inflamed: 11 terminal ileum + 6 ascending colon; 16 non-inflamed: 8 + 8). 45 analysed immune states (51 states minus 6 'exclude' states, 680 cells). Smallest biopsy has 517 analysed immune cells.
The biopsy is the unit; all p-values are exploratory; BH q-values are computed across the 45 states. Adjacent-to-inflamed and non-IBD biopsies are not used here.

## Two independent views of abundance
1. **Per-biopsy proportions** (script 20): each state as a share of the biopsy's immune cells. Primary test = site-stratified Wilcoxon (van Elteren) with 10,000 within-site permutations of status. Effect = log2 ratio of site-standardised mean proportions.
2. **MrVI local density** (scripts 19, 21): for each cell, log ratio of its density under the inflamed vs the non-inflamed biopsies' aggregated posteriors (own biopsy omitted; reproduces MrVI's own enrichment output exactly). Primary state statistic = net fraction of cells at least 2-fold denser in one group (bounded, because raw log ratios are heavy-tailed). Null = 10,000 within-site permutations of biopsy labels, redone including own-biopsy omission; the null mean (not zero) is subtracted to give a site-adjusted score.

## Lineage level (share of immune cells, site-standardised means)
| Lineage | inflamed | non-inflamed | BH q (proportions) | BH q (MrVI) |
|---|---|---|---|---|
| T/NK | 34.2 % | 50.6 % | 0.0035 | 0.031 |
| Mast | 4.8 % | 2.0 % | 0.018 | 0.031 |
| Myeloid | 14.0 % | 8.2 % | 0.031 | 0.47 |
| B | 12.0 % | 9.2 % | 0.83 | 0.73 |
| Plasma | 35.0 % | 30.1 % | 0.45 | 0.73 |

## State level: 8 of 45 states have q < 0.05 in the per-biopsy test; 3 are also q < 0.05 in MrVI
| Tier | State (confidence) | inflamed vs non-inflamed (% of immune cells) | log2 ratio (TI / AC) | q proportions | q MrVI |
|---|---|---|---|---|---|
| A | Mast, activated (IEG-high) (medium) | 1.99 vs 0.56 | +1.81 (+1.60 / +2.05) | 0.007 | 0.041 |
| A | CD8 tissue-resident memory (high) | 5.99 vs 15.35 | -1.36 (-1.23 / -1.64) | 0.009 | 0.037 |
| A | Mast, core (KIT+ MS4A2+) (high) | 1.46 vs 0.50 | +1.53 (+1.14 / +1.83) | 0.028 | 0.037 |
| B | pDC (high) | 0.15 vs 0.01 | +2.91 | 0.004 | 0.28 |
| B | gamma-delta / innate-like IEL (high) | 1.59 vs 4.65 | -1.54 | 0.009 | 0.089 |
| B | Mast, low-information (low) | 0.57 vs 0.14 | +1.91 | 0.009 | 0.12 |
| C | Plasmablast (B-lineage cluster) (high) | 0.27 vs 0.05 | +2.30 | 0.009 | 0.84 |
| C | Plasma low-UMI (ribosomal-high) (low) | 2.96 vs 1.74 | +0.76 | 0.030 | 0.31 |

Tier A = q < 0.05 in both views; B = q < 0.05 in one and p < 0.05 in the other, same direction; C = q < 0.05 in one view only. Direction is the same in terminal ileum and ascending colon for every tier A-C state.
Tier D (nominal p < 0.05, q > 0.05): mature DC (LAMP3+), neutrophil, inflammatory macrophage/monocyte, classical monocyte, IL7R+ memory CD4, ILC3, MAIT/Th17-like.

## Robustness (13 analysis variants, fig07c)
Variants: OLS-logit permutation, quasi-binomial GLM, adding biologic therapy and age, 300 random one-biopsy-per-patient draws, leave-one-biopsy-out, high+medium states only, three denominators (all immune cells including excluded, all cells in the biopsy, own lineage), stricter 10 % MT cap, MrVI cell-weighted and biopsy-balanced.
Mast activated is p < 0.05 in 13/13 variants; gamma-delta IEL 13/13; CD8 Trm 12/13 (it loses significance when the denominator is all cells in the biopsy, i.e. its fall is relative to other immune cells and is not shown against total tissue cells); mast core 12/13; pDC 12/12 tested.
Within their own lineage, pDC (q 0.016), plasmablast (q 0.020), NK (q 0.030), Treg-like CD4 (q 0.038, up) and resting mast (q 0.016, down) stay significant; CD8 Trm and gamma-delta IEL are weaker within T/NK (q 0.091 and 0.051), so part of their fall is a lineage-composition effect (T/NK share falls while mast and myeloid shares rise). The two views agree across states (Spearman 0.77, same direction in 76 % of states).

## Plasma isotype use (fig07d)
IgG share of plasma-cell heavy-chain UMI is higher in inflamed biopsies (12.8 % vs 2.1 %; site-standardised +10.7 percentage points; p = 0.0025, q = 0.018 over 7 tests); cell-level calls point the same way (+15.2 pp, p = 0.018, q = 0.061). IgA falls correspondingly but is not significant (biopsy level q = 0.39), IgM is unchanged, and IgA2 within IgA is unchanged. The IgG signal is heterogeneous: 8 of 17 inflamed biopsies exceed the highest non-inflamed value (6.3 %), the other 9 do not. Ambient immunoglobulin RNA affects both measures.

## Caveats
- The signals for **inflammatory macrophage/monocyte and neutrophil** are concentrated: three biopsies hold 88 % of the macrophage/monocyte cells (HA27TI_inf, HA98AC_inf, HA82TI_inf) and 86 % of the neutrophils (HA27TI_inf, HA60TI_inf, HA98AC_inf); three non-inflamed biopsies (HA26TI_non, HA35TI_non, HA36TI_non) have at least 1 % neutrophils. Biopsy-level q = 0.15 and 0.11; treat as suggestive. HA98AC_inf is the epithelium-poor outlier and HA27 is the only patient with both states.
- **Mast, activated (IEG-high)** is defined by immediate-early genes, which can arise from tissue dissociation; the core mast state (KIT+ MS4A2+, high confidence) shows the same increase.
- **pDC** is rare (63 cells in 17 biopsies, up to 0.6 % of immune cells in a biopsy) and the MrVI view is not significant at q < 0.05; pDC, plasmablast and the two low-confidence states (mast low-information, plasma low-UMI) should not be interpreted biologically on their own.
- Proportions are compositional: a fall in one state can follow a rise in another; the own-lineage denominator addresses part of this.
- Choices made after seeing data (declared): the van Elteren test became the primary test (instead of OLS-logit) because many states are sparse; the bounded net-fraction became the primary MrVI statistic because the mean log ratio is heavy-tailed (the mean-score permutation p-values are in the table). Both alternatives are reported next to the primary result.
- 33 biopsies from 28 patients; 5 patients contribute two biopsies (HA27 contributes both states). One-biopsy-per-patient draws agree with the main result. The statistical power for small effects is limited.
- States are labelled by markers only; the authors' 79-type labels are not available yet (email sent).
