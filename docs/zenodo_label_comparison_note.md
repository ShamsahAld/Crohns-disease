# Comparison of our immune states with the authors' Zenodo labels

**Data.** Zenodo record 10.5281/zenodo.21508422, `scRNA_cell_metadata.csv` (CC-BY-4.0, md5 verified): 213,275 cells, labels in `lineage`, `cluster_label1`, `cluster_label2`. The authors state that their Zenodo annotations differ slightly from their final annotated dataset (still under embargo), so this is a comparison with a related labelling, not with the 79-type annotation.

**Matching.** Cells were matched on patient + site + status + barcode. 28,226 of our 61,538 immune cells lie in biopsies present in Zenodo (18 of our 33 primary biopsies); 25,531 of them (90.5%) share a barcode. 188 matched cells are in our excluded clusters and 44 are called epithelial/stromal by the authors; 25,343 cells (44 of our 45 analysed states have matched cells) are used for the state-level comparison.

## Results
- **Lineage level:** 99.2% of matched cells fall in the same broad lineage (T/NK 100%, plasma 99.9%, mast 99.4%, B 98.2%, myeloid 96.6%; 44 of our myeloid cells are epithelial in the authors' labels, both in our epithelial-contaminated clusters).
- **State level (our 44 states vs their 21 immune labels):** ARI 0.446, NMI 0.728 (0.555 / 0.775 with high+medium confidence states only). 83% of cells sit in the majority authors' label of their own state; 40 of 41 states with >=20 matched cells have >=50% of their cells in one authors' label. Our states are finer, so the reverse purity is lower (55%).
- **Clear matches:** neutrophil 100%; NK 93%; ILC3 -> ILCs 96%; gd/innate-like IEL -> gdT 92%; CD8 tissue-resident memory -> CD103+ CD8 79%; CD8 GZMK+ effector-memory -> KLRG1+ CD8 79%; CD4 naive/central memory -> CD4 Tcm 93%; naive and memory B states -> Naive/Memory B 88-99%; germinal-centre and cycling B -> GC B 97%; resident and C1Q+ macrophages -> Macrophage 92-94%; classical monocyte and inflammatory macrophage/monocyte -> Monocyte 98-100%; all mast states -> Mast 98-100%; cDC1 -> DC 91%.
- **Partial matches:** Treg-like CD4 -> Treg 57% (CD4 Tcm 40%); the authors' Treg label goes 85% to our Treg-like state. cDC2 -> DC 57% (Monocyte 40%).
- **Naming differences to review:** our 'MAIT/Th17-like' state (977 matched cells) is 81% 'CD103+ CD4' in the authors' labels and only 66% of their 55 MAIT cells; our 'Inflammatory macrophage/monocyte' is 100% 'Monocyte' in theirs.
- **Where the labelings are orthogonal:** plasma cells. The authors split IgA vs IgG plasma; our plasma states split by quality/stress/secretion and are IgA-dominant. Their 492 IgG plasma cells are spread across our plasma states (largest: IgA plasma core 27%, plasma low-UMI 26%). Likewise the authors have one Mast label, so our four mast sub-states cannot be validated by it.
- **Not testable with these labels:** pDC (1 of 37 cells in covered biopsies is in the Zenodo set), mature DC (0 of 29), plasmablast B-lineage cluster (1 of 45). Other states are >=98% present (cDC1: 32%). The reason these cells are missing from the authors' set is unknown (their filtering or their clustering); the pDC increase is therefore not cross-checked by this comparison.

## Caveats
- Matched cells are 41% of our immune cells, from 18 of 33 biopsies (7 of 17 inflamed, 11 of 16 non-inflamed).
- The two labelings come from the same cohort and largely the same cells, so agreement is not independent replication.
- ARI/NMI depend on the number of clusters; ours are finer, so they understate agreement in subtype terms.
- No labels, states or results in our analysis were changed by this comparison.

Files: `zenodo_label_agreement_metrics.csv`, `zenodo_label_confusion_counts.csv`, `zenodo_label_best_match_our_states.csv`, `zenodo_label_best_match_authors_labels.csv`, `zenodo_lineage_confusion_counts.csv`, `fig11_label_comparison_zenodo.png`.
