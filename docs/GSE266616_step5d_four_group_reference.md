# GSE266616 - step 5d: four-group reference (non-IBD control / non-inflamed / adjacent-to-inflamed / inflamed)

Scope: all 85 biopsies, same QC and normalisation as the primary analysis (>=200 genes, <=25% mitochondrial reads, Scrublet per sample; 3,000 HVGs per sample).
This is a **reference** analysis. The primary comparison (steps 6-10) stays inflamed vs non-inflamed CD (33 biopsies).
All p-values are exploratory. Group is confounded with treatment (controls are all untreated), age and site mix, and the design is between-patient.

## Data check
- All 85 per-sample archives: SHA-256 identical to the copies in the user's HPC folder (85/85), byte size equal to the size GEO lists (85/85), archives read end to end, matrix dimensions consistent (85/85). Table: [download_audit_all85.csv]({{artifact:f3db87f8-3342-406b-982f-f15bcb7ffa89}}). 55 files were already in my scratch folder from earlier GEO downloads; 30 were fetched in this run.
- Cells: 497,181 raw -> 337,984 after QC (68.0%); 6,333 doublets removed. 45 of 85 biopsies keep <70% of their cells (flag HIGH_LOSS). Median percent kept by group: control 67.1, non-inflamed 66.7, adjacent 68.9, inflamed 72.8 (Kruskal-Wallis p = 0.11). Tables: [qc_summary_per_sample_all85.csv]({{artifact:e41c74a1-7bc4-4af1-94c7-3e02139361bc}}), [qc_thresholds_all85.json]({{artifact:a42d4e84-a578-4320-a4ff-03e1712ed512}}).

## 1. Cell level: lineage, not tissue state, organises the embedding
- PC1-20 eta-squared maxima: cell lineage 0.84, sample 0.17, site 0.09, tissue state 0.04, biologic 0.02, platform 0.00 ([pca_eta2_all_groups_cells.csv]({{artifact:a641bff8-a280-4690-8efc-e4626281ebbc}})).
- Every group occupies every lineage island on the UMAP; the four groups do not form separate clusters. Differences are in how many cells fall in each island.
- Lineage share of QC-passed cells (%), control / non-inflamed / adjacent / inflamed: epithelial 50.6 / 60.4 / 38.4 / 28.1; myeloid 2.6 / 2.7 / 5.3 / 8.9; stromal 3.2 / 3.0 / 3.9 / 7.0; plasma 14.4 / 11.1 / 13.9 / 18.3 ([lineage_share_by_group_cells.csv]({{artifact:abad7e08-4736-4ee2-af5a-2431bdeb2850}})).

![UMAP of all QC-passed cells]({{artifact:e7bb12e7-15fe-4dd3-9d9d-7c04c3f9eac8}})

![eta-squared heatmap]({{artifact:fd089d8c-739b-4836-8a51-176912521420}})

## 2. Sample level (pseudobulk PCA, one point per biopsy)
- All cells: PC1 (51%) is tissue site (eta2 0.75). PC2 (16%) is tissue state (eta2 0.48); biologic therapy also loads on PC2 (0.22) because it is collinear with group (all controls untreated).
- On PC2, relative to control (OLS adjusted for site): non-inflamed -0.04 SD (p = 0.85), adjacent +1.00 SD (p = 1.3e-5), inflamed +1.60 SD (p = 2.4e-10). Tables: [vs control]({{artifact:fccf322a-50ef-4134-82cb-87e989c2462a}}), [eta2]({{artifact:e0bcd8b7-c462-47d3-bfb8-3985e6e41dc0}}).
- Within lineages the same kind of axis separates states (plasma PC1 eta2 0.43; myeloid PC2 0.50; epithelial PC2 0.55; T/NK PC2 0.25).

![pseudobulk PCA]({{artifact:f9d8dfe3-a626-410f-8abb-096f85c9144f}})

## 3. Held-out inflammation axis
Method: a PC is fitted on inflamed and non-inflamed biopsies only (chosen by site-adjusted association with state). Control and adjacent biopsies are then projected onto it and z-scored to the non-inflamed training mean and SD. Inflamed vs non-inflamed on this axis is circular by construction. The evidence is where control and adjacent fall.
- Control vs non-inflamed (SD units): all cells +0.09, immune -0.00, myeloid -0.13, plasma +0.21, epithelial -0.25 (all p > 0.2); T/NK +0.46 (BH q = 0.014; patient-clustered q = 0.020).
- Adjacent vs non-inflamed: +0.75 all cells, +0.56 immune, +0.96 T/NK, +1.06 plasma, +0.76 epithelial (q < 0.03 each); myeloid +0.38 (q = 0.07).
- Adjacent vs inflamed: adjacent is lower than inflamed in all cells (-0.47, q = 0.026), immune (-0.62), myeloid (-0.73), epithelial (-0.86), plasma (-0.46); T/NK -0.19 (q = 0.34). Adjacent is intermediate and closest to inflamed in T/NK.
- Patient-clustered SEs (several patients contribute more than one biopsy) leave the pattern unchanged. BH q for adjacent vs inflamed: all cells 0.046, immune 0.0015, myeloid 0.0015, plasma 0.0091, epithelial < 0.0001, T/NK 0.10 ([axis_sensitivity_patient_clustered.csv]({{artifact:282101fc-e9aa-427e-86c0-b7559c28b95b}})).
- Paired check, same patient and same site (17 pairs; 16 for epithelial): adjacent is lower than inflamed in 13/16 epithelial, 12/17 immune, 13/17 myeloid, 13/17 plasma, 10/17 T/NK and 10/17 all-cell pairs (Wilcoxon p = 0.0006, 0.0038, 0.0056, 0.020, 0.089, 0.098) ([axis_adjacent_vs_inflamed_paired.csv]({{artifact:43bc1526-046d-4218-be24-e6a407adbd7c}})).
- Tables: [scores]({{artifact:304b5097-815a-4b53-a373-d8b119411727}}), [group tests]({{artifact:97512b8f-ca15-4996-adfd-11e606c67a8b}}), [medians]({{artifact:6ce9c85b-b81d-4a0b-881e-0232b1b4a69a}}).

![inflammation axis by group and site]({{artifact:dafe2e86-f96d-4bea-b709-8186246e47c8}})

## Reading the result
Tissue state separates at sample level but not at cell level. Control and non-inflamed biopsies sit together, inflamed biopsies are at the other end, and adjacent biopsies lie between them. This supports (i) using non-inflamed CD as the comparator for the primary analysis and (ii) keeping adjacent as its own group instead of pooling it with non-inflamed. It does not show that adjacent tissue is inflamed. All-cell axis scores partly reflect composition (epithelial fraction); the lineage-specific rows do not.

## Provenance
- Job c7eddce2 ran the download audit, QC, normalisation and 05d on all 85 biopsies. It crashed at the last figure of 05d (undefined random generator) after every table had been written.
- The axis figure was redrawn from inflammation_axis_scores.csv by [05d_plot_axis_figure.py]({{artifact:2c8ba373-8823-417a-945a-9a91f89939f3}}). The bug is fixed in [05d_all_groups_reference.py]({{artifact:8246d8da-e645-44de-9cb8-7f06cf87844b}}), together with a shared y-scale across the site panels. The fixed 05d was NOT re-run end to end.
- [05e_axis_sensitivity.py]({{artifact:1d4983e4-59b0-4e25-a712-5d763ea634b4}}) produced the clustered and paired tables from the saved scores (run locally).
