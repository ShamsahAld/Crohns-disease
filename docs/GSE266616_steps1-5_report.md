# GSE266616 — Crohn's disease, inflamed vs non-inflamed: roadmap steps 1–5

Scope: download → metadata checks → QC → normalisation → PCA/UMAP, for the **33 Crohn's disease biopsies labelled inflamed (17) or non-inflamed (16)**, from 28 patients (195,633 deposited cells). Adjacent-to-inflamed and non-IBD control biopsies are not part of this comparison.
All numbers below are read from the result files of one clean end-to-end run of the final scripts on the HPC (the run reproduced an earlier run exactly: same cells kept, identical tables).

## 1. Data and verification

* GEO lists **85** samples for GSE266616 (the roadmap says 79); your earlier cell-count table sums to 497,181 cells over the 85 matrices, which is the roadmap's figure (I re-counted the 33 analysed matrices directly: 195,633 cells, identical to that table). By sample ID there are 17 control and 29 Crohn's patients (the roadmap says 32).
* Design of the 33 analysed biopsies: inflamed 11 ileum (TI) + 6 colon (AC); non-inflamed 8 TI + 8 AC. **Only patient HA27 contributes both states, and from different sites**; four patients (HA26, HA31, HA36, HA37) contribute two non-inflamed biopsies. 24 of 33 biopsies are from patients on a biologic. One sample (HA48TI_inf) was sequenced on a different platform (GPL18573).
  The comparison is therefore **between patients**, with site and treatment as covariates ([metadata_checks.txt]({{artifact:a4a27d23-2c48-4dfb-8761-f86333e9992b}})).
* **Your files on the HPC are legitimate.** A fresh download from NCBI GEO was compared with `~/crohns_project/data/extracted_geo`: all 85 per-sample files and the RAW tar have exactly GEO's byte sizes; for the 33 analysed samples the SHA-256 of the fresh download is identical to your copy (33/33); every archive reads end-to-end and every matrix's dimensions agree with its barcode and gene files; your earlier metadata table agrees with GEO on status, site, patient, age, platform and disease for all 85 samples ([data_verification_GSE266616.csv]({{artifact:c41f8b29-57c6-4845-a35c-1f77f86d2c6f}}), [download_audit_fresh_geo_download.csv]({{artifact:4ed8df21-0f19-4f54-8e2e-b3fff8f44193}})). The analysis uses the fresh copy; your derived `.h5ad` files were not used and nothing in `~/crohns_project` was modified.

## 2. Methods and the biological reasoning behind each choice

| Step | Script | What was done | Why |
|---|---|---|---|
| 1 | [01_download_geo.py]({{artifact:d3c63a97-39f0-4e80-8f0e-d192edb9b327}}) | Download the 33 per-sample 10x archives, check size, gzip/tar integrity, matrix dimensions, SHA-256 | Establish that the data are intact before any analysis |
| 2 | [02_load_and_check.py]({{artifact:aa05742e-7a14-4965-8f18-4d0cae99a299}}) | Load raw counts (36,601 genes, identical in all samples), design and label checks, cells per sample | Counts are integer UMIs, no duplicate barcodes, no zero-UMI cells; cells per sample 2,140–16,089 |
| 3 | [03_qc.py]({{artifact:1a460261-7c77-402a-a6f6-91e190f6ed98}}), [qc_figures.py]({{artifact:27d73827-225c-464e-8305-b4090734ef18}}) | Remove cells with < 200 genes or > 25% mitochondrial reads; Scrublet doublets per sample | Absorptive epithelium is mitochondria-rich, so 25% (not the usual 10–20%) is used; no adaptive lower bound on UMI/genes because resting lymphocytes are naturally RNA-poor and are the cells of interest; no upper UMI cut-off because plasma/epithelial cells are legitimately large |
| 4 | [04_normalize.py]({{artifact:28ed33e3-19f8-4768-b905-5030d0db1ae1}}) | Counts kept in `layers['counts']`; X = log1p(counts / library size × 10,000); 3,000 HVGs (Seurat-v3, selected per sample then combined) | Per-sample HVG selection keeps genes variable across many biopsies; MT, ribosomal, haemoglobin, Ig/TCR V(D)J and sex-linked genes are excluded from HVG selection because they carry technical, clone- or sex-specific signal. No regression of library size or MT fraction (it would remove real epithelial biology) |
| 5a | [05a_pca_umap_all_cells.py]({{artifact:88e93a7f-3150-4254-9a61-483d82e02abc}}) | PCA (50 computed, 30 used; scree elbow at PC10) → 30-neighbour graph → UMAP (min_dist 0.3) → Leiden (resolution 0.5); provisional lineage from marker scores per cluster | Whole-biopsy view |
| 5b | [05b_immune_pca_umap.py]({{artifact:a30ad5b8-8e74-4215-87e9-e22427e63d7e}}) | Same, restricted to immune cells (65,191 cells, 47.8% of cells), HVGs re-selected inside the immune compartment | The research question is about immune cells; in the whole biopsy the epithelial/immune split uses up the leading PCs |
| 5c | [05c_pseudobulk_pca.py]({{artifact:c969cc28-117b-4a73-a710-5b6ce5c0c2ce}}) | Sum counts per biopsy (all cells, immune, T/NK, myeloid, plasma) → log2 CPM → top 2,000 variable genes → PCA of the 33 samples | Cells from one biopsy are not independent replicates and the design is between patients, so the sample is the correct unit for asking whether inflammation groups biopsies |

Support files: [common.py]({{artifact:8498c67e-e5fa-4eea-916d-1fa1bb2723a2}}), [00_setup_env.sh]({{artifact:4e9875d3-aa69-49ff-8cd6-9e6ed61cfe78}}), [00_env_witness.py]({{artifact:bd5dfda8-c826-48fd-b1f9-af965a64b439}}), [run_pipeline.sh]({{artifact:cc6cb212-02e8-4693-8d8c-d75ee1a47b9b}}) (re-runs everything), [environment_pip_freeze.txt]({{artifact:2bb045c9-d7e0-442a-9b6f-69219a34287d}}) (package versions: scanpy 1.12.4, anndata 0.13.4, scrublet 0.2.3, etc.).

## 3. Quality control

195,633 cells → **136,284 kept (69.7%)**: 57,078 removed for > 25% mitochondrial reads and/or < 200 genes, plus 2,271 predicted doublets among the rest.

![Per-biopsy QC distributions]({{artifact:d2aa1d96-1531-4169-be15-7b7cebb0a480}})
![Why cells were removed]({{artifact:d84d9902-462e-4e1d-a740-411f588184c9}})

Points to know before trusting later analyses:

* **Mitochondrial filtering is the main source of loss** (55,569 cells above 25%); 14 of 33 biopsies keep < 70% of cells. The loss is larger in non-inflamed biopsies (median 66.7% kept vs 72.8% in inflamed; Mann–Whitney p = 0.069, a trend only; [qc_loss_by_status.csv]({{artifact:a089a57a-e885-43de-aaaa-e49cbde38df9}})). The cells most at risk are epithelial (median mitochondrial fraction of kept cells: epithelial 15.9% vs 2.1–8.1% in the other lineages), so the loss is probably concentrated there, but lineage was not assigned to the removed cells. **Proportions should therefore be compared within the immune compartment** and the 25% threshold should be tested at 20%/30% before any composition conclusion (16,800 kept cells, 12.3%, lie between 20 and 25%).
* **Scrublet's automatic threshold is unreliable here** (0.25–0.67; it called almost no doublets in large epithelial samples). The applied threshold is the automatic one capped at 0.25 (with the automatic threshold the median predicted rate would have been 0.56% per sample), giving 0.7–3.4% predicted doublets per sample, which is still modest (droplet multiplet rates rise by roughly 0.4–0.8% per 1,000 recovered cells depending on chemistry, i.e. several percent for the larger samples); a lineage-mixing screen in step 5a found no cluster with mixed lineage signatures ([fig03d_doublet_diagnostics.png]({{artifact:58ee443e-fcbd-451e-ae99-f1db9a145f42}})).
* Ambient RNA could not be corrected (only called cells are deposited). IGHM, IGLC1, IGHG1 are among the four top-ranked HVGs; plasma cells and ambient immunoglobulin RNA are the likely sources, which were not separated here.
* Per-sample tables: [qc_summary_per_sample.csv]({{artifact:f0e73472-c83a-48c9-9010-8858d12f4620}}), [qc_thresholds.json]({{artifact:32226f4f-dd73-458e-bf9a-77085d052a23}}).

## 4. Normalisation

![Normalisation and HVG recurrence]({{artifact:0808e344-3252-4773-838c-7acc3be89232}})

Library sizes span more than one order of magnitude and have two modes, so cells cannot be compared on raw counts; differences in cell size between cell types are a likely reason (median genes per cell: T/NK 1106, epithelial 2890), although I did not test this directly. After normalisation every cell sums to 10,000 before the log (checked). HVG list: [hvg_list.csv]({{artifact:7018e001-76fd-49f7-9e38-244e53206aa5}}).

## 5. PCA and UMAP — what should it show, and what it showed

**Biologically expected for a whole intestinal biopsy** (and how it was checked): (i) cell identity dominates the leading PCs; (ii) technical covariates do not define structure; (iii) tissue site mostly separates epithelial cells; (iv) samples mix within a lineage but keep some donor structure; (v) inflamed vs non-inflamed is a subtle, cell-type-specific shift, not a global axis — so it must be judged within cell types and at sample level.

### 5.1 Whole biopsy (136,284 cells)

![Lineage markers on the UMAP]({{artifact:5e79c1ac-9062-431b-8e34-f20fbb79c6ed}})

Lineage was assigned per Leiden cluster (28 clusters) from marker scores — provisional, not the final annotation (roadmap step 7): 61,507 epithelial (45.1%), T/NK 19.4%, plasma 14.5%, B 6.4%, myeloid 5.9%, stromal 4.9%, endothelial 2.0%, mast 1.7%, glial 0.2%. As a check, PTPRC is detected in 91% of T/NK, 85% of B and 72% of myeloid cells but 3% of epithelial cells, and EPCAM in 89% of epithelial vs ≤ 9% of any other lineage ([lineage_summary.csv]({{artifact:840178e6-688a-44c3-8c2d-c0d4c75d3b9c}}), [leiden_lineage_scores.csv]({{artifact:ed7413c8-e5d4-4904-bfa0-4cc3edb7595e}})).

![PCA overview]({{artifact:6cb12627-ac05-472a-ac6e-84aacb18e86a}})

* **Cell identity dominates:** lineage explains 0.84 and 0.81 of the variance of PC1 and PC2 (up to 0.84 across PC1–20). The best any other factor does on any of PC1–20: sample 0.20, tissue site 0.08, inflamed vs non-inflamed 0.07, biologic therapy 0.02, platform 0.00 ([pca_eta2_all_cells.csv]({{artifact:4db21aec-891b-4593-817a-a43f1f9d8df2}})).
* **No technical structure:** UMI count, genes and mitochondrial fraction differ between cell types as expected (epithelium has the highest mitochondrial fraction) but no cluster is defined by them; the doublet score is diffuse; and the cells of the single GPL18573 sample are spread across all clusters ([fig05a_all_umap_technical.png]({{artifact:83162574-ec1b-4feb-8381-44f2d0d54a6b}})).
* **Donor structure is real:** a cell's nearest neighbours are 2.5–5.5× more often from the same biopsy than if cells mixed freely, in every lineage (platform: 1.0×). This is the structure MrVI (step 6) is designed to model.
* **Site and inflammation, tested against a null that keeps donor structure** (status or site re-assigned at random across biopsies, 200 permutations; observed / null mean):

| Lineage | Cells | Same-biopsy neighbours (× expected) | Same-platform | Site (obs / null) | Inflamed vs non-inflamed (obs / null) |
|---|---|---|---|---|---|
| epithelial | 61,507 | 5.5x | 1.01x | 1.65 / 1.26 (p=0.005) | 1.36 / 1.27 (p=0.070) |
| T_NK | 26,378 | 4.1x | 1.00x | 1.21 / 1.12 (p=0.005) | 1.21 / 1.12 (p=0.005) |
| plasma | 19,764 | 3.3x | 1.00x | 1.13 / 1.10 (p=0.174) | 1.24 / 1.11 (p=0.005) |
| B | 8,668 | 2.5x | 1.00x | 1.18 / 1.14 (p=0.169) | 1.10 / 1.15 (p=0.896) |
| myeloid | 8,084 | 4.7x | 1.00x | 1.21 / 1.22 (p=0.667) | 1.16 / 1.24 (p=0.990) |
| stromal | 6,640 | 4.9x | 1.00x | 1.40 / 1.27 (p=0.035) | 1.31 / 1.28 (p=0.299) |
| endothelial | 2,702 | 3.1x | 1.01x | 1.15 / 1.11 (p=0.035) | 1.16 / 1.11 (p=0.015) |
| mast | 2,297 | 2.8x | 1.00x | 1.15 / 1.10 (p=0.035) | 1.11 / 1.10 (p=0.373) |

  Site separates **epithelial** cells beyond donor structure (ileal vs colonic epithelium), as expected; immune lineages largely do not. Inflamed vs non-inflamed shows structure beyond donors only in **T/NK and plasma cells** (and marginally endothelial); epithelial, stromal, B, myeloid and mast cells are within the null ([knn_mixing_by_lineage.csv]({{artifact:1756e139-8d6c-410d-9c1d-2e8e11312fdb}})).

![Inflamed vs non-inflamed, site and sample on the UMAP]({{artifact:6fe2eb87-65da-4731-9f0f-675dea007c38}})

### 5.2 Immune cells only (65,191 cells)

![Immune UMAP with lineage and subset markers]({{artifact:37541575-7e02-484c-b22e-a081fead2515}})

Re-clustering the immune compartment resolves T cells (CD4⁺, CD8⁺ and FOXP3⁺ regions, plus an NKG7-high NK-like region), B, plasma, myeloid and mast cells. Immune lineage explains up to 0.89 of PC variance; sample 0.17; site 0.04; inflamed vs non-inflamed 0.02; platform 0.00. Inside a single lineage (PC1–5), inflamed vs non-inflamed explains at most 0.03–0.12 of variance, versus up to 0.49 for sample identity (myeloid) ([fig05b_immune_pca_eta2.png]({{artifact:99efa68b-516a-42e4-933c-6a6918243e06}}), [immune_eta2_within_lineage.csv]({{artifact:88812915-f44a-45c4-bbe2-79f18cff8988}})). **Conclusion at cell level: inflammation is not a dominant axis in single-cell PCA/UMAP; it is subtle and confounded with donor.**

![Immune UMAP: inflammation, site, sample]({{artifact:0f2c6d08-dad1-4ee6-9e56-2641f5dec649}})

### 5.3 Sample level (pseudobulk, n = 33 biopsies)

![Pseudobulk PCA]({{artifact:40569fb5-8c00-41e4-9738-6cb29c4151d0}})

This is where inflammation becomes visible. In the all-cell pseudobulk, **PC1 (40% of variance) is tissue site** (eta² = 0.73) and **PC2 (28%) separates inflamed from non-inflamed biopsies within each site** (eta² = 0.41). Within immune cells alone, PC1 is again site (0.53) and PC2 inflammation (0.32). Within single lineages the strongest associations with inflammation are **plasma cells PC1** (eta² 0.52) and **T/NK cells PC2** (0.40); myeloid cells show a weaker one. Platform and biologic therapy explain little (≤ 0.13 of any of PC1–3).

Linear models of each PC score on status, adjusted for site (n = 33; 9 of the 15 PC × set tests have BH q < 0.05; PC signs are arbitrary, so only the size of the difference matters):

| Sample-level PC | Inflamed − non-inflamed (SD of PC score) | p (unadjusted) | p (site-adjusted) | q (BH, 15 tests) |
|---|---|---|---|---|
| all cells PC2 | -1.29 | 5.4e-05 | 9.0e-05 | 0.000 |
| immune cells PC1 | +0.51 | 1.5e-01 | 2.3e-01 | 0.284 |
| immune cells PC2 | -1.12 | 6.7e-04 | 5.1e-04 | 0.002 |
| T/NK cells PC1 | -0.93 | 6.3e-03 | 5.7e-05 | 0.000 |
| T/NK cells PC2 | +1.27 | 7.9e-05 | 9.5e-06 | 0.000 |
| myeloid cells PC1 | +1.02 | 2.4e-03 | 3.6e-03 | 0.009 |
| myeloid cells PC2 | -0.90 | 8.8e-03 | 4.4e-03 | 0.009 |
| plasma cells PC1 | -1.45 | 2.1e-06 | 1.3e-07 | 0.000 |

![Variance of sample-level PCs explained by each factor]({{artifact:e9d059d5-95b5-40e5-8d9a-436de658e5c2}})

Exploratory composition (median share of QC-passed cells per biopsy; descriptive only — differential abundance is roadmap step 8):

| Lineage | Non-inflamed, all cells | Inflamed, all cells | Non-inflamed, immune cells | Inflamed, immune cells |
|---|---|---|---|---|
| epithelial | 60.2% | 31.8% | – | – |
| T_NK | 21.2% | 21.4% | 51.7% | 35.9% |
| plasma | 9.2% | 16.8% | 29.8% | 32.9% |
| B | 1.3% | 3.3% | 5.8% | 5.2% |
| myeloid | 3.3% | 7.6% | 7.7% | 13.5% |
| stromal | 2.7% | 4.4% | – | – |
| endothelial | 1.5% | 2.1% | – | – |
| mast | 0.7% | 2.8% | 1.6% | 4.8% |
| glial | 0.1% | 0.1% | – | – |

Inflamed biopsies have a smaller epithelial share and larger myeloid, mast and plasma-cell shares. The epithelial difference may partly reflect the higher QC loss in non-inflamed biopsies and how much epithelium each biopsy captured. HA98AC_inf is an extreme sample (epithelial share 0.0%) and is the outlier at the bottom left of the all-cell pseudobulk PCA ([fig05c_composition_per_sample.png]({{artifact:525cd480-224c-4d10-b1cc-66b81256f64f}}), [composition_median_by_status.csv]({{artifact:6ba3524b-9e76-4f9f-80d4-f7af23dc3630}})).

## 6. Caveats

1. **Between-patient design, partial confounding with site**: 28 patients in 33 biopsies, one patient with both states from different sites; four patients with two non-inflamed biopsies. Site and (to a lesser degree) treatment are covariates, not controlled factors; p-values treat biopsies as independent and are exploratory (PCs are chosen post hoc).
2. QC removes more cells from non-inflamed biopsies (trend, p = 0.07) and from epithelium; check the 20% / 30% mitochondrial thresholds before composition analyses.
3. Doublet detection is under-sensitive; ambient RNA is uncorrected; the lineage labels are provisional marker-based calls.
4. HVG recurrence is bimodal: 334 of the 3,000 HVGs are variable in only one biopsy and 348 in all 33 ([hvg_list.csv]({{artifact:7018e001-76fd-49f7-9e38-244e53206aa5}})). Sensitivity of the embeddings to the HVG set was not tested.
5. Directional statements about cell types (e.g. what drives plasma-cell PC1) are hypotheses for steps 6–9, not findings.

## 7. Where things are

* Figures and tables are attached to this session (links above). The remaining figures: [fig02_cells_per_sample.png]({{artifact:3b6b0d15-5dea-42cc-9471-866650d3f467}}), [fig03b_qc_thresholds.png]({{artifact:3fe41fed-5885-4cc9-9650-6ab6565c318f}}), [fig05a_all_umap_technical.png]({{artifact:83162574-ec1b-4feb-8381-44f2d0d54a6b}}).
* Large AnnData files stay on the cluster in `~/.claude-science-scratch/gse266616/` (the local disk is full): `processed/01_raw.h5ad` (2.0 GB, raw counts), `02_qc.h5ad` (1.5 GB), `03_norm.h5ad` (2.9 GB, counts layer + log-normalised + HVGs), `04_embedded.h5ad` (2.9 GB, PCA/UMAP/lineage), `05_immune.h5ad` (0.85 GB); `data/` holds the verified raw matrices. `03_norm.h5ad` (raw counts in `layers['counts']`, sample ID and covariates in `obs`) is the input for MrVI.
* Sample tables: [sample_metadata_inf_vs_non.csv]({{artifact:4ad0842b-1a25-40d5-8eb7-8830f82f9c1a}}) (analysed), [sample_metadata_GSE266616_all85.csv]({{artifact:27527cd9-33e2-483e-8aff-687dcf0a640c}}) (all 85).
