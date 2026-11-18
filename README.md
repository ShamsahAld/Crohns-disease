# Immune cells in inflamed vs non-inflamed Crohn's disease biopsies (GSE266616)

Course project (CBIO7102). Single-cell RNA-seq of intestinal biopsies from the public dataset **GEO GSE266616**: we compare **inflamed vs non-inflamed Crohn's disease (CD)** biopsies, focusing on immune cell composition and expression within cell types. Adjacent-to-inflamed biopsies and non-IBD controls are used only as reference groups. All results are **exploratory** (33 biopsies from 28 patients; the biopsy is the unit of analysis).

## Data
- GSE266616: 85 biopsies from 46 donors (terminal ileum and ascending colon), 497,181 deposited barcodes. 
- Primary comparison: 33 CD biopsies (17 inflamed, 16 non-inflamed) from 28 patients. The immune-only analysis uses 61,538 immune cells (mitochondrial reads <= 15%).
- **No data are stored in this repository.** `scripts/01_download_geo.py` downloads and verifies the GEO files.

## Folder layout
```
scripts/          numbered analysis scripts (Python and R), see the table below
result_tables/    small result tables, grouped by analysis step (01 ... 09)
figures/          figures fig02 ... fig11
docs/             written notes per step, literature comparison, label comparison
environment/      package versions of the environments used (cluster venvs, local Python, R)
```

## Pipeline (run in this order)
| Step | Script(s) | What it does |
|---|---|---|
| 0 | `00_setup_env.sh`, `00_env_witness.py`, `run_pipeline.sh` | Build the scanpy environment; `run_pipeline.sh` re-runs steps 1-5 |
| 1-4 | `01_download_geo.py`, `02_load_and_check.py`, `03_qc.py`, `qc_figures.py`, `04_normalize.py` | Download GEO, load into one AnnData of raw counts, cell QC (genes, mitochondrial reads, Scrublet doublets), normalisation and variable genes |
| 5 | `05a_` ... `05e_*.py`, `11_composition_tests.py` | PCA/UMAP of all and immune cells, pseudobulk (biopsy-level) PCA, four-group reference axis, biopsy-level composition tests |
| 6 | `13_immune_mt_cutoff.py`, `14_build_immune_mrvi_input.py` | Choose the immune mitochondrial cutoff (15%), build the immune-only MrVI input with an ambient-epithelial-RNA safeguard |
| 6 | `15a/15b_install_*.sh`, `15c_mrvi_witness.py`, `16_mrvi_pilot.py`, `17_mrvi_train_full.py`, `17b_mrvi_mixing_diagnostics.py` | MrVI environment, pilot on 15,000 cells, full training (CPU, 8 cores) and diagnostics |
| 7 | `18_annotate_immune_states.py` | Per-lineage Leiden clusters on the MrVI embedding, marker evidence, 51 named immune states with confidence classes |
| 7 | `19_mrvi_da_persample.py`, `20_da_proportions.py`, `21_mrvi_da_permutation.py`, `22_export_cell_qc_isotypes.py`, `23_plasma_isotypes.py` | Differential abundance (per-biopsy proportions and MrVI densities, permutation tests) and plasma-cell isotype use |
| 8 | `24a/24b/24c_*.py`, `25_pseudobulk_de.R`, `25a_prepare_pseudobulk.py`, `26_de_calls.py` | Within-cell-type differential expression: pseudobulk edgeR / limma-voom combined with MrVI per-cell fold changes |
| 9 | `27_pathway_enrichment.R` | Hallmark and Reactome enrichment (fgsea, over-representation) |
| 10 | `25b_*.R`, `28_mrvi_hvg_sensitivity.py`, `28b_mrvi_hvg_rerun.py`, `29_pca_robustness.py`, `30_da_robustness.py` | Robustness checks (drop extreme biopsy, covariates, patient clustering, MT cap, MrVI with 2,000 / 5,000 genes) |
| 11 | `31_zenodo_label_comparison.py` | Compare our states with the authors' cell labels in a public Zenodo record (annotation check only) |

## Main results (exploratory)
- **Plasma cells shift toward IgG** in inflamed biopsies: IgG share of heavy-chain counts 12.8% vs 2.1% (q = 0.018).
- **CD8 tissue-resident memory T cells are lower** (6.0% vs 15.4% of immune cells, q = 0.009), but only relative to other immune cells. **Mast cells are higher** (activated and core states). **pDCs are higher and gamma-delta/innate-like T cells lower**, mainly in the per-biopsy proportion tests. Neutrophils and inflammatory macrophage/monocytes are higher but not significant after correction (a few biopsies carry them).
- **Expression:** 621 genes called (pseudobulk and MrVI agree), 389 of them also pass stricter checks; most are in plasma and T/NK cells. Pathways: interferon, TNF-alpha/NF-kB and inflammatory-response programmes in myeloid, B and cDC2 cells; cell-cycle programmes in T/NK cells. These are generic inflammatory signals.
- **Robustness:** 91 checks. Per-biopsy results are stable to dropping the most extreme biopsy and to patient-level resampling; mast cells weaken when therapy, age and disease location are added. MrVI-based support depends on the gene set (q < 0.05 with 3,000 and 5,000 genes, not with 2,000).
- **Annotation check:** 99.2% of matched cells fall in the same broad lineage as the authors' Zenodo labels; state-level adjusted Rand index 0.45, normalised mutual information 0.73. See `docs/zenodo_label_comparison_note.md`.

## Reproducing
- Data location: `common.py` and `run_pipeline.sh` read the environment variable `GSE_BASE` (default `~/.claude-science-scratch/gse266616`; raw data in `data/`, AnnData files in `processed/`). Most later scripts accept `GSE_PROC` and `GSE_OUT`. Six scripts have the default path written in the code (`17b_`, `24a_`, `24c_`, `28_`, `28b_`, `29_`): edit the `PROC = ...` line near the top, or link your data folder to that path.
- `run_pipeline.sh` also needs `filelist.csv` (the GEO /suppl/ listing, not included) and `sample_metadata_inf_vs_non.csv` (included in `result_tables/01_qc_normalisation_pca/`).
- MrVI steps were run on CPU (8 cores, 56 GB RAM). Versions are in `environment/`: scvi-tools 1.5.2, torch 2.15.0 (CPU), scanpy 1.12.5 (MrVI environment) / 1.12.4 (QC environment), R 4.5.3 with edgeR 4.10.7, limma 3.68.5, fgsea 1.36.2; the pathway step used msigdbr 7.5.1 in a separate R 4.4.3 environment.
- Random seeds are fixed (seed 0 / script-specific) but CPU training is not bit-for-bit reproducible: two runs of the same 5,000-gene MrVI model stopped at different epochs.

## What is and is not included
- **Included:** all scripts, small result tables, figures, notes.
- **Not included (large or regenerable):** raw and processed data, `.h5ad` files, trained model weights, MrVI per-sample density arrays (`.npz`), per-cell score tables, the merged DE table with MrVI columns, and the pseudobulk robustness table (re-created by the scripts).
- **Produced interactively (no script):** the combined table `da_state_results_combined.csv`, `robustness_summary.csv`, the IBDverse name comparison table, and figures `fig06a-f`, `fig07a-d`, `fig08a-c`, `fig09a-b`, `fig10`, `fig11`. Tables from the annotation step were renamed when copied (`step5_*`, `cluster_to_state_labels.csv`, `state_summary.csv`, `cell_state_labels.csv.gz`).
- Two log files had the cluster home path shortened to `~`. Slide-building scripts are not included.

## External resources and credits
- Dataset: GEO GSE266616. Source study: Zhou R, Xie B, Weber C, et al. bioRxiv 10.1101/2025.03.26.645552 (preprint; later titled "Mucosal ecosystems of the ileum and colon shape the therapeutic response in Crohn's Disease").
- Cell labels from Zenodo record 10.5281/zenodo.21508422 (Zhao Y and Pott S, CC-BY-4.0) were used **only to check our annotation**. The authors state that these labels differ slightly from their final annotation, which is not yet public.
- The IBDverse annotation table (EBI BioStudies S-BSST2922, CC0) was used only for a name-level comparison of cell-type names (`docs/ibdverse_vocabulary_check_note.md`).
- Methods: MrVI (Boyeau et al., Nat Methods 2025;22(11):2264-2274) in scvi-tools; scanpy; edgeR; limma-voom; fgsea; MSigDB Hallmark and Reactome via msigdbr.

## Limitations
33 biopsies from 28 patients and mostly treated patients (24 of 33 on biologics); immune states are named from marker genes by us and are not validated against the authors' full annotation; shares are compositional; MrVI gene-level estimates have no calibrated error and were used only as an agreement filter. No treatment-response analysis was done.
