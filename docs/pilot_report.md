# MrVI pilot report (plan steps 2-3)  -  GSE266616 immune-only input

**Environment.** New venv `~/.claude-science-scratch/envs/mrvi` (scvi-tools 1.5.2, Python 3.12). Your own conda env and `~/crohns_project` were not touched.
Two findings about the software that differ from the older MrVI documentation:
1. In scvi-tools 1.5.2 MrVI is a **PyTorch** model (`scvi.external.mrvi.MRVAE`), not JAX. The JAX/flax packages that were installed are unused.
2. The cluster's NVIDIA driver supports CUDA <= 12.8, but the default PyPI torch wheel is built for CUDA 13, which crashes Lightning ("driver too old"). The venv therefore uses the **CPU-only torch wheel**, so MrVI trains on CPU. JAX itself does see the GPU (4096x4096 matrix multiply 183x faster than CPU), but MrVI does not use it.
A GPU is not needed: see timings below.

**Pilot design.** 15,015 immune cells subsampled from the 61,538-cell immune-only input (cap 455 cells per biopsy so that small biopsies are not lost; all 33 biopsies present), 3,000 vote-first HVGs, raw counts, `sample_key` = biopsy, no `batch_key`. Pilot lineages: T/NK 6,538; plasma 4,745; myeloid 1,735; B 1,515; mast 482.

| Setting | Seconds/epoch (15,015 cells, 8 CPUs) | Validation ELBO at epoch 12 | Projected full-set time, 130 epochs | Projected, 300 epochs |
|---|---|---|---|---|
| batch 128 (library default) | 2.31 | 960 | 20 min | 47 min |
| batch 512 | 1.53 | 992 | 14 min | 31 min |

(Full-set projection = pilot time x 4.10, assuming time is linear in cells.)

**Convergence.** Neither run is converged after 12 epochs: validation ELBO still fell by 0.8% (batch 128) and 1.5% (batch 512) over the last three epochs. Batch 128 reaches a lower loss at every epoch, so it is the better choice despite being slower per epoch.

**Downstream timings at pilot scale** (15,015 cells, 3,000 genes, covariates status + site + biologic; model trained 12 epochs, so the numbers show cost, not results): latent space 0.3 s; local sample distances 2.1 s; differential abundance 2.7 s; differential expression 24 s (store_lfc off, mc_samples 20). Projected for the full set: DE about 1.7 min, DA about 0.2 min.

**Decision for the full run (plan step 4).** Train on all 61,538 immune cells with `batch_size=128`, `max_epochs=300`, early stopping on validation loss, fixed seed, on a CPU node (partition ws-ia, 8 CPUs, 48 GB); expected wall time 20-50 min. The covariates for the later tests are `status` (reference non_inflamed), `site` (reference TI) and `biologic` (reference none); each categorical covariate is expanded to level columns (e.g. `status_inflamed`, `site_AC`, `biologic_anti_TNF`).

**Caveats.** The pilot used 12 epochs, so no biological conclusion is drawn from it. MrVI treats each biopsy as independent: the five patients with two biopsies are not modelled as a random effect (the pseudobulk cross-check in plan step 7 handles that). The timings are for 8 CPUs on a shared node and may vary.
