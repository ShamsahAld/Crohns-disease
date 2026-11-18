#!/bin/bash
# Re-run the whole analysis (roadmap steps 1-5) for the 33 Crohn's disease inflamed / non-inflamed biopsies of GSE266616.
#
#   one-time:   bash 00_setup_env.sh                       # builds the Python environment (CPU only)
#   then:       sbatch run_pipeline.sh                     # or:  bash run_pipeline.sh
#
# Needs in the working directory: all *.py scripts, sample_metadata_inf_vs_non.csv, filelist.csv (GEO /suppl/ listing).
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=03:00:00
set -eo pipefail
export GSE_BASE=${GSE_BASE:-$HOME/.claude-science-scratch/gse266616}      # data/ (raw + unpacked matrices) and processed/ (AnnData files) live here
export GSE_OUT=${GSE_OUT:-$PWD/out}                                       # figures and tables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} MPLBACKEND=Agg
PY=${PY:-$HOME/.claude-science-scratch/envs/scrna_py312/bin/python}
mkdir -p "$GSE_OUT"

# Step 1  download + integrity checks (add --audit-dir <folder> to compare byte-for-byte with an existing copy)
$PY 01_download_geo.py --outdir "$GSE_BASE/data"
# Steps 2-5
for step in 02_load_and_check.py 03_qc.py 04_normalize.py 05a_pca_umap_all_cells.py 05b_immune_pca_umap.py 05c_pseudobulk_pca.py; do
  echo "################ $step"; $PY -u $step
done
