#!/bin/bash
# Step 0 - build a private Python venv for the scRNA-seq pipeline (CPU only) and validate it.
# Lives under the app's scratch area so nothing in the user's own conda envs / project folders is touched.
# Host note: python3-venv/ensurepip is NOT installed on this cluster, so the venv is created with --without-pip
# and packages are installed with the system pip's --python option (pip >= 22.3).
set -euo pipefail
ENVDIR="$HOME/.claude-science-scratch/envs/scrna_py312"
mkdir -p "$(dirname "$ENVDIR")"
rm -rf "$ENVDIR"                                   # always rebuild cleanly (this dir is created only by this script)
python3 -m venv --without-pip "$ENVDIR"
PIP="python3 -m pip --python $ENVDIR/bin/python"
$PIP install --quiet \
    "scanpy>=1.10" anndata igraph leidenalg scikit-misc scrublet \
    h5py pyarrow matplotlib seaborn umap-learn pynndescent scikit-learn statsmodels
$PIP freeze > pip_freeze.txt
"$ENVDIR/bin/python" 00_env_witness.py
