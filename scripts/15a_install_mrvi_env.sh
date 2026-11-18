#!/bin/bash
# 15a_install_mrvi_env.sh - create a NEW venv for MrVI (scvi-tools >= 1.3, JAX/Flax backend) under scratch. Does not touch ~/.conda or ~/crohns_project.
set -eo pipefail
ENV=$HOME/.claude-science-scratch/envs/mrvi
PIP="python3 -m pip --python $ENV/bin/python"
echo "python: $(python3 --version)"
rm -rf "$ENV"                                  # my own env from an earlier attempt, never the user's
python3 -m venv --without-pip "$ENV"
$PIP install --quiet "scvi-tools" scanpy anndata statsmodels scikit-misc leidenalg igraph
$PIP install --quiet -U "jax[cuda12]"          # CUDA JAX build (the MrVI model runs on JAX); falls back to CPU if no GPU is visible
$PIP list 2>/dev/null | grep -iE "^(scvi-tools|jax|jaxlib|jax-cuda12|flax|optax|numpyro|torch|scanpy|anndata|numpy|pandas|scipy|lightning|pytorch-lightning|statsmodels) " || true
$ENV/bin/python - << 'PY'
import scvi, jax, flax, numpy
print("scvi", scvi.__version__, "| jax", jax.__version__, "| flax", flax.__version__, "| numpy", numpy.__version__)
print("jax devices on this (CPU) node:", jax.devices())
from scvi.external import MRVI
print("MRVI import OK:", MRVI)
PY
echo INSTALL_DONE
