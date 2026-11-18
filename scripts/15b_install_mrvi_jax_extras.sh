#!/bin/bash
# 15b - add the JAX modelling libraries MrVI needs to the venv created by 15a (scvi-tools 1.5 treats them as optional)
set -eo pipefail
ENV=$HOME/.claude-science-scratch/envs/mrvi
PIP="python3 -m pip --python $ENV/bin/python"
$PIP install --quiet flax optax numpyro
$PIP list 2>/dev/null | grep -iE "^(scvi-tools|jax|jaxlib|jax-cuda[0-9a-z-]*|flax|optax|numpyro|chex|nvidia-[a-z0-9-]*|torch|numpy) " || true
$ENV/bin/python - << 'PY'
import scvi, jax, flax, optax, numpyro, numpy
print("scvi", scvi.__version__, "| jax", jax.__version__, "| flax", flax.__version__, "| optax", optax.__version__, "| numpyro", numpyro.__version__)
print("jax backend on this node:", jax.default_backend(), jax.devices())
from scvi.external import MRVI
import inspect
print("MRVI import OK")
for n in ["setup_anndata", "train", "get_latent_representation", "get_local_sample_distances", "differential_expression", "differential_abundance"]:
    try: print(n, inspect.signature(getattr(MRVI, n)))
    except Exception as e: print(n, "ERR", e)
PY
echo INSTALL_B_DONE
