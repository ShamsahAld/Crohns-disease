"""15c_mrvi_witness.py - kernel-dispatch witness for the MrVI environment (toy data, ~minutes).
Exercises: setup_anndata(sample_key), train (device/time per epoch), latent u and z, local sample distances, differential_expression and
differential_abundance with several sample-level covariates (status 2 levels, site 2 levels, biologic 3 levels). Prints MRVI WITNESS OK at the end."""
import os, time, subprocess, threading, traceback, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, jax, scvi
from scvi.external import MRVI
print("scvi", scvi.__version__, "| jax", jax.__version__, "| devices:", jax.devices(), "| backend:", jax.default_backend(), flush=True)
print(subprocess.run("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader", shell=True, capture_output=True, text=True).stdout.strip(), flush=True)
rng = np.random.default_rng(0); n_s, n_cells, n_genes = 14, 7000, 400
samp = rng.integers(0, n_s, n_cells); status = np.array(["non_inflamed", "inflamed"])[(np.arange(n_s) % 2)[samp]]
site = np.array(["TI", "AC"])[((np.arange(n_s) // 2) % 2)[samp]]; bio = np.array(["none", "anti_TNF", "anti_IL12_23"])[(np.arange(n_s) % 3)[samp]]
ct = rng.integers(0, 2, n_cells); base = np.exp(rng.normal(0, 1, (2, n_genes))); lib = rng.lognormal(0, 0.3, n_cells)
eff = np.ones((n_cells, n_genes)); eff[np.ix_((ct == 0) & (status == "inflamed"), np.arange(30))] = 3.0
counts = rng.poisson(base[ct] * lib[:, None] * eff * 3).astype(np.float32)
obs = pd.DataFrame({"sample_id": pd.Categorical([f"S{i:02d}" for i in samp]), "status": pd.Categorical(status, categories=["non_inflamed", "inflamed"]), "site": pd.Categorical(site, categories=["TI", "AC"]),
                    "biologic": pd.Categorical(bio, categories=["none", "anti_TNF", "anti_IL12_23"]), "ctype": ct.astype(str)}, index=[f"c{i}" for i in range(n_cells)])
a = ad.AnnData(sp.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])); a.layers["counts"] = a.X.copy()
ok = {}
def stage(name, fn):
    t = time.time()
    try: r = fn(); ok[name] = True; print(f"[OK] {name} ({time.time()-t:.1f}s)", flush=True); return r
    except Exception: ok[name] = False; print(f"[FAIL] {name}"); traceback.print_exc(); return None
MRVI.setup_anndata(a, layer="counts", sample_key="sample_id"); model = MRVI(a)
gpu_idx = os.environ.get("CUDA_VISIBLE_DEVICES", "0").split(",")[0]; util, stop = [], threading.Event()
def poll():
    while not stop.is_set():
        o = subprocess.run(f"nvidia-smi -i {gpu_idx} --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits", shell=True, capture_output=True, text=True).stdout.strip().split(",")
        try: util.append((float(o[0]), float(o[1])))
        except Exception: pass
        time.sleep(0.5)
th = threading.Thread(target=poll, daemon=True); th.start()
# Lightning's accelerator only governs the (unused) PyTorch side; the torch wheel here is built for CUDA 13 but the driver is 12.8, so keep Lightning on CPU - the JAX model still runs on the GPU
t0 = time.time(); stage("train", lambda: model.train(max_epochs=6, batch_size=128, accelerator="cpu")); el = time.time() - t0
stop.set(); th.join(timeout=2)
print(f"train 6 epochs: {el:.1f}s ({el/6:.1f}s/epoch incl. startup/compile) | GPU {gpu_idx} utilisation during training: max {max([u for u, _ in util], default=float('nan')):.0f}% mean {np.mean([u for u, _ in util]) if util else float('nan'):.0f}% | memory used max {max([m for _, m in util], default=float('nan')):.0f} MiB", flush=True)
try:
    leaf = jax.tree_util.tree_leaves(model.module.train_state.params)[0]; print("model parameter device:", leaf.devices(), flush=True)
except Exception as e: print("param device check skipped:", repr(e)[:120], flush=True)
u = stage("latent_u", lambda: model.get_latent_representation(give_z=False)); z = stage("latent_z", lambda: model.get_latent_representation(give_z=True))
print("u", None if u is None else u.shape, "| z", None if z is None else z.shape, flush=True)
d = stage("local_sample_distances", lambda: model.get_local_sample_distances(batch_size=256)); print(d if d is None else {k: tuple(v.shape) for k, v in d.data_vars.items()}, flush=True)
print("sample_info:", list(model.sample_info.columns), "| n samples:", len(model.sample_info), flush=True)
cov = ["status", "site", "biologic"]
de = stage("differential_expression", lambda: model.differential_expression(sample_cov_keys=cov, store_lfc=True, use_vmap=False, mc_samples=20)); 
if de is not None: print("DE vars:", {k: tuple(v.shape) for k, v in de.data_vars.items()}, "| dims:", dict(de.sizes), flush=True)
da = stage("differential_abundance", lambda: model.differential_abundance(sample_cov_keys=cov)); 
if da is not None: print("DA vars:", {k: tuple(v.shape) for k, v in da.data_vars.items()}, flush=True)
if de is not None:
    try:
        print("DE coords:", {k: list(map(str, v.values[:6])) for k, v in de.coords.items() if v.ndim == 1 and k not in ("cell_name",)}, flush=True)
    except Exception: traceback.print_exc()
print("stages:", ok, flush=True); print("MRVI WITNESS OK" if all(ok.values()) else "MRVI WITNESS INCOMPLETE", flush=True)
