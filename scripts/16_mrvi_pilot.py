"""16_mrvi_pilot.py - roadmap step 6c: pilot MrVI on a ~15,000-cell subsample of the immune-only input; time everything; extrapolate to the full set.

Input : $GSE_PROC/06_immune_mrvi_input_hvg.h5ad (61,538 immune cells x 3,000 vote-first HVGs, X = raw counts)
Stages: (0) GPU-vs-CPU matmul benchmark + where the model arrays live; (1) stratified subsample (<= cap cells per biopsy, total ~15k);
        (2) MRVI.setup_anndata(sample_key = biopsy, no batch_key) and train at batch 128 and 512 (JAX on the GPU; Lightning kept on CPU because the torch wheel is CPU-only);
        (3) timings of get_local_sample_distances, differential_abundance and differential_expression (status, site, biologic) at the real gene count.
Outputs: out/pilot_summary.json, out/pilot_history_b128.csv, out/pilot_history_b512.csv, out/pilot_log.txt (stdout).  Nothing outside scratch is written."""
import os, sys, time, json, subprocess, threading, traceback, numpy as np, pandas as pd, anndata as ad, jax, jax.numpy as jnp, scvi
from scvi.external import MRVI
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
PROC = os.environ.get("GSE_PROC", os.path.expanduser("~/.claude-science-scratch/gse266616/processed")); OUT = os.environ.get("GSE_OUT", "out"); os.makedirs(OUT, exist_ok=True)
TARGET = int(os.environ.get("PILOT_CELLS", 15000)); EPOCHS = int(os.environ.get("PILOT_EPOCHS", 12)); SEED = 0
S = {"versions": dict(scvi=scvi.__version__, jax=jax.__version__), "devices": [str(d) for d in jax.devices()]}
print("scvi", scvi.__version__, "| jax", jax.__version__, "| devices", jax.devices(), flush=True)
gpu_idx = os.environ.get("CUDA_VISIBLE_DEVICES", "0").split(",")[0]
def gpu_query():
    o = subprocess.run(f"nvidia-smi -i {gpu_idx} --query-gpu=name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits", shell=True, capture_output=True, text=True).stdout.strip().split(",")
    return [x.strip() for x in o]
S["gpu"] = gpu_query(); print("GPU:", S["gpu"], flush=True)

# ---- stage 0: benchmark
def bench(dev, n=4096, reps=5):
    x = jax.device_put(jnp.ones((n, n), dtype=jnp.float32), dev); f = jax.jit(lambda m: m @ m); y = f(x).block_until_ready(); t = time.time()
    for _ in range(reps): f(x).block_until_ready()
    return (time.time() - t) / reps, str(y.devices())
try:
    tg, dg = bench(jax.devices("gpu")[0]); tc, dc = bench(jax.devices("cpu")[0], n=2048, reps=3); tc *= 8   # 2048^3 -> 4096^3 is 8x the flops
    S["benchmark"] = dict(gpu_s=tg, gpu_dev=dg, cpu_s_scaled=tc, cpu_dev=dc, speedup=tc / tg); print(f"matmul 4096^2: GPU {tg*1e3:.1f} ms on {dg} | CPU (scaled) {tc*1e3:.0f} ms -> {tc/tg:.0f}x", flush=True)
except Exception: traceback.print_exc(); S["benchmark"] = "failed"

# ---- stage 1: subsample
a = ad.read_h5ad(f"{PROC}/06_immune_mrvi_input_hvg.h5ad"); print("input:", a.shape, flush=True)
for c_ in ["sample_id", "status", "site", "biologic"]: a.obs[c_] = a.obs[c_].astype(str).astype("category")
a.obs["status"] = a.obs["status"].cat.reorder_categories(["non_inflamed", "inflamed"]); a.obs["site"] = a.obs["site"].cat.reorder_categories(["TI", "AC"])
bio_levels = [x for x in ["none", "anti_TNF", "anti_IL12_IL23"] if x in a.obs["biologic"].cat.categories] + [x for x in a.obs["biologic"].cat.categories if x not in ["none", "anti_TNF", "anti_IL12_IL23"]]
a.obs["biologic"] = a.obs["biologic"].cat.reorder_categories(bio_levels); print("biologic levels:", bio_levels, flush=True)
n_b = a.obs["sample_id"].value_counts()
lo, hi = 50, 5000
for _ in range(30):
    cap = (lo + hi) // 2; tot = int(np.minimum(n_b.values, cap).sum()); (lo, hi) = (cap, hi) if tot < TARGET else (lo, cap)
cap = hi; rng = np.random.default_rng(SEED)
idx = np.concatenate([rng.choice(np.where(a.obs["sample_id"].values == s)[0], size=min(int(n), cap), replace=False) for s, n in n_b.items()]); idx.sort()
p = a[idx].copy(); p.layers["counts"] = p.X.copy(); print(f"pilot subsample: {p.n_obs} cells, cap {cap}/biopsy, {p.obs.sample_id.nunique()} biopsies, lineages {p.obs.immune_lineage.value_counts().to_dict()}", flush=True)
S.update(pilot_cells=int(p.n_obs), cap_per_biopsy=int(cap), full_cells=int(a.n_obs), n_genes=int(a.n_vars))

# ---- stage 2: train at two batch sizes
def run_train(bs):
    q = p.copy(); MRVI.setup_anndata(q, layer="counts", sample_key="sample_id"); m = MRVI(q); util, stop = [], threading.Event()
    def poll():
        while not stop.is_set():
            try: g = gpu_query(); util.append((float(g[1]), float(g[2])))
            except Exception: pass
            time.sleep(1.0)
    th = threading.Thread(target=poll, daemon=True); th.start(); t = time.time()
    m.train(max_epochs=EPOCHS, batch_size=bs, accelerator="cpu", check_val_every_n_epoch=1); el = time.time() - t; stop.set(); th.join(timeout=3)
    try:
        h = m.history; hist = pd.concat([pd.Series(v.values.ravel() if hasattr(v, "values") else np.asarray(v).ravel(), name=k) for k, v in h.items()], axis=1)
    except Exception: traceback.print_exc(); hist = pd.DataFrame()
    hist.to_csv(f"{OUT}/pilot_history_b{bs}.csv", index=False)
    try: plat = pd.Series([d.platform for arr in jax.live_arrays() for d in arr.devices()]).value_counts().to_dict()
    except Exception: plat = {}
    return m, dict(batch_size=bs, epochs=EPOCHS, seconds=el, sec_per_epoch=el / EPOCHS, gpu_util_max=max([u for u, _ in util], default=np.nan), gpu_util_mean=float(np.mean([u for u, _ in util])) if util else np.nan,
                   gpu_mem_max_mib=max([m_ for _, m_ in util], default=np.nan), live_array_platforms=plat, history_cols=list(hist.columns))
res_train = {}
for bs in (128, 512):
    try: m_, info = run_train(bs); res_train[bs] = (m_, info); print("TRAIN", json.dumps(info, default=str), flush=True)
    except Exception: traceback.print_exc(); print(f"[FAIL] train batch {bs}", flush=True)
S["train"] = {str(k): v[1] for k, v in res_train.items()}
if res_train:
    bs_use = min(res_train, key=lambda k: res_train[k][1]["sec_per_epoch"]); model = res_train[bs_use][0]; S["batch_used_for_downstream"] = int(bs_use)
    # ---- stage 3: downstream timings
    cov = ["status", "site", "biologic"]
    def timed(name, fn):
        t = time.time()
        try: r = fn(); S[f"time_{name}_s"] = time.time() - t; print(f"[OK] {name}: {time.time()-t:.1f}s", flush=True); return r
        except Exception: traceback.print_exc(); S[f"time_{name}_s"] = None; print(f"[FAIL] {name}", flush=True)
    u = timed("latent_u", lambda: model.get_latent_representation(give_z=False)); print("u", None if u is None else u.shape, flush=True)
    timed("local_sample_distances", lambda: model.get_local_sample_distances(batch_size=256, keep_cell=False, groupby="immune_lineage"))
    timed("differential_abundance", lambda: model.differential_abundance(sample_cov_keys=cov))
    de = timed("differential_expression", lambda: model.differential_expression(sample_cov_keys=cov, store_lfc=False, use_vmap=False, mc_samples=20))
    if de is not None: print("DE vars:", {k: tuple(v.shape) for k, v in de.data_vars.items()}, flush=True)
    # extrapolation to the full immune set (time ~ linear in cells)
    f = S["full_cells"] / S["pilot_cells"]; spe = res_train[bs_use][1]["sec_per_epoch"]; S["extrapolation"] = dict(factor=f, sec_per_epoch_full=spe * f, hours_for_130_epochs=spe * f * 130 / 3600, hours_for_400_epochs=spe * f * 400 / 3600,
        default_epochs_full=int(min(round(20000 / S["full_cells"] * 400), 400)), de_minutes_full=(S.get("time_differential_expression_s") or np.nan) * f / 60, da_minutes_full=(S.get("time_differential_abundance_s") or np.nan) * f / 60)
json.dump(S, open(f"{OUT}/pilot_summary.json", "w"), indent=1, default=str); print(json.dumps(S.get("extrapolation", {}), default=str), flush=True); print("PILOT DONE", flush=True)
