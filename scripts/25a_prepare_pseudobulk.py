"""25a_prepare_pseudobulk.py - plan step 8c helper: split pseudobulk_counts.npz (24b) into per-unit count tables and metadata for 25_pseudobulk_de.R.
Units: the 5 lineages and every state with confidence high or medium. Biopsies enter a unit only with >= 20 cells of it; genes with < 10 summed UMI in the unit are dropped before edgeR's filterByExpr.
Outputs (folder pbde/): pb_units.csv, pb_<i>_counts.csv.gz (genes x biopsies), pb_<i>_meta.csv (sample_id, status, site, patient, biologic, age, n_cells)."""
import os, numpy as np, pandas as pd
z = np.load("pseudobulk_counts.npz"); genes = z["genes"].astype(str); samples = z["samples"].astype(str); units = pd.read_csv("pseudobulk_units.csv"); bmeta = pd.read_csv("biopsy_metadata_mrvi.csv").set_index("sample_id").loc[samples]
os.makedirs("pbde", exist_ok=True); sel = units[(units.kind == "lineage") | (units.confidence.isin(["high", "medium"]))].reset_index(drop=True)
for i, u in sel.iterrows():
    U_ = z["U::" + u.unit]; N_ = z["N::" + u.unit]; ok = N_ >= 20; sub = U_[ok]; keepg = sub.sum(0) >= 10
    pd.DataFrame(sub[:, keepg].T, index=genes[keepg], columns=samples[ok]).to_csv(f"pbde/pb_{i + 1:02d}_counts.csv.gz")
    m = bmeta.loc[samples[ok], ["status", "site", "patient_id", "biologic", "age"]].rename(columns={"patient_id": "patient"}).copy(); m.insert(0, "sample_id", samples[ok]); m["n_cells"] = N_[ok]; m.to_csv(f"pbde/pb_{i + 1:02d}_meta.csv", index=False)
sel[["unit"]].to_csv("pbde/pb_units.csv", index=False); print(len(sel), "units written")
