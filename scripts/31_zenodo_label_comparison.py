"""31_zenodo_label_comparison.py - compare our immune cell states with the authors' cell labels in the Zenodo record 10.5281/zenodo.21508422.

Inputs : cell_state_labels.csv.gz (18_annotate_immune_states.py), qc_summary_per_sample_all85.csv (03_qc.py),
         scRNA_cell_metadata.csv from https://zenodo.org/records/21508422 (CC-BY-4.0; downloaded here, md5 checked; ~43 MB).
Method : cells are matched on patient + site + tissue status + barcode; our excluded clusters and cells the authors call non-immune are dropped
         for the state-level comparison. Outputs: confusion tables, best matches in both directions, ARI / NMI / homogeneity / completeness / purity.
Output : zenodo_sample_overlap.csv, zenodo_immune_label_inventory.csv, zenodo_lineage_confusion_counts.csv, zenodo_label_confusion_counts.csv,
         zenodo_label_best_match_our_states.csv, zenodo_label_best_match_authors_labels.csv, zenodo_label_agreement_metrics.csv
Use    : labels are used ONLY to check our annotation; the authors state their Zenodo annotations differ slightly from their final (embargoed) ones.
"""
import os, re, hashlib, requests, numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_score, completeness_score
OUT = os.environ.get("OUT", "."); REC = "21508422"
if not os.path.exists("scRNA_cell_metadata.csv"):
    j = requests.get(f"https://zenodo.org/api/records/{REC}", timeout=60).json(); f = [x for x in j["files"] if x["key"] == "scRNA_cell_metadata.csv"][0]
    b = requests.get(f["links"]["self"], timeout=300).content; assert hashlib.md5(b).hexdigest() == f["checksum"].split(":")[1]; open("scRNA_cell_metadata.csv", "wb").write(b)
md = pd.read_csv("scRNA_cell_metadata.csv", low_memory=False)
norm = lambda s: re.sub(r"[-_]", "", re.match(r"^((?:HA|OR)\d+)[-_]?(AC|TI)", s).group(0))
zmap = {"CD|inf": "inflamed", "CD|nonInf": "non_inflamed", "CD|adjInf": "adjacent", "Control|NA": "control"}
md["key"] = md.Sample_ID_scRNA.map(norm) + "|" + (md.disease_status.astype(str) + "|" + md.inflammation_status.fillna("NA").astype(str)).map(zmap)
md = md.rename(columns={"barcode": "bc"})[["key", "bc", "lineage", "cluster_label1", "cluster_label2"]]; assert not md.duplicated(["key", "bc"]).any()
q85 = pd.read_csv("qc_summary_per_sample_all85.csv"); q85["key"] = q85.sample_id.map(norm) + "|" + q85.status
zk, ok85 = set(md.key), set(q85.key); cnt_z = md.groupby("key").size().rename("n_zenodo")
# 1. sample overlap
ov_tab = q85[["sample_id", "key", "status", "site", "n_kept"]].copy(); ov_tab["in_primary_33"] = ov_tab.status.isin(["inflamed", "non_inflamed"]); ov_tab["in_zenodo_scRNA"] = ov_tab.key.isin(zk)
ov_tab = ov_tab[["sample_id", "key", "status", "site", "n_kept", "in_primary_33", "in_zenodo_scRNA"]].merge(cnt_z.rename("zenodo_cells"), left_on="key", right_index=True, how="left")
extra = pd.DataFrame([dict(sample_id="(not in GEO-85 list)", key=k, status=k.split("|")[1], site=k.split("|")[0][-2:], n_kept=np.nan, in_primary_33=False, in_zenodo_scRNA=True, zenodo_cells=int(cnt_z[k])) for k in sorted(zk - ok85)])
ov_tab = pd.concat([ov_tab, extra], ignore_index=True).rename(columns={"sample_id": "geo_sample_id", "n_kept": "our_cells_after_QC", "key": "patient_site_status"}); ov_tab.to_csv(f"{OUT}/zenodo_sample_overlap.csv", index=False)
# 2. match cells
cs = pd.read_csv("cell_state_labels.csv.gz"); cs["key"] = cs.sample_id.map(norm) + "|" + cs.status; cs["bc"] = cs.barcode.str.split("_").str[-1]
ov = cs.merge(md, on=["key", "bc"], how="left", indicator=True); mt = ov[ov._merge == "both"].drop(columns="_merge").copy()
print("our immune cells:", len(cs), "| in covered biopsies:", int(cs.key.isin(zk).sum()), "| matched by barcode:", len(mt))
def grp(r):
    if r.lineage == "tcell": return "T_NK"
    if r.lineage == "bcell": return "plasma" if "plasma" in r.cluster_label2 else "B"
    if r.lineage == "myeloid": return "mast" if r.cluster_label2 == "Mast" else "myeloid"
    return "non-immune"
mt["their_group"] = mt.apply(grp, axis=1); lin_order = ["T_NK", "B", "plasma", "myeloid", "mast"]
lc = pd.crosstab(mt.immune_lineage, mt.their_group).reindex(index=lin_order, columns=lin_order + ["non-immune"]).fillna(0).astype(int); lc.to_csv(f"{OUT}/zenodo_lineage_confusion_counts.csv")
lab = md[md.lineage.isin(["tcell", "bcell", "myeloid"])].groupby(["lineage", "cluster_label2"]).size().rename("n_cells_zenodo").reset_index()
m_ = mt.groupby("cluster_label2").size().rename("n_cells_also_in_our_immune_set"); lab = lab.merge(m_, left_on="cluster_label2", right_index=True, how="left").fillna({"n_cells_also_in_our_immune_set": 0}).astype({"n_cells_also_in_our_immune_set": int})
lab.sort_values(["lineage", "n_cells_zenodo"], ascending=[True, False]).to_csv(f"{OUT}/zenodo_immune_label_inventory.csv", index=False)
# 3. state-level comparison
m = mt[(mt.confidence != "exclude") & (mt.their_group != "non-immune")].copy(); ct = pd.crosstab(m.state, m.cluster_label2); ct.to_csv(f"{OUT}/zenodo_label_confusion_counts.csv")
conf = m.groupby("state").confidence.first(); lin_of = m.groupby("state").immune_lineage.first(); rows = []
for s_ in ct.index:
    v = ct.loc[s_].sort_values(ascending=False); n = v.sum()
    rows.append(dict(our_state=s_, lineage=lin_of[s_], confidence=conf[s_], n_matched_cells=int(n), top_authors_label=v.index[0], top_share=round(v.iloc[0] / n, 3), second_label=(v.index[1] if v.iloc[1] > 0 else ""), second_share=round(v.iloc[1] / n, 3)))
pd.DataFrame(rows).sort_values(["lineage", "n_matched_cells"], ascending=[True, False]).to_csv(f"{OUT}/zenodo_label_best_match_our_states.csv", index=False)
rows = []
for l_ in ct.columns:
    v = ct[l_].sort_values(ascending=False); n = v.sum()
    rows.append(dict(authors_label=l_, n_matched_cells=int(n), top_our_state=v.index[0], top_share=round(v.iloc[0] / n, 3), second_state=v.index[1], second_share=round(v.iloc[1] / n, 3)))
pd.DataFrame(rows).sort_values("n_matched_cells", ascending=False).to_csv(f"{OUT}/zenodo_label_best_match_authors_labels.csv", index=False)
rows = []
def met(d, name):
    rows.append(dict(subset=name, n_cells=len(d), our_states=d.state.nunique(), their_labels=d.cluster_label2.nunique(), ARI=adjusted_rand_score(d.cluster_label2, d.state), NMI=normalized_mutual_info_score(d.cluster_label2, d.state),
        homogeneity_their_labels_within_our_states=homogeneity_score(d.cluster_label2, d.state), completeness_our_states_within_their_labels=completeness_score(d.cluster_label2, d.state),
        purity_of_our_states=d.groupby("state").cluster_label2.agg(lambda s: s.value_counts().iloc[0]).sum() / len(d), purity_of_their_labels=d.groupby("cluster_label2").state.agg(lambda s: s.value_counts().iloc[0]).sum() / len(d)))
met(m, "all immune (excluded clusters removed)"); met(m[m.confidence.isin(["high", "medium"])], "high+medium confidence states only")
for L in lin_order: met(m[m.immune_lineage == L], f"our lineage {L}")
pd.DataFrame(rows).round(3).to_csv(f"{OUT}/zenodo_label_agreement_metrics.csv", index=False)
print("lineage-level agreement:", round((mt.immune_lineage == mt.their_group).mean() * 100, 1), "% | state-level ARI/NMI:", round(rows[0]["ARI"], 3), round(rows[0]["NMI"], 3))
