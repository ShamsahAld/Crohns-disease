"""11_composition_tests.py  -  biopsy-level tests behind the slide-7 dot plots (EXPLORATORY).

For each panel on slide 7 (epithelial = % of all cells; T/NK, plasma, myeloid, mast, B = % of immune cells):
  * Mann-Whitney U, inflamed (n=17) vs non-inflamed (n=16) biopsies, two-sided
  * OLS of the percentage on status + site (TI/AC), p-value for the status term
  * Benjamini-Hochberg over the six panels, separately for each p-value family
The biopsy (not the cell) is the unit. Between-patient design; 5 patients contribute two biopsies; site and therapy are
not balanced, so treat as supporting evidence only. Immune total = T_NK + plasma + myeloid + mast + B.

Inputs : pseudobulk_composition.csv (fraction of QC-passed cells per lineage per biopsy), sample_metadata_inf_vs_non.csv
Output : composition_tests_by_biopsy.csv
Usage  : python 11_composition_tests.py --comp pseudobulk_composition.csv --meta sample_metadata_inf_vs_non.csv --out composition_tests_by_biopsy.csv
"""
import argparse, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
ap = argparse.ArgumentParser(); ap.add_argument("--comp", required=True); ap.add_argument("--meta", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
comp, meta = pd.read_csv(a.comp), pd.read_csv(a.meta)
meta["sample_id"] = meta["sample_folder"].str.replace(r"\.tar\.gz$", "", regex=True).str.replace(r"^GSM\d+_", "", regex=True).str.replace("-", "_")
d = comp.merge(meta[["sample_id", "status", "site", "patient_id"]], on="sample_id", how="left")
assert d.status.notna().all() and len(d) == 33
IMM = ["T_NK", "plasma", "myeloid", "mast", "B"]
d["immune_total"] = d[IMM].sum(axis=1)
rows = []
for lin in ["epithelial"] + IMM:
    d["y"] = 100 * (d[lin] if lin == "epithelial" else d[lin] / d["immune_total"])
    non, inf = d.loc[d.status == "non_inflamed", "y"], d.loc[d.status == "inflamed", "y"]
    fit = smf.ols("y ~ C(status, Treatment('non_inflamed')) + C(site)", data=d).fit()
    rows.append(dict(panel=lin, denominator="all cells" if lin == "epithelial" else "immune cells", n_non=len(non), n_inf=len(inf),
                     median_non=non.median(), median_inf=inf.median(), mwu_p=mannwhitneyu(non, inf, alternative="two-sided").pvalue,
                     site_adj_diff_pct_points=fit.params.filter(like="status").iloc[0], site_adj_p=fit.pvalues.filter(like="status").iloc[0]))
r = pd.DataFrame(rows)
r["mwu_q_BH"] = multipletests(r.mwu_p, method="fdr_bh")[1]; r["site_adj_q_BH"] = multipletests(r.site_adj_p, method="fdr_bh")[1]
r.round(5).to_csv(a.out, index=False); print(r.round(4).to_string())
