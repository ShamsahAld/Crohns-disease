#!/usr/bin/env python
"""Step 5e - robustness of the held-out inflammation-axis group tests (output of 05d).

05d tests each group on the axis with OLS (score ~ group + site) that treats biopsies as independent. They are not:
the 33 control biopsies come from 17 patients (AC + TI), several non-inflamed patients contribute two biopsies, and 17 of 19
adjacent biopsies come from a patient who also gave an inflamed biopsy at the same site. This script therefore
  (1) refits the same models with patient-clustered (cluster-robust) standard errors, BH over the three non-circular contrasts per set;
  (2) runs a paired comparison adjacent vs inflamed on the same patient AND the same site (Wilcoxon signed-rank).
Input : inflammation_axis_scores.csv (from 05d), sample metadata.   Output: axis_sensitivity_patient_clustered.csv, axis_adjacent_vs_inflamed_paired.csv
Usage : python 05e_axis_sensitivity.py --scores out/inflammation_axis_scores.csv --outdir out
"""
import argparse, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser(); ap.add_argument("--scores", required=True); ap.add_argument("--outdir", default="."); a = ap.parse_args()
sc = pd.read_csv(a.scores)
LEVELS = ["non_inflamed", "control", "adjacent", "inflamed"]
rows = []
for s, dd in sc.groupby("set"):
    grp = pd.factorize(dd.patient_id)[0]
    for ref, tests in [("non_inflamed", ["control", "adjacent", "inflamed"]), ("inflamed", ["adjacent"])]:
        d2 = dd.assign(status=pd.Categorical(dd.status, [ref] + [x for x in LEVELS if x != ref]))
        f_ind = smf.ols("score ~ C(status) + C(site)", d2).fit()
        f_cl = smf.ols("score ~ C(status) + C(site)", d2).fit(cov_type="cluster", cov_kwds={"groups": grp})
        for g in tests:
            k = f"C(status)[T.{g}]"
            rows.append(dict(set=s, contrast=f"{g} - {ref}", diff_SD=f_cl.params[k], p_independent=f_ind.pvalues[k], p_patient_clustered=f_cl.pvalues[k],
                             n_biopsies=len(dd), n_patients=dd.patient_id.nunique(), circular=(g == "inflamed" and ref == "non_inflamed")))
res = pd.DataFrame(rows); res["q_BH_noncircular"] = np.nan
m = ~res.circular
res.loc[m, "q_BH_noncircular"] = multipletests(res.loc[m, "p_patient_clustered"], method="fdr_bh")[1]
res.round(5).to_csv(f"{a.outdir}/axis_sensitivity_patient_clustered.csv", index=False)

pr = []
for s, dd in sc.groupby("set"):
    A = dd[dd.status == "adjacent"].groupby(["patient_id", "site"]).score.mean().rename("adjacent")
    I = dd[dd.status == "inflamed"].groupby(["patient_id", "site"]).score.mean().rename("inflamed")
    j = pd.concat([A, I], axis=1, join="inner")
    if len(j) >= 6:
        pr.append(dict(set=s, n_pairs=len(j), median_adjacent_minus_inflamed=(j.adjacent - j.inflamed).median(),
                       n_adjacent_lower=int((j.adjacent < j.inflamed).sum()), wilcoxon_p=wilcoxon(j.adjacent, j.inflamed).pvalue))
pd.DataFrame(pr).round(5).to_csv(f"{a.outdir}/axis_adjacent_vs_inflamed_paired.csv", index=False)
print(res.round(4).to_string()); print(pd.DataFrame(pr).round(4).to_string())
