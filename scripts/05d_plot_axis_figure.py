#!/usr/bin/env python
"""Step 5d (figure only) - redraw fig05d_inflammation_axis_by_group.png from inflammation_axis_scores.csv.

The same drawing code lives at the end of 05d_all_groups_reference.py. This stand-alone copy exists because the first full
run of 05d crashed at that last figure (an undefined random generator, since fixed) AFTER all tables had been written; the
figure is therefore redrawn from the saved scores table instead of recomputing the 85-biopsy embedding. Numbers are not recomputed here.
Usage: GSE_STATUS_ORDER=control,non_inflamed,adjacent,inflamed GSE_OUT=out python 05d_plot_axis_figure.py --scores out/inflammation_axis_scores.csv
"""
import argparse, numpy as np, pandas as pd
from common import *
style()
ap = argparse.ArgumentParser(); ap.add_argument("--scores", required=True); a = ap.parse_args()
ax_df = pd.read_csv(a.scores); G = list(STATUS_ORDER); assert len(G) == 4
sets_ax = [s_ for s_ in ["all cells", "immune cells", "T/NK cells", "plasma cells", "myeloid cells", "epithelial cells"] if s_ in set(ax_df["set"])]
fig, axes = plt.subplots(len(sets_ax), 2, figsize=(7.2, 1.55 * len(sets_ax) + 0.6), sharex=True, sharey="row", squeeze=False)   # same y-scale for ileum and colon within a row
rng = np.random.default_rng(0)
for r_, nm in enumerate(sets_ax):
    d_ = ax_df[ax_df["set"] == nm]
    for c_, site in enumerate(["TI", "AC"]):
        ax_ = axes[r_, c_]; dd = d_[d_.site == site]
        for k, g in enumerate(G):
            v = dd[dd.status == g].score.values; held = g in ("control", "adjacent")
            ax_.scatter(k + rng.uniform(-0.18, 0.18, len(v)), v, s=16, marker="o", facecolors="none" if held else STATUS_COL[g], edgecolors=STATUS_COL[g], linewidths=0.9 if held else 0.3, zorder=3)
            if len(v): ax_.hlines(np.median(v), k - 0.3, k + 0.3, color="k", lw=1.2, zorder=4)
        ax_.axhline(0, color="0.85", lw=0.5, zorder=0); ax_.margins(x=0.08)
        if c_ == 0: ax_.set_ylabel(f"{nm}\n(SD units)", fontsize=6.5)
        if r_ == 0: ax_.set_title({"TI": "Terminal ileum", "AC": "Ascending colon"}[site], loc="left")
WRAP = {"control": "non-IBD\ncontrol", "non_inflamed": "non-\ninflamed", "adjacent": "adjacent\nto inflamed", "inflamed": "inflamed"}
for ax_ in axes[-1]: ax_.set_xticks(range(4)); ax_.set_xticklabels([WRAP[g] for g in G], fontsize=5.5)
fig.suptitle("Position of each biopsy on the inflamed-vs-non-inflamed axis (open circles: groups that did not define the axis)", x=0.01, ha="left", fontsize=7.5, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965)); savefig(fig, "fig05d_inflammation_axis_by_group.png")
