"""QC figures for step 3 (kept separate so they can be regenerated from the saved per-cell QC table without recomputing QC)."""
import numpy as np, pandas as pd
from common import *

def make_figures(o, summ, min_genes, max_pct_mt, diag=None):
    style()
    summ = summ.copy(); stat = summ.set_index("sample_id")["status"]; site = summ.set_index("sample_id")["site"]
    order = summ.assign(status=pd.Categorical(summ.status, STATUS_ORDER)).sort_values(["status", "site", "sample_id"])["sample_id"].tolist()
    status_legend = [mpl.patches.Patch(color=STATUS_COL[k], label=DISPLAY[k]) for k in STATUS_ORDER]
    def color_site_labels(ax):
        for t in ax.get_xticklabels(): t.set_color(SITE_COL[site[t.get_text()]])

    # ---- fig03a: per-sample distributions. Violins are drawn explicitly so each body's colour is its sample's status (no hue mapping)
    grp = {sid: d.sample(min(len(d), 2500), random_state=0) for sid, d in o.groupby("sample_id")}
    specs = [("total_counts", "UMI per cell", True, [1e3, 3e3, 1e4, 3e4, 1e5], ["1k", "3k", "10k", "30k", "100k"], None),
             ("n_genes_by_counts", "Genes per cell", True, [200, 500, 1e3, 3e3, 1e4], ["200", "500", "1k", "3k", "10k"], min_genes),
             ("pct_counts_mt", "Mitochondrial reads (%)", False, [0, 20, 40, 60, 80, 100], ["0", "20", "40", "60", "80", "100"], max_pct_mt)]
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 5.4), sharex=True)
    for ax, (col, yl, logy, tv, tl, line) in zip(axes, specs):
        data = [np.log10(np.clip(grp[s][col].values, 1, None)) if logy else grp[s][col].values for s in order]
        parts = ax.violinplot(data, positions=np.arange(len(order)), widths=0.85, showextrema=False)
        for body, sid in zip(parts["bodies"], order):
            body.set_facecolor(STATUS_COL[stat[sid]]); body.set_edgecolor("k"); body.set_linewidth(0.3); body.set_alpha(1.0)
        if line is not None: ax.axhline(np.log10(line) if logy else line, color="k", lw=0.6, ls="--")
        ax.set_yticks(np.log10(tv) if logy else tv); ax.set_yticklabels(tl); ax.set_ylabel(yl); ax.margins(x=0.01)
    axes[-1].set_xticks(range(len(order))); axes[-1].set_xticklabels(order, rotation=90, fontsize=5.5); color_site_labels(axes[-1])
    axes[0].legend(handles=status_legend + [mpl.lines.Line2D([], [], color="k", ls="--", lw=0.6, label="QC threshold")], ncol=3, loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=6, borderaxespad=0.1)
    axes[0].set_title("QC metrics per biopsy", loc="left")
    savefig(fig, "fig03a_qc_violins.png")

    # ---- fig03b: thresholds in context
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6)); d = o.sample(min(len(o), 150000), random_state=0)
    axes[0].hexbin(d["n_genes_by_counts"], d["pct_counts_mt"], gridsize=60, xscale="log", bins="log", cmap="viridis", mincnt=1, linewidths=0)
    axes[0].axvline(min_genes, color="k", lw=0.6, ls="--"); axes[0].axhline(max_pct_mt, color="k", lw=0.6, ls="--")
    axes[0].set_xlabel("Genes per cell"); axes[0].set_ylabel("Mitochondrial reads (%)"); axes[0].set_title("Genes vs mitochondrial reads", loc="left")
    axes[1].hist(o.loc[o.doublet_score.notna(), "doublet_score"], bins=60, color="#777777", linewidth=0)
    axes[1].set_yscale("log"); axes[1].set_xlabel("Scrublet doublet score"); axes[1].set_ylabel("Cells"); axes[1].set_title("Doublet scores, all samples", loc="left")
    d2 = o.sample(min(len(o), 60000), random_state=1)
    axes[2].scatter(d2["total_counts"], d2["n_genes_by_counts"], s=0.4, c="#555555", linewidths=0, rasterized=True)
    axes[2].set_xscale("log"); axes[2].set_yscale("log"); axes[2].set_xlabel("UMI per cell"); axes[2].set_ylabel("Genes per cell"); axes[2].set_title("UMI vs genes", loc="left")
    for ax_, l in zip(axes, "abc"): panel_letter(ax_, l)
    fig.tight_layout(); savefig(fig, "fig03b_qc_thresholds.png")

    # ---- fig03c: why cells were removed
    s2 = summ.set_index("sample_id").loc[order]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    parts_ = [("n_kept", "kept", "#4d4d4d"), ("n_doublet", "predicted doublet", "#E69F00"), ("n_high_mt", f"high MT (>{max_pct_mt:g}%)", "#CC79A7"), ("n_low_genes", f"<{min_genes} genes", "#56B4E9")]
    bottom = np.zeros(len(s2))
    for col, lab, colr in parts_:
        v = 100 * s2[col].values / s2["n_raw"].values; ax.bar(range(len(s2)), v, bottom=bottom, color=colr, label=lab, width=0.8); bottom += v
    ax.set_xticks(range(len(s2))); ax.set_xticklabels(order, rotation=90, fontsize=5.5); color_site_labels(ax)
    ax.set_ylabel("Cells (% of deposited)"); ax.set_ylim(0, 100); ax.margins(x=0.01)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6)
    savefig(fig, "fig03c_qc_removed_per_sample.png")

    # ---- fig03d: per-sample doublet score histograms (only when the Scrublet score arrays are available)
    if diag:
        ncol = 6; nrow = int(np.ceil(len(order) / ncol)); fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.15 * nrow), sharex=True)
        for ax_, sid in zip(axes.ravel(), order):
            sc_obs, sc_sim, auto_, used_ = diag[sid]; bins = np.linspace(0, 1, 41)
            ax_.hist(sc_obs, bins=bins, density=True, color="#777777", linewidth=0)
            if len(sc_sim): ax_.hist(sc_sim, bins=bins, density=True, histtype="step", color="#E69F00", lw=0.6)
            if not np.isnan(auto_): ax_.axvline(auto_, color="0.5", lw=0.6, ls=":")
            ax_.axvline(used_, color="k", lw=0.6); ax_.set_yscale("log"); ax_.set_title(sid, fontsize=5, loc="left", pad=1); ax_.tick_params(labelsize=4.5, length=2)
        for ax_ in axes.ravel()[len(order):]: ax_.axis("off")
        fig.supxlabel("Doublet score (grey: cells; orange: simulated doublets; dotted: automatic threshold; solid: applied threshold)", fontsize=6)
        fig.tight_layout(); savefig(fig, "fig03d_doublet_diagnostics.png")
