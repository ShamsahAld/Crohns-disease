"""Shared helpers for the GSE266616 pipeline: paths, palettes, lineage marker sets, figure style, embedding plots."""
import os, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap

BASE = os.environ.get("GSE_BASE", os.path.expanduser("~/.claude-science-scratch/gse266616"))
DATA = f"{BASE}/data"
PROC = os.environ.get("GSE_PROC", f"{BASE}/processed")
OUT = os.environ.get("GSE_OUT", "out")
for d in (PROC, OUT): os.makedirs(d, exist_ok=True)

# ---- colours (Okabe-Ito based; colour is bound to an entity and reused across every figure)
STATUS_ORDER = os.environ.get("GSE_STATUS_ORDER", "non_inflamed,inflamed").split(",")      # groups present, in display order (set to control,non_inflamed,adjacent,inflamed for the all-groups reference)
STATUS_COL = {"control": "#7F7F7F", "non_inflamed": "#0072B2", "adjacent": "#E69F00", "inflamed": "#D55E00"}
SITE_COL = {"TI": "#009E73", "AC": "#CC79A7"}
BIOLOGIC_COL = {"none": "#999999", "anti_TNF": "#56B4E9", "anti_IL12_IL23": "#E69F00"}
PLATFORM_COL = {"GPL24676": "#999999", "GPL18573": "#D55E00"}
LINEAGE_COL = {"epithelial": "#1f77b4", "T_NK": "#d62728", "B": "#2ca02c", "plasma": "#9467bd", "myeloid": "#ff7f0e",
               "mast": "#8c564b", "stromal": "#e377c2", "endothelial": "#17becf", "glial": "#bcbd22",
               "mixed": "#444444", "unassigned": "#bbbbbb"}
IMMUNE = ["T_NK", "B", "plasma", "myeloid", "mast"]
COMPARTMENT = {**{k: "immune" for k in IMMUNE}, "epithelial": "epithelial", "stromal": "stromal/other",
               "endothelial": "stromal/other", "glial": "stromal/other", "mixed": "mixed", "unassigned": "unassigned"}

# Marker sets for PROVISIONAL coarse lineage calls (not the final cell-type annotation; that is roadmap step 7)
LINEAGE_MARKERS = {
    "epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "VIL1", "CDX2", "FABP2", "FABP1", "LGALS4"],
    "T_NK":        ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "CCL5", "NKG7", "GNLY", "TRAC"],
    "B":           ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22"],
    "plasma":      ["MZB1", "JCHAIN", "DERL3", "TNFRSF17", "PRDM1"],   # no IGHA1/IGKC/SDC1: ambient-RNA prone or also epithelial
    "myeloid":     ["LYZ", "CD68", "CD14", "CSF1R", "C1QA", "C1QB", "S100A8", "S100A9", "CD163"],
    "mast":        ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "stromal":     ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "ACTA2", "MYH11", "RGS5"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP", "ACKR1", "RAMP2"],
    "glial":       ["S100B", "PLP1", "GFAP", "SOX10"],
}

# Genes excluded from highly-variable-gene selection (technical / clone- or sex-specific, not cell-state):
HVG_EXCLUDE_REGEX = r"^(MT-|RPS\d|RPL\d|RPSA$|IG[HKL][VJD]\d|TR[ABDG][VJD]\d|HB[ABDEGMQZ]\d?$)"
SEX_GENES = ["XIST", "TSIX", "RPS4Y1", "RPS4Y2", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "PRKY", "TBL1Y"]

DISPLAY = {"non_inflamed": "non-inflamed", "inflamed": "inflamed", "control": "non-IBD control", "adjacent": "adjacent to inflamed", "TI": "terminal ileum", "AC": "ascending colon",
           "anti_TNF": "anti-TNF", "anti_IL12_IL23": "anti-IL-12/23", "none": "no biologic", "GPL24676": "GPL24676", "GPL18573": "GPL18573 (1 sample)"}

def style():
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out", "legend.frameon": False,
        "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42, "axes.titleweight": "normal",
        "axes.titlelocation": "left", "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    })

EXPR_CMAP = LinearSegmentedColormap.from_list("expr", ["#e6e6e6", "#08306b"])

def panel_letter(ax, letter, dx=-0.02, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")

def embed_frame(ax, xlab="UMAP1", ylab="UMAP2"):
    """Embedding axes: no ticks/spines, small corner arrow pair naming the axes (figure-style 6.6)."""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    kw = dict(xycoords="axes fraction", textcoords="axes fraction",
              arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.2", shrinkA=0, shrinkB=0, mutation_scale=5))
    ax.annotate("", xy=(0.16, 0.0), xytext=(0.0, 0.0), **kw)
    ax.annotate("", xy=(0.0, 0.16), xytext=(0.0, 0.0), **kw)
    ax.text(0.0, -0.02, xlab, transform=ax.transAxes, ha="left", va="top", fontsize=6)        # runs right from the origin
    ax.text(-0.02, 0.0, ylab, transform=ax.transAxes, ha="right", va="bottom", rotation=90, fontsize=6)   # runs up from the origin

def embed_cat(ax, xy, labels, palette, s=0.3, label_groups=False, legend=False, title=None, xlab="UMAP1", ylab="UMAP2",
              min_label_cells=200, legend_title=None, seed=0, legend_ncol=3):
    labels = np.asarray(labels).astype(str)
    order = np.random.default_rng(seed).permutation(len(labels))       # shuffled draw order: no category hides another
    cols = np.array([palette.get(l, "#999999") for l in labels])
    ax.scatter(xy[order, 0], xy[order, 1], c=cols[order], s=s, linewidths=0, rasterized=True)
    ax.margins(0.03); embed_frame(ax, xlab, ylab)
    if label_groups:
        for g in pd.unique(labels):
            m = labels == g
            if m.sum() < min_label_cells: continue
            ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), DISPLAY.get(g, g), fontsize=6, ha="center", va="center",
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])
    if legend:                                                          # legend sits below the panel, right-aligned (never over data / neighbours)
        present = set(labels)
        hs = [mpl.lines.Line2D([], [], marker="o", ls="", ms=3.5, color=palette[g], label=DISPLAY.get(g, g)) for g in palette if g in present]
        ax.legend(handles=hs, loc="upper right", bbox_to_anchor=(1.0, -0.09), frameon=False, fontsize=6, handletextpad=0.1,
                  columnspacing=0.8, borderaxespad=0, ncol=min(len(hs), legend_ncol), title=legend_title, title_fontsize=6)
    if title: ax.set_title(title, loc="left")

def embed_cont(ax, xy, v, title=None, cmap=EXPR_CMAP, vmin=None, vmax=None, s=0.3, cbar=True, xlab="UMAP1", ylab="UMAP2",
               italic=False, cbar_label=None):
    v = np.asarray(v, float); o = np.argsort(v)                           # highest values drawn last
    if vmax is None: vmax = np.nanpercentile(v, 99.5)
    if vmin is None: vmin = np.nanmin(v)
    if vmax <= vmin: vmax = vmin + 1e-6
    h = ax.scatter(xy[o, 0], xy[o, 1], c=v[o], s=s, linewidths=0, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    ax.margins(0.03); embed_frame(ax, xlab, ylab)
    if cbar:                                                            # inset colour bar: does not shrink the panel
        cax = ax.inset_axes([1.02, 0.30, 0.04, 0.40]); cb = ax.figure.colorbar(h, cax=cax)
        cb.ax.tick_params(labelsize=5, length=2, width=0.4); cb.outline.set_visible(False)
        if cbar_label: cb.set_label(cbar_label, fontsize=5)
    if title: ax.set_title(title, loc="left", fontstyle="italic" if italic else "normal")

def embed_highlight(ax, xy, mask, color, title, s=0.3, **kw):
    """Grey background of all cells, one subset in colour (split view)."""
    ax.scatter(xy[~mask, 0], xy[~mask, 1], c="#dddddd", s=s, linewidths=0, rasterized=True)
    idx = np.where(mask)[0]; idx = np.random.default_rng(0).permutation(idx)
    ax.scatter(xy[idx, 0], xy[idx, 1], c=color, s=s, linewidths=0, rasterized=True, alpha=0.6)
    ax.margins(0.03); embed_frame(ax, **kw); ax.set_title(title, loc="left")

def eta2(y, g):
    """Fraction of variance of y explained by categorical grouping g (one-way ANOVA effect size)."""
    codes = pd.Categorical(np.asarray(g).astype(str)).codes
    cnt = np.bincount(codes); mu = np.bincount(codes, weights=y) / cnt
    ss_tot = ((y - y.mean()) ** 2).sum()
    return float((cnt * (mu - y.mean()) ** 2).sum() / ss_tot) if ss_tot > 0 else np.nan

def eta2_table(scores, obs, factors, n_pcs=20):
    out = {f: [eta2(scores[:, i], obs[f].values) for i in range(min(n_pcs, scores.shape[1]))] for f in factors}
    return pd.DataFrame(out, index=[f"PC{i+1}" for i in range(min(n_pcs, scores.shape[1]))])

def heat_with_values(ax, df, vmax=None, cmap="Blues", fs=5, title=None):
    v = df.values.astype(float); vmax = vmax or max(0.05, np.nanmax(v))
    im = ax.imshow(v, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns, rotation=40, ha="right")
    ax.set_yticks(range(df.shape[0])); ax.set_yticklabels(df.index)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            if np.isfinite(v[i, j]):
                ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center", fontsize=fs, color="white" if v[i, j] > 0.6 * vmax else "black")
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if title: ax.set_title(title, loc="left")
    return im

def heat_grouped(ax, df, groups, vmax=None, fs=5.5, title=None):
    """Heatmap whose columns are in blocks; groups = [(block label, n_columns), ...] drawn above the blocks."""
    im = heat_with_values(ax, df, vmax=vmax, fs=fs)
    x0 = 0
    for lab, n in groups:
        ax.text(x0 + (n - 1) / 2, -0.85, lab, ha="center", va="bottom", fontsize=6)
        if x0 > 0: ax.axvline(x0 - 0.5, color="white", lw=2.5)
        x0 += n
    if title: ax.set_title(title, loc="left", pad=32)
    return im

def savefig(fig, name):
    p = os.path.join(OUT, name); fig.savefig(p); plt.close(fig); print("saved", p, flush=True); return p


def assign_lineage(adata, cluster_key="leiden"):
    """Provisional lineage per Leiden cluster from marker-gene scores. Adds score_<lineage>, coarse_lineage and compartment to adata.obs;
    returns the cluster x lineage mean-score table with the assigned label."""
    import scanpy as sc
    for name, genes in LINEAGE_MARKERS.items():
        sc.tl.score_genes(adata, [g for g in genes if g in adata.var_names], score_name=f"score_{name}", random_state=0)
    S = adata.obs[[f"score_{n}" for n in LINEAGE_MARKERS]].copy(); S.columns = list(LINEAGE_MARKERS)
    cm = S.groupby(adata.obs[cluster_key].values, observed=True).mean()
    z = (cm - cm.mean()) / cm.std(ddof=0).replace(0, 1)
    assign = {}
    for kk in cm.index:
        best = z.loc[kk].idxmax(); rawbest = cm.loc[kk, best]; other = cm.loc[kk].drop(best)
        assign[kk] = "unassigned" if rawbest < 0.1 else ("mixed" if (other.max() >= 0.8 * rawbest and other.max() >= 0.5) else best)
    adata.obs["coarse_lineage"] = adata.obs[cluster_key].map(assign).astype(str)
    adata.obs["compartment"] = adata.obs["coarse_lineage"].map(COMPARTMENT).astype(str)
    return cm.assign(assigned=pd.Series(assign), n_cells=adata.obs[cluster_key].value_counts().reindex(cm.index).values)
