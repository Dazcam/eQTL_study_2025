#--------------------------------------------------------------------------------------
#
#    Generate Figure 4: NPC-to-Glu-A pseudotime trajectory
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  14MANUSCRIPT_PLOTS | Rule: fig4_plt
#            Upstream:   17PSEUDOTIME (NPC-to-Glu-UL palantir object and per-bin
#                        tensorQTL perm results); 05TENSORQTL perm (L1 + L2 eGenes)
#            Downstream: None (manuscript figure)
#
# Purpose:   A: Glu-A trajectory -- pseudotime UMAP | PAX6 | SATB2 feature plots
#            B: Glu-A trajectory -- per Q4 bin, proportion of eGenes (qval < 0.05)
#               shared with vs novel to the L1 + L2 eGene set
#            C: Glu-A trajectory -- distance-to-TSS density of significant eQTL
#               (qval < 0.05), one line per Q4 bin. Previously Supplementary
#               Figure 9 (manuscript_plot_tss_density.R, now retired).
#            The Glu-B (NPC-to-Glu-DL) trajectory is no longer shown.
#
# Inputs:    glu_a_h5ad      NPC-to-Glu-UL pseudotime object (palantir_pseudotime,
#                            X_umap)
#            perm_dir        L1 / L2 tensorQTL perm results
#            pseudotime_dir  Per-bin tensorQTL perm results per trajectory
#            cell_types      config['cell_types']; L1 + L2 define the reference
#                            eGene set (trajectory bins excluded)
#            exp_pc_map      Final expression PC count per cell type / bin
#
# Outputs:   out_file        Figure 4 (TIFF, 600 dpi, LZW)
#
# Notes:     geno_pc (4) and norm_method (quantile) are fixed to the final eQTL
#            settings used throughout the manuscript.
#            Panel C reproduces ggplot2::geom_density() defaults: Gaussian kernel,
#            bw.nrd0 bandwidth, 512 points across the x-axis limits, and eQTL
#            beyond +/-600 kb dropped before the density is estimated (as
#            scale_x_continuous(limits = ...) did in the R version).
#
#--------------------------------------------------------------------------------------

##  Load packages, functions and variables  -------------------------------------------
import logging
import os
import sys
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import LinearSegmentedColormap, to_hex
from matplotlib.ticker import PercentFormatter

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Set up logging for Snakemake (uncaught errors also go to the log)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Generating Figure 4 (panels A/B/C) ...")

# Input and output paths
h5ad_path      = str(snakemake.input.glu_a_h5ad)
perm_dir       = str(snakemake.params.perm_dir)
pseudotime_dir = str(snakemake.params.pseudotime_dir)
cell_types     = list(snakemake.params.cell_types)
exp_pc_map     = dict(snakemake.params.exp_pc_map)
out_file       = str(snakemake.output[0])

trajectory    = "NPC-to-Glu-UL"   # pipeline label; shown as "NPC to Glu-A"
terminal_gene = "SATB2"
geno_pc       = 4
norm_method   = "quantile"
quantile_n    = 4
tss_limit_kb  = 600               # panel C x-axis limits (+/- kb)

# Report what each variable is set to
variables = {
    "h5ad_path":      h5ad_path,
    "perm_dir":       perm_dir,
    "pseudotime_dir": pseudotime_dir,
    "cell_types":     ", ".join(cell_types),
    "out_file":       out_file,
    "trajectory":     trajectory,
    "terminal_gene":  terminal_gene,
    "geno_pc":        geno_pc,
    "norm_method":    norm_method,
    "quantile_n":     quantile_n,
    "tss_limit_kb":   tss_limit_kb,
}
width = max(map(len, variables))
log.info("Variables\n============================\n" +
         "\n".join(f"{k:<{width}}  {v}" for k, v in variables.items()) +
         "\n============================")

# Font conventions: match the rest of the manuscript and R/ggplot's default
# device font (falls back to DejaVu Sans if Arial/Helvetica aren't installed)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 12,
})

# Blue feature-plot colormap (same as umaps_manuscript.png)
blue_cmap = LinearSegmentedColormap.from_list("light_dark_blue", ["#f7fbff", "#08306b"])

# Panel letters share size and offset from each panel's top-left corner
PANEL_LABEL_X    = -0.06
PANEL_LABEL_Y    = 1.08
PANEL_LABEL_SIZE = 22

# Panel B styling
BIN_ORDER     = [4, 3, 2, 1]   # bottom -> top, so Bin 1 is at the top
COLOURS       = {"Novel": "#B3B3B3", "Shared": "#B2DF8A"}   # gray70 / light green
BAR_LINEWIDTH = 1.8

# Panel C styling: same bin palette as SF9, i.e. R's
# colorRampPalette(c("#FDE725", "#21908C", "#440154"))(6), bins 1..4 use the first 4
_ramp = LinearSegmentedColormap.from_list("bin_ramp", ["#FDE725", "#21908C", "#440154"])
BIN_PALETTE = {b: to_hex(_ramp(x)) for b, x in zip(range(1, 7), np.linspace(0, 1, 6))}
DENSITY_LINEWIDTH = 1.4        # ~ ggplot linewidth 0.5 mm
DENSITY_N         = 512        # ggplot geom_density default


def relabel(x):
    """Pipeline label -> figure label, e.g. NPC-to-Glu-UL -> NPC to Glu-A."""
    return x.replace("Glu-UL", "Glu-A").replace("Glu-DL", "Glu-B").replace("-to-", " to ")


def read_sig(perm_file):
    """Significant (qval < 0.05) rows of a tensorQTL perm file, or None if missing."""
    if not os.path.exists(perm_file):
        log.warning(f"Not found: {perm_file}")
        return None
    log.info(f"Reading: {perm_file}")
    df = pd.read_csv(perm_file, sep="\t")
    return df.loc[df["qval"].notna() & (df["qval"] < 0.05)]


def read_egenes_ct(ct):
    """eGenes (qval < 0.05) for one L1 / L2 cell type, or None if the file is missing."""
    pc = exp_pc_map[ct]
    perm_file = os.path.join(
        perm_dir,
        f"{ct}_{norm_method}_genPC_{geno_pc}_expPC_{pc}",
        f"{ct}_{norm_method}_perm.cis_qtl.txt.gz")
    sig = read_sig(perm_file)
    return None if sig is None else set(sig["phenotype_id"])


def read_bin_sig(bin_n):
    """Significant perm rows for one pseudotime bin, or None if the file is missing."""
    pc = exp_pc_map[f"{trajectory}-Q{quantile_n}-Bin{bin_n}"]
    perm_file = os.path.join(
        pseudotime_dir, trajectory, "tensorqtl", "perm",
        f"Q{quantile_n}_bin{bin_n}_genPC{geno_pc}_expPC{pc}",
        f"Q{quantile_n}_bin{bin_n}_perm.cis_qtl.txt.gz")
    return read_sig(perm_file)


def bw_nrd0(x):
    """R's bw.nrd0 (Silverman's rule of thumb), the geom_density default."""
    hi = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))   # R quantile type 7
    lo = min(hi, iqr / 1.34)
    if lo == 0:
        lo = hi or abs(x[0]) or 1.0
    return 0.9 * lo * len(x) ** -0.2


def gaussian_density(x, grid):
    """Gaussian KDE of x evaluated on grid, bandwidth bw.nrd0 (as geom_density)."""
    bw = bw_nrd0(x)
    z = (grid[:, None] - x[None, :]) / bw
    return np.exp(-0.5 * z ** 2).sum(axis=1) / (len(x) * bw * np.sqrt(2 * np.pi))


## Load per-bin data (used by panels B and C)  ------------------------------------------
log.info(f"Loading {trajectory} Q{quantile_n} bin perm results ...")
bin_sig = {}   # {bin: significant perm rows}
for b in range(1, quantile_n + 1):
    sig = read_bin_sig(b)
    if sig is not None:
        bin_sig[b] = sig

if not bin_sig:
    raise FileNotFoundError(f"No bin perm files found for {trajectory} in {pseudotime_dir}")


## Panel B data: Q4 bin overlap with L1 + L2 eGene set  --------------------------------
log.info("Loading L1 + L2 reference eGenes ...")
paper_cell_types = [ct for ct in cell_types if "-to-" not in ct]
paper_egene_set = set()
for ct in paper_cell_types:
    genes = read_egenes_ct(ct)
    if genes is not None:
        paper_egene_set |= genes
log.info(f"Unique L1 + L2 eGenes: {len(paper_egene_set)} from {len(paper_cell_types)} cell types")

overlap_counts = {}   # {bin: {"Novel": n, "Shared": n}}
for b, sig in bin_sig.items():
    bin_genes = set(sig["phenotype_id"])
    overlap_counts[b] = {"Shared": len(bin_genes & paper_egene_set),
                         "Novel":  len(bin_genes - paper_egene_set)}
    log.info(f"Bin {b}: {overlap_counts[b]}")


## Panel C data: TSS distances of significant eQTL  ------------------------------------
tss_kb = {}   # {bin: distances (kb) within the axis limits}
n_outside = 0
for b, sig in bin_sig.items():
    d = sig["start_distance"].dropna().to_numpy(dtype=float) / 1000
    n_outside += int((np.abs(d) > tss_limit_kb).sum())
    tss_kb[b] = d[np.abs(d) <= tss_limit_kb]
    log.info(f"Bin {b}: {len(d)} significant eQTL with a TSS distance")

if n_outside > 0:
    log.info(f"Note: {n_outside} eQTL beyond +/-{tss_limit_kb} kb are outside the "
             "panel C x-axis limits and excluded from the density")


## Build figure  ------------------------------------------------------------------------
# Row 1: A (three UMAP panels, 2 cols each).
# Row 2: B (bar plot) | B legend | C (density) | C legend
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(nrows=2, ncols=6, height_ratios=[1, 0.8], wspace=0.4, hspace=0.5)
gs_row2 = gs[1, :].subgridspec(1, 4, width_ratios=[1, 0.55, 1, 0.12], wspace=0)


## Panel A: pseudotime UMAP + feature plots  -------------------------------------------
log.info(f"Panel A: loading {h5ad_path}")
adata = ad.read_h5ad(h5ad_path)
log.info(f"Panel A: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

umap1, umap2 = adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1]
pt_vals = adata.obs["palantir_pseudotime"].values
pt_size = max(0.2, 20000 / adata.n_obs)


def get_expr(gene):
    """Dense expression vector for one gene, or None if absent."""
    if gene not in adata.var_names:
        log.warning(f"Panel A: {gene} not found in adata.var_names")
        return None
    expr = adata.X[:, adata.var_names.get_loc(gene)]
    return np.asarray(expr.todense() if sp.issparse(expr) else expr).flatten()


# Pseudotime UMAP (viridis), no title
ax0 = fig.add_subplot(gs[0, 0:2])
sc0 = ax0.scatter(umap1, umap2, c=pt_vals, cmap="viridis",
                  s=pt_size, alpha=0.5, rasterized=True, linewidths=0)
plt.colorbar(sc0, ax=ax0, shrink=0.6, label="Pseudotime")
ax0.axis("off")
ax0.text(PANEL_LABEL_X, PANEL_LABEL_Y, "A", transform=ax0.transAxes,
         fontsize=PANEL_LABEL_SIZE, fontweight="bold")

# Feature plots: progenitor (PAX6) and terminal marker
for col_start, gene_name in zip([2, 4], ["PAX6", terminal_gene]):
    ax = fig.add_subplot(gs[0, col_start:col_start + 2])
    expr = get_expr(gene_name)
    if expr is not None:
        vmax = float(np.percentile(expr, 99))
        scg = ax.scatter(umap1, umap2, c=expr, cmap=blue_cmap,
                         s=pt_size, alpha=0.5, rasterized=True,
                         linewidths=0, vmin=0, vmax=max(vmax, 1e-6))
        plt.colorbar(scg, ax=ax, shrink=0.6, label="Expression")
    else:
        ax.text(0.5, 0.5, f"{gene_name}\nnot found", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="grey")
    ax.set_title(gene_name, fontsize=12, fontweight="bold")
    ax.axis("off")

del adata


## Panel B: Q4 overlap bar plot  -------------------------------------------------------
ax = fig.add_subplot(gs_row2[0, 0])

y_pos = np.arange(len(BIN_ORDER))
novel_vals  = np.array([overlap_counts.get(b, {}).get("Novel", 0) for b in BIN_ORDER], dtype=float)
shared_vals = np.array([overlap_counts.get(b, {}).get("Shared", 0) for b in BIN_ORDER], dtype=float)
totals = novel_vals + shared_vals
totals[totals == 0] = 1   # avoid divide-by-zero on missing bins
novel_prop  = novel_vals / totals
shared_prop = shared_vals / totals

# Shared first (left), Novel after (right)
ax.barh(y_pos, shared_prop, color=COLOURS["Shared"], edgecolor="black",
        linewidth=BAR_LINEWIDTH, height=0.7)
ax.barh(y_pos, novel_prop, left=shared_prop, color=COLOURS["Novel"],
        edgecolor="black", linewidth=BAR_LINEWIDTH, height=0.7)

# Count labels centred within each stacked segment
for yi, (sv, sp_, nv, np_) in enumerate(zip(shared_vals, shared_prop, novel_vals, novel_prop)):
    if sv > 0:
        ax.text(sp_ / 2, yi, int(sv), ha="center", va="center", fontsize=10, fontweight="bold")
    if nv > 0:
        ax.text(sp_ + np_ / 2, yi, int(nv), ha="center", va="center", fontsize=10, fontweight="bold")

ax.set_yticks(y_pos)
ax.set_yticklabels([f"Bin {b}" for b in BIN_ORDER], fontweight="bold", fontsize=11)
# Pad beyond 0/1 so the bars' left/right edges aren't clipped to half width
ax.set_xlim(-0.005, 1.005)
ax.xaxis.set_major_formatter(PercentFormatter(xmax=1))
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.tick_params(axis="both", length=0, labelsize=11)
ax.set_xlabel("Proportion of eGenes", fontsize=13)
ax.set_ylabel("Pseudotime bin", fontsize=13)
ax.set_title(relabel(trajectory), fontweight="bold")

# No axis lines; light grey vertical gridlines behind the bars
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="lightgrey", linewidth=0.8)
ax.yaxis.grid(False)

ax.text(PANEL_LABEL_X, PANEL_LABEL_Y, "B", transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE, fontweight="bold")

# Legend to the right of panel B: square patches with black borders to match the bars
legend_handles = [
    mpatches.Patch(facecolor=COLOURS["Shared"], edgecolor="black", linewidth=BAR_LINEWIDTH, label="Shared"),
    mpatches.Patch(facecolor=COLOURS["Novel"],  edgecolor="black", linewidth=BAR_LINEWIDTH, label="Novel"),
]
ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.03, 0.5),
          frameon=False, fontsize=12, handlelength=1.2, handleheight=1.2)


## Panel C: TSS-distance density per bin  ----------------------------------------------
ax = fig.add_subplot(gs_row2[0, 2])

grid = np.linspace(-tss_limit_kb, tss_limit_kb, DENSITY_N)
density_handles = []
for b in sorted(tss_kb):
    d = tss_kb[b]
    if len(d) < 2:
        log.warning(f"Panel C: bin {b} has {len(d)} eQTL in range; density not drawn")
        continue
    ax.plot(grid, gaussian_density(d, grid), color=BIN_PALETTE[b],
            linewidth=DENSITY_LINEWIDTH)
    density_handles.append(mlines.Line2D([], [], color=BIN_PALETTE[b],
                                         linewidth=DENSITY_LINEWIDTH * 2, label=str(b)))

ax.axvline(0, linestyle="--", color="black", linewidth=1)
ax.set_xlim(-tss_limit_kb, tss_limit_kb)
ax.set_xticks(np.arange(-tss_limit_kb, tss_limit_kb + 1, 200))
ax.set_ylim(bottom=0)
ax.tick_params(axis="both", length=0, labelsize=11)
ax.set_xlabel("Distance to TSS (kb)", fontsize=13)
ax.set_ylabel("Density", fontsize=13)
ax.set_title(relabel(trajectory), fontweight="bold")

# Match panel B: no axis lines, light grey gridlines behind the data
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_axisbelow(True)
ax.grid(True, color="lightgrey", linewidth=0.8)

ax.text(PANEL_LABEL_X, PANEL_LABEL_Y, "C", transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE, fontweight="bold")

ax.legend(handles=density_handles, title="Bin", loc="center left",
          bbox_to_anchor=(1.03, 0.5), frameon=False, fontsize=12,
          title_fontproperties={"size": 12, "weight": "bold"}, handlelength=1.5)


## Save  --------------------------------------------------------------------------------
log.info(f"Writing: {out_file}")
fig.savefig(out_file, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
