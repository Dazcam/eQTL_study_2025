# --------------------------------------------------------------------------------------
#
#    manuscript_plots_fig4.py
#
#    Full Figure 4, generated entirely in Python (no R/raster-embedding step):
#
#      A: Glu-A (NPC-to-Glu-UL) -- pseudotime UMAP | PAX6 | SATB2
#      B: Glu-B (NPC-to-Glu-DL) -- pseudotime UMAP | PAX6 | TLE4
#      C: Glu-A -- Q4 eGene overlap with full paper eGene set (Novel/Shared only)
#      D: Glu-B -- Q4 eGene overlap with full paper eGene set (Novel/Shared only)
#
#    Usage: called via Snakemake rule fig4_plt
#
# --------------------------------------------------------------------------------------

import sys
import os
import logging
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import anndata as ad

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Row definitions for panels A/B: (row label, h5ad path, trajectory, terminal gene)
ROWS = [
    ("A", str(snakemake.input.glu_a_h5ad), "NPC-to-Glu-UL", "SATB2"),
    ("B", str(snakemake.input.glu_b_h5ad), "NPC-to-Glu-DL", "TLE4"),
]

perm_dir       = str(snakemake.params.perm_dir)
pseudotime_dir = str(snakemake.params.pseudotime_dir)
cell_types     = list(snakemake.params.cell_types)
exp_pc_map     = dict(snakemake.params.exp_pc_map)

geno_pc = 4
norm_method = "quantile"

out_file = str(snakemake.output[0])

# ── Logging ───────────────────────────────────────────────────────────────────

log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)

logger.info("Generating full Figure 4 (panels A/B/C/D) ...")

# ── Font conventions -- match the rest of the manuscript (umaps_manuscript.png),
# and match R/ggplot's default device font (falls back to DejaVu Sans if
# Arial/Helvetica aren't installed in this container) ------------------------

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 12,
})

# ── Blue feature-plot colormap -- same as umaps_manuscript.png's blue_cmap ------

blue_cmap = LinearSegmentedColormap.from_list("light_dark_blue", ["#f7fbff", "#08306b"])

# Panel letters (A-D) all share the same size and the same relative offset
# from their own panel's top-left corner, so they read as one consistent set
# rather than four differently-sized labels.
PANEL_LABEL_X = -0.06
PANEL_LABEL_Y = 1.08
PANEL_LABEL_SIZE = 22


def relabel(x):
    return x.replace("Glu-UL", "Glu-A").replace("Glu-DL", "Glu-B").replace("-to-", " to ")

# ── Panels C/D data: Q4 overlap with full paper eGene set (Novel/Shared) ---------

def read_egenes_ct(ct, exp_pc_map, perm_dir, geno_pc=4, norm_method="quantile"):
    pc = exp_pc_map[ct]
    perm_file = os.path.join(
        perm_dir,
        f"{ct}_{norm_method}_genPC_{geno_pc}_expPC_{pc}",
        f"{ct}_{norm_method}_perm.cis_qtl.txt.gz"
    )
    if not os.path.exists(perm_file):
        logger.warning(f"Not found: {perm_file}")
        return None
    df = pd.read_csv(perm_file, sep="\t")
    df = df[df['qval'].notna() & (df['qval'] < 0.05)]
    return df['phenotype_id'].unique()


def read_perm_egenes(pseudotime_dir, trajectory, bin_n, exp_pc_map, geno_pc=4, quantile_n=4):
    pc_key = f"{trajectory}-Q{quantile_n}-Bin{bin_n}"
    pc = exp_pc_map[pc_key]
    perm_file = os.path.join(
        pseudotime_dir, trajectory, "tensorqtl", "perm",
        f"Q{quantile_n}_bin{bin_n}_genPC{geno_pc}_expPC{pc}",
        f"Q{quantile_n}_bin{bin_n}_perm.cis_qtl.txt.gz"
    )
    if not os.path.exists(perm_file):
        logger.warning(f"Not found: {perm_file}")
        return None
    df = pd.read_csv(perm_file, sep="\t")
    df = df[df['qval'].notna() & (df['qval'] < 0.05)]
    return set(df['phenotype_id'].unique())


logger.info("Loading L1+L2 paper reference eGenes...")
paper_cell_types = [ct for ct in cell_types if "-to-" not in ct]
paper_egene_lists = []
for ct in paper_cell_types:
    genes = read_egenes_ct(ct, exp_pc_map, perm_dir, geno_pc=geno_pc, norm_method=norm_method)
    if genes is not None:
        paper_egene_lists.append(genes)
paper_egene_set = set(np.concatenate(paper_egene_lists)) if paper_egene_lists else set()
logger.info(f"Total unique paper eGenes (L1 + L2): {len(paper_egene_set)}")

logger.info("Loading Q4 pseudotime bin eGenes and computing overlap...")
overlap_counts = {}  # {trajectory: {bin: {"Novel": n, "Shared": n}}}
for _, _, trajectory, _ in ROWS:
    overlap_counts[trajectory] = {}
    for b in range(1, 5):
        bin_genes = read_perm_egenes(pseudotime_dir, trajectory, b, exp_pc_map, geno_pc=geno_pc)
        if bin_genes is None:
            continue
        n_shared = len(bin_genes & paper_egene_set)
        n_novel = len(bin_genes - paper_egene_set)
        overlap_counts[trajectory][b] = {"Novel": n_novel, "Shared": n_shared}

# ── Build figure: 3 rows x 6 cols (A/B each span 2 cols x 3 panels; C/D each span 3 cols) --

fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(nrows=3, ncols=6, height_ratios=[1, 1, 0.8], wspace=0.4, hspace=0.5)

# --- Panels A/B: pseudotime UMAP + feature plots
for row_idx, (row_label, h5ad_path, trajectory, terminal_gene) in enumerate(ROWS):
    logger.info(f"Row {row_label}: loading {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    logger.info(f"Row {row_label}: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    umap1 = adata.obsm['X_umap'][:, 0]
    umap2 = adata.obsm['X_umap'][:, 1]
    pt_vals = adata.obs['palantir_pseudotime'].values
    pt_size = max(0.2, 20000 / adata.n_obs)

    def get_expr(gene, adata=adata):
        if gene not in adata.var_names:
            logger.warning(f"Row {row_label}: {gene} not found in adata.var_names")
            return None
        idx = adata.var_names.get_loc(gene)
        expr = adata.X[:, idx]
        if sp.issparse(expr):
            expr = np.array(expr.todense()).flatten()
        else:
            expr = np.array(expr).flatten()
        return expr

    # Column 1: pseudotime UMAP (viridis) -- no title, per manuscript request
    ax0 = fig.add_subplot(gs[row_idx, 0:2])
    sc0 = ax0.scatter(umap1, umap2, c=pt_vals, cmap='viridis',
                       s=pt_size, alpha=0.5, rasterized=True, linewidths=0)
    plt.colorbar(sc0, ax=ax0, shrink=0.6, label='Pseudotime')
    ax0.axis('off')
    ax0.text(PANEL_LABEL_X, PANEL_LABEL_Y, row_label, transform=ax0.transAxes,
              fontsize=PANEL_LABEL_SIZE, fontweight='bold')

    # Columns 2-3: gene feature plots (PAX6, terminal marker) -- blue_cmap
    for col_start, gene_name in zip([2, 4], ['PAX6', terminal_gene]):
        ax = fig.add_subplot(gs[row_idx, col_start:col_start + 2])
        expr = get_expr(gene_name)
        if expr is not None:
            vmax = float(np.percentile(expr, 99))
            scg = ax.scatter(umap1, umap2, c=expr, cmap=blue_cmap,
                              s=pt_size, alpha=0.5, rasterized=True,
                              linewidths=0, vmin=0, vmax=max(vmax, 1e-6))
            plt.colorbar(scg, ax=ax, shrink=0.6, label='Expression')
        else:
            ax.text(0.5, 0.5, f'{gene_name}\nnot found', ha='center', va='center',
                    transform=ax.transAxes, fontsize=10, color='grey')
        ax.set_title(gene_name, fontsize=12, fontweight='bold')
        ax.axis('off')

    del adata

# --- Panels C/D: Q4 overlap barplots
bin_order = [4, 3, 2, 1]  # bottom -> top, so Bin 1 ends up at the top (matches prior R version)
colours = {"Novel": "#B3B3B3", "Shared": "#B2DF8A"}  # gray70 / light green, matching the reference figure
BAR_LINEWIDTH = 1.8

# Full-width span (0:6) so C's left edge lines up with A/B's -- narrowing to
# a margin (1:5) fixed D's y-label overlap but knocked C out of alignment and
# left unused whitespace either side. wspace alone is enough to keep D's
# y-axis label clear of C's bars without needing that margin.
cd_gs = gs[2, 0:6].subgridspec(1, 2, wspace=0.45)

for cd_idx, (row_label, _, trajectory, _) in enumerate(ROWS):
    ax = fig.add_subplot(cd_gs[0, cd_idx])
    counts = overlap_counts[trajectory]

    y_pos = np.arange(len(bin_order))
    novel_vals = np.array([counts.get(b, {}).get("Novel", 0) for b in bin_order], dtype=float)
    shared_vals = np.array([counts.get(b, {}).get("Shared", 0) for b in bin_order], dtype=float)
    totals = novel_vals + shared_vals
    totals[totals == 0] = 1  # avoid divide-by-zero on missing bins
    novel_prop = novel_vals / totals
    shared_prop = shared_vals / totals

    # Shared first (left), Novel after (right) -- matches the reference figure's ordering
    ax.barh(y_pos, shared_prop, color=colours["Shared"], edgecolor='black',
            linewidth=BAR_LINEWIDTH, height=0.7, label='Shared')
    ax.barh(y_pos, novel_prop, left=shared_prop, color=colours["Novel"],
            edgecolor='black', linewidth=BAR_LINEWIDTH, height=0.7, label='Novel')

    # Count labels centred within each stacked segment
    for yi, (sv, sp_, nv, np_) in enumerate(zip(shared_vals, shared_prop, novel_vals, novel_prop)):
        if sv > 0:
            ax.text(sp_ / 2, yi, int(sv), ha='center', va='center',
                    fontsize=10, fontweight='bold')
        if nv > 0:
            ax.text(sp_ + np_ / 2, yi, int(nv), ha='center', va='center',
                    fontsize=10, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Bin {b}" for b in bin_order], fontweight='bold', fontsize=11)
    # Slight padding beyond 0/1: the bars' left/right edge lines sit exactly
    # at 0 and 1, and coinciding with the axes clip boundary was clipping off
    # the outward half of the stroke, rendering those edges at half the
    # thickness of the top/bottom edges. Padding moves the clip boundary
    # clear of the bars so the full line width renders on every side.
    ax.set_xlim(-0.005, 1.005)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis='both', length=0, labelsize=11)  # no tick marks, labels kept
    ax.set_xlabel("Proportion of eGenes", fontsize=13)
    ax.set_ylabel("Pseudotime bin", fontsize=13)
    ax.set_title(relabel(trajectory), fontweight='bold')

    # No visible axis lines -- light grey vertical gridlines at the percentage
    # breaks instead, drawn behind the bars
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color='lightgrey', linewidth=0.8)
    ax.yaxis.grid(False)

    panel_label = "C" if cd_idx == 0 else "D"
    ax.text(PANEL_LABEL_X, PANEL_LABEL_Y, panel_label, transform=ax.transAxes,
            fontsize=PANEL_LABEL_SIZE, fontweight='bold')

# Single shared legend for panels C/D -- explicit square Patch handles (not
# the default barh swatches) so they render as solid boxes with a black
# border, matching the bars themselves; positioned just above row C/D, not
# between rows A/B.
legend_handles = [
    mpatches.Patch(facecolor=colours["Shared"], edgecolor='black', linewidth=BAR_LINEWIDTH, label='Shared'),
    mpatches.Patch(facecolor=colours["Novel"], edgecolor='black', linewidth=BAR_LINEWIDTH, label='Novel'),
]
fig.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, 0.30),
           ncol=2, frameon=False, fontsize=12, handlelength=1.2, handleheight=1.2)

fig.savefig(out_file, dpi=600, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
logger.info(f"Figure 4 saved: {out_file}")
logger.info("Done.")
