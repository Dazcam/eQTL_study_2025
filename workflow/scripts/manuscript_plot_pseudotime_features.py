# --------------------------------------------------------------------------------------
#
#    manuscript_plots_pseudotime_features.py
#
#    Figure 4 panels A/B: pseudotime UMAP + gene feature plots, one row per
#    trajectory 
#
#    Usage: called via Snakemake rule fig4_pseudotime_features_plt
#
# --------------------------------------------------------------------------------------

import sys
import logging
import warnings
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import anndata as ad

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Row definitions: (h5ad path, row label, PAX6-partner terminal gene)
ROWS = [
    ("A", str(snakemake.input.glu_a_h5ad), "SATB2"),  # NPC-to-Glu-UL -> Glu-A
    ("B", str(snakemake.input.glu_b_h5ad), "TLE4"),   # NPC-to-Glu-DL -> Glu-B
]

output_png = str(snakemake.output[0])

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

logger.info("Generating Figure 4 panels A/B (pseudotime UMAP + feature plots) ...")

# ── Font conventions -- match the rest of the manuscript (umaps_manuscript.png) ---

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.titlesize": 12,
})

# ── Build figure: 2 rows x 3 cols -----------------------------------------------

fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(nrows=2, ncols=3, figure=fig, wspace=0.35, hspace=0.4)

for row_idx, (row_label, h5ad_path, terminal_gene) in enumerate(ROWS):
    logger.info(f"Row {row_label}: loading {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    logger.info(f"Row {row_label}: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    umap1 = adata.obsm['X_umap'][:, 0]
    umap2 = adata.obsm['X_umap'][:, 1]
    pt_vals = adata.obs['palantir_pseudotime'].values
    pt_size = max(0.2, 20000 / adata.n_obs)

    def get_expr(gene):
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

    # Column 1: pseudotime UMAP
    ax0 = fig.add_subplot(gs[row_idx, 0])
    sc0 = ax0.scatter(umap1, umap2, c=pt_vals, cmap='viridis',
                       s=pt_size, alpha=0.5, rasterized=True, linewidths=0)
    plt.colorbar(sc0, ax=ax0, shrink=0.6, label='Pseudotime')
    ax0.set_title('Pseudotime', fontsize=12)
    ax0.axis('off')
    ax0.set_title(row_label, loc='left', fontweight='bold', fontsize=14)

    # Columns 2-3: gene feature plots (PAX6, terminal marker)
    for col_idx, gene in enumerate(['PAX6', terminal_gene], start=1):
        ax = fig.add_subplot(gs[row_idx, col_idx])
        expr = get_expr(gene)
        if expr is not None:
            vmax = float(np.percentile(expr, 99))
            scg = ax.scatter(umap1, umap2, c=expr, cmap='Reds',
                              s=pt_size, alpha=0.5, rasterized=True,
                              linewidths=0, vmin=0, vmax=max(vmax, 1e-6))
            plt.colorbar(scg, ax=ax, shrink=0.6, label='Expression')
        else:
            ax.text(0.5, 0.5, f'{gene}\nnot found', ha='center', va='center',
                    transform=ax.transAxes, fontsize=10, color='grey')
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.axis('off')

    del adata

fig.savefig(output_png, dpi=600, bbox_inches='tight')
plt.close(fig)
logger.info(f"Figure saved: {output_png}")
logger.info("Done.")
