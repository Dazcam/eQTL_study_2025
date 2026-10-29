# --------------------------------------------------------------------------------------
#
#    Palantir pseudotime for one neurogenic trajectory
#
# --------------------------------------------------------------------------------------
#
# Pipeline:  17PSEUDOTIME | Rule: palantair_trajectory (one job per trajectory)
#            Upstream:   subset_trajectory
#            Downstream: pseudobulk_trajectory
#
# Purpose:   1. Load the trajectory subset (NPC + one Glu population)
#            2. Validate inputs (Harmony PCA embedding, counts layer)
#            3. Recompute the neighbour graph on the Harmony PCA embedding
#            4. Palantir diffusion map and multi-scale space
#            5. Root cell: max root_gene expression in the root type, centroid fallback
#            6. Terminal cell: max TERMINAL_MARKERS gene in the terminal type, same
#               fallback
#            7. Run Palantir pseudotime with the explicit root and terminal cells
#            8. Recompute UMAP on the subset
#            9. Six-panel figure (cell type, pseudotime, PCW, PAX6, terminal marker,
#               violin); 9b. marker gene feature plots + terminal overlay panel
#            10. Export per-cell summary CSV
#            11. Save annotated h5ad
#
# Inputs:    h5ad             Trajectory subset (adata_subset.h5ad)
#
# Outputs:   h5ad             Annotated h5ad with palantir_pseudotime / _entropy
#            (undeclared, written next to the h5ad)
#            pseudotime_umap.png, pseudotime_markers.png, pseudotime_summary.csv
#
# Notes:     TERMINAL_MARKERS, MARKER_GENES and the overlay genes are set below and
#            can be changed without touching the rest of the script. The Harmony
#            embedding from the full atlas is reused to keep the batch correction.
#
# --------------------------------------------------------------------------------------

import sys
import logging
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import palantir

from pathlib import Path
from scanpy_gene_lists import custom_palette
from pseudotime_functions import get_gene_expr, select_anchor_cell

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

RANDOM_SEED = 1234
np.random.seed(RANDOM_SEED)

# Terminal marker genes for each Glu population — used to select the terminal
# cell for Palantir. Modify here if markers change.
TERMINAL_MARKERS = {
    'Glu-DL': 'TLE4',
    'Glu-UL': 'SATB2',
}

# Marker genes for feature plot — upper and deep layer cortical markers.
# Genes absent from adata.var_names are shown as blank panels with a warning.
MARKER_GENES = [
    'ADRA2A', 'B3GALT2', 'CBLN2',  'CUX1',   'CUX2',
    'ETV1',   'FOXP2',   'GLRA3',  'NECTIN3', 'NR4A2',
    'PDE1A',  'SYT6',    'TLE4',   'TOX',     'EOMES',
    'GAS7',   'GRIK3',   'HS3ST4', 'UNC5D',   'NKAIN2',
    'PLXNA4', 'BCL11B',  'FEZF2',  'RORB',    'SATB2',
    'PDGFRA', 'GAD1'
]

# Genes to mark terminal cell positions on the overlay panel.
# Shows where the max-expression cell for each gene sits on the UMAP.
# Does not affect Palantir terminal selection — TLE4/SATB2 are still used for that.
TERMINAL_OVERLAY_GENES   = ['TLE4', 'FEZF2', 'BCL11B', 'FOXP2', 'HS3ST4']
TERMINAL_OVERLAY_COLOURS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#fce205']

## Set up logging for Snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

sc.settings.verbosity = 1

log.info("Running Palantir pseudotime ...")

##  Variables  -------------------------------------------------------------------------
input_h5ad    = str(snakemake.input.h5ad)
output_h5ad   = str(snakemake.output.h5ad)
output_dir    = Path(output_h5ad).parent

n_neighbors   = snakemake.params.n_neighbors
n_pcs         = snakemake.params.n_pcs
n_dcs         = snakemake.params.n_dcs
num_waypoints = snakemake.params.num_waypoints
root_gene     = snakemake.params.root_gene

output_figure  = str(output_dir / "pseudotime_umap.png")
output_markers = str(output_dir / "pseudotime_markers.png")
output_csv     = str(output_dir / "pseudotime_summary.csv")

variables = {
    "input_h5ad":     input_h5ad,
    "output_h5ad":    output_h5ad,
    "output_figure":  output_figure,
    "output_markers": output_markers,
    "output_csv":     output_csv,
    "n_neighbors":    n_neighbors,
    "n_pcs":          n_pcs,
    "n_dcs":          n_dcs,
    "num_waypoints":  num_waypoints,
    "root_gene":      root_gene,
}
width = max(map(len, variables))
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  1. Load data  ----------------------------------------------------------------------
log.info(f"Loading subset adata from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
log.info(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# Read trajectory metadata stored by subset script
traj       = adata.uns['trajectory']
trajectory = traj['name']
root_type  = traj['root_type']
term_type  = traj['terminal_type']
cell_types = traj['cell_types']

log.info(f"Trajectory:    {trajectory}")
log.info(f"Root:          {root_type}")
log.info(f"Terminal:      {term_type}")

cell_dist = adata.obs['cell_type'].value_counts().sort_index()
log.info(f"\nCell type distribution:\n{cell_dist.to_string()}")

##  2. Validate inputs  ----------------------------------------------------------------
if 'X_pca_harmony' not in adata.obsm:
    raise KeyError(
        "adata.obsm['X_pca_harmony'] not found. "
        "Harmony PCA embedding is required for neighbour graph recomputation."
    )
log.info(f"X_pca_harmony shape: {adata.obsm['X_pca_harmony'].shape}")

if 'counts' not in adata.layers:
    raise KeyError(
        "adata.layers['counts'] not found — required for downstream pseudobulk."
    )
log.info("counts layer present")

##  3. Recompute neighbour graph  ------------------------------------------------------
#
#  The stored graph from the full atlas has edges to all 7 populations.
#  Recompute on the subset using Harmony PCA to preserve batch correction,
#  with a configurable n_pcs appropriate for the two-population subset.

log.info(f"Recomputing neighbour graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
sc.pp.neighbors(
    adata,
    n_neighbors  = n_neighbors,
    n_pcs        = n_pcs,
    use_rep      = 'X_pca_harmony',
    random_state = RANDOM_SEED
)
log.info("Neighbour graph recomputed.")

##  4. Palantir diffusion map  ---------------------------------------------------------
#
#  Palantir computes multi-scale diffusion distances using an adaptive Gaussian
#  kernel. This is more robust than Scanpy DPT for trajectories with weak
#  inter-population connectivity because it weights diffusion components by
#  their eigenvalues, reducing the influence of noisy high-order components.
#
#  Input: Harmony PCA as a DataFrame — Palantir requires a DataFrame with
#  cell barcodes as the index.

log.info(f"Computing Palantir diffusion map (n_components={n_dcs})...")
pca_df = pd.DataFrame(
    adata.obsm['X_pca_harmony'],
    index=adata.obs_names
)

dm_res = palantir.utils.run_diffusion_maps(
    pca_df,
    n_components = n_dcs,
    seed         = RANDOM_SEED
)

log.info("Computing multi-scale space...")
ms_data = palantir.utils.determine_multiscale_space(dm_res)

# Store diffusion components in adata for reference
adata.obsm['X_palantir_dm'] = dm_res['EigenVectors'].values
adata.obsm['X_palantir_ms'] = ms_data.values
log.info(f"Diffusion map shape: {adata.obsm['X_palantir_dm'].shape}")

##  5. Root cell selection  ------------------------------------------------------------
#
#  Primary: cell with highest root_gene expression within root cell type.
#  PAX6 is highest in apical radial glia — the earliest, most uncommitted
#  progenitor state — making PAX6-max the most biologically principled root.
#
#  Fallback: NPC cell closest to centroid in Palantir multi-scale space.

root_cell_id, root_method = select_anchor_cell(adata, root_type, root_gene)

log.info(f"Root cell:             {root_cell_id}")
log.info(f"Root selection method: {root_method}")

##  6. Terminal cell selection  --------------------------------------------------------
#
#  Select the cell with highest terminal marker gene expression within the
#  terminal cell type. Mirrors the root cell selection strategy — deterministic
#  and biologically motivated.
#
#  Fallback: terminal cell closest to centroid in Palantir multi-scale space.

terminal_gene = TERMINAL_MARKERS.get(term_type)
if terminal_gene is None:
    raise ValueError(
        f"No terminal marker defined for '{term_type}' in TERMINAL_MARKERS. "
        f"Available: {list(TERMINAL_MARKERS.keys())}"
    )

terminal_cell_id, terminal_method = select_anchor_cell(adata, term_type, terminal_gene)

log.info(f"Terminal cell:             {terminal_cell_id}")
log.info(f"Terminal selection method: {terminal_method}")

##  7. Run Palantir pseudotime  --------------------------------------------------------
log.info(f"Running Palantir (num_waypoints={num_waypoints})...")
try:
    pr_res = palantir.core.run_palantir(
        ms_data,
        early_cell      = root_cell_id,
        terminal_states = [terminal_cell_id],
        num_waypoints   = num_waypoints,
        seed            = RANDOM_SEED
    )

    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obs['palantir_entropy']    = pr_res.entropy

    pt_summary = adata.obs.groupby('cell_type')['palantir_pseudotime'].agg(
        ['mean', 'median', 'std', 'min', 'max']
    ).round(4)
    log.info(f"\nPalantir pseudotime summary by cell type:\n{pt_summary.to_string()}")

    palantir_success = True

except Exception as e:
    log.error(f"Palantir run failed: {e}")
    palantir_success = False

if not palantir_success:
    raise RuntimeError(
        "Palantir pseudotime computation failed — see log for details."
    )

##  8. Recompute UMAP on subset  -------------------------------------------------------
log.info("Computing UMAP on subset...")
sc.tl.umap(adata, random_state=RANDOM_SEED)
log.info("UMAP computed.")

##  Common plot variables  -------------------------------------------------------------
palette  = {ct: custom_palette[ct] for ct in cell_types if ct in custom_palette}
umap1    = adata.obsm['X_umap'][:, 0]
umap2    = adata.obsm['X_umap'][:, 1]
pt_vals  = adata.obs['palantir_pseudotime'].values
pt_size  = max(0.2, 20000 / adata.n_obs)

##  9. Six-panel figure  ---------------------------------------------------------------
#
#  A: cell type                D: PAX6 expression
#  B: Palantir pseudotime      E: terminal marker gene (TLE4 or SATB2)
#  C: PCW                      F: pseudotime violin by cell type
#
#  All UMAP panels use direct matplotlib scatter to avoid scanpy rendering
#  issues at high cell counts.

log.info("Generating six-panel figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

# Panel A — cell type
for ct in adata.obs['cell_type'].cat.categories:
    mask  = adata.obs['cell_type'] == ct
    color = palette.get(ct, '#AAAAAA')
    axes[0].scatter(
        umap1[mask], umap2[mask],
        c=color, s=pt_size, alpha=0.5,
        label=ct, rasterized=True, linewidths=0
    )
axes[0].legend(markerscale=6, fontsize=8, loc='best', framealpha=0.8)
axes[0].set_title('Cell type', fontsize=11)
axes[0].axis('off')

# Panel B — Palantir pseudotime
sc_b = axes[1].scatter(
    umap1, umap2, c=pt_vals, cmap='viridis',
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_b, ax=axes[1], shrink=0.6, label='Pseudotime')
axes[1].set_title('Palantir pseudotime', fontsize=11)
axes[1].axis('off')

# Panel C — PCW
pcw_vals = adata.obs['PCW'].values.astype(float)
pcw_cmap = plt.get_cmap('RdYlBu_r', len(np.unique(pcw_vals)))
sc_c = axes[2].scatter(
    umap1, umap2, c=pcw_vals, cmap=pcw_cmap,
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_c, ax=axes[2], shrink=0.6, label='PCW')
axes[2].set_title('Post-conception week (PCW)', fontsize=11)
axes[2].axis('off')

# Panel D — root gene (PAX6) / Panel E — terminal marker gene
plot_genes = ['PAX6', TERMINAL_MARKERS[term_type]]

for ax, gene in zip(axes[3:5], plot_genes):
    expr = get_gene_expr(adata, gene)
    if expr is not None:
        vmax  = float(np.percentile(expr, 99))
        sc_g  = ax.scatter(
            umap1, umap2, c=expr, cmap='Reds',
            s=pt_size, alpha=0.5, rasterized=True,
            linewidths=0, vmin=0, vmax=vmax
        )
        plt.colorbar(sc_g, ax=ax, shrink=0.6, label='Expression')
    else:
        ax.text(0.5, 0.5, f'{gene}\nnot found',
                ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
        ax.axis('off')
        log.warning(f"{gene} not found in adata.var_names — panel left blank")
    ax.set_title(gene, fontsize=13, fontweight='bold')
    ax.axis('off')

# Panel F — pseudotime violin by cell type
ct_order = (adata.obs.groupby('cell_type')['palantir_pseudotime']
            .median()
            .sort_values()
            .index.tolist())

sc.pl.violin(
    adata,
    keys    = 'palantir_pseudotime',
    groupby = 'cell_type',
    order   = ct_order,
    palette = palette,
    ax      = axes[5],
    show    = False,
    rotation= 0,
    use_raw = False
)
axes[5].set_title('Palantir pseudotime by cell type', fontsize=11)
axes[5].set_ylabel('Pseudotime')
axes[5].set_xlabel('')

fig.suptitle(
    f'Pseudotime trajectory: {trajectory}\n'
    f'(n={adata.shape[0]:,} cells, root={root_type}, '
    f'method={root_method})',
    fontsize=12, fontweight='bold', y=1.01
)
fig.tight_layout()
fig.savefig(output_figure, dpi=200, bbox_inches='tight')
plt.close(fig)
log.info(f"Figure saved: {output_figure}")

##  9b. Marker gene feature plots  -----------------------------------------------------
#
#  6x5 grid of UMAPs coloured by normalised expression for each gene in
#  MARKER_GENES. Genes absent from adata.var_names shown as blank panels.
#  Percentage of cells expressing each gene shown in the panel title.
#
#  Final panel: cell type UMAP with TERMINAL_OVERLAY_GENES max-expression
#  cell positions marked. Shows where each gene's most committed cell sits
#  on the UMAP — purely for visualisation, does not affect Palantir.

log.info("Generating marker gene feature plots...")

n_genes  = len(MARKER_GENES)
n_cols   = 6
n_panels = n_genes + 1                          # +1 for terminal overlay
n_rows   = int(np.ceil(n_panels / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
axes      = axes.flatten()

for i, gene in enumerate(MARKER_GENES):
    ax = axes[i]
    expr = get_gene_expr(adata, gene)
    if expr is not None:
        pct_expr = (expr > 0).mean() * 100
        vmax     = float(np.percentile(expr, 99))
        sc_g     = ax.scatter(
            umap1, umap2, c=expr, cmap='Reds',
            s=pt_size, alpha=0.5, rasterized=True,
            linewidths=0, vmin=0, vmax=max(vmax, 1e-6)
        )
        plt.colorbar(sc_g, ax=ax, shrink=0.6, label='Expression')
        ax.set_title(f'{gene}\n({pct_expr:.1f}% expressed)',
                     fontsize=10, fontweight='bold')
    else:
        ax.text(0.5, 0.5, f'{gene}\nnot found',
                ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='grey')
        ax.set_title(gene, fontsize=10, fontweight='bold')
        log.warning(f"Marker gene not found in adata.var_names: {gene}")
    ax.axis('off')

# Terminal overlay panel — cell type background + max-expression cell positions
overlay_ax = axes[n_genes]
for ct in adata.obs['cell_type'].cat.categories:
    mask  = adata.obs['cell_type'] == ct
    color = palette.get(ct, '#AAAAAA')
    overlay_ax.scatter(
        umap1[mask], umap2[mask],
        c=color, s=pt_size, alpha=0.3,
        rasterized=True, linewidths=0
    )

for gene, colour in zip(TERMINAL_OVERLAY_GENES, TERMINAL_OVERLAY_COLOURS):
    expr = get_gene_expr(adata, gene)
    if expr is None:
        log.warning(f"Terminal overlay gene not found: {gene}")
        continue
    max_idx   = np.argmax(expr)
    max_u1    = umap1[max_idx]
    max_u2    = umap2[max_idx]
    overlay_ax.scatter(
        max_u1, max_u2,
        c=colour, s=150, zorder=5,
        edgecolors='black', linewidths=0.8,
        label=f'{gene}-max'
    )
    overlay_ax.annotate(
        gene,
        xy=(max_u1, max_u2),
        xytext=(5, 5), textcoords='offset points',
        fontsize=9, fontweight='bold', color=colour
    )

overlay_ax.legend(markerscale=1.2, fontsize=8, loc='best', framealpha=0.8)
overlay_ax.set_title('Deep/upper layer terminal gene\nmax-expression positions',
                     fontsize=10, fontweight='bold')
overlay_ax.axis('off')

# Hide remaining unused axes
for j in range(n_panels, len(axes)):
    axes[j].axis('off')

fig.suptitle(
    f'Marker gene expression: {trajectory}',
    fontsize=12, fontweight='bold', y=1.01
)
fig.tight_layout()
fig.savefig(output_markers, dpi=200, bbox_inches='tight')
plt.close(fig)
log.info(f"Marker figure saved: {output_markers}")

##  10. Export summary CSV  ------------------------------------------------------------
log.info("Exporting pseudotime summary CSV...")
summary_cols = [
    'cell_type', 'leiden_harmony_0.2', 'sample', 'PCW',
    'palantir_pseudotime', 'palantir_entropy'
]
summary_df = adata.obs[summary_cols].copy()
summary_df.index.name = 'cell_id'
summary_df.to_csv(output_csv)
log.info(f"Summary CSV written: {output_csv} ({summary_df.shape})")

##  11. Save annotated h5ad  -----------------------------------------------------------
adata.uns['pseudotime_params'] = {
    'trajectory'      : trajectory,
    'root_type'       : root_type,
    'terminal_type'   : term_type,
    'root_gene'       : root_gene,
    'root_cell_id'    : root_cell_id,
    'root_method'     : root_method,
    'terminal_gene'   : terminal_gene,
    'terminal_cell_id': terminal_cell_id,
    'terminal_method' : terminal_method,
    'n_neighbors'     : n_neighbors,
    'n_pcs'           : n_pcs,
    'n_dcs'           : n_dcs,
    'num_waypoints'   : num_waypoints,
    'random_seed'     : RANDOM_SEED,
}

log.info(f"Writing annotated adata to: {output_h5ad}")
adata.write_h5ad(output_h5ad)
log.info("Done.")

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
