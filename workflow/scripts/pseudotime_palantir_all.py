# --------------------------------------------------------------------------------------
#
#    pseudotime_palantir_all.py
#
#    Palantir pseudotime and fate probability estimation across all neurogenic
#    trajectories of the prenatal human cerebral cortex snRNA-seq atlas.
#
#    Pipeline:
#      1.  Load full atlas h5ad
#      2.  Assign L1 cell type labels from leiden_harmony_0.2
#      3.  Filter to neurogenic populations (exclude MG, GABA, Endo-Peri)
#      4.  Recompute neighbour graph on Harmony PCA embedding
#      5.  Palantir diffusion map and multi-scale space
#      6.  Root cell selection — max root_gene expression in NPC
#      7.  Terminal cell selection — max marker gene per terminal population
#      8.  Run Palantir with multiple terminal states
#      9.  Figure 1 — six-panel UMAP (cell type, pseudotime, PCW, PAX6, entropy, violin)
#      10. Figure 2 — fate probability UMAPs (one per terminal + dominant fate)
#      11. Figure 3 — NPC fate assignment (UMAP, bar, violin, scatter)
#      12. Figure 4 — fate probability by PCW
#      13. Figure 5 — marker gene expression UMAPs (3x3 grid)
#      14. Export fate probability CSV
#      15. Save annotated h5ad
#
#    Memory management:
#      - Full atlas loaded with backed='r' then subsetted before copy
#      - No duplicate adata objects created during filtering
#      - gc.collect() called after figures
#
#    Usage: called via Snakemake rule pseudotime_palantir_all
#
# --------------------------------------------------------------------------------------

import sys
import gc
import logging
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import anndata as ad
import palantir

from pathlib import Path
from scanpy_gene_lists import cluster_anns, custom_palette

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

RANDOM_SEED = 1234
np.random.seed(RANDOM_SEED)

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

sc.settings.verbosity = 1

# ── Parameters ────────────────────────────────────────────────────────────────

input_h5ad    = str(snakemake.input.h5ad)
output_h5ad   = str(snakemake.output.h5ad)
output_dir    = Path(output_h5ad).parent
output_dir.mkdir(parents=True, exist_ok=True)

exclude_types = list(snakemake.params.exclude_types)
root_gene     = str(snakemake.params.root_gene)
terminals     = dict(snakemake.params.terminals)
marker_genes  = list(snakemake.params.marker_genes)
n_neighbors   = int(snakemake.params.n_neighbors)
n_pcs         = int(snakemake.params.n_pcs)
n_dcs         = int(snakemake.params.n_dcs)
num_waypoints = int(snakemake.params.num_waypoints)

# Derive output paths from output_h5ad directory
fate_csv   = str(output_dir / "npc_fate_probabilities.csv")
fig_umap   = str(output_dir / "palantir_all_umap.png")
fig_fate   = str(output_dir / "palantir_all_fate_probs.png")
fig_npc    = str(output_dir / "palantir_all_npc_fates.png")
fig_pcw    = str(output_dir / "palantir_all_fate_by_pcw.png")
fig_markers= str(output_dir / "palantir_all_markers.png")

terminal_names = list(terminals.keys())
terminal_genes = list(terminals.values())

logger.info(f"Input h5ad:       {input_h5ad}")
logger.info(f"Output directory: {output_dir}")
logger.info(f"Exclude types:    {exclude_types}")
logger.info(f"Root gene:        {root_gene}")
logger.info(f"Terminals:        {terminals}")
logger.info(f"Marker genes:     {marker_genes}")
logger.info(f"n_neighbors:      {n_neighbors}")
logger.info(f"n_pcs:            {n_pcs}")
logger.info(f"n_dcs:            {n_dcs}")
logger.info(f"num_waypoints:    {num_waypoints}")

# ── 1. Load data ──────────────────────────────────────────────────────────────

logger.info(f"\nLoading full atlas from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad, backed='r')
logger.info(f"Backed load: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# ── 2. Assign L1 cell type labels ────────────────────────────────────────────

logger.info("Assigning L1 cell type labels from leiden_harmony_0.2...")
cell_type_series = adata.obs['leiden_harmony_0.2'].map(cluster_anns)

n_unmapped = cell_type_series.isna().sum()
if n_unmapped > 0:
    unmapped = adata.obs.loc[
        cell_type_series.isna(), 'leiden_harmony_0.2'
    ].unique().tolist()
    raise ValueError(
        f"{n_unmapped:,} cells could not be mapped to L1 labels. "
        f"Unmapped cluster IDs: {unmapped}. Check cluster_anns in scanpy_gene_lists."
    )

full_dist = cell_type_series.value_counts().sort_index()
logger.info(f"\nFull atlas L1 distribution:\n{full_dist.to_string()}")

# ── 3. Filter to neurogenic populations ──────────────────────────────────────
#
#  Exclude non-neurogenic cell types (MG, GABA, Endo-Peri) in backed mode
#  before copying to memory — avoids loading the full atlas into RAM.
#  Single copy operation only; no intermediate adata objects created.

logger.info(f"\nFiltering out: {exclude_types}")
keep_mask = ~cell_type_series.isin(exclude_types)
n_keep    = keep_mask.sum()
n_remove  = (~keep_mask).sum()
logger.info(f"  Cells retained: {n_keep:,}")
logger.info(f"  Cells removed:  {n_remove:,}")

adata = adata[keep_mask].to_memory()
adata.obs['cell_type'] = cell_type_series[keep_mask].values

filtered_dist = adata.obs['cell_type'].value_counts().sort_index()
logger.info(f"\nFiltered cell type distribution:\n{filtered_dist.to_string()}")

# Validate required fields
for field in ['X_pca_harmony', 'X_umap']:
    if field not in adata.obsm:
        raise KeyError(f"adata.obsm['{field}'] not found")

if 'counts' not in adata.layers:
    raise KeyError("adata.layers['counts'] not found")

logger.info(f"X_pca_harmony shape: {adata.obsm['X_pca_harmony'].shape}")
logger.info(f"X_umap shape:        {adata.obsm['X_umap'].shape}")

# ── 4. Recompute neighbour graph ──────────────────────────────────────────────

logger.info(f"\nRecomputing neighbour graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
sc.pp.neighbors(
    adata,
    n_neighbors  = n_neighbors,
    n_pcs        = n_pcs,
    use_rep      = 'X_pca_harmony',
    random_state = RANDOM_SEED
)
logger.info("Neighbour graph recomputed.")

# ── 5. Palantir diffusion map ─────────────────────────────────────────────────

logger.info(f"\nComputing Palantir diffusion map (n_components={n_dcs})...")
pca_df = pd.DataFrame(
    adata.obsm['X_pca_harmony'],
    index=adata.obs_names
)

dm_res = palantir.utils.run_diffusion_maps(
    pca_df,
    n_components = n_dcs,
    seed         = RANDOM_SEED
)

logger.info("Computing multi-scale space...")
ms_data = palantir.utils.determine_multiscale_space(dm_res)

adata.obsm['X_palantir_dm'] = dm_res['EigenVectors'].values
adata.obsm['X_palantir_ms'] = ms_data.values
logger.info(f"Diffusion map shape: {adata.obsm['X_palantir_dm'].shape}")

# ── 6. Root cell selection ────────────────────────────────────────────────────
#
#  Select NPC cell with maximum root_gene expression.
#  Fallback: NPC cell closest to centroid in multi-scale space.

logger.info(f"\nSelecting root cell (gene: {root_gene})...")
npc_mask = adata.obs['cell_type'] == 'NPC'
if npc_mask.sum() == 0:
    raise ValueError("No NPC cells found after filtering")

npc_idx      = np.where(npc_mask)[0]
root_cell_id = None
root_method  = None

if root_gene in adata.var_names:
    gene_col  = adata.var_names.get_loc(root_gene)
    gene_expr = adata.X[npc_idx, gene_col]
    if sp.issparse(gene_expr):
        gene_expr = np.array(gene_expr.todense()).flatten()
    else:
        gene_expr = np.array(gene_expr).flatten()

    max_expr = gene_expr.max()
    if max_expr > 0:
        root_local_idx  = np.argmax(gene_expr)
        root_global_idx = npc_idx[root_local_idx]
        root_cell_id    = adata.obs_names[root_global_idx]
        root_method     = f"{root_gene} max expression ({max_expr:.3f})"
    else:
        logger.warning(f"{root_gene} expression is zero in all NPCs — falling back to centroid")
else:
    logger.warning(f"{root_gene} not in adata.var_names — falling back to centroid")

if root_cell_id is None:
    npc_ms          = adata.obsm['X_palantir_ms'][npc_idx]
    centroid        = npc_ms.mean(axis=0)
    dists           = np.linalg.norm(npc_ms - centroid, axis=1)
    root_local_idx  = np.argmin(dists)
    root_global_idx = npc_idx[root_local_idx]
    root_cell_id    = adata.obs_names[root_global_idx]
    root_method     = f"NPC centroid (dist={dists[root_local_idx]:.4f})"

logger.info(f"Root cell:   {root_cell_id}")
logger.info(f"Root method: {root_method}")

# ── 7. Terminal cell selection ────────────────────────────────────────────────
#
#  For each terminal population, select the cell with maximum marker gene
#  expression within that cell type.
#  Fallback: cell closest to centroid in multi-scale space.

logger.info("\nSelecting terminal cells...")
terminal_cell_ids = {}
terminal_methods  = {}

for traj_name, term_gene in terminals.items():
    logger.info(f"\n  Terminal: {traj_name} (gene: {term_gene})")

    term_mask = adata.obs['cell_type'] == traj_name
    if term_mask.sum() == 0:
        raise ValueError(f"No cells found for terminal population '{traj_name}'")

    term_idx         = np.where(term_mask)[0]
    terminal_cell_id = None
    terminal_method  = None

    if term_gene in adata.var_names:
        gene_col  = adata.var_names.get_loc(term_gene)
        gene_expr = adata.X[term_idx, gene_col]
        if sp.issparse(gene_expr):
            gene_expr = np.array(gene_expr.todense()).flatten()
        else:
            gene_expr = np.array(gene_expr).flatten()

        max_expr = gene_expr.max()
        if max_expr > 0:
            term_local_idx   = np.argmax(gene_expr)
            term_global_idx  = term_idx[term_local_idx]
            terminal_cell_id = adata.obs_names[term_global_idx]
            terminal_method  = f"{term_gene} max expression ({max_expr:.3f})"
        else:
            logger.warning(f"{term_gene} zero in all {traj_name} — falling back to centroid")
    else:
        logger.warning(f"{term_gene} not in adata.var_names — falling back to centroid")

    if terminal_cell_id is None:
        term_ms          = adata.obsm['X_palantir_ms'][term_idx]
        centroid         = term_ms.mean(axis=0)
        dists            = np.linalg.norm(term_ms - centroid, axis=1)
        term_local_idx   = np.argmin(dists)
        term_global_idx  = term_idx[term_local_idx]
        terminal_cell_id = adata.obs_names[term_global_idx]
        terminal_method  = f"{traj_name} centroid (dist={dists[term_local_idx]:.4f})"

    terminal_cell_ids[traj_name] = terminal_cell_id
    terminal_methods[traj_name]  = terminal_method
    logger.info(f"  Cell:   {terminal_cell_id}")
    logger.info(f"  Method: {terminal_method}")

# ── 8. Run Palantir ───────────────────────────────────────────────────────────

logger.info(f"\nRunning Palantir (num_waypoints={num_waypoints})...")
logger.info(f"  Root:      {root_cell_id}")
logger.info(f"  Terminals: {terminal_cell_ids}")

try:
    pr_res = palantir.core.run_palantir(
        ms_data,
        early_cell      = root_cell_id,
        terminal_states = list(terminal_cell_ids.values()),
        num_waypoints   = num_waypoints,
        seed            = RANDOM_SEED
    )

    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obs['palantir_entropy']    = pr_res.entropy

    # Store fate probabilities — columns named by terminal cell IDs
    # Rename to trajectory names for interpretability
    fate_df = pr_res.branch_probs.copy()
    cell_id_to_name = {v: k for k, v in terminal_cell_ids.items()}
    fate_df.columns = [cell_id_to_name.get(c, c) for c in fate_df.columns]

    for traj_name in terminal_names:
        col = f'fate_prob_{traj_name}'
        if traj_name in fate_df.columns:
            adata.obs[col] = fate_df[traj_name]
        else:
            logger.warning(f"Fate probability column missing for {traj_name}")
            adata.obs[col] = np.nan

    # Dominant fate per cell
    fate_cols = [f'fate_prob_{t}' for t in terminal_names]
    adata.obs['dominant_fate'] = adata.obs[fate_cols].idxmax(axis=1).str.replace(
        'fate_prob_', ''
    )

    pt_summary = adata.obs.groupby('cell_type')['palantir_pseudotime'].agg(
        ['mean', 'median', 'std', 'min', 'max']
    ).round(4)
    logger.info(f"\nPseudotime summary by cell type:\n{pt_summary.to_string()}")

    # Fate probability summary per cell type
    for traj_name in terminal_names:
        col = f'fate_prob_{traj_name}'
        fp_summary = adata.obs.groupby('cell_type')[col].agg(
            ['mean', 'median', 'std']
        ).round(4)
        logger.info(f"\nFate probability ({traj_name}) by cell type:\n{fp_summary.to_string()}")

    # NPC fate assignment summary
    npc_obs = adata.obs[adata.obs['cell_type'] == 'NPC']
    npc_fate_counts = npc_obs['dominant_fate'].value_counts()
    logger.info(f"\nNPC dominant fate assignment:\n{npc_fate_counts.to_string()}")

    palantir_success = True

except Exception as e:
    logger.error(f"Palantir run failed: {e}")
    palantir_success = False

if not palantir_success:
    raise RuntimeError("Palantir failed — see log for details.")

# ── Common plot variables ─────────────────────────────────────────────────────

cell_types_present = adata.obs['cell_type'].unique().tolist()
palette    = {ct: custom_palette.get(ct, '#AAAAAA') for ct in cell_types_present}
umap1      = adata.obsm['X_umap'][:, 0]
umap2      = adata.obsm['X_umap'][:, 1]
pt_vals    = adata.obs['palantir_pseudotime'].values
ent_vals   = adata.obs['palantir_entropy'].values
pt_size    = max(0.2, 20000 / adata.n_obs)

fate_colours = {
    'Glu-DL': '#00B6EB',
    'Glu-UL': '#4363d8',
    'OPC':    '#f58231',
    'ambiguous': '#AAAAAA'
}

# ── 9. Figure 1 — Six-panel UMAP ─────────────────────────────────────────────

logger.info("\nGenerating Figure 1 — six-panel UMAP...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

# Panel A — cell type
for ct in sorted(cell_types_present):
    mask  = adata.obs['cell_type'] == ct
    color = palette.get(ct, '#AAAAAA')
    axes[0].scatter(
        umap1[mask], umap2[mask],
        c=color, s=pt_size, alpha=0.5,
        label=ct, rasterized=True, linewidths=0
    )
axes[0].legend(markerscale=6, fontsize=7, loc='best', framealpha=0.8)
axes[0].set_title('Cell type', fontsize=11)
axes[0].axis('off')

# Panel B — pseudotime
sc_b = axes[1].scatter(
    umap1, umap2, c=pt_vals, cmap='viridis',
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_b, ax=axes[1], shrink=0.6, label='Pseudotime')
axes[1].set_title('Palantir pseudotime', fontsize=11)
axes[1].axis('off')

# Panel C — PCW
pcw_vals = adata.obs['PCW'].values.astype(float)
pcw_cmap = plt.cm.get_cmap('RdYlBu_r', len(np.unique(pcw_vals)))
sc_c = axes[2].scatter(
    umap1, umap2, c=pcw_vals, cmap=pcw_cmap,
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_c, ax=axes[2], shrink=0.6, label='PCW')
axes[2].set_title('Post-conceptional week (PCW)', fontsize=11)
axes[2].axis('off')

# Panel D — root gene (PAX6)
if root_gene in adata.var_names:
    gene_idx = adata.var_names.get_loc(root_gene)
    expr     = adata.X[:, gene_idx]
    if sp.issparse(expr):
        expr = np.array(expr.todense()).flatten()
    else:
        expr = np.array(expr).flatten()
    vmax  = float(np.percentile(expr, 99))
    sc_d  = axes[3].scatter(
        umap1, umap2, c=expr, cmap='Reds',
        s=pt_size, alpha=0.5, rasterized=True,
        linewidths=0, vmin=0, vmax=vmax
    )
    plt.colorbar(sc_d, ax=axes[3], shrink=0.6, label='Expression')
axes[3].set_title(root_gene, fontsize=13, fontweight='bold')
axes[3].axis('off')

# Panel E — entropy
sc_e = axes[4].scatter(
    umap1, umap2, c=ent_vals, cmap='magma',
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_e, ax=axes[4], shrink=0.6, label='Entropy')
axes[4].set_title('Palantir entropy', fontsize=11)
axes[4].axis('off')

# Panel F — pseudotime violin by cell type
ct_order = (adata.obs.groupby('cell_type')['palantir_pseudotime']
            .median().sort_values().index.tolist())
sc.pl.violin(
    adata, keys='palantir_pseudotime', groupby='cell_type',
    order=ct_order, palette=palette, ax=axes[5],
    show=False, rotation=0, use_raw=False
)
axes[5].set_title('Pseudotime by cell type', fontsize=11)
axes[5].set_ylabel('Pseudotime')
axes[5].set_xlabel('')

fig.suptitle(
    f'Palantir pseudotime — filtered atlas\n'
    f'(n={adata.shape[0]:,} cells, root={root_gene}, '
    f'terminals={", ".join(terminal_names)})',
    fontsize=12, fontweight='bold', y=1.01
)
fig.tight_layout()
fig.savefig(fig_umap, dpi=200, bbox_inches='tight')
plt.close(fig)
gc.collect()
logger.info(f"Figure 1 saved: {fig_umap}")

# ── 10. Figure 2 — Fate probability UMAPs ────────────────────────────────────

logger.info("Generating Figure 2 — fate probability UMAPs...")

n_terms  = len(terminal_names)
fig, axes = plt.subplots(1, n_terms + 1, figsize=(6 * (n_terms + 1), 5))

for i, traj_name in enumerate(terminal_names):
    col    = f'fate_prob_{traj_name}'
    fp_vals = adata.obs[col].values
    sc_f   = axes[i].scatter(
        umap1, umap2, c=fp_vals, cmap='viridis',
        s=pt_size, alpha=0.6, rasterized=True,
        linewidths=0, vmin=0, vmax=1
    )
    plt.colorbar(sc_f, ax=axes[i], shrink=0.7, label='Fate probability')
    axes[i].set_title(f'Fate → {traj_name}', fontsize=11, fontweight='bold')
    axes[i].axis('off')

# Final panel — dominant fate
dom_fate = adata.obs['dominant_fate'].values
for traj_name in terminal_names:
    mask  = dom_fate == traj_name
    color = fate_colours.get(traj_name, '#AAAAAA')
    axes[-1].scatter(
        umap1[mask], umap2[mask],
        c=color, s=pt_size, alpha=0.6,
        label=traj_name, rasterized=True, linewidths=0
    )
axes[-1].legend(markerscale=6, fontsize=8, loc='best', framealpha=0.8)
axes[-1].set_title('Dominant fate', fontsize=11, fontweight='bold')
axes[-1].axis('off')

fig.suptitle('Palantir fate probabilities', fontsize=12, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(fig_fate, dpi=200, bbox_inches='tight')
plt.close(fig)
gc.collect()
logger.info(f"Figure 2 saved: {fig_fate}")

# ── 11. Figure 3 — NPC fate assignment ───────────────────────────────────────

logger.info("Generating Figure 3 — NPC fate assignment...")

npc_mask_bool = adata.obs['cell_type'] == 'NPC'
npc_umap1     = umap1[npc_mask_bool]
npc_umap2     = umap2[npc_mask_bool]
npc_obs       = adata.obs[npc_mask_bool]
npc_pt_size   = max(0.5, 10000 / npc_mask_bool.sum())

fig, axes = plt.subplots(1, 4, figsize=(24, 5))

# Panel A — NPC UMAP coloured by dominant fate
for traj_name in terminal_names:
    mask  = npc_obs['dominant_fate'] == traj_name
    color = fate_colours.get(traj_name, '#AAAAAA')
    axes[0].scatter(
        npc_umap1[mask], npc_umap2[mask],
        c=color, s=npc_pt_size, alpha=0.6,
        label=traj_name, rasterized=True, linewidths=0
    )
axes[0].legend(markerscale=6, fontsize=8, loc='best', framealpha=0.8)
axes[0].set_title('NPC dominant fate', fontsize=11, fontweight='bold')
axes[0].axis('off')

# Panel B — bar chart of NPC fate counts
fate_counts = npc_obs['dominant_fate'].value_counts()
bar_colours = [fate_colours.get(f, '#AAAAAA') for f in fate_counts.index]
axes[1].bar(fate_counts.index, fate_counts.values,
            color=bar_colours, edgecolor='black', linewidth=0.8)
for j, (fate, count) in enumerate(fate_counts.items()):
    axes[1].text(j, count + 50, str(count),
                 ha='center', va='bottom', fontsize=10, fontweight='bold')
axes[1].set_title('NPC fate assignment counts', fontsize=11, fontweight='bold')
axes[1].set_ylabel('Number of NPC cells')
axes[1].set_xlabel('Dominant fate')
axes[1].tick_params(axis='x', labelsize=10)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Panel C — violin of fate probability per terminal split by dominant fate
fate_long = npc_obs[['dominant_fate'] + [f'fate_prob_{t}' for t in terminal_names]].melt(
    id_vars='dominant_fate',
    var_name='fate_terminal',
    value_name='probability'
)
fate_long['fate_terminal'] = fate_long['fate_terminal'].str.replace('fate_prob_', '')

positions = np.arange(len(terminal_names))
width     = 0.25
offsets   = np.linspace(-width, width, len(terminal_names))

for k, dom in enumerate(terminal_names):
    subset = fate_long[fate_long['dominant_fate'] == dom]
    for j, term in enumerate(terminal_names):
        vals = subset[subset['fate_terminal'] == term]['probability'].dropna().values
        if len(vals) == 0:
            continue
        vp = axes[2].violinplot(
            vals, positions=[j + offsets[k]],
            widths=width * 0.8, showmedians=True
        )
        color = fate_colours.get(dom, '#AAAAAA')
        for pc in vp['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        vp['cmedians'].set_color('black')
        vp['cbars'].set_color('black')
        vp['cmins'].set_color('black')
        vp['cmaxes'].set_color('black')

axes[2].set_xticks(positions)
axes[2].set_xticklabels(terminal_names, fontsize=9)
axes[2].set_title('Fate probability by dominant fate', fontsize=11, fontweight='bold')
axes[2].set_ylabel('Fate probability')
axes[2].set_xlabel('Terminal')
legend_patches = [mpatches.Patch(color=fate_colours.get(t, '#AAAAAA'), label=t)
                  for t in terminal_names]
axes[2].legend(handles=legend_patches, fontsize=8, loc='upper right')
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)

# Panel D — scatter of Glu-DL vs Glu-UL fate probability coloured by dominant fate
fp_dl = npc_obs['fate_prob_Glu-DL'].values if 'fate_prob_Glu-DL' in npc_obs.columns else np.zeros(len(npc_obs))
fp_ul = npc_obs['fate_prob_Glu-UL'].values if 'fate_prob_Glu-UL' in npc_obs.columns else np.zeros(len(npc_obs))

for traj_name in terminal_names:
    mask  = npc_obs['dominant_fate'].values == traj_name
    color = fate_colours.get(traj_name, '#AAAAAA')
    axes[3].scatter(
        fp_dl[mask], fp_ul[mask],
        c=color, s=npc_pt_size * 2, alpha=0.5,
        label=traj_name, rasterized=True, linewidths=0
    )
axes[3].set_xlabel('Fate prob → Glu-DL', fontsize=10)
axes[3].set_ylabel('Fate prob → Glu-UL', fontsize=10)
axes[3].set_title('Glu-DL vs Glu-UL fate probability\n(NPC cells only)',
                  fontsize=11, fontweight='bold')
axes[3].legend(markerscale=4, fontsize=8, loc='best', framealpha=0.8)
axes[3].plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=0.8)
axes[3].set_xlim(-0.05, 1.05)
axes[3].set_ylim(-0.05, 1.05)
axes[3].spines['top'].set_visible(False)
axes[3].spines['right'].set_visible(False)

fig.suptitle('NPC fate assignment', fontsize=12, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(fig_npc, dpi=200, bbox_inches='tight')
plt.close(fig)
gc.collect()
logger.info(f"Figure 3 saved: {fig_npc}")

# ── 12. Figure 4 — Fate probability by PCW ───────────────────────────────────

logger.info("Generating Figure 4 — fate probability by PCW...")

pcw_vals_str = adata.obs['PCW'].astype(str)
pcw_unique   = sorted(adata.obs['PCW'].unique())

fig, axes = plt.subplots(1, n_terms, figsize=(6 * n_terms, 5))
if n_terms == 1:
    axes = [axes]

for i, traj_name in enumerate(terminal_names):
    col = f'fate_prob_{traj_name}'
    data_by_pcw = [
        adata.obs.loc[adata.obs['PCW'] == pcw, col].dropna().values
        for pcw in pcw_unique
    ]
    bp = axes[i].boxplot(
        data_by_pcw,
        labels     = [str(p) for p in pcw_unique],
        patch_artist = True,
        medianprops  = dict(color='black', linewidth=1.5),
        flierprops   = dict(marker='.', markersize=2, alpha=0.3)
    )
    color = fate_colours.get(traj_name, '#AAAAAA')
    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    axes[i].set_title(f'Fate prob → {traj_name} by PCW',
                      fontsize=11, fontweight='bold')
    axes[i].set_xlabel('PCW', fontsize=10)
    axes[i].set_ylabel('Fate probability', fontsize=10)
    axes[i].tick_params(axis='x', labelsize=8)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)

fig.suptitle('Fate probability by post-conceptional week',
             fontsize=12, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(fig_pcw, dpi=200, bbox_inches='tight')
plt.close(fig)
gc.collect()
logger.info(f"Figure 4 saved: {fig_pcw}")

# ── 13. Figure 5 — Marker gene expression UMAPs ──────────────────────────────

logger.info("Generating Figure 5 — marker gene expression UMAPs...")

n_genes  = len(marker_genes)
n_cols   = 3
n_rows   = int(np.ceil(n_genes / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
axes      = axes.flatten()

for i, gene in enumerate(marker_genes):
    ax = axes[i]
    if gene in adata.var_names:
        gene_idx = adata.var_names.get_loc(gene)
        expr     = adata.X[:, gene_idx]
        if sp.issparse(expr):
            expr = np.array(expr.todense()).flatten()
        else:
            expr = np.array(expr).flatten()
        vmax  = float(np.percentile(expr, 99))
        sc_g  = ax.scatter(
            umap1, umap2, c=expr, cmap='Reds',
            s=pt_size, alpha=0.6, rasterized=True,
            linewidths=0, vmin=0, vmax=max(vmax, 1e-6)
        )
        plt.colorbar(sc_g, ax=ax, shrink=0.7, label='Expression')
        pct_expr = (expr > 0).mean() * 100
        ax.set_title(f'{gene}\n({pct_expr:.1f}% expressed)',
                     fontsize=11, fontweight='bold')
    else:
        ax.text(0.5, 0.5, f'{gene}\nnot found',
                ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
        ax.set_title(gene, fontsize=11, fontweight='bold')
        logger.warning(f"{gene} not found in adata.var_names")
    ax.axis('off')

# Hide unused axes
for j in range(n_genes, len(axes)):
    axes[j].axis('off')

fig.suptitle('Marker gene expression', fontsize=12, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(fig_markers, dpi=200, bbox_inches='tight')
plt.close(fig)
gc.collect()
logger.info(f"Figure 5 saved: {fig_markers}")

# ── 14. Export fate probability CSV ──────────────────────────────────────────

logger.info("\nExporting fate probability CSV...")
fate_export_cols = (
    ['cell_type', 'palantir_pseudotime', 'palantir_entropy', 'dominant_fate'] +
    [f'fate_prob_{t}' for t in terminal_names]
)
fate_export_cols = [c for c in fate_export_cols if c in adata.obs.columns]

fate_export_df             = adata.obs[fate_export_cols].copy()
fate_export_df.index.name  = 'cell_id'
fate_export_df.to_csv(fate_csv)
logger.info(f"Fate CSV written: {fate_csv} ({fate_export_df.shape})")

# NPC summary for log
npc_fate_df = fate_export_df[fate_export_df['cell_type'] == 'NPC']
logger.info(f"\nNPC fate probability summary:")
logger.info(f"  Total NPC cells: {len(npc_fate_df):,}")
for traj_name in terminal_names:
    col = f'fate_prob_{traj_name}'
    if col in npc_fate_df.columns:
        logger.info(
            f"  {traj_name}: mean={npc_fate_df[col].mean():.3f}, "
            f"median={npc_fate_df[col].median():.3f}, "
            f"dominant={( npc_fate_df['dominant_fate'] == traj_name).sum():,} cells"
        )

# ── 15. Save annotated h5ad ───────────────────────────────────────────────────

logger.info(f"\nSaving annotated h5ad to: {output_h5ad}")
adata.uns['palantir_all_params'] = {
    'exclude_types'    : exclude_types,
    'root_gene'        : root_gene,
    'root_cell_id'     : root_cell_id,
    'root_method'      : root_method,
    'terminals'        : terminals,
    'terminal_cell_ids': terminal_cell_ids,
    'terminal_methods' : terminal_methods,
    'n_neighbors'      : n_neighbors,
    'n_pcs'            : n_pcs,
    'n_dcs'            : n_dcs,
    'num_waypoints'    : num_waypoints,
    'random_seed'      : RANDOM_SEED,
}

adata.write_h5ad(output_h5ad)
logger.info("Done.")
