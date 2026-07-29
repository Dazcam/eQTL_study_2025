# --------------------------------------------------------------------------------------
#
#    pseudotime_subset.py
#
#    Subset the full snRNA-seq atlas to the two cell populations defining a
#    single neurogenic trajectory, in preparation for pseudotime analysis.
#
#    Trajectories are defined in the config as 'NPC-to-Glu-DL' etc. The two
#    cell type labels are extracted by splitting on '-to-', so the wildcard
#    name is the sole input needed at runtime.
#
#    Pipeline position:
#      Input:  adata_clusters_v3.h5ad  (full 679,738-cell atlas)
#      Output: adata_subset.h5ad       (trajectory-specific subset)
#      Next:   pseudotime_paga.py      (neighbour graph, PAGA, DPT)
#
#    Usage: called via Snakemake rule subset_trajectory
#
# --------------------------------------------------------------------------------------

import sys
import logging
import anndata as ad

from scanpy_gene_lists import cluster_anns

# ── Logging ───────────────────────────────────────────────────────────────────────────

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

# ── Parameters ────────────────────────────────────────────────────────────────────────

input_h5ad  = str(snakemake.input.h5ad)
output_h5ad = str(snakemake.output.h5ad)
trajectory  = snakemake.wildcards.trajectory

# Derive the two cell type labels directly from the wildcard name
cell_types = trajectory.split('-to-')
if len(cell_types) != 2:
    raise ValueError(
        f"Trajectory wildcard '{trajectory}' could not be split into two cell types. "
        "Expected format: 'CellTypeA-to-CellTypeB'"
    )
root_type, terminal_type = cell_types
logger.info(f"Trajectory:    {trajectory}")
logger.info(f"Root:          {root_type}")
logger.info(f"Terminal:      {terminal_type}")

# ── 1. Load data ──────────────────────────────────────────────────────────────────────

logger.info(f"Loading adata from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad, backed='r').to_memory()
logger.info(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# ── 2. Assign L1 cell type labels ────────────────────────────────────────────────────

logger.info("Assigning L1 cell type labels from leiden_harmony_0.2...")
adata.obs['cell_type'] = adata.obs['leiden_harmony_0.2'].map(cluster_anns)

n_unmapped = adata.obs['cell_type'].isna().sum()
if n_unmapped > 0:
    unmapped = adata.obs.loc[
        adata.obs['cell_type'].isna(), 'leiden_harmony_0.2'
    ].unique().tolist()
    raise ValueError(
        f"{n_unmapped:,} cells could not be mapped to L1 labels. "
        f"Unmapped cluster IDs: {unmapped}. Check cluster_anns in scanpy_gene_lists."
    )
logger.info("All cells successfully mapped to L1 labels")

cell_type_dist = adata.obs['cell_type'].value_counts().sort_index()
logger.info(f"\nFull atlas L1 cell type distribution:\n{cell_type_dist.to_string()}")

# ── 3. Validate requested cell types are present ──────────────────────────────────────

for ct in cell_types:
    if ct not in adata.obs['cell_type'].values:
        raise ValueError(
            f"Cell type '{ct}' not found in adata.obs['cell_type']. "
            f"Available types: {sorted(adata.obs['cell_type'].unique().tolist())}"
        )

# ── 4. Subset ─────────────────────────────────────────────────────────────────────────

mask   = adata.obs['cell_type'].isin(cell_types)
adata  = adata[mask].copy()
logger.info(f"\nSubset to {cell_types}:")
logger.info(f"  Cells retained: {adata.shape[0]:,} of {mask.shape[0]:,}")

subset_dist = adata.obs['cell_type'].value_counts().sort_index()
logger.info(f"\nSubset cell type distribution:\n{subset_dist.to_string()}")

# ── 5. Validate counts layer ──────────────────────────────────────────────────────────

if 'counts' not in adata.layers:
    raise KeyError(
        "adata.layers['counts'] not found in subset. "
        "Raw counts are required for downstream pseudobulk aggregation."
    )
logger.info(f"counts layer present: shape {adata.layers['counts'].shape}")

# ── 6. Store trajectory metadata ──────────────────────────────────────────────────────

adata.uns['trajectory'] = {
    'name'          : trajectory,
    'root_type'     : root_type,
    'terminal_type' : terminal_type,
    'cell_types'    : cell_types,
}
logger.info(f"Trajectory metadata stored in adata.uns['trajectory']")

# ── 7. Save subset h5ad ───────────────────────────────────────────────────────────────

logger.info(f"Writing subset adata to: {output_h5ad}")
adata.write_h5ad(output_h5ad)
logger.info(f"Done.")
