# --------------------------------------------------------------------------------------
#
#    Subset the atlas to one neurogenic trajectory
#
# --------------------------------------------------------------------------------------
#
# Pipeline:  17PSEUDOTIME | Rule: subset_trajectory (one job per trajectory)
#            Upstream:   02SCANPY clustering (adata_clusters_v3.h5ad)
#            Downstream: palantair_trajectory
#
# Purpose:   - Assign L1 cell type labels from leiden_harmony_0.2
#            - Keep only the two cell types named in the trajectory wildcard,
#              split on '-to-' (e.g. NPC-to-Glu-DL -> NPC, Glu-DL)
#            - Check the raw counts layer is present for pseudobulk
#            - Store trajectory metadata in adata.uns['trajectory']
#
# Inputs:    h5ad   Full clustered atlas
#
# Outputs:   h5ad   Trajectory subset (adata_subset.h5ad)
#
# --------------------------------------------------------------------------------------

import sys
import logging
import anndata as ad

from scanpy_gene_lists import cluster_anns

## Set up logging for Snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Subsetting atlas to trajectory ...")

##  Variables  -------------------------------------------------------------------------
input_h5ad  = str(snakemake.input.h5ad)
output_h5ad = str(snakemake.output.h5ad)
trajectory  = snakemake.wildcards.trajectory

variables = {
    "input_h5ad":  input_h5ad,
    "output_h5ad": output_h5ad,
    "trajectory":  trajectory,
}
width = max(map(len, variables))
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

# Derive the two cell type labels directly from the wildcard name
cell_types = trajectory.split('-to-')
if len(cell_types) != 2:
    raise ValueError(
        f"Trajectory wildcard '{trajectory}' could not be split into two cell types. "
        "Expected format: 'CellTypeA-to-CellTypeB'"
    )
root_type, terminal_type = cell_types
log.info(f"Trajectory:    {trajectory}")
log.info(f"Root:          {root_type}")
log.info(f"Terminal:      {terminal_type}")

##  1. Load data  ----------------------------------------------------------------------
log.info(f"Loading adata from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad, backed='r').to_memory()
log.info(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

##  2. Assign L1 cell type labels  -----------------------------------------------------
log.info("Assigning L1 cell type labels from leiden_harmony_0.2...")
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
log.info("All cells successfully mapped to L1 labels")

cell_type_dist = adata.obs['cell_type'].value_counts().sort_index()
log.info(f"\nFull atlas L1 cell type distribution:\n{cell_type_dist.to_string()}")

##  3. Validate requested cell types are present  --------------------------------------
for ct in cell_types:
    if ct not in adata.obs['cell_type'].values:
        raise ValueError(
            f"Cell type '{ct}' not found in adata.obs['cell_type']. "
            f"Available types: {sorted(adata.obs['cell_type'].unique().tolist())}"
        )

##  4. Subset  -------------------------------------------------------------------------
mask   = adata.obs['cell_type'].isin(cell_types)
adata  = adata[mask].copy()
log.info(f"\nSubset to {cell_types}:")
log.info(f"  Cells retained: {adata.shape[0]:,} of {mask.shape[0]:,}")

subset_dist = adata.obs['cell_type'].value_counts().sort_index()
log.info(f"\nSubset cell type distribution:\n{subset_dist.to_string()}")

##  5. Validate counts layer  ----------------------------------------------------------
if 'counts' not in adata.layers:
    raise KeyError(
        "adata.layers['counts'] not found in subset. "
        "Raw counts are required for downstream pseudobulk aggregation."
    )
log.info(f"counts layer present: shape {adata.layers['counts'].shape}")

##  6. Store trajectory metadata  ------------------------------------------------------
adata.uns['trajectory'] = {
    'name'          : trajectory,
    'root_type'     : root_type,
    'terminal_type' : terminal_type,
    'cell_types'    : cell_types,
}
log.info("Trajectory metadata stored in adata.uns['trajectory']")

##  7. Save subset h5ad  ---------------------------------------------------------------
log.info(f"Writing subset adata to: {output_h5ad}")
adata.write_h5ad(output_h5ad)
log.info("Done.")

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
