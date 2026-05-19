#!/usr/bin/env python3
# --------------------------------------------------------------------------------------
#
#    prepare_singler_query.py
#
#    Preprocesses the snRNA-seq query data for SingleR classification against
#    the Qian et al. 2025 MERFISH reference.
#    - Applies cluster_anns mapping to assign L1 cell type labels
#    - Sets adata.X to raw counts from adata.layers['counts']
#    - Subsets features to genes shared with the Qian MERFISH 300-gene panel
#    - Reports gene overlap statistics in log
#    - Writes a query h5ad ready for SingleR in R
#
#    Usage: called via Snakemake rule prep_singler_query
#
# --------------------------------------------------------------------------------------

import sys
import logging
import numpy as np
import anndata as ad
import pandas as pd
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

input_query_h5ad = str(snakemake.input.query_h5ad)
input_ref_h5ad   = str(snakemake.input.ref_h5ad)
output_h5ad      = str(snakemake.output.h5ad)
output_umap_csv  = str(snakemake.output.umap_csv)

# L1 cell type mapping from leiden_harmony_0.2 clusters
logger.info("Cluster annotations mapping:")
for k, v in cluster_anns.items():
    logger.info(f"  leiden_harmony_0.2 cluster {k} -> {v}")

# ── Load query ────────────────────────────────────────────────────────────────────────

logger.info(f"Loading query snRNA-seq data from: {input_query_h5ad}")
adata = ad.read_h5ad(input_query_h5ad, backed='r').to_memory()
logger.info(f"Query shape: {adata.shape} "
            f"({adata.shape[0]:,} cells x {adata.shape[1]:,} genes)")

# ── Assign L1 cell type labels ────────────────────────────────────────────────────────

logger.info("Assigning L1 cell type labels from leiden_harmony_0.2...")
adata.obs['cell_type'] = adata.obs['leiden_harmony_0.2'].map(cluster_anns)

# Report any unmapped clusters
n_unmapped = adata.obs['cell_type'].isna().sum()
if n_unmapped > 0:
    unmapped = adata.obs.loc[adata.obs['cell_type'].isna(), 
                              'leiden_harmony_0.2'].unique().tolist()
    logger.warning(f"WARNING: {n_unmapped:,} cells have unmapped cluster IDs: {unmapped}")
else:
    logger.info("All cells successfully mapped to L1 labels")

cell_type_dist = adata.obs['cell_type'].value_counts().sort_index()
logger.info(f"\nL1 cell type distribution:\n{cell_type_dist.to_string()}")

# ── Set X to raw counts ───────────────────────────────────────────────────────────────

logger.info("Setting adata.X to raw counts from adata.layers['counts']...")
adata.X = adata.layers['counts'].copy()

# Sanity check
import scipy.sparse as sp
x_sample = adata.X[0, :]
if sp.issparse(x_sample):
    x_dense = np.array(x_sample.todense()).flatten()
else:
    x_dense = np.array(x_sample).flatten()

all_whole = np.all(x_dense == np.floor(x_dense))
x_max = adata.X.max()
logger.info(f"X dtype: {adata.X.dtype}")
logger.info(f"X max value: {x_max:.1f}")
logger.info(f"X values are whole numbers: {all_whole}")

if not all_whole:
    logger.warning("WARNING: X values are not all whole numbers - "
                   "check that layers['counts'] contains raw counts")

# ── Load Qian reference to get MERFISH gene panel ─────────────────────────────────────

logger.info(f"\nLoading Qian reference to extract MERFISH gene panel: {input_ref_h5ad}")
qian = ad.read_h5ad(input_ref_h5ad, backed='r')
merfish_genes = qian.var.index.tolist()
logger.info(f"MERFISH panel size: {len(merfish_genes)} genes")

# ── Gene overlap ──────────────────────────────────────────────────────────────────────

query_genes  = set(adata.var.index.tolist())
merfish_set  = set(merfish_genes)

shared_genes  = query_genes & merfish_set
missing_genes = merfish_set - query_genes

logger.info(f"\n--- Gene overlap report ---")
logger.info(f"Genes in MERFISH panel:          {len(merfish_set)}")
logger.info(f"Genes in snRNA-seq query:        {len(query_genes):,}")
logger.info(f"Shared genes:                    {len(shared_genes)} "
            f"({100 * len(shared_genes) / len(merfish_set):.1f}% of MERFISH panel)")
logger.info(f"MERFISH genes absent from query: {len(missing_genes)}")
logger.info(f"Absent genes (likely unexpressed in snRNA-seq): {sorted(missing_genes)}")

if len(shared_genes) < 200:
    logger.warning(f"WARNING: Only {len(shared_genes)} shared genes — "
                   f"label transfer quality may be reduced")

# ── Subset query to shared genes ──────────────────────────────────────────────────────

logger.info(f"\nSubsetting query to {len(shared_genes)} shared genes...")
shared_genes_ordered = [g for g in merfish_genes if g in query_genes]
adata = adata[:, shared_genes_ordered].copy()
logger.info(f"Query shape after subsetting: {adata.shape}")

# ── Export UMAP coordinates and metadata ──────────────────────────────────────────────

logger.info("Extracting UMAP coordinates and cell metadata...")

# Read UMAP from original h5ad - need to reload as we subsetted adata
adata_full = ad.read_h5ad(input_query_h5ad, backed='r')

umap_df = pd.DataFrame(
    adata_full.obsm['X_umap'],
    index   = adata_full.obs_names,
    columns = ['UMAP1', 'UMAP2']
)

# Add cell metadata
umap_df['cell_type']          = adata_full.obs['leiden_harmony_0.2'].map(cluster_anns)
umap_df['leiden_harmony_0.2'] = adata_full.obs['leiden_harmony_0.2'].values
umap_df['sample']             = adata_full.obs['sample'].values
umap_df['PCW']                = adata_full.obs['PCW'].values

umap_df.index.name = 'cell_id'

logger.info(f"UMAP dataframe shape: {umap_df.shape}")
logger.info(f"Cell type distribution in UMAP export:\n{umap_df['cell_type'].value_counts().sort_index().to_string()}")

logger.info(f"Writing UMAP coordinates to: {output_umap_csv}")
umap_df.to_csv(output_umap_csv)
logger.info("UMAP export done.")

# ── Write ─────────────────────────────────────────────────────────────────────────────

logger.info(f"\nWriting query to: {output_h5ad}")
adata.write_h5ad(output_h5ad)
logger.info("Done.")
