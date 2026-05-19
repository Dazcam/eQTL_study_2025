# --------------------------------------------------------------------------------------
#
#    prepare_qian_reference.py
#
#    Preprocesses the Qian et al. 2025 MERFISH dataset for use as a SingleR reference.
#    - Excludes V1/V2 cortical areas (transcriptionally discontinuous per Qian et al.)
#    - Retains all gestational weeks (GW15, GW20, GW22, GW34)
#    - Subsamples 500 cells per H2 cluster for computational tractability
#    - Writes a subsampled h5ad for use as SingleR reference in R
#
#    Usage: called via Snakemake rule prepare_qian_reference
#
# --------------------------------------------------------------------------------------

import sys
import logging
import numpy as np
import anndata as ad
import scipy.sparse as sp

# ── Logging ───────────────────────────────────────────────────────────────────────────

log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger()

# Also log to stdout for cluster job output
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# ── Parameters ────────────────────────────────────────────────────────────────────────

input_h5ad   = str(snakemake.input)
output_h5ad  = str(snakemake.output)
n_cells      = int(snakemake.params.n_cells_per_cluster)
random_seed  = int(snakemake.params.random_seed)

AREA_EXCLUDE = ['A-V1', 'B-V1', 'A-V2', 'B-V2', 'C-V2']

# ── Load ──────────────────────────────────────────────────────────────────────────────

logger.info(f"Loading Qian MERFISH data in backed mode from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad, backed='r')
adata.obs_names_make_unique()
logger.info(f"Full dataset shape: {adata.shape} "
            f"({adata.shape[0]:,} cells x {adata.shape[1]:,} genes)")

# ── Filter ────────────────────────────────────────────────────────────────────────────

# Retain all GWs - H2 cluster labels were assigned across all timepoints jointly
gw_counts = adata.obs['gw'].value_counts().sort_index()
logger.info(f"Gestational week distribution (all retained):\n{gw_counts.to_string()}")

# Exclude V1/V2 - transcriptionally discontinuous per Qian et al.
area_mask = ~adata.obs['area'].isin(AREA_EXCLUDE) | adata.obs['area'].isna()
n_excluded = (~area_mask).sum()
logger.info(f"Excluding V1/V2 areas: {AREA_EXCLUDE}")
logger.info(f"Cells excluded (V1/V2): {n_excluded:,}")
logger.info(f"Cells retained after area filter: {area_mask.sum():,}")

filtered_obs = adata.obs[area_mask].copy()

# ── Subsample ─────────────────────────────────────────────────────────────────────────

logger.info(f"Subsampling up to {n_cells} cells per H2 cluster...")
np.random.seed(random_seed)

sampled_indices = []
cluster_summary = {}

for cluster in sorted(filtered_obs['H2_annotation'].unique()):
    cluster_idx = filtered_obs.index[filtered_obs['H2_annotation'] == cluster]
    n = min(n_cells, len(cluster_idx))
    sampled = np.random.choice(cluster_idx, size=n, replace=False)
    sampled_indices.extend(sampled)
    cluster_summary[cluster] = {'available': len(cluster_idx), 'sampled': n}
    logger.info(f"  {cluster}: {len(cluster_idx):,} available, {n} sampled")

logger.info(f"\nTotal reference cells after subsampling: {len(sampled_indices):,}")

# ── Load subsampled matrix into memory ────────────────────────────────────────────────

logger.info("Loading subsampled matrix into memory...")

# Get integer positions for sampled indices
# adata.raw.X requires positional integer indexing
sampled_positions = adata.obs_names.get_indexer(sampled_indices)

logger.info("Extracting raw counts from adata.raw...")

raw_X    = adata.raw.X[sampled_positions, :]
raw_var  = adata.raw.var
obs_meta = adata.obs.iloc[sampled_positions].copy()

# Build a clean AnnData from scratch using raw counts
adata_ref = ad.AnnData(
    X   = sp.csr_matrix(raw_X.astype('float32')),
    obs = obs_meta,
    var = raw_var
)

logger.info(f"Reference shape: {adata_ref.shape}")

# Sanity check
x_max    = adata_ref.X.max()
all_whole = all(v == int(v) for v in adata_ref.X.data[:100])
logger.info(f"Max value in X: {x_max:.4f}")
logger.info(f"Values are whole numbers (first 100 non-zero): {all_whole}")
if not all_whole:
    logger.warning("WARNING: X values are not whole numbers - check adata.raw contains raw counts")


# ── GW distribution in final reference ───────────────────────────────────────────────

gw_dist = adata_ref.obs['gw'].value_counts().sort_index()
logger.info(f"\nGW distribution in final reference:\n{gw_dist.to_string()}")

h2_dist = adata_ref.obs['H2_annotation'].value_counts().sort_index()
logger.info(f"\nH2 annotation distribution in final reference:\n{h2_dist.to_string()}")

# ── Write ─────────────────────────────────────────────────────────────────────────────

logger.info(f"\nWriting reference to: {output_h5ad}")
adata_ref.write_h5ad(output_h5ad)
logger.info("Done.")
