# --------------------------------------------------------------------------------------
#
#    review_kbet.py
#
#    Assesses Harmony batch correction quality using kBET (k-nearest neighbour
#    Batch Effect Test) on subsampled cells from the snRNA-seq object.
#
#    Runs kBET on both pre- and post-Harmony PCA embeddings and compares
#    rejection rates. A drop in rejection rate post-Harmony indicates effective
#    batch mixing. Results are stratified by cell type.
#
#    Inputs:
#        - adata_clusters_v3.h5ad  (contains X_pca_pre_harmony, X_pca_harmony,
#                                   obs['sample'], obs['leiden_harmony_0.2'])
#
#    Outputs:
#        - CSV of per-cell-type and overall kBET rejection rates (pre/post)
#        - Bar plot comparing rejection rates
#
#    Usage: called via Snakemake rule batch_correction_kbet
#
# --------------------------------------------------------------------------------------

import sys
import logging
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sklearn.neighbors import NearestNeighbors
from scipy.stats import chi2
from scipy.sparse import issparse

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

input_h5ad    = str(snakemake.input.h5ad)
output_csv    = str(snakemake.output.csv)
output_plot   = str(snakemake.output.plot)

batch_col     = snakemake.params.batch_col        # 'sample'
cell_type_col = snakemake.params.cell_type_col    # 'cell_type'
n_subsample   = int(snakemake.params.n_subsample) # default 10000
random_seed   = int(snakemake.params.random_seed) # default 1234
k0            = int(snakemake.params.k0)          # neighbourhood size, default 30

rng = np.random.default_rng(random_seed)

# ── kBET implementation ───────────────────────────────────────────────────────────────

def run_kbet(embedding, batch_labels, k=30):
    """
    Pure-Python kBET implementation.

    For each cell in the embedding, finds its k nearest neighbours and performs
    a chi-squared goodness-of-fit test comparing the observed batch label
    distribution in that neighbourhood to the global batch label distribution.

    A high rejection rate indicates poor batch mixing (batches are spatially
    clustered). A low rejection rate indicates good mixing.

    Parameters
    ----------
    embedding    : np.ndarray, shape (n_cells, n_pcs)
    batch_labels : np.ndarray of str/int, shape (n_cells,)
    k            : int, neighbourhood size

    Returns
    -------
    rejection_rate : float, proportion of cells where H0 (random mixing) is rejected
    accept_rate    : float, 1 - rejection_rate
    test_stats     : np.ndarray of per-cell chi-squared statistics
    p_values       : np.ndarray of per-cell p-values
    """
    batches        = np.unique(batch_labels)
    n_batches      = len(batches)
    n_cells        = embedding.shape[0]
    batch_to_idx   = {b: i for i, b in enumerate(batches)}
    batch_idx      = np.array([batch_to_idx[b] for b in batch_labels])

    # Global batch frequency (expected distribution)
    global_counts  = np.bincount(batch_idx, minlength=n_batches)
    global_freq    = global_counts / global_counts.sum()
    expected       = global_freq * k  # expected counts per neighbourhood

    # Build kNN graph
    logger.info(f"    Building kNN graph (k={k}, n_cells={n_cells:,})...")
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree', n_jobs=-1)
    nn.fit(embedding)
    indices = nn.kneighbors(embedding, return_distance=False)
    # Exclude self (first neighbour)
    indices = indices[:, 1:]

    # Chi-squared test per cell
    test_stats = np.zeros(n_cells)
    p_values   = np.zeros(n_cells)

    for i in range(n_cells):
        neighbour_batches  = batch_idx[indices[i]]
        observed           = np.bincount(neighbour_batches, minlength=n_batches).astype(float)

        # Avoid division by zero for batches with zero expected count
        mask               = expected > 0
        chi2_stat          = np.sum(
            (observed[mask] - expected[mask]) ** 2 / expected[mask]
        )
        df                 = mask.sum() - 1
        p_val              = 1 - chi2.cdf(chi2_stat, df=max(df, 1))

        test_stats[i]      = chi2_stat
        p_values[i]        = p_val

    rejection_rate = np.mean(p_values < 0.05)
    accept_rate    = 1 - rejection_rate

    return rejection_rate, accept_rate, test_stats, p_values

# ── Load data ─────────────────────────────────────────────────────────────────────────

logger.info(f"Loading h5ad from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
logger.info(f"Object shape: {adata.shape} ({adata.shape[0]:,} cells x {adata.shape[1]:,} genes)")

# ── Assign cell type labels ───────────────────────────────────────────────────────────

logger.info("Assigning L1 cell type labels from leiden_harmony_0.2...")
adata.obs[cell_type_col] = adata.obs['leiden_harmony_0.2'].map(cluster_anns)

n_unmapped = adata.obs[cell_type_col].isna().sum()
if n_unmapped > 0:
    logger.warning(f"WARNING: {n_unmapped:,} cells could not be mapped to a cell type label")
else:
    logger.info("All cells successfully mapped to L1 labels")

logger.info(f"Cell type distribution:\n{adata.obs[cell_type_col].value_counts().sort_index().to_string()}")

# ── Validate embeddings ───────────────────────────────────────────────────────────────

for key in ['X_pca_pre_harmony', 'X_pca_harmony']:
    if key not in adata.obsm:
        raise KeyError(f"Required embedding '{key}' not found in adata.obsm. "
                       f"Available keys: {list(adata.obsm.keys())}")
logger.info("Both X_pca_pre_harmony and X_pca_harmony found in adata.obsm")

# ── Stratified subsample ──────────────────────────────────────────────────────────────

logger.info(f"Stratified subsampling {n_subsample:,} cells (seed={random_seed})...")

cell_types   = adata.obs[cell_type_col].dropna().unique()
n_per_type   = max(1, n_subsample // len(cell_types))
sampled_idx  = []

for ct in cell_types:
    ct_idx   = np.where(adata.obs[cell_type_col] == ct)[0]
    n_sample = min(n_per_type, len(ct_idx))
    chosen   = rng.choice(ct_idx, size=n_sample, replace=False)
    sampled_idx.extend(chosen.tolist())
    logger.info(f"  {ct}: sampled {n_sample:,} / {len(ct_idx):,} cells")

sampled_idx = np.array(sampled_idx)
logger.info(f"Total subsampled cells: {len(sampled_idx):,}")

emb_pre  = adata.obsm['X_pca_pre_harmony'][sampled_idx]
emb_post = adata.obsm['X_pca_harmony'][sampled_idx]
batches  = adata.obs[batch_col].values[sampled_idx]
ctypes   = adata.obs[cell_type_col].values[sampled_idx]

# ── Run kBET overall ─────────────────────────────────────────────────────────────────

logger.info("Running kBET on full subsample (pre-Harmony)...")
rr_pre, ar_pre, stats_pre, pvals_pre = run_kbet(emb_pre, batches, k=k0)
logger.info(f"  Pre-Harmony  rejection rate: {rr_pre:.4f} | acceptance rate: {ar_pre:.4f}")

logger.info("Running kBET on full subsample (post-Harmony)...")
rr_post, ar_post, stats_post, pvals_post = run_kbet(emb_post, batches, k=k0)
logger.info(f"  Post-Harmony rejection rate: {rr_post:.4f} | acceptance rate: {ar_post:.4f}")

# ── Run kBET per cell type ────────────────────────────────────────────────────────────

logger.info("Running kBET per cell type...")
records = []

for ct in sorted(cell_types):
    ct_mask = ctypes == ct
    n_ct    = ct_mask.sum()

    if n_ct < k0 + 1:
        logger.warning(f"  Skipping {ct}: only {n_ct} cells after subsample (need > {k0})")
        continue

    ct_batches  = batches[ct_mask]
    n_unique    = len(np.unique(ct_batches))

    if n_unique < 2:
        logger.warning(f"  Skipping {ct}: only {n_unique} unique batch(es) in subsample")
        continue

    logger.info(f"  {ct} (n={n_ct:,}, batches={n_unique})...")

    rr_pre_ct,  _, _, _ = run_kbet(emb_pre[ct_mask],  ct_batches, k=min(k0, n_ct - 1))
    rr_post_ct, _, _, _ = run_kbet(emb_post[ct_mask], ct_batches, k=min(k0, n_ct - 1))

    logger.info(f"    Pre:  {rr_pre_ct:.4f} | Post: {rr_post_ct:.4f}")

    records.append({
        'cell_type'          : ct,
        'n_cells_subsampled' : int(n_ct),
        'n_batches'          : int(n_unique),
        'rejection_rate_pre' : round(rr_pre_ct,  4),
        'rejection_rate_post': round(rr_post_ct, 4),
        'delta_rejection'    : round(rr_pre_ct - rr_post_ct, 4)
    })

# Add overall row
records.insert(0, {
    'cell_type'          : 'ALL',
    'n_cells_subsampled' : len(sampled_idx),
    'n_batches'          : len(np.unique(batches)),
    'rejection_rate_pre' : round(rr_pre,  4),
    'rejection_rate_post': round(rr_post, 4),
    'delta_rejection'    : round(rr_pre - rr_post, 4)
})

results_df = pd.DataFrame(records)
logger.info(f"\nkBET results:\n{results_df.to_string(index=False)}")

# ── Save CSV ──────────────────────────────────────────────────────────────────────────

results_df.to_csv(output_csv, index=False)
logger.info(f"Results saved to: {output_csv}")

# ── Plot ──────────────────────────────────────────────────────────────────────────────

logger.info("Generating plot...")

plot_df  = results_df.copy()
labels   = plot_df['cell_type'].tolist()
x        = np.arange(len(labels))
width    = 0.35

fig, ax  = plt.subplots(figsize=(max(8, len(labels) * 1.2), 5))

bars_pre  = ax.bar(x - width / 2, plot_df['rejection_rate_pre'],  width,
                   color='#d62728', alpha=0.85, label='Pre-Harmony')
bars_post = ax.bar(x + width / 2, plot_df['rejection_rate_post'], width,
                   color='#1f77b4', alpha=0.85, label='Post-Harmony')

ax.set_xlabel('Cell type', fontsize=12)
ax.set_ylabel('kBET rejection rate', fontsize=12)
ax.set_title(f'kBET batch correction assessment\n'
             f'(k={k0}, n_subsampled={len(sampled_idx):,}, seed={random_seed})',
             fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=10)
ax.set_ylim(0, 1.05)
ax.axhline(0.05, color='grey', linestyle='--', linewidth=0.8,
           label='α = 0.05 (random mixing threshold)')
ax.legend(fontsize=10)

# Annotate delta
for i, row in enumerate(records):
    delta = row['delta_rejection']
    ypos  = max(row['rejection_rate_pre'], row['rejection_rate_post']) + 0.02
    ax.text(i, ypos, f'Δ{delta:+.2f}', ha='center', va='bottom',
            fontsize=8, color='black')

plt.tight_layout()
fig.savefig(output_plot, dpi=300, bbox_inches='tight')
logger.info(f"Plot saved to: {output_plot}")
logger.info("kBET assessment complete.")
