# --------------------------------------------------------------------------------------
#
#    review_lisi.py
#
#    Assesses Harmony batch correction quality using LISI (Local Inverse
#    Simpson's Index) on both pre- and post-Harmony PCA embeddings.
#
#    Two metrics are computed:
#      - iLISI (integration LISI): measures batch mixing in the embedding.
#            Score near 1 = poor mixing (one batch dominates each neighbourhood).
#            Score near N_batches = good mixing.
#            Computed on BOTH pre- and post-Harmony embeddings for comparison.
#
#      - cLISI (cell-type LISI): measures cell type separation in the embedding.
#            Score near 1 = cell types are well separated (desired).
#            Score near N_cell_types = cell types are blended (over-correction).
#            Computed on post-Harmony embedding only.
#
#    kNN is built using pynndescent (approximate NN) for scalability at ~680K cells.
#    This is the same backend Scanpy uses for UMAP and gives equivalent accuracy
#    for LISI purposes.
#
#    Inputs:
#        - adata_clusters_v3.h5ad  (contains X_pca_pre_harmony, X_pca_harmony,
#                                   obs['sample'], obs['leiden_harmony_0.2'])
#
#    Outputs:
#        - CSV of per-cell iLISI (pre/post) and cLISI (post) scores
#        - CSV summary per cell type (_summary suffix, same directory)
#        - Violin plot comparing distributions
#
#    Usage: called via Snakemake rule batch_correction_lisi
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

from pynndescent import NNDescent
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
random_seed   = int(snakemake.params.random_seed) # default 1234
perplexity    = int(snakemake.params.perplexity)  # neighbourhood size, default 30
n_threads     = int(snakemake.resources.threads)

# ── LISI implementation ───────────────────────────────────────────────────────────────

def compute_lisi(embedding, labels, perplexity=30, random_seed=1234, n_threads=1,
                 label_name='labels'):
    """
    Compute the Local Inverse Simpson's Index (LISI) for each cell.

    Uses pynndescent for approximate kNN construction, which scales to hundreds
    of thousands of cells in high dimensions. The per-cell LISI computation is
    fully vectorised using numpy.

    For each cell, builds a perplexity-scaled Gaussian kernel over its k nearest
    neighbours to define a soft neighbourhood, then computes the Inverse
    Simpson's Index over the label distribution in that neighbourhood.

    Score near 1              = one label dominates (poor mixing / good separation)
    Score near N_unique_labels = labels evenly distributed (good mixing / poor separation)

    Parameters
    ----------
    embedding   : np.ndarray, shape (n_cells, n_pcs)
    labels      : np.ndarray of str/int, shape (n_cells,)
    perplexity  : int, controls neighbourhood size (analogous to tSNE perplexity)
    random_seed : int
    n_threads   : int, passed to NNDescent
    label_name  : str, used in log messages only

    Returns
    -------
    lisi_scores : np.ndarray of float, shape (n_cells,)
    """
    n_cells       = embedding.shape[0]
    unique_labels = np.unique(labels)
    label_to_idx  = {l: i for i, l in enumerate(unique_labels)}
    label_idx     = np.array([label_to_idx[l] for l in labels])
    n_labels      = len(unique_labels)

    # Standard LISI practice: k = 3 * perplexity
    k = min(3 * perplexity, n_cells - 1)

    logger.info(f"    [{label_name}] Building approximate kNN graph "
                f"(k={k}, perplexity={perplexity}, "
                f"n_cells={n_cells:,}, n_labels={n_labels}, "
                f"n_threads={n_threads})...")

    index = NNDescent(
        embedding,
        n_neighbors=k + 1,
        random_state=random_seed,
        n_jobs=n_threads,
        verbose=False
    )
    indices, distances = index.neighbor_graph

    # NNDescent includes self as first neighbour — drop it
    indices   = indices[:, 1:].astype(np.int32)
    distances = distances[:, 1:].astype(np.float32)

    logger.info(f"    [{label_name}] kNN graph built. Computing LISI scores...")

    # ── Vectorised per-cell LISI ──────────────────────────────────────────────────────
    lisi_scores = np.zeros(n_cells, dtype=np.float32)
    target_h    = np.log(float(perplexity))

    # Process in chunks to avoid peak memory blow-up at 680K x 90
    chunk_size = 50000
    for start in range(0, n_cells, chunk_size):
        end   = min(start + chunk_size, n_cells)
        d     = distances[start:end].astype(np.float64)   # (chunk, k)
        idx   = indices[start:end]                         # (chunk, k)
        chunk = end - start

        # Binary search for per-cell bandwidth beta = 1/(2*sigma^2)
        beta     = np.ones(chunk,  dtype=np.float64)
        beta_min = np.full(chunk, -np.inf, dtype=np.float64)
        beta_max = np.full(chunk,  np.inf, dtype=np.float64)

        for _ in range(50):
            p     = np.exp(-d * beta[:, None])               # (chunk, k)
            p_sum = p.sum(axis=1, keepdims=True)
            p_sum = np.where(p_sum == 0, 1e-10, p_sum)
            p     = p / p_sum
            h     = -np.sum(p * np.log(p + 1e-10), axis=1)  # (chunk,)

            h_diff    = h - target_h
            converged = np.abs(h_diff) < 1e-5
            increase  = (h_diff > 0) & ~converged
            decrease  = (h_diff <= 0) & ~converged

            beta_min = np.where(increase, beta, beta_min)
            beta_max = np.where(decrease, beta, beta_max)
            beta     = np.where(
                increase,
                np.where(beta_max == np.inf,  beta * 2,       (beta + beta_max) / 2),
                np.where(
                    decrease,
                    np.where(beta_min == -np.inf, beta / 2,   (beta + beta_min) / 2),
                    beta
                )
            )
            if converged.all():
                break

        # Soft label distribution: weight neighbour labels by kernel probability p
        # idx shape: (chunk, k),  label_idx[idx] shape: (chunk, k)
        neighbour_labels = label_idx[idx]                    # (chunk, k)
        one_hot = (neighbour_labels[:, :, None] ==
                   np.arange(n_labels)[None, None, :])        # (chunk, k, n_labels)
        label_probs = (p[:, :, None] * one_hot).sum(axis=1)  # (chunk, n_labels)

        # Inverse Simpson's Index
        simpson = (label_probs ** 2).sum(axis=1)              # (chunk,)
        lisi_scores[start:end] = np.where(simpson > 0, 1.0 / simpson, 1.0)

    logger.info(f"    [{label_name}] Done. "
                f"mean={lisi_scores.mean():.3f}  median={np.median(lisi_scores):.3f}")

    return lisi_scores

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

logger.info(f"Cell type distribution:\n"
            f"{adata.obs[cell_type_col].value_counts().sort_index().to_string()}")

# ── Validate embeddings ───────────────────────────────────────────────────────────────

for key in ['X_pca_pre_harmony', 'X_pca_harmony']:
    if key not in adata.obsm:
        raise KeyError(f"Required embedding '{key}' not found in adata.obsm. "
                       f"Available keys: {list(adata.obsm.keys())}")
logger.info("Both X_pca_pre_harmony and X_pca_harmony found in adata.obsm")

# ── Extract arrays ────────────────────────────────────────────────────────────────────

emb_pre  = np.array(adata.obsm['X_pca_pre_harmony'], dtype=np.float32)
emb_post = np.array(adata.obsm['X_pca_harmony'],     dtype=np.float32)
batches  = adata.obs[batch_col].values
ctypes   = adata.obs[cell_type_col].values

# Drop cells with unmapped cell types
valid_mask = ~pd.isnull(ctypes)
if valid_mask.sum() < len(ctypes):
    logger.warning(f"Dropping {(~valid_mask).sum():,} cells with unmapped cell types")
    emb_pre  = emb_pre[valid_mask]
    emb_post = emb_post[valid_mask]
    batches  = batches[valid_mask]
    ctypes   = ctypes[valid_mask]

n_cells = emb_pre.shape[0]
logger.info(f"Running LISI on {n_cells:,} cells with {n_threads} threads")

# ── Compute iLISI pre-Harmony ─────────────────────────────────────────────────────────

logger.info("Computing iLISI on pre-Harmony embedding...")
ilisi_pre = compute_lisi(emb_pre, batches, perplexity=perplexity,
                         random_seed=random_seed, n_threads=n_threads,
                         label_name='iLISI-pre')
logger.info(f"  iLISI pre-Harmony:  mean={ilisi_pre.mean():.3f}  "
            f"median={np.median(ilisi_pre):.3f}")

# ── Compute iLISI post-Harmony ────────────────────────────────────────────────────────

logger.info("Computing iLISI on post-Harmony embedding...")
ilisi_post = compute_lisi(emb_post, batches, perplexity=perplexity,
                          random_seed=random_seed, n_threads=n_threads,
                          label_name='iLISI-post')
logger.info(f"  iLISI post-Harmony: mean={ilisi_post.mean():.3f}  "
            f"median={np.median(ilisi_post):.3f}")

# ── Compute cLISI post-Harmony ────────────────────────────────────────────────────────

logger.info("Computing cLISI on post-Harmony embedding...")
clisi_post = compute_lisi(emb_post, ctypes, perplexity=perplexity,
                          random_seed=random_seed, n_threads=n_threads,
                          label_name='cLISI-post')
logger.info(f"  cLISI post-Harmony: mean={clisi_post.mean():.3f}  "
            f"median={np.median(clisi_post):.3f}  "
            f"(ideal ~1.0; n_cell_types={len(np.unique(ctypes))})")

# ── Assemble per-cell results ─────────────────────────────────────────────────────────

results_df = pd.DataFrame({
    'cell_id'     : adata.obs_names[valid_mask],
    batch_col     : batches,
    cell_type_col : ctypes,
    'ilisi_pre'   : ilisi_pre,
    'ilisi_post'  : ilisi_post,
    'clisi_post'  : clisi_post
})

# ── Per-cell-type summary ─────────────────────────────────────────────────────────────

logger.info("\nPer-cell-type LISI summary:")
summary_records = []
for ct in sorted(np.unique(ctypes)):
    mask = results_df[cell_type_col] == ct
    rec  = {
        'cell_type'         : ct,
        'n_cells'           : int(mask.sum()),
        'ilisi_pre_mean'    : round(float(results_df.loc[mask, 'ilisi_pre'].mean()),    3),
        'ilisi_pre_median'  : round(float(results_df.loc[mask, 'ilisi_pre'].median()),  3),
        'ilisi_post_mean'   : round(float(results_df.loc[mask, 'ilisi_post'].mean()),   3),
        'ilisi_post_median' : round(float(results_df.loc[mask, 'ilisi_post'].median()), 3),
        'clisi_post_mean'   : round(float(results_df.loc[mask, 'clisi_post'].mean()),   3),
        'clisi_post_median' : round(float(results_df.loc[mask, 'clisi_post'].median()), 3),
    }
    summary_records.append(rec)
    logger.info(f"  {ct} (n={rec['n_cells']:,}): "
                f"iLISI pre={rec['ilisi_pre_mean']:.3f} "
                f"-> post={rec['ilisi_post_mean']:.3f} | "
                f"cLISI post={rec['clisi_post_mean']:.3f}")

summary_df = pd.DataFrame(summary_records)

# ── Save CSVs ─────────────────────────────────────────────────────────────────────────

results_df.to_csv(output_csv, index=False)
logger.info(f"Per-cell results saved to: {output_csv}")

summary_path = output_csv.replace('.csv', '_summary.csv')
summary_df.to_csv(summary_path, index=False)
logger.info(f"Per-cell-type summary saved to: {summary_path}")

# ── Plot ──────────────────────────────────────────────────────────────────────────────

logger.info("Generating violin plot...")

cell_types_ordered = sorted(np.unique(ctypes))
colours = {
    'iLISI pre'  : '#d62728',
    'iLISI post' : '#1f77b4',
    'cLISI post' : '#2ca02c'
}

fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

def violin_by_celltype(ax, data_dict, title, ylabel, cell_types, df, ct_col):
    """Draw side-by-side violins per cell type, one series per key in data_dict."""
    labels_seen   = set()
    label_handles = []
    n_series = len(data_dict)
    spacing  = 1.0
    width    = 0.25

    for ci, ct in enumerate(cell_types):
        mask = df[ct_col] == ct
        base = ci * (n_series * width + spacing)
        for si, (label, col_name) in enumerate(data_dict.items()):
            pos    = base + si * width
            values = df.loc[mask, col_name].values
            vp     = ax.violinplot(values, positions=[pos], widths=width * 0.9,
                                   showmedians=True, showextrema=False)
            for pc in vp['bodies']:
                pc.set_facecolor(colours[label])
                pc.set_alpha(0.75)
            vp['cmedians'].set_color('black')
            vp['cmedians'].set_linewidth(1.5)
            if label not in labels_seen:
                label_handles.append(
                    mpatches.Patch(color=colours[label], alpha=0.75, label=label)
                )
                labels_seen.add(label)

    tick_positions = [
        ci * (n_series * width + spacing) + (n_series - 1) * width / 2
        for ci in range(len(cell_types))
    ]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(cell_types, rotation=35, ha='right', fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(handles=label_handles, fontsize=9)

# Panel 1: iLISI pre vs post
violin_by_celltype(
    axes[0],
    {'iLISI pre': 'ilisi_pre', 'iLISI post': 'ilisi_post'},
    'iLISI: batch mixing\n(higher = better mixing)',
    'iLISI score',
    cell_types_ordered, results_df, cell_type_col
)

# Panel 2: cLISI post
violin_by_celltype(
    axes[1],
    {'cLISI post': 'clisi_post'},
    'cLISI: cell type separation post-Harmony\n(lower = better separation)',
    'cLISI score',
    cell_types_ordered, results_df, cell_type_col
)

# Panel 3: Overall median summary bars
ax3     = axes[2]
metrics = ['ilisi_pre', 'ilisi_post', 'clisi_post']
labels3 = ['iLISI\npre-Harmony', 'iLISI\npost-Harmony', 'cLISI\npost-Harmony']
cols3   = [colours['iLISI pre'], colours['iLISI post'], colours['cLISI post']]
medians = [float(np.median(results_df[m])) for m in metrics]
bars    = ax3.bar(labels3, medians, color=cols3, alpha=0.85,
                  edgecolor='black', linewidth=0.7)
ax3.set_title('Overall median LISI scores', fontsize=11)
ax3.set_ylabel('Median LISI score', fontsize=10)
for bar, val in zip(bars, medians):
    ax3.text(bar.get_x() + bar.get_width() / 2, val + 0.01, f'{val:.3f}',
             ha='center', va='bottom', fontsize=10)

n_batches_total = len(np.unique(batches))
n_ct_total      = len(np.unique(ctypes))
ax3.set_ylim(0, max(medians) * 1.25)
ax3.text(0.5, 0.97,
         f'N cells = {n_cells:,} | N batches = {n_batches_total} | '
         f'N cell types = {n_ct_total}\nperplexity = {perplexity}',
         transform=ax3.transAxes, ha='center', va='top', fontsize=8, color='grey')

fig.suptitle('LISI batch correction assessment (pre- and post-Harmony)',
             fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(output_plot, dpi=300, bbox_inches='tight')
logger.info(f"Plot saved to: {output_plot}")
logger.info("LISI assessment complete.")
