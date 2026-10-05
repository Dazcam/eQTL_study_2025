# --------------------------------------------------------------------------------------
#
#    manuscript_plot_snrna_qc.py
#
#    Reviewer request (line 473): snRNA-seq QC summaries.
#
#    Figure 1 (per sample), rows top -> bottom:
#      - Number of nuclei per sample, stacked by cell type (bar)   [A: nuclei per sample]
#      - UMIs per nucleus (boxplot per sample)                     [B: UMI distributions]
#      - Genes detected per nucleus (boxplot per sample)           [C: genes per sample]
#
#    Figure 2 (per cell type):
#      - Genes detected per nucleus (boxplot per cell type)        [D: genes per cell type]
#
#    Table (Excel, not tracked by Snakemake): nuclei, median UMIs and median genes per
#    nucleus, per sample and per cell type.
#
#    All metrics are recomputed from layers['counts'] of the final object (adata.X holds
#    scaled log-normalised values, and obs['n_counts'/'n_genes'] may pre-date later QC).
#    Only obs and layers['counts'] are read from the h5ad, so the dense X is never loaded.
#
#    Usage: called via Snakemake rule snrna_qc_plt
#
#    Layout: everything between the EXPORTABLE BLOCK markers has no Snakemake dependency
#    and takes plain pandas/numpy/scipy objects, so it can be moved into scanpy_utils.py
#    unchanged (needs numpy, pandas, scipy.sparse as sp, matplotlib.pyplot as plt, seaborn
#    as sns, plus the CELLTYPE_ORDER constant).
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
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

logger = logging.getLogger("snrna_qc")

# Hard genes-per-nucleus QC thresholds (fixed; only used to log a sanity check)
MIN_GENES_PER_CELL = 1000
MAX_GENES_PER_CELL = 5000

# Summary tables are written here by hand, not tracked by Snakemake
TABLE_DIR = "../results/13MANUSCRIPT_TABLES/"
TABLE_FILE = "snrna_qc_summary.xlsx"

# Glutamatergic populations were renamed after annotation; colours are retained
CELLTYPE_RENAME = {"Glu-UL": "Glu-A", "Glu-DL": "Glu-B"}

# Display order of cell types (legend, x-axis, table)
CELLTYPE_ORDER = ["Glu-A", "Glu-B", "GABA", "NPC", "OPC", "MG", "Endo-Peri"]


# ======================================================================================
#  EXPORTABLE BLOCK (start)
# ======================================================================================

# ── Loading and metric calculation ────────────────────────────────────────────────────

def load_obs_and_counts(h5ad_path, layer="counts"):
    """
    Read only `obs` and one counts layer from an h5ad (skips the dense scaled `X`).

    Returns
    -------
    obs : pd.DataFrame
    counts : scipy.sparse.csr_matrix   (nuclei x genes, raw counts)
    var_names : pd.Index
    """
    import h5py
    try:
        from anndata.io import read_elem            # anndata >= 0.11
    except ImportError:
        from anndata.experimental import read_elem  # older anndata

    with h5py.File(h5ad_path, "r") as f:
        if "layers" not in f or layer not in f["layers"]:
            raise KeyError(f"Layer '{layer}' not found in {h5ad_path}")
        obs = read_elem(f["obs"])
        var_names = read_elem(f["var"]).index
        counts = read_elem(f["layers"][layer])

    counts = sp.csr_matrix(counts)
    counts.eliminate_zeros()  # so nnz per row == genes detected

    # Sanity check: these must be integer counts
    head = counts.data[:100000]
    if not np.all(np.mod(head, 1) == 0):
        raise ValueError(f"layers['{layer}'] does not look like integer counts.")
    return obs, counts, var_names


def compute_cell_qc(counts):
    """UMIs and detected genes per nucleus from a raw-count csr matrix."""
    n_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.diff(counts.indptr)
    return n_counts, n_genes


# ── Summary tables ────────────────────────────────────────────────────────────────────

def _order_celltypes(present, order=None):
    """Cell types in `order` (default CELLTYPE_ORDER), any others sorted after."""
    order = CELLTYPE_ORDER if order is None else order
    present = set(present)
    first = [c for c in order if c in present]
    return first + sorted(present - set(first))


def build_qc_tables(qc, sample_col="sample", celltype_col="cell_type"):
    """
    Per-sample and per-cell-type summaries of the numbers behind the QC figures.
    Medians are taken over nuclei.

    Returns {"per_sample": DataFrame, "per_cell_type": DataFrame}
    """
    def summarise(df, col, order, label):
        g = df.groupby(col, observed=True)
        tbl = pd.DataFrame({
            "n_nuclei": g.size(),
            "median_umis_per_nucleus": g["n_counts"].median(),
            "median_genes_per_nucleus": g["n_genes"].median(),
        }).reindex(order)
        tbl.index.name = label
        return tbl.reset_index()

    qc = qc.assign(**{sample_col: qc[sample_col].astype(str)})
    ct_df = qc.dropna(subset=[celltype_col])
    ct_df = ct_df.assign(**{celltype_col: ct_df[celltype_col].astype(str)})
    return {
        "per_sample": summarise(qc, sample_col, sorted(qc[sample_col].unique()), "sample"),
        "per_cell_type": summarise(ct_df, celltype_col,
                                   _order_celltypes(ct_df[celltype_col].unique()), "cell_type"),
    }


# ── Plot helpers ──────────────────────────────────────────────────────────────────────

def _resolve_palette(categories, palette=None):
    """Return {category: colour}. `palette` may be None, a dict, or a list (in category order)."""
    if isinstance(palette, dict):
        fallback = plt.cm.tab10.colors
        return {c: palette.get(c, fallback[i % len(fallback)]) for i, c in enumerate(categories)}
    if palette is None or len(palette) == 0:
        palette = plt.cm.tab10.colors
    return {c: palette[i % len(palette)] for i, c in enumerate(categories)}


def _black_boxplot(ax, df, x, y, order, linewidth=0.8, fliersize=1, whis=1.5):
    """
    Boxplot in the style of plot_stacked_figure (white boxes, small fliers) but with
    black lines. Seaborn's own fliers are switched off and the outliers (beyond
    whis x IQR from the quartiles, i.e. the same rule) are drawn as ONE rasterised
    scatter, because one rasterised artist per group makes the PDF writer run out of
    memory with ~1e5 outlier points.
    """
    sns.boxplot(data=df, x=x, y=y, order=order, ax=ax, color="white",
                showfliers=False, whis=whis, linewidth=linewidth)
    for patch in ax.patches:
        patch.set_facecolor("white")
        patch.set_edgecolor("black")
    for line in ax.lines:
        line.set_color("black")

    pos = {g: i for i, g in enumerate(order)}
    q = df.groupby(x, observed=True)[y].quantile([0.25, 0.75]).unstack()
    iqr = q[0.75] - q[0.25]
    lo = df[x].map(q[0.25] - whis * iqr)
    hi = df[x].map(q[0.75] + whis * iqr)
    out = df[(df[y] < lo) | (df[y] > hi)]
    ax.scatter(out[x].map(pos).values, out[y].values, s=fliersize, c="black",
               linewidths=0, rasterized=True, zorder=3)


def _median_label(ax, text, fontsize=10):
    """Median annotation above the axes (left-aligned), so it never overlaps the boxes."""
    ax.text(0.0, 1.015, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=fontsize, color="black")


def _style_axis(ax, ylabel):
    ax.set_ylabel(ylabel, fontsize=10, color="black")
    ax.yaxis.grid(True, color="lightgrey", linewidth=0.8)
    ax.set_axisbelow(True)


def _set_group_ticks(ax, labels, rotation=90, fontsize=6):
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=rotation, fontsize=fontsize)
    ax.set_xlim(-0.5, len(labels) - 0.5)


# ── Figure 1: per sample ──────────────────────────────────────────────────────────────

def plot_qc_per_sample(qc, sample_col="sample", celltype_col="cell_type",
                       palette=None, figsize=(24, 11), stack_by_celltype=True):
    """
    Per-sample QC figure: nuclei per sample (top), UMIs per nucleus, genes per nucleus.

    Parameters
    ----------
    qc : pd.DataFrame with one row per nucleus and columns
         [sample_col, celltype_col, 'n_counts', 'n_genes']
    stack_by_celltype : colour the nuclei-count bars by cell type (False = single grey bar)

    Returns
    -------
    fig : matplotlib Figure
    """
    qc = qc.assign(**{sample_col: qc[sample_col].astype(str)})
    sample_order = sorted(qc[sample_col].unique())          # sorted by sample label
    n = len(sample_order)

    fig, axs = plt.subplots(3, 1, figsize=figsize, sharex=True,
                            gridspec_kw={"height_ratios": [1.3, 1, 1], "hspace": 0.18})

    # Row 1: number of nuclei per sample
    counts_tbl = (qc.groupby([sample_col, celltype_col], observed=True).size()
                    .unstack(fill_value=0).reindex(sample_order).fillna(0))
    if stack_by_celltype:
        cts = _order_celltypes(counts_tbl.columns)
        colours = _resolve_palette(cts, palette)
        bottom = np.zeros(n)
        for ct in cts:
            vals = counts_tbl[ct].values
            axs[0].bar(range(n), vals, bottom=bottom, width=1, color=colours[ct],
                       edgecolor="white", linewidth=0.3, label=ct)
            bottom += vals
        axs[0].legend(bbox_to_anchor=(1.005, 1), loc="upper left", title="Cell type",
                      frameon=False, fontsize=9, title_fontsize=10)
    else:
        axs[0].bar(range(n), counts_tbl.sum(axis=1).values, width=0.8,
                   color="lightgrey", edgecolor="black", linewidth=0.5)
    _style_axis(axs[0], "Number of nuclei")
    _median_label(axs[0], f"Median: {int(counts_tbl.sum(axis=1).median()):,} nuclei per sample")

    # Row 2: UMIs per nucleus
    _black_boxplot(axs[1], qc, sample_col, "n_counts", sample_order)
    _style_axis(axs[1], "Number of UMIs per nucleus")
    _median_label(axs[1], f"Median: {int(np.median(qc['n_counts'])):,} UMIs")

    # Row 3: genes per nucleus
    _black_boxplot(axs[2], qc, sample_col, "n_genes", sample_order)
    _style_axis(axs[2], "Number of genes per nucleus")
    _median_label(axs[2], f"Median: {int(np.median(qc['n_genes'])):,} genes")

    # x labels only on the bottom panel, all samples shown
    for ax in axs[:-1]:
        ax.tick_params(axis="x", labelbottom=False)
    _set_group_ticks(axs[2], sample_order)
    axs[2].set_xlabel("Sample", fontsize=11)
    return fig


# ── Figure 2: per cell type ───────────────────────────────────────────────────────────

def plot_qc_per_celltype(qc, celltype_col="cell_type", figsize=(6.5, 4.8)):
    """
    Per-cell-type QC figure: genes per nucleus (boxplot), in CELLTYPE_ORDER.

    Parameters
    ----------
    qc : pd.DataFrame with one row per nucleus and columns [celltype_col, 'n_genes']
    """
    qc = qc.dropna(subset=[celltype_col])
    qc = qc.assign(**{celltype_col: qc[celltype_col].astype(str)})
    order = _order_celltypes(qc[celltype_col].unique())
    n_cells = qc[celltype_col].value_counts()
    labels = [f"{ct}\n(n = {n_cells[ct]:,})" for ct in order]

    fig, ax = plt.subplots(figsize=figsize)
    _black_boxplot(ax, qc, celltype_col, "n_genes", order, linewidth=1)
    _style_axis(ax, "Number of genes per nucleus")
    _median_label(ax, f"Median: {int(np.median(qc['n_genes'])):,} genes")
    ax.set_xlabel("")
    _set_group_ticks(ax, labels, rotation=45, fontsize=8)
    for lab in ax.get_xticklabels():
        lab.set_ha("right")
    return fig

# ======================================================================================
#  EXPORTABLE BLOCK (end)
# ======================================================================================


def main(h5ad_path, out_sample, out_celltype, scripts_dir,
         celltype_col="cell_type", leiden_col="leiden_harmony_0.2", sample_col="sample"):

    logger.info(f"Reading obs and layers['counts'] from {h5ad_path} ...")
    obs, counts, var_names = load_obs_and_counts(h5ad_path)
    logger.info(f"Loaded {counts.shape[0]:,} nuclei x {counts.shape[1]:,} genes")

    # ── Cell types: use existing column, else map clusters with cluster_anns ─────────
    palette = None
    sys.path.insert(0, scripts_dir)
    try:
        import scanpy_gene_lists as gl
        palette = getattr(gl, "custom_palette", None)
        if isinstance(palette, dict):   # same colours under the new names
            palette = {CELLTYPE_RENAME.get(k, k): v for k, v in palette.items()}
    except ImportError:
        gl = None
        logger.warning("scanpy_gene_lists not importable; default palette used")

    if celltype_col not in obs.columns:
        if gl is None or leiden_col not in obs.columns:
            raise KeyError(f"'{celltype_col}' missing and cannot be rebuilt from '{leiden_col}'")
        logger.info(f"Building '{celltype_col}' from {leiden_col} via cluster_anns ...")
        obs[celltype_col] = obs[leiden_col].astype(str).map(gl.cluster_anns)
    obs[celltype_col] = obs[celltype_col].astype(object).replace(CELLTYPE_RENAME)
    n_unmapped = int(obs[celltype_col].isna().sum())
    logger.info(f"Nuclei without a cell type: {n_unmapped}")

    # ── Metrics from raw counts ──────────────────────────────────────────────────────
    n_counts, n_genes = compute_cell_qc(counts)
    del counts
    qc = pd.DataFrame({
        sample_col: obs[sample_col].astype(str).values,
        celltype_col: obs[celltype_col].astype(object).values,
        "n_counts": n_counts,
        "n_genes": n_genes,
    })

    logger.info(f"Samples: {qc[sample_col].nunique()}  Nuclei: {len(qc):,}")
    logger.info(f"UMIs per nucleus: median {np.median(n_counts):,.0f}  range {n_counts.min():,.0f}-{n_counts.max():,.0f}")
    logger.info(f"Genes per nucleus: median {np.median(n_genes):,.0f}  range {n_genes.min():,}-{n_genes.max():,}")
    outside = int(((n_genes < MIN_GENES_PER_CELL) | (n_genes > MAX_GENES_PER_CELL)).sum())
    logger.info(f"Nuclei outside the {MIN_GENES_PER_CELL}-{MAX_GENES_PER_CELL} gene threshold (counting final gene set only): {outside:,}")

    for col, recomputed in (("n_counts", n_counts), ("n_genes", n_genes)):
        if col in obs.columns:
            agree = np.mean(np.isclose(obs[col].values.astype(float), recomputed))
            logger.info(f"obs['{col}'] matches recomputed values for {agree:.1%} of nuclei")

    # ── Figures ──────────────────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "pdf.fonttype": 42,   # editable text in the PDF
    })

    fig = plot_qc_per_sample(qc, sample_col=sample_col, celltype_col=celltype_col,
                             palette=palette)
    fig.savefig(out_sample, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved {out_sample}")

    fig = plot_qc_per_celltype(qc, celltype_col=celltype_col)
    fig.savefig(out_celltype, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved {out_celltype}")

    # ── Summary tables (Excel; not a Snakemake output) ───────────────────────────────
    tables = build_qc_tables(qc, sample_col=sample_col, celltype_col=celltype_col)
    os.makedirs(TABLE_DIR, exist_ok=True)
    table_path = os.path.join(TABLE_DIR, TABLE_FILE)
    with pd.ExcelWriter(table_path) as writer:
        for sheet, tbl in tables.items():
            tbl.to_excel(writer, sheet_name=sheet, index=False)
    logger.info(f"Saved summary tables (sheets: {', '.join(tables)}) to "
                f"{os.path.abspath(table_path)}")


if "snakemake" in globals():
    log_file = snakemake.log[0]
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.info("Generating snRNA-seq QC figures ...")

    main(
        h5ad_path=str(snakemake.input.h5ad),
        out_sample=str(snakemake.output.per_sample),
        out_celltype=str(snakemake.output.per_celltype),
        scripts_dir=str(snakemake.scriptdir),
    )
    logger.info("Done.")
