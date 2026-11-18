#--------------------------------------------------------------------------------------
#
#    L1 cell type UMAP, marker feature plots and cell type proportions (Fig 1B-D)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: umap_fig
#            Upstream:   02 scanpy final clustered object (adata_clusters_v3.h5ad)
#            Downstream: Fig 1B-D (Fig 1A is made outside the pipeline)
#
# Purpose:   Plotting code from scanpy_annotation.ipynb:
#            - Fig 1B-C: UMAP coloured by L1 cell type, with feature plots of
#              the marker genes in scanpy_gene_lists.final_genes
#              (plot_celltype_and_gene_features, no panel labels)
#            - Fig 1D: proportion of nuclei per L1 cell type in each sample
#            Cell types from leiden_harmony_0.2 with cluster_anns; manuscript
#            labels (config manuscript_labels), ordered as config fetal_l1_types;
#            colours from scanpy_gene_lists.custom_palette
#
# Inputs:    h5ad        Final clustered object
#            labels      config['manuscript_labels']
#            l1_types    config['fetal_l1_types'] (cell type order)
#
# Outputs:   bc_pdf, bc_svg   Fig 1B-C (PDF and SVG)
#            d_pdf, d_svg     Fig 1D (PDF and SVG)
#
#--------------------------------------------------------------------------------------

##  Load packages, functions and variables  -------------------------------------------
import logging
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.sparse import issparse

# Set up logging for Snakemake (uncaught errors also go to the log)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Plotting Fig 1B-D ...")

# Load functions and gene lists
sys.path.append(snakemake.scriptdir)
from scanpy_utils import plot_celltype_and_gene_features  # noqa: E402
from scanpy_gene_lists import cluster_anns, custom_palette, final_genes  # noqa: E402

LABELS   = dict(snakemake.params.labels)
L1_TYPES = list(snakemake.params.l1_types)

variables = {
    "h5ad":        snakemake.input.h5ad,
    "labels":      ", ".join(f"{k} -> {v}" for k, v in LABELS.items()),
    "l1_types":    ", ".join(L1_TYPES),
    "final_genes": ", ".join(final_genes),
    "bc_pdf":      snakemake.output.bc_pdf,
    "bc_svg":      snakemake.output.bc_svg,
    "d_pdf":       snakemake.output.d_pdf,
    "d_svg":       snakemake.output.d_svg,
}
width = max(map(len, variables))
log.info("Variables\n============================\n" +
         "\n".join(f"{k:<{width}}  {v}" for k, v in variables.items()) +
         "\n============================")


def relabel(x):
    """Pipeline -> manuscript label, applying config manuscript_labels in order."""
    for old, new in LABELS.items():
        x = x.replace(old, new)
    return x


def plot_stacked_figure(adata, sample_column, color_column=None, barplot=True, boxplots=False, recalculate_columns=True):
    """
    Create a stacked figure. Supports boxplots for QC and/or a stacked barplot for proportions.
    Returns the figure object for saving.
    """

    # 1. Data Processing & Recalculation logic
    adata.obs[sample_column] = adata.obs[sample_column].astype(str)

    if recalculate_columns:
        # Simple check for raw counts to update n_counts/n_genes
        first_values = adata.X[:5, :5].toarray() if issparse(adata.X) else adata.X[:5, :5]
        if np.all(first_values >= 0) and np.all(np.isclose(first_values, first_values.astype(int))):
            adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if issparse(adata.X) else adata.X.sum(axis=1)
            adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten() if issparse(adata.X) else (adata.X > 0).sum(axis=1)

    # 2. Setup Plotting Data
    sample_order = adata.obs[sample_column].unique()
    cluster_colors = adata.uns.get(f'{color_column}_colors', None) if color_column else None

    # 3. Define Grid Logic
    # Dynamically set rows: 2 for boxplots, 1 for barplot
    rows_needed = []
    if boxplots: rows_needed.extend([1, 1]) # UMI and Gene rows
    if barplot: rows_needed.append(2)        # Barplot row (weighted higher)

    if not rows_needed:
        raise ValueError("At least one of 'barplot' or 'boxplots' must be True.")

    fig, axs = plt.subplots(
        len(rows_needed), 1,
        figsize=(24, 4 * len(rows_needed)),
        sharex=True,
        gridspec_kw={'height_ratios': rows_needed}
    )

    # Handle axs indexing if only one plot exists
    if len(rows_needed) == 1:
        axs = [axs]

    ax_idx = 0

    # 4. Boxplots (Only if requested)
    if boxplots:
        # UMI Plot
        umi_data = pd.DataFrame({'Sample': np.repeat(sample_order, [len(adata.obs[adata.obs[sample_column]==s]) for s in sample_order]),
                                 'UMIs': adata.obs.groupby(sample_column)["n_counts"].apply(list).loc[sample_order].explode()})
        sns.boxplot(data=umi_data, x="Sample", y="UMIs", ax=axs[ax_idx], color="white", fliersize=1, linewidth=1)
        axs[ax_idx].set_ylabel("UMIs per cell")
        ax_idx += 1

        # Gene Plot
        gene_data = pd.DataFrame({'Sample': np.repeat(sample_order, [len(adata.obs[adata.obs[sample_column]==s]) for s in sample_order]),
                                  'Genes': adata.obs.groupby(sample_column)["n_genes"].apply(list).loc[sample_order].explode()})
        sns.boxplot(data=gene_data, x="Sample", y="Genes", ax=axs[ax_idx], color="white", fliersize=1, linewidth=1)
        axs[ax_idx].set_ylabel("Genes per cell")
        ax_idx += 1

    # 5. Barplot (Only if requested)
    if barplot:
        proportions = adata.obs.groupby([sample_column, color_column]).size().unstack(fill_value=0)
        proportions = proportions.div(proportions.sum(axis=1), axis=0).reindex(sample_order)

        proportions.plot(kind="bar", stacked=True, ax=axs[ax_idx], width=1, legend=False, color=cluster_colors)
        axs[ax_idx].set_ylabel("Cell class proportions")
        axs[ax_idx].set_xlabel("Sample")

        if cluster_colors is not None:
            axs[ax_idx].legend(bbox_to_anchor=(1.01, 1), loc="upper left", title=color_column)

    # Clean up X-axis labels for the bottom-most plot
    axs[-1].tick_params(axis='x', rotation=90, labelsize=7)

    plt.tight_layout()
    return fig


##  Load data and assign cell types  --------------------------------------------------
log.info(f"Loading: {snakemake.input.h5ad}")
adata = sc.read(snakemake.input.h5ad)
log.info(f"Object dimensions (cells x genes): {adata.shape}")

adata.obs['cell_type'] = adata.obs['leiden_harmony_0.2'].map(cluster_anns)
if adata.obs['cell_type'].isna().any():
    raise ValueError("Some leiden_harmony_0.2 clusters have no entry in cluster_anns")

# Manuscript labels, in config fetal_l1_types order
order = [relabel(ct) for ct in L1_TYPES]
adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'].astype(str).map(relabel), categories=order)
if adata.obs['cell_type'].isna().any():
    raise ValueError("Some cell types are not in config fetal_l1_types")
palette = {relabel(k): v for k, v in custom_palette.items()}
log.info(f"Nuclei per cell type:\n{adata.obs['cell_type'].value_counts(sort=False).to_string()}")


##  Fig 1B-C: cell type UMAP and marker feature plots  --------------------------------
log.info("Plotting UMAP and feature plots ...")
fig = plot_celltype_and_gene_features(
    adata=adata,
    final_genes=final_genes,
    umap_palette=palette,
    cell_type_column='cell_type',
    ncols=4,
    figsize_width=10.5,
    base_height_per_row=1.0,
    panel_labels=False
)
fig.tight_layout()
fig.savefig(snakemake.output.bc_pdf, dpi=300, bbox_inches="tight")
fig.savefig(snakemake.output.bc_svg, bbox_inches="tight")
plt.close(fig)
log.info(f"Saved {snakemake.output.bc_pdf} and {snakemake.output.bc_svg}")


##  Fig 1D: cell type proportions per sample  -----------------------------------------
log.info("Plotting cell type proportions per sample ...")
fig = plot_stacked_figure(
    adata,
    sample_column="sample",
    color_column="cell_type",
    barplot=True,
    boxplots=False,
    recalculate_columns=True
)
fig.savefig(snakemake.output.d_pdf, dpi=300, bbox_inches="tight")
fig.savefig(snakemake.output.d_svg, bbox_inches="tight")
plt.close(fig)
log.info(f"Saved {snakemake.output.d_pdf} and {snakemake.output.d_svg}")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
