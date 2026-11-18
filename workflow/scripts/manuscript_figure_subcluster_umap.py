#--------------------------------------------------------------------------------------
#
#    L2 subcluster UMAPs (Supp Fig 1)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: subcluster_umap_fig
#            Upstream:   02 scanpy subclustering (adata_{l1_pop}_subclusters.h5ad)
#            Downstream: Supp Fig 1
#
# Purpose:   Plotting code from scanpy_extra2.ipynb: a 2 x 2 grid with one UMAP
#            per L2 population, coloured by subcluster with on-data labels; panels
#            A-D in panel_order (the order given in the Supp Fig 1 legend).
#            Every subcluster in each object is shown. Manuscript labels (config
#            manuscript_labels).
#
# Inputs:    h5ad         One subcluster object per population, in panel_order
#            panel_order  L1 populations, one per panel (pipeline labels)
#            labels       config['manuscript_labels']
#
# Outputs:   pdf          Supp Fig 1 (PDF)
#
#--------------------------------------------------------------------------------------

##  Load packages, functions and variables  -------------------------------------------
import logging
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scanpy as sc

# Set up logging for Snakemake (uncaught errors also go to the log)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Plotting Supp Fig 1 ...")

LABELS      = dict(snakemake.params.labels)
PANEL_ORDER = list(snakemake.params.panel_order)
H5AD_FILES  = list(snakemake.input.h5ad)

variables = {
    "h5ad":        ", ".join(H5AD_FILES),
    "panel_order": ", ".join(PANEL_ORDER),
    "labels":      ", ".join(f"{k} -> {v}" for k, v in LABELS.items()),
    "pdf":         snakemake.output.pdf,
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


##  Load subcluster objects  ----------------------------------------------------------
adata_dict = {}
for cell_type, path in zip(PANEL_ORDER, H5AD_FILES):
    log.info(f"Loading L2 data for {cell_type}: {path}")
    adata_sub = sc.read(path)

    if 'subcluster' not in adata_sub.obs.columns:
        raise KeyError(f"'subcluster' column missing in L2 object for {cell_type}")

    # Manuscript labels for the on-data legend
    adata_sub.obs['subcluster'] = (
        adata_sub.obs['subcluster']
        .astype(str)
        .map(relabel)
        .astype('category')
    )
    log.info(f"{cell_type} subclusters: {', '.join(adata_sub.obs['subcluster'].cat.categories)}")

    adata_dict[cell_type] = adata_sub


##  Supp Fig 1: 2 x 2 UMAP grid  ------------------------------------------------------
fig = plt.figure(figsize=(10, 9.5))
gs = gridspec.GridSpec(nrows=2, ncols=2,
                       figure=fig,
                       wspace=0.1, hspace=0.1)

# Titles / panel labels
panel_labels = ["A", "B", "C", "D"]

for i, cell_type in enumerate(PANEL_ORDER):
    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col])

    adata = adata_dict[cell_type]

    sc.pl.umap(
        adata,
        color="subcluster",
        ax=ax,
        palette=None,
        legend_loc="on data",
        legend_fontsize=9,
        legend_fontoutline=2,
        frameon=False,
        title="",
        show=False
    )

    # Panel label (A/B/C/D) in top-left of each UMAP
    ax.text(-0.03, 1.03, panel_labels[i],
            transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

plt.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig(snakemake.output.pdf, dpi=300, bbox_inches="tight")
log.info(f"Saved 2 x 2 UMAP figure to: {snakemake.output.pdf}")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
