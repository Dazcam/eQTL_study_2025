#--------------------------------------------------------------------------------------
#
#    Generate nuclei counts per L1 cluster, L2 subcluster and pseudotime bin
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  13MANUSCRIPT_TABLES | Rule: cell_counts_tbl
#            Upstream:   03SCANPY L1 / L2 objects; 17PSEUDOTIME pseudotime objects
#            Downstream: None (supplementary table / methods numbers)
#
# Purpose:   Count nuclei per L1 cluster, L2 subcluster and Q4 pseudotime bin,
#            using manuscript labels (Glu-UL -> Glu-A, Glu-DL -> Glu-B, "-Q4-"
#            dropped). Bins are re-derived with pd.qcut, as for eQTL and DE.
#
# Inputs:    l1   adata_clusters_v3.h5ad
#            l2   adata_{l1_pop}_subclusters.h5ad, one per L2 population
#            pt   {trajectory}/adata_pseudotime.h5ad, one per trajectory
#
# Outputs:   xlsx Sheets: L1, L2, Pseudotime bins
#
# Notes:     Only .obs is read from each object. L2 population and trajectory
#            names are taken from the file paths.
#
#--------------------------------------------------------------------------------------

##  Load packages, functions and variables  -------------------------------------------
import logging
import os
import sys

import h5py
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

# Set up logging for Snakemake (uncaught errors also go to the log)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Counting nuclei per L1 cluster, L2 subcluster and pseudotime bin ...")

# Load functions
sys.path.append(snakemake.scriptdir)
from scanpy_gene_lists import cluster_anns  # noqa: E402  (L1 cluster ID -> cell type)

L1_CLUSTER_COL = "leiden_harmony_0.2"
N_BINS = 4

log.info("Variables\n============================\n"
         f"l1_h5ad   {snakemake.input.l1}\n"
         f"l2_h5ad   {', '.join(snakemake.input.l2)}\n"
         f"pt_h5ad   {', '.join(snakemake.input.pt)}\n"
         f"out_file  {snakemake.output[0]}\n"
         "============================")


def relabel(x):
    """Pipeline label -> manuscript label."""
    return x.replace("Glu-UL", "Glu-A").replace("Glu-DL", "Glu-B").replace("-Q4-", "-")


def read_obs(path):
    """Read only .obs from an h5ad."""
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def count(labels, name):
    """Nuclei per label and % of the total."""
    tbl = labels.map(relabel).value_counts().sort_index().rename_axis(name).reset_index(name="Nuclei")
    tbl["%"] = (100 * tbl["Nuclei"] / tbl["Nuclei"].sum()).round(2)
    return tbl


## L1  ----------------------------------------------------------------------------------
l1 = read_obs(snakemake.input.l1)[L1_CLUSTER_COL].astype(str).map(cluster_anns)
if l1.isna().any():
    raise ValueError("Some L1 clusters have no entry in cluster_anns")
l1_tbl = count(l1, "L1 cluster")
log.info(f"L1:\n{l1_tbl.to_string(index=False)}")


## L2  ----------------------------------------------------------------------------------
l2_tbls = []
for path in snakemake.input.l2:
    pop = os.path.basename(path).removeprefix("adata_").removesuffix("_subclusters.h5ad")
    tbl = count(read_obs(path)["subcluster"].astype(str), "L2 subcluster")
    tbl.insert(0, "L1 cluster", relabel(pop))
    l2_tbls.append(tbl)
l2_tbl = pd.concat(l2_tbls, ignore_index=True).rename(columns={"%": "% of L1 cluster"})
log.info(f"L2:\n{l2_tbl.to_string(index=False)}")


## Pseudotime bins  ---------------------------------------------------------------------
pt_tbls = []
for path in snakemake.input.pt:
    traj = os.path.basename(os.path.dirname(path))
    pt = read_obs(path)["palantir_pseudotime"].astype(float)
    bins = pd.qcut(pt, q=N_BINS, labels=list(range(1, N_BINS + 1))).astype(int)
    tbl = count(f"{traj}-Q{N_BINS}-Bin" + bins.astype(str), "Bin")
    tbl.insert(0, "Trajectory", relabel(traj))
    pt_tbls.append(tbl)
pt_tbl = pd.concat(pt_tbls, ignore_index=True).rename(columns={"%": "% of trajectory"})
log.info(f"Pseudotime bins:\n{pt_tbl.to_string(index=False)}")


## Write output  ------------------------------------------------------------------------
log.info(f"Writing: {snakemake.output[0]}")
with pd.ExcelWriter(snakemake.output[0]) as writer:
    l1_tbl.to_excel(writer, sheet_name="L1", index=False)
    l2_tbl.to_excel(writer, sheet_name="L2", index=False)
    pt_tbl.to_excel(writer, sheet_name="Pseudotime bins", index=False)

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
