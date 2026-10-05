#--------------------------------------------------------------------------------------
#
#    Build the shareable Scanpy object and cell-to-label mapping table
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  15DATA_SHARING | Rule: final_scanpy_obj
#            Upstream:   03SCANPY L1 clustering, L2 subclustering and palantir
#                        pseudotime objects
#            Downstream: None (final data-sharing output; in rule all)
#
# Purpose:   Attach manuscript cell type labels (Glu-A / Glu-B naming) to every cell:
#            - L1: leiden_harmony_0.2 clusters mapped through cluster_anns
#            - L2: per-population subclusters (4 of 7 L1 populations; rest NA)
#            - PT: palantir pseudotime + Q4 bin per trajectory (NA if not on it);
#                  bins re-derived with pd.qcut as in manuscript_table_diff_exp.ipynb
#            Then write a slim h5ad (raw counts as X, selected obs/obsm only).
#
# Inputs:    l1     L1 clustered object (adata_clusters_v3.h5ad)
#            l2     L2 subclustering objects, one per population in l2_populations
#            pt     Pseudotime objects, one per trajectory in pt_trajectories
#
# Outputs:   h5ad     Slim shareable AnnData object
#            mapping  Cell ID -> L1 / L2 / pseudotime label table (TSV, gzipped)
#
# Memory:    Only .obs is read from the L2 / PT objects. From L1, only obs, var,
#            layers['counts'] and the requested obsm keys are read; the dense
#            scaled X is never loaded.
#
#--------------------------------------------------------------------------------------

##  Load packages, functions and variables  -------------------------------------------
import gc
import logging
import os
import sys

import anndata as ad
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

log.info("Building shareable Scanpy object and cell label table ...")

p = snakemake.params

# Report what each variable is set to
variables = {
    "l1_h5ad":         snakemake.input.l1,
    "l2_h5ad":         ", ".join(snakemake.input.l2),
    "pt_h5ad":         ", ".join(snakemake.input.pt),
    "out_h5ad":        snakemake.output.h5ad,
    "out_mapping":     snakemake.output.mapping,
    "l1_cluster_col":  p.l1_cluster_col,
    "l2_populations":  ", ".join(p.l2_populations),
    "l1_without_l2":   ", ".join(p.l1_without_l2),
    "pt_trajectories": ", ".join(p.pt_trajectories),
    "pt_n_bins":       p.pt_n_bins,
    "keep_obs":        ", ".join(p.keep_obs),
    "keep_obsm":       ", ".join(p.keep_obsm),
    "x_from_counts":   p.x_from_counts,
}
width = max(map(len, variables))
log.info("Variables\n============================\n" +
         "\n".join(f"{k:<{width}}  {v}" for k, v in variables.items()) +
         "\n============================")

# Load functions
sys.path.append(snakemake.scriptdir)


def relabel_cell_type(x):
    """Pipeline label -> manuscript label. KEEP IN SYNC with relabel_cell_type() in the Snakefile."""
    return x.replace("Glu-UL", "Glu-A").replace("Glu-DL", "Glu-B").replace("-Q4-", "-")


def read_obs(path):
    """Read only .obs from an h5ad."""
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


## L1 labels  --------------------------------------------------------------------------
log.info(f"Reading L1 obs from {snakemake.input.l1}")
l1_obs = read_obs(snakemake.input.l1)
assert l1_obs.index.is_unique, "L1 obs_names are not unique"
cells = l1_obs.index
log.info(f"L1: {len(cells)} cells")

l1_id = l1_obs[p.l1_cluster_col].astype(str)
l1_old = l1_id.map(cluster_anns)
if l1_old.isna().any():
    raise ValueError(
        f"cluster_anns has no entry for L1 clusters {sorted(l1_id[l1_old.isna()].unique())}. "
        "Check the key type in cluster_anns (str vs int).")

# Populations without L2 should be exactly those expected (e.g. MG, Endo-Peri, OPC)
no_l2 = set(l1_old.unique()) - set(p.l2_populations)
assert no_l2 == set(p.l1_without_l2), (
    f"L1 populations without an L2 object are {sorted(no_l2)}, expected {sorted(p.l1_without_l2)}")


## L2 labels  --------------------------------------------------------------------------
l2_old = pd.Series(pd.NA, index=cells, dtype="object")   # e.g. Glu-UL-3
l2_id = pd.Series(pd.NA, index=cells, dtype="object")    # leiden ID within the L1 population

for pop, path in zip(p.l2_populations, snakemake.input.l2):
    log.info(f"Reading L2 obs for {pop} from {path}")
    obs = read_obs(path)
    assert obs.index.is_unique, f"{pop}: L2 obs_names not unique"
    assert len(obs.index.difference(cells)) == 0, f"{pop}: L2 cells absent from the L1 object"

    l1_cells = set(cells[l1_old.values == pop])
    assert l1_cells == set(obs.index), (
        f"{pop}: L2 object has {obs.shape[0]} cells but L1 population has {len(l1_cells)}")

    leiden_cols = [c for c in obs.columns if c.startswith(f"leiden_{pop}_L2_harmony_")]
    assert len(leiden_cols) == 1, f"{pop}: expected one post-Harmony L2 leiden column, found {leiden_cols}"
    ids = obs[leiden_cols[0]].astype(str)
    labels = ids.map(lambda x: f"{pop}-{x}")
    if "subcluster" in obs.columns:
        assert (obs["subcluster"].astype(str) == labels).all(), \
            f"{pop}: 'subcluster' disagrees with {leiden_cols[0]}"

    l2_old.loc[obs.index] = labels.values
    l2_id.loc[obs.index] = ids.values
    log.info(f"{pop}: {obs.shape[0]} cells, {ids.nunique()} L2 clusters")
    del obs
    gc.collect()


## Pseudotime labels  ------------------------------------------------------------------
pt_cols = {}
observed_bins = set()

for traj, path in zip(p.pt_trajectories, snakemake.input.pt):
    log.info(f"Reading pseudotime obs for {traj} from {path}")
    obs = read_obs(path)
    assert obs.index.is_unique, f"{traj}: obs_names not unique"
    assert len(obs.index.difference(cells)) == 0, \
        f"{traj}: {len(obs.index.difference(cells))} cells absent from the L1 object"
    if "palantir_pseudotime" not in obs.columns:
        raise KeyError(f"'palantir_pseudotime' missing in pseudotime object for {traj}")

    # Re-derive quantile bins and check against any stored bins
    n = p.pt_n_bins
    pt = obs["palantir_pseudotime"].astype(float)
    q = pd.qcut(pt, q=n, labels=list(range(1, n + 1)))
    if f"quantile_Q{n}" in obs.columns:
        assert (obs[f"quantile_Q{n}"].astype(int) == q.astype(int)).all(), \
            f"{traj}: re-derived bins differ from stored quantile_Q{n}"
    bin_old = q.astype(object).map(lambda v: f"{traj}-Q{n}-Bin{v}" if pd.notna(v) else pd.NA)  # e.g. NPC-to-Glu-DL-Q4-Bin1
    bin_new = bin_old.map(lambda x: relabel_cell_type(x) if isinstance(x, str) else pd.NA)      # e.g. NPC-to-Glu-B-Bin1

    new_traj = relabel_cell_type(traj)
    pt_cols[f"{new_traj}_pseudotime"] = pt.reindex(cells)
    pt_cols[f"{new_traj}_bin"] = bin_new.reindex(cells)
    observed_bins |= set(bin_new.dropna())

    log.info(f"{traj}: {obs.shape[0]} cells; L1 composition:\n"
             f"{l1_old.loc[obs.index].value_counts().to_string()}")
    log.info(f"{traj}: bin counts:\n{bin_new.value_counts().sort_index().to_string()}")
    del obs, pt, q
    gc.collect()

bad_bins = sorted(observed_bins - set(p.manuscript_cell_types))
assert not bad_bins, f"Pseudotime bin labels not in config['cell_types']: {bad_bins}"


## Build annotation table (.obs of final object)  -------------------------------------
keep = [c for c in p.keep_obs if c in l1_obs.columns]
absent = [c for c in p.keep_obs if c not in l1_obs.columns]
if absent:
    log.warning(f"keep_obs columns not found and skipped: {absent}")
log.info(f"Dropping {l1_obs.shape[1] - len(keep)} .obs columns: "
         f"{[c for c in l1_obs.columns if c not in keep]}")

obs_final = l1_obs[keep].copy()
obs_final["L1_cluster_id"] = l1_id.values
obs_final["L1_cell_type"] = l1_old.map(relabel_cell_type).values
obs_final["L2_cluster_id"] = l2_id.values
obs_final["L2_cell_type"] = l2_old.map(lambda x: relabel_cell_type(x) if isinstance(x, str) else pd.NA).values
for name, col in pt_cols.items():
    obs_final[name] = col.values
for c in obs_final.columns:
    if obs_final[c].dtype == object:
        obs_final[c] = obs_final[c].astype("category")   # missing values stay as real NA

del l1_obs, l1_id, l1_old, l2_id, l2_old, pt_cols
gc.collect()

observed_l2 = set(obs_final["L2_cell_type"].dropna().astype(str))
not_in_config = sorted(observed_l2 - set(p.manuscript_cell_types))
if not_in_config:
    log.info(f"L2 labels not in config['cell_types'] (subclusters not carried forward): {not_in_config}")

# Write the lightweight table first: no large object is in memory at this point
log.info("Writing cell label mapping table ...")
os.makedirs(os.path.dirname(snakemake.output.mapping) or ".", exist_ok=True)
obs_final.rename_axis("cell_id").to_csv(snakemake.output.mapping, sep="\t", na_rep="NA")   # .gz extension -> gzip
log.info(f"Wrote {snakemake.output.mapping}")
log.info("L1 x L2 counts:\n" +
         str(pd.crosstab(obs_final["L1_cell_type"],
                         obs_final["L2_cell_type"].astype(object).fillna("NA"))))


## Write slim h5ad  --------------------------------------------------------------------
with h5py.File(snakemake.input.l1, "r") as f:
    if p.x_from_counts:
        if "counts" not in f["layers"]:
            raise KeyError("layers['counts'] not found in the L1 object")
        X = read_elem(f["layers/counts"])
        log.info("Using raw counts from layers['counts'] as X")
    else:
        X = read_elem(f["X"])
    var = read_elem(f["var"])
    obsm = {k: read_elem(f["obsm"][k]) for k in p.keep_obsm if "obsm" in f and k in f["obsm"]}
log.info(f"Keeping obsm: {list(obsm)}")

final = ad.AnnData(X=X, obs=obs_final, var=var, obsm=obsm)
del X, var, obsm
gc.collect()

final.uns["L1_to_L2"] = {
    l1: sorted(obs_final.loc[obs_final["L1_cell_type"] == l1, "L2_cell_type"]
               .dropna().astype(str).unique().tolist())
    for l1 in sorted(obs_final["L1_cell_type"].astype(str).unique())
}
log.info(f"Final object: {final.n_obs} x {final.n_vars}; obs columns: {list(final.obs.columns)}")

log.info("Writing slim h5ad ...")
os.makedirs(os.path.dirname(snakemake.output.h5ad) or ".", exist_ok=True)
final.write_h5ad(snakemake.output.h5ad, compression="gzip")
log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
