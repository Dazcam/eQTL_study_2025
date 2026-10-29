# --------------------------------------------------------------------------------------
#
#    Pseudotime quantile bins and pseudobulk counts for one trajectory
#
# --------------------------------------------------------------------------------------
#
# Pipeline:  17PSEUDOTIME | Rule: pseudobulk_trajectory (one job per trajectory)
#            Upstream:   palantair_trajectory
#            Downstream: pseudotime_prep_tensorqtl (one job per Q value x bin)
#
# Purpose:   1. Load the pseudotime h5ad
#            2. Load the genotype sample list and find overlapping samples
#            3. Restore raw counts from layers['counts'] into adata.X
#            4. Assign equal-cell-number quantile bins for each Q value
#            5. QC figure: cell type, pseudotime and one panel per Q value
#            6. Pseudobulk CSV per Q value x bin: raw counts summed per sample,
#               genotyped samples only; samples as rows, genes as columns
#            7. Per-bin summary CSV
#            8. Sentinel file
#
# Inputs:    h5ad      Annotated pseudotime h5ad (adata_pseudotime.h5ad)
#            geno      Genotyped sample list
#
# Outputs:   sentinel  pseudobulk_complete.txt; its folder holds the outputs below
#            csvs      Q{n}/Q{n}_bin{b}_{trajectory}_pseudobulk.csv
#            (undeclared) pseudobulk_qc.png, pseudobulk_summary.csv
#
# Notes:     Q values come from PT_Q_VALUES in 17pseudotime.smk. The QC figure has
#            five panels, so it assumes three Q values. Each bin's subset is deleted
#            and gc.collect() called to limit memory use.
#
# --------------------------------------------------------------------------------------

import sys
import gc
import logging
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad

from pathlib import Path
from datetime import datetime
from scanpy_gene_lists import custom_palette
from pseudotime_functions import is_integer_valued

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

RANDOM_SEED  = 1234

## Set up logging for Snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

sc.settings.verbosity = 1

log.info("Binning pseudotime and generating pseudobulk counts ...")

##  Variables  -------------------------------------------------------------------------
input_h5ad       = str(snakemake.input.h5ad)
geno_sample_file = str(snakemake.input.geno)
sentinel_file    = str(snakemake.output.sentinel)
output_dir       = Path(sentinel_file).parent
output_figure    = str(output_dir / "pseudobulk_qc.png")
output_summary   = str(output_dir / "pseudobulk_summary.csv")
Q_VALUES         = [int(n) for n in snakemake.params.q_values]  # PT_Q_VALUES in 17pseudotime.smk

variables = {
    "input_h5ad":       input_h5ad,
    "geno_sample_file": geno_sample_file,
    "sentinel_file":    sentinel_file,
    "output_dir":       output_dir,
    "output_figure":    output_figure,
    "output_summary":   output_summary,
    "Q_VALUES":         ", ".join(map(str, Q_VALUES)),
}
width = max(map(len, variables))
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  1. Load data  ----------------------------------------------------------------------
log.info(f"Loading pseudotime adata from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
log.info(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# Read trajectory metadata
traj       = adata.uns['trajectory']
trajectory = traj['name']
root_type  = traj['root_type']
term_type  = traj['terminal_type']
cell_types = traj['cell_types']

log.info(f"Trajectory:    {trajectory}")
log.info(f"Root:          {root_type}")
log.info(f"Terminal:      {term_type}")

# Validate required columns
for col in ['palantir_pseudotime', 'sample', 'cell_type']:
    if col not in adata.obs.columns:
        raise KeyError(f"Required column '{col}' not found in adata.obs")

if 'counts' not in adata.layers:
    raise KeyError("adata.layers['counts'] not found — required for pseudobulk")

if 'X_umap' not in adata.obsm:
    raise KeyError("adata.obsm['X_umap'] not found — required for QC figure")

##  2. Load genotype sample list  ------------------------------------------------------
log.info(f"Loading genotype sample list from: {geno_sample_file}")
with open(geno_sample_file, 'r') as f:
    sample_set = set(f.read().splitlines())

all_samples      = set(adata.obs['sample'].unique())
overlapping      = all_samples & sample_set
missing_from_geno = all_samples - sample_set

log.info(f"Samples in adata:          {len(all_samples)}")
log.info(f"Samples in genotype list:  {len(sample_set)}")
log.info(f"Overlapping samples:       {len(overlapping)}")
log.info(f"Samples without genotype:  {len(missing_from_geno)}")
if missing_from_geno:
    log.info(f"Excluded samples: {sorted(missing_from_geno)}")

##  3. Restore raw counts  -------------------------------------------------------------
#
#  Replace normalised adata.X with raw integer counts from layers['counts'].
#  Pseudobulk is generated by summing these counts across cells per sample.

log.info("Restoring raw counts from layers['counts']...")
adata.X = adata.layers['counts'].copy()

# Sanity check — confirm integer values
is_integer = is_integer_valued(adata.X[:100, :100])
log.info(f"Raw counts are integers: {is_integer}")
if not is_integer:
    raise ValueError(
        "adata.layers['counts'] does not contain integer values. "
        "Check that raw counts were correctly stored during preprocessing."
    )

##  4. Assign quantile bins  -----------------------------------------------------------
#
#  For each Q value, assign cells to N equal-cell-number bins based on
#  palantir_pseudotime using pd.qcut. Labels are integers 1 to N.
#  Equal cell numbers (not equal pseudotime intervals) ensures each bin
#  has sufficient cells for pseudobulk aggregation.

log.info("Assigning pseudotime quantile bins...")
pt = adata.obs['palantir_pseudotime']

for n in Q_VALUES:
    col = f'quantile_Q{n}'
    adata.obs[col] = pd.qcut(
        pt,
        q     = n,
        labels= list(range(1, n + 1))
    ).astype(int)

    bin_counts = adata.obs[col].value_counts().sort_index()
    log.info(f"\nQ{n} bin cell counts:\n{bin_counts.to_string()}")

    # Log pseudotime boundaries for each bin
    boundaries = pd.qcut(pt, q=n, retbins=True)[1]
    log.info(f"Q{n} pseudotime boundaries: {np.round(boundaries, 4).tolist()}")

##  5. QC Figure  ----------------------------------------------------------------------
#
#  Five-panel UMAP figure plus summary table.
#  All UMAP panels use direct matplotlib scatter to avoid scanpy rendering
#  issues at high cell counts.

log.info("Generating QC figure...")

palette  = {ct: custom_palette[ct] for ct in cell_types if ct in custom_palette}
umap1    = adata.obsm['X_umap'][:, 0]
umap2    = adata.obsm['X_umap'][:, 1]
pt_vals  = adata.obs['palantir_pseudotime'].values
pt_size  = max(0.2, 20000 / adata.n_obs)

# Quantile colours matching report Panel E scheme
QUANTILE_COLOURS = {
    1: '#e6a023', 2: '#61b6e5', 3: '#efe245',
    4: '#0173b2', 5: '#2aa47e', 6: '#d4632e'
}

def get_bin_colours(n):
    return [QUANTILE_COLOURS.get(i + 1, '#AAAAAA') for i in range(n)]

fig = plt.figure(figsize=(24, 6))
gs  = fig.add_gridspec(1, 5, hspace=0.35, wspace=0.3)

umap_axes = [fig.add_subplot(gs[0, i]) for i in range(5)]

# Panel A — cell type
for ct in adata.obs['cell_type'].cat.categories:
    mask  = adata.obs['cell_type'] == ct
    color = palette.get(ct, '#AAAAAA')
    umap_axes[0].scatter(
        umap1[mask], umap2[mask],
        c=color, s=pt_size, alpha=0.5,
        label=ct, rasterized=True, linewidths=0
    )
umap_axes[0].legend(markerscale=6, fontsize=7, loc='best', framealpha=0.8)
umap_axes[0].set_title('Cell type', fontsize=10)
umap_axes[0].axis('off')

# Panel B — pseudotime
sc_b = umap_axes[1].scatter(
    umap1, umap2, c=pt_vals, cmap='viridis',
    s=pt_size, alpha=0.5, rasterized=True, linewidths=0
)
plt.colorbar(sc_b, ax=umap_axes[1], shrink=0.6, label='Pseudotime')
umap_axes[1].set_title('Palantir pseudotime', fontsize=10)
umap_axes[1].axis('off')

# Panels C-E — Q4, Q5, Q6 bins
for ax_idx, n in enumerate(Q_VALUES):
    col     = f'quantile_Q{n}'
    bins    = adata.obs[col].values
    colours = get_bin_colours(n)
    ax      = umap_axes[ax_idx + 2]

    for b in range(1, n + 1):
        mask = bins == b
        ax.scatter(
            umap1[mask], umap2[mask],
            c=[colours[b - 1]], s=pt_size, alpha=0.5,
            label=f'Q{b}', rasterized=True, linewidths=0
        )
    ax.legend(markerscale=6, fontsize=7, loc='best',
              framealpha=0.8, title=f'Q{n}', title_fontsize=7)
    ax.set_title(f'Q{n} bins', fontsize=10)
    ax.axis('off')

fig.suptitle(
    f'Pseudotime quantile bins: {trajectory}\n'
    f'(n={adata.shape[0]:,} cells)',
    fontsize=12, fontweight='bold', y=1.01
)
fig.savefig(output_figure, dpi=200, bbox_inches='tight')
plt.close(fig)
log.info(f"QC figure saved: {output_figure}")

##  6. Generate pseudobulk CSVs  -------------------------------------------------------
#
#  For each Q value and each bin:
#    - Subset adata to cells in that bin
#    - Further subset to genotype-matched samples
#    - Sum raw counts per sample
#    - Validate integer values
#    - Write CSV (samples as rows, genes as columns)
#    - Delete subset and collect garbage
#
#  One CSV per bin per Q value, written to Q{n}/ subdirectory.

summary_rows = []

for n in Q_VALUES:
    col     = f'quantile_Q{n}'
    q_dir   = output_dir / f'Q{n}'
    q_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"\n── Generating Q{n} pseudobulk ───────────────────────────")

    for b in range(1, n + 1):
        log.info(f"  Q{n} bin {b}...")

        # Subset to bin
        bin_mask = adata.obs[col] == b
        bin_adata = adata[bin_mask].copy()

        # Filter to genotype-matched samples
        geno_mask  = bin_adata.obs['sample'].isin(sample_set)
        bin_adata  = bin_adata[geno_mask].copy()

        n_cells   = bin_adata.n_obs
        samples   = bin_adata.obs['sample'].unique()
        n_samples = len(samples)
        log.info(f"    Cells: {n_cells:,} | Samples: {n_samples}")

        if n_samples == 0:
            log.warning(f"    No genotype-matched samples in Q{n} bin {b} — skipping")
            del bin_adata
            gc.collect()
            continue

        # Generate pseudobulk — sum counts per sample
        pseudobulk_rows = []
        for sample in samples:
            samp_mask   = bin_adata.obs['sample'] == sample
            samp_adata  = bin_adata[samp_mask]
            counts      = samp_adata.X

            if sp.issparse(counts):
                pseudobulk_counts = np.array(counts.sum(axis=0)).flatten()
            else:
                pseudobulk_counts = np.array(counts.sum(axis=0)).flatten()

            pseudobulk_rows.append(pseudobulk_counts)

        pseudobulk_df = pd.DataFrame(
            pseudobulk_rows,
            index   = samples,
            columns = bin_adata.var_names
        )

        # Validate integer values in pseudobulk
        is_integer  = is_integer_valued(pseudobulk_df.values)
        log.info(f"    Pseudobulk values are integers: {is_integer}")
        if not is_integer:
            raise ValueError(
                f"Non-integer values found in pseudobulk for Q{n} bin {b}. "
                "Check that raw counts were correctly restored into adata.X."
            )

        # Write CSV
        out_csv = str(q_dir / f"Q{n}_bin{b}_{trajectory}_pseudobulk.csv")
        pseudobulk_df.to_csv(out_csv)
        log.info(f"    Written: {out_csv} ({pseudobulk_df.shape})")

        # Collect summary stats
        cells_per_sample = bin_adata.obs.groupby('sample').size()
        summary_rows.append({
            'trajectory'          : trajectory,
            'quantile_n'          : n,
            'bin'                 : b,
            'n_cells_total'       : n_cells,
            'n_samples'           : n_samples,
            'min_cells_per_sample': int(cells_per_sample.min()),
            'median_cells_per_sample': int(cells_per_sample.median()),
            'max_cells_per_sample': int(cells_per_sample.max()),
            'n_samples_below_10'  : int((cells_per_sample < 10).sum()),
        })

        # Clean up
        del bin_adata, pseudobulk_df, pseudobulk_rows
        gc.collect()

##  7. Write summary CSV  --------------------------------------------------------------
log.info("Writing pseudobulk summary CSV...")
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(output_summary, index=False)
log.info(f"Summary CSV written: {output_summary} ({summary_df.shape})")
log.info(f"\nPseudobulk summary:\n{summary_df.to_string(index=False)}")

##  8. Write sentinel file  ------------------------------------------------------------
log.info("Writing sentinel file...")
with open(sentinel_file, 'w') as f:
    f.write(f"pseudobulk_trajectory completed\n")
    f.write(f"trajectory:  {trajectory}\n")
    f.write(f"timestamp:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"cells:       {adata.shape[0]:,}\n")
    f.write(f"samples:     {len(overlapping)}\n")
    f.write(f"Q values:    {Q_VALUES}\n")

log.info(f"Sentinel written: {sentinel_file}")
log.info("Done.")

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
