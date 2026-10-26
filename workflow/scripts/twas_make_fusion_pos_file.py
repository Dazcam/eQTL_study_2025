#--------------------------------------------------------------------------------------
#
#    Generate FUSION .pos file
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  11TWAS | Rule: make_pos_file (one job per cell type)
#            Upstream:   compute_weights (one .wgt.RDat per gene),
#                        prep_exp_data (gene coordinate file)
#            Downstream: twas_weights_report; 12CTWAS copy_fusion_weights
#
# Purpose:   - List every non-empty .wgt.RDat file for the cell type in a FUSION
#              .pos file, with gene ID, chromosome, P0, P1 and sample size
#
# Inputs:    weight_dir  Folder of per-gene weight files (folder of the
#                        compute_weights sentinel file)
#            coord_file  Gene coordinate file: chr, start, end, gene_id
#            cis_window  Window (bp) added to gene start/end for P0/P1
#
# Outputs:   output_pos  FUSION .pos file: PANEL, WGT, ID, CHR, P0, P1, N
#
# Notes:     sample_size (N) is fixed at 133.
#            Empty .wgt.RDat files (genes with no cis SNPs or a failed FUSION
#            run) are skipped.
#
#--------------------------------------------------------------------------------------

import logging
import os
import sys

import pandas as pd

## Set up logging for Snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Generating FUSION .pos file ...")

##  Load variables  -------------------------------------------------------------------
cis_window = int(snakemake.params.cis_window)
cell_type  = snakemake.wildcards.cell_type
coord_file = snakemake.params.coord_file
weight_dir = os.path.dirname(snakemake.input[0])
output_pos = snakemake.output[0]
sample_size = 133

variables = {
    "cell_type":   cell_type,
    "weight_dir":  weight_dir,
    "coord_file":  coord_file,
    "cis_window":  cis_window,
    "sample_size": sample_size,
    "output_pos":  output_pos,
}
width = max(len(k) for k in variables)
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  Find weight files  ----------------------------------------------------------------
weight_files = [
    os.path.abspath(os.path.join(weight_dir, f))
    for f in os.listdir(weight_dir)
    if f.endswith('.wgt.RDat') and os.path.getsize(os.path.join(weight_dir, f)) > 0
]
log.info(f"Found {len(weight_files)} valid .wgt.RDat files")

##  Load coordinates  -----------------------------------------------------------------
log.info(f"Loading coordinates: {coord_file}")
coords = pd.read_csv(coord_file, sep='\t')

# Standardise column names (case-insensitive)
cols = {c.strip().lower(): c for c in coords.columns}
required = ['chr', 'start', 'end', 'gene_id']
missing = [r for r in required if r not in cols]
if missing:
    log.error(f"Missing columns in coord file: {missing}")
    log.error(f"Available: {list(coords.columns)}")
    raise KeyError(f"Missing columns: {missing}")

# Map to expected names
coords = coords.rename(columns={
    cols['chr']: 'CHR',
    cols['start']: 'P0',
    cols['end']: 'P1',
    cols['gene_id']: 'ID'
})

# Ensure types
coords['CHR'] = coords['CHR'].astype(str).str.replace('^chr', '', regex=True)  # remove 'chr' if present
coords['P0']  = pd.to_numeric(coords['P0'], errors='coerce')
coords['P1']  = pd.to_numeric(coords['P1'], errors='coerce')
coords['ID']  = coords['ID'].astype(str)

# Add 'chr' prefix option
#coords['CHR'] = 'chr' + coords['CHR']

##  Write .pos file  ------------------------------------------------------------------
written = 0
skipped = 0

with open(output_pos, 'w') as f:
    f.write("PANEL\tWGT\tID\tCHR\tP0\tP1\tN\n")  # FUSION expects this header

    for i, wgt_path in enumerate(weight_files, 1):
        gene = os.path.basename(wgt_path).replace('.wgt.RDat', '')
        row = coords[coords['ID'] == gene]

        if row.empty:
            skipped += 1
            continue

        chr_fusion = row['CHR'].iloc[0]
        p0 = max(0, int(row['P0'].iloc[0]) - cis_window)
        p1 = int(row['P1'].iloc[0]) + cis_window

        f.write(f"{cell_type}\t{os.path.basename(wgt_path)}\t{gene}\t{chr_fusion}\t{p0}\t{p1}\t{sample_size}\n")
        written += 1

        if i % 200 == 0 or i == len(weight_files):
            log.info(f"Processed {i}/{len(weight_files)} | Written: {written} | Skipped: {skipped}")

log.info(f"Output: {output_pos}")
log.info(f"Total written: {written}, skipped: {skipped}")
log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
