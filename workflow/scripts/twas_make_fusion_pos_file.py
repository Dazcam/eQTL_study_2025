import pandas as pd
import os
import sys

# ------------------------------------------------------------------
# Input parameters
# ------------------------------------------------------------------
cis_window = int(snakemake.params.cis_window)
cell_type  = snakemake.wildcards.cell_type
coord_file = snakemake.params.coord_file
weight_dir = os.path.dirname(snakemake.input[0])  # directory containing weights
weight_files = [
    os.path.abspath(os.path.join(weight_dir, f))
    for f in os.listdir(weight_dir)
    if f.endswith('.wgt.RDat') and os.path.getsize(os.path.join(weight_dir, f)) > 0
]
output_pos = snakemake.output[0]
sample_size = 133

# ------------------------------------------------------------------
# Helper: log to stdout (Snakemake captures stdout → log file)
# ------------------------------------------------------------------
def log(msg):
    print(f"[make_pos:{cell_type}] {msg}", flush=True)

log("Starting")
log(f"Loading coordinates: {coord_file}")
log(f"Found {len(weight_files)} valid .wgt.RDat files")

# ------------------------------------------------------------------
# Load coordinates 
# ------------------------------------------------------------------
try:
    coords = pd.read_csv(coord_file, sep='\t')
except Exception as e:
    log(f"ERROR reading coord file: {e}")
    raise

# Standardize column names (case-insensitive)
cols = {c.strip().lower(): c for c in coords.columns}
required = ['chr', 'start', 'end', 'gene_id']
missing = [r for r in required if r not in cols]
if missing:
    log(f"ERROR: Missing columns in coord file: {missing}")
    log(f"Available: {list(coords.columns)}")
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

# ------------------------------------------------------------------
# Write .pos file
# ------------------------------------------------------------------
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
            log(f"Processed {i}/{len(weight_files)} | Written: {written} | Skipped: {skipped}")

log(f"Done. Output: {output_pos}")
log(f"Total written: {written}, skipped: {skipped}")
