import pandas as pd

# Access Snakemake input/output
lifted_file = snakemake.input.lifted
sumstats_file = snakemake.input.sumstats
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Optional: redirect stdout/stderr to log file (if needed)
import sys
sys.stderr = open(log_file, "w")
sys.stdout = sys.stderr

# Process files
lifted = pd.read_csv(lifted_file, sep='\t', header=None, comment='#')
lifted.columns = ['CHR', 'START', 'END', 'SNP']
lifted['CHR'] = lifted['CHR'].str.replace('chr', '', regex=False)
lifted['BP'] = lifted['START'] + 1  # BED is 0-based
lifted = lifted[['SNP', 'CHR', 'BP']]

sumstats = pd.read_csv(sumstats_file, sep='\t')
merged = pd.merge(lifted, sumstats, on='SNP', how='inner')
merged['CHR'] = merged['CHR_x']
merged['BP'] = merged['BP_x']
merged = merged.drop(columns=['CHR_x', 'CHR_y', 'BP_x', 'BP_y'])

columns_names = ['SNP', 'CHR', 'BP', 'PVAL', 'A1', 'A2', 'BETA', 'SE', 'Z', 'N']
merged = merged[columns_names]

merged.to_csv(output_file, sep='\t', index=False)