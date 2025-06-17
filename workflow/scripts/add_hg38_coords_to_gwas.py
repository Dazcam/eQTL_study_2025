import pandas as pd
import sys

# Redirect all output and error messages to the Snakemake log
log_file = snakemake.log[0]
sys.stdout = open(log_file, "w")
sys.stderr = sys.stdout

print("=== Starting script: add_hg38_coords_to_sumstats.py ===")

# Get file paths from Snakemake
lifted_file = str(snakemake.input["lifted"])
sumstats_file = str(snakemake.input["sumstats"])
output_file = str(snakemake.output[0])

print(f"Input liftover BED file: {lifted_file}")
print(f"Input summary stats file: {sumstats_file}")
print(f"Output file will be written to: {output_file}")

# Step 1: Read lifted BED file
print("Reading lifted BED file...")
lifted = pd.read_csv(lifted_file, sep='\t', header=None, comment='#')

# Add column names and convert to 1-based coordinates
print("Processing BED coordinates and SNP IDs...")
lifted.columns = ['CHR', 'START', 'END', 'SNP']
lifted['CHR'] = lifted['CHR'].str.replace('chr', '', regex=False)
lifted['BP'] = lifted['START'] + 1  # Convert from 0-based to 1-based
lifted = lifted[['SNP', 'CHR', 'BP']]
print(f"Lifted file processed: {lifted.shape[0]} SNPs")

# Step 2: Read original summary stats file
print("Reading summary statistics file...")
sumstats = pd.read_csv(sumstats_file, sep='\t')
print(f"Summary stats read: {sumstats.shape[0]} SNPs")

# Step 3: Merge lifted BED with summary stats by SNP ID
print("Merging lifted coordinates with summary statistics...")
merged = pd.merge(lifted, sumstats, on='SNP', how='inner')
print(f"Merged dataset: {merged.shape[0]} SNPs")

# Step 4: Replace CHR and BP with hg38 versions
print("Replacing CHR and BP columns with hg38 values...")
merged['CHR'] = merged['CHR_x']
merged['BP'] = merged['BP_x']
merged = merged.drop(columns=['CHR_x', 'CHR_y', 'BP_x', 'BP_y'])

# Step 5: Reorder columns
print("Reordering columns...")
expected_columns = ['SNP', 'CHR', 'BP', 'PVAL', 'A1', 'A2', 'BETA', 'SE', 'Z', 'N']
missing_cols = [col for col in expected_columns if col not in merged.columns]
if missing_cols:
    raise ValueError(f"Missing expected columns in merged DataFrame: {missing_cols}")

merged = merged[expected_columns]

# Step 6: Write output
print("Writing output file...")
merged.to_csv(output_file, sep='\t', index=False)
print("Output successfully written.")

print("=== Script completed successfully ===")
