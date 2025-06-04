# Script to extract genotype, expression, and covariate data for SuSiE fine-mapping
# This script processes cis-windows for eQTLs, extracts relevant data, and prepares it for SuSiE analysis.  
# It uses PLINK2 for genotype extraction and parallel processing to handle multiple genes efficiently.
# Author: Dazcam   
# Date: 04-06-2025
# Requirements: pandas, numpy, os, Pool, PLINK2 installed and accessible in PATH
    
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool  

# File paths
cis_windows_file = snakemake.input[['sig_eGenes']]  
expr_file = snakemake.input[['pseudobulk']]  
geno_file = snakemake.input[['geno_file']] 
geno_prefix = os.path.splitext(geno_file)[0]  
covar_file = snakemake.input[['covar_file']] 
output_file = snakemake.output[0]  
output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# Load cis-windows
print(f"Loading gene list from {cis_windows_file}")
cis_windows = pd.read_csv(cis_windows_file, sep="\t")

# Load covariate data
print(f"Loading covariate data from {covar_file}")
covar_data = pd.read_csv(covar_file, sep="\t")
covar_samples = covar_data.columns[1:].tolist()  # Skip 'Unnamed: 0' column
covar_matrix = covar_data.iloc[:, 1:].to_numpy().T  # Transpose so samples are rows

# Function to extract data for one gene
def extract_data_for_gene(row):
    gene, chrom, cis_start, cis_end = row
    output_prefix = f"{output_dir}/{gene}"
    print(f"\n\nProcessing gene: {gene}, chromosome: {chrom}, start: {cis_start}, end: {cis_end}")
    print(f"Output dir: {output_prefix}")
    
    # Step 1: Generate .bed, .bim, .fam files for SNP metadata
    plink_bed_cmd = f"plink2 --pfile {geno_prefix} --chr {chrom} --from-bp {cis_start} --to-bp {cis_end} --make-bed --out {output_prefix}"
    os.system(plink_bed_cmd)
    
    # Step 2: Export genotype matrix as .raw file
    plink_raw_cmd = f"plink2 --pfile {geno_prefix} --chr {chrom} --from-bp {cis_start} --to-bp {cis_end} --export A --out {output_prefix}"
    os.system(plink_raw_cmd)
    
    # Read genotype data from .raw file
    try:
        raw_data = pd.read_csv(f"{output_prefix}.raw", sep=r"\s+")
    except Exception as e:
        print(f"Error reading {output_prefix}.raw: {e}")
        return None
    
    # Extract sample IDs and genotype matrix
    geno_samples = raw_data['IID'].tolist()
    geno_mat = raw_data.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']).to_numpy()
    
    # Read SNP information from .bim file
    try:
        bim = pd.read_csv(f"{output_prefix}.bim", sep="\t", header=None,
                          names=["chr", "snp", "cm", "pos", "a1", "a2"])
    except Exception as e:
        print(f"Error reading {output_prefix}.bim: {e}")
        return None
    
    # Load expression data
    print(f"Loading Expression data from {expr_file}")
    expr_data = pd.read_csv(expr_file, sep="\t")
    gene_expr = expr_data[expr_data["TargetID"] == gene].drop(columns=["#Chr", "start", "end", "TargetID"])
    expr_samples = gene_expr.columns.tolist()
    expr_vector = gene_expr.iloc[0].to_numpy()
    
    # Align samples
    print(f"Checking if samples align ...")
    common_samples = list(set(geno_samples) & set(expr_samples) & set(covar_samples))
    if not common_samples:
        print(f"No common samples for gene {gene}")
        return None
    
    geno_idx = [geno_samples.index(s) for s in common_samples]
    expr_idx = [expr_samples.index(s) for s in common_samples]
    covar_idx = [covar_samples.index(s) for s in common_samples]
    
    geno_mat = geno_mat[geno_idx, :]
    expr_vector = expr_vector[expr_idx]
    covar_subset = covar_matrix[covar_idx, :]
    
    # Save data for SuSiE
    np.savez(f"{output_prefix}_data.npz", X=geno_mat, y=expr_vector, Z=covar_subset, variants=bim)
    
    # Clean up temporary files
    for ext in [".raw", ".bim", ".fam", ".log"]:
        try:
            os.remove(f"{output_prefix}{ext}")
        except FileNotFoundError:
            pass
    
    return gene


# Uncomment the following lines to test with a single gene and check npz file contents
#single_gene_row = cis_windows[cis_windows['ensembl_gene_id'] == 'ENSG00000001460'].values[0]
#result = extract_data_for_gene(single_gene_row)

# Path to the .npz file
# npz_file = "../results/05SLDSR/cis_data/ENSG00000001460_data.npz"

# # Load the .npz file with allow_pickle=True
# data = np.load(npz_file, allow_pickle=True)

# # List all arrays in the .npz file
# print("Arrays in .npz file:", data.files)

# # Inspect each array
# print("\nGenotype matrix (X):")
# print("Shape:", data['X'].shape)
# print("First few rows and columns:\n", data['X'][:5, :5])
# print("Data type:", data['X'].dtype)
# print("Missing values (NaN count):", np.isnan(data['X']).sum())

# print("\nExpression vector (y):")
# print("Shape:", data['y'].shape)
# print("First few values:", data['y'][:5])
# print("Data type:", data['y'].dtype)

# print("\nCovariate matrix (Z):")
# print("Shape:", data['Z'].shape)
# print("First few rows and columns:\n", data['Z'][:5, :5])
# print("Data type:", data['Z'].dtype)

# print("\nVariants (SNP metadata):")
# variants = pd.DataFrame(data['variants'], columns=["chr", "snp", "cm", "pos", "a1", "a2"])
# print("Shape:", variants.shape)
# print("First few rows:\n", variants.head())
  
# Parallel processing
with Pool() as pool:
    results = pool.map(extract_data_for_gene, cis_windows.values)

# Filter successful extractions
successful_genes = [r for r in results if r is not None]
print(f"Processed {len(successful_genes)} genes")

# Create touch file for snakemake
touch_file = os.path.join(output_file)
with open(touch_file, 'w') as f:
    f.write(f"Processed {len(successful_genes)} genes for SuSiE fine-mapping.\n")
    f.write("Genes processed:\n")
    for gene in successful_genes:
        f.write(f"{gene}\n")
