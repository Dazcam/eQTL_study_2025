# Script to extract genotype, expression, and covariate data for SuSiE fine-mapping
# This script processes cis-windows for eQTLs, extracts relevant data, and prepares it for SuSiE analysis.
# It uses PLINK2 for genotype extraction and parallel processing to handle multiple genes efficiently.
# Author: Dazcam
# Date: 04-06-2025
# Requirements: pandas, numpy, os, multiprocessing, PLINK2 installed and accessible in PATH

import pandas as pd
import numpy as np
import os
import logging
from multiprocessing import Pool
import subprocess

# Configure logging to write only to Snakemake's log file
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force=True
)

# File paths
cis_windows_file = snakemake.input['sig_eGenes']
expr_file = snakemake.input['pseudobulk']
geno_file = snakemake.input['genotypes']
covar_file = snakemake.input['covariates']

# For testing
cis_windows_file = "../results/05SLDSR/eGenes_fdr_0.05/ExN-1/ExN-1_eGenes_fdr_0.05.tsv"
expr_file = "../results/03SCANPY/pseudobulk/ExN-1_tmm.bed"
geno_file = "../results/04TENSORQTL/chrALL_final.filt.pgen"
covar_file = "../results/03SCANPY/pseudobulk/ExN-1_covariates.txt"
output_file = "../results/05SLDSR/eGenes_for_SuSiE/ExN-1/ExN-1_prep_susie_input.done"

output_file = snakemake.output[0]
geno_prefix = os.path.splitext(geno_file)[0]
output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)



# Log inputs for debugging
logging.info(f"Inputs: sig_eGenes={cis_windows_file}, pseudobulk={expr_file}, genotypes={geno_file}, covariates={covar_file}")
logging.info(f"Output: {output_file}")

# Load cis-windows
logging.info(f"Loading gene list from {cis_windows_file}")
cis_windows = pd.read_csv(cis_windows_file, sep="\t")
cis_windows = cis_windows.head(1)

# Load covariate data
logging.info(f"Loading covariate data from {covar_file}")
covar_data = pd.read_csv(covar_file, sep="\t")
covar_samples = covar_data.columns[1:].tolist()  # Skip 'Unnamed: 0' column
covar_matrix = covar_data.iloc[:, 1:].to_numpy().T  # Transpose so samples are rows

# Function to extract data for one gene
def extract_data_for_gene(row):

    logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force=True)

    try:
        # Extract gene information (line below for testing)
        # gene, chrom, cis_start, cis_end = cis_windows.iloc[0][['ensembl_gene_id', 'chromosome_name', 'cis_start', 'cis_end']].values
        gene, chrom, cis_start, cis_end = row
        output_prefix = f"{output_dir}/{gene}"
        logging.info(f"Processing gene: {gene}, chromosome: {chrom}, start: {cis_start}, end: {cis_end}")
        logging.info(f"Output dir: {output_prefix}")
        
        # Step 1: Generate .bed, .bim, .fam files for SNP metadata
        plink_bed_cmd = f"plink2 --pfile {geno_prefix} --chr {chrom} --from-bp {cis_start} --to-bp {cis_end} --make-bed --out {output_prefix}"
        os.system(plink_bed_cmd)
        
        # Step 2: Export genotype matrix as .raw file
        plink_raw_cmd = f"plink2 --pfile {geno_prefix} --chr {chrom} --from-bp {cis_start} --to-bp {cis_end} --export A-transpose --out {output_prefix}"
        os.system(plink_raw_cmd)
        
        # Read genotype data from .raw file
        try:
            raw_data = pd.read_csv(f"{output_prefix}.traw", sep="\t")
        except Exception as e:
            logging.error(f"Error reading {output_prefix}.raw: {e}")
            return None
        
        # Extract sample IDs and genotype matrix (notes: .traw has samples in columns and adds '0_' prefix to sample IDs)
        logging.info(f"Extracting genotype data for gene {gene}...")
        geno_samples = [s.replace('0_', '') for s in raw_data.columns[6:].tolist()]  # Skip CHR, SNP, POS, COUNTED, ALT, etc.
        geno_mat = raw_data.iloc[:, 6:].to_numpy()  # Variants x samples
        
        # Read SNP information from .bim file
        try:
            bim = pd.read_csv(f"{output_prefix}.bim", sep="\t", header=None,
                              names=["chr", "snp", "cm", "pos", "a1", "a2"])
        except Exception as e:
            logging.error(f"Error reading {output_prefix}.bim: {e}")
            return None
        
        # Load expression data
        logging.info(f"Loading expression data from {expr_file}")
        expr_data = pd.read_csv(expr_file, sep="\t")
        gene_expr = expr_data[expr_data["TargetID"] == gene].drop(columns=["#Chr", "start", "end", "TargetID"])
        expr_samples = gene_expr.columns.tolist()
        expr_vector = gene_expr.iloc[0].to_numpy()
        
        # Align samples
        logging.info(f"Checking if samples align for gene {gene}...")
        common_samples = list(set(geno_samples) & set(expr_samples) & set(covar_samples))
        if not common_samples:
            logging.warning(f"No common samples for gene {gene}")
            return None
        
        geno_idx = [geno_samples.index(s) for s in common_samples]
        expr_idx = [expr_samples.index(s) for s in common_samples]
        covar_idx = [covar_samples.index(s) for s in common_samples]
        
        geno_mat = geno_mat[geno_idx, :]
        expr_vector = expr_vector[expr_idx]
        covar_subset = covar_matrix[covar_idx, :]
        
        # Save data as RDS for SuSiE
        rds_file = f"{output_prefix}_data.rds"
        temp_geno = f"temp_geno_{gene}.tsv"
        temp_expr = f"temp_expr_{gene}.tsv"
        temp_covar = f"temp_covar_{gene}.tsv"
        temp_variants = f"temp_variants_{gene}.tsv"
        temp_r_script = f"temp_r_{gene}.R"

        # Save temporary TSVs
        pd.DataFrame(geno_mat).to_csv(temp_geno, sep="\t", index=False, header=False)
        pd.Series(expr_vector).to_csv(temp_expr, sep="\t", index=False, header=False)
        pd.DataFrame(covar_subset).to_csv(temp_covar, sep="\t", index=False, header=False)
        bim.to_csv(temp_variants, sep="\t", index=False, header=True)

        # Generate RDS with base R
        r_script = f"""
        data <- list(
        X = as.matrix(read.table('{temp_geno}', sep='\t')),
        y = as.numeric(read.table('{temp_expr}', sep='\t')[,1]),
        Z = as.matrix(read.table('{temp_covar}', sep='\t')),
        variants = read.table('{temp_variants}', sep='\t', header=TRUE)
        )
        saveRDS(data, '{rds_file}')
        """
        with open(temp_r_script, "w") as f:
            f.write(r_script)
        try:
            subprocess.run(["Rscript", temp_r_script], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Rscript failed for gene {gene}: {e.stderr}")
            raise

        # Clean up temporary files
        for temp_file in [temp_geno, temp_expr, temp_covar, temp_variants, temp_r_script]:
            try:
                os.remove(temp_file)
            except FileNotFoundError:
                pass
        for ext in [".traw", ".bed", ".bim", ".fam", ".log"]:
            try:
                os.remove(f"{output_prefix}{ext}")
            except FileNotFoundError:
                pass
    
    except Exception as e:
        logging.error(f"Error processing gene {gene}: {e}")
        return None

# Parallel processing with specified threads
logging.info(f"Starting parallel processing with {snakemake.threads} threads")
with Pool(processes=snakemake.threads) as pool:
    results = pool.map(extract_data_for_gene, cis_windows.values)

# Filter successful extractions
successful_genes = [r for r in results if r is not None]
logging.info(f"Processed {len(successful_genes)} genes")

# Create touch file for Snakemake
with open(output_file, 'w') as f:
    f.write(f"Processed {len(successful_genes)} genes for SuSiE fine-mapping.\n")
    f.write("Genes processed:\n")
    for gene in successful_genes:
        f.write(f"{gene}\n")