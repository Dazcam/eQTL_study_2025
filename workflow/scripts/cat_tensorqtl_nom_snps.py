#!/usr/bin/env python3

"""
Concatenate and Summarize TensorQTL Nominal Output

This script reads TensorQTL nominal pass output files (Parquet format, one per chromosome)
for a specified cell type, concatenates them into a single TSV file for input to SMR
preprocessing, and summarizes the number of SNPs per gene, total SNPs, and total genes.
The summary verifies that the nominal output contains multiple SNPs per gene for SMR's
HEIDI test.

Arguments:
    nom_dir (str): Directory containing TensorQTL nominal output Parquet files
                   (e.g., "../results/04TENSORQTL/ExN-1").
    cell_type (str): Cell type name (e.g., "ExN-1").
    concat_out (str): Path to output concatenated TSV file.
    summary_out (str): Path to output summary TSV file.

Output:
    - concat_out: TSV with columns: phenotype_id, variant_id, slope, slope_se, pval_nominal, af.
    - summary_out: TSV with columns: Cell_Type, Total_SNPs, Total_Genes, SNPs_Per_Gene_Mean,
      SNPs_Per_Gene_Median, SNPs_Per_Gene_Min, SNPs_Per_Gene_Max.
    - summary_out_detailed.tsv: Detailed TSV with phenotype_id and n_snps.

Details:
    Assumes Parquet files are named as `<cell_type>_nom.cis_qtl_pairs.<chr>.parquet`
    for chromosomes 1-22. Uses pandas and pyarrow for efficient Parquet processing.
    The concatenated TSV is compatible with prep_eQTL_for_smr.R for SMR input preparation.

Author: Darren Cameron
Date: 2025-07-29

Note:
    Run within a Snakemake pipeline or standalone. Requires pandas and pyarrow.
    Ensure Parquet files cover all chromosomes (1-22).

Example:
    python concat_tensorqtl_nom_snps.py \
        --nom_dir ../results/04TENSORQTL/ExN-1 \
        --cell_type ExN-1 \
        --concat_out ../results/04TENSORQTL/ExN-1/ExN-1_nom.cis_qtl_pairs.tsv \
        --summary_out ../results/04TENSORQTL/ExN-1/snp_summary.tsv
"""

import pandas as pd
import glob
import os
import argparse
import logging

def concat_tensorqtl_nom_snps(nom_dir, cell_type, concat_out, summary_out):
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # List Parquet files for chromosomes 1-22
    parquet_files = sorted(glob.glob(os.path.join(nom_dir, f"{cell_type}_nom.cis_qtl_pairs.[1-9]*.parquet")))
    logger.info(f"Found {len(parquet_files)} Parquet files for {cell_type}")
    
    if not parquet_files:
        raise FileNotFoundError(f"No Parquet files found in {nom_dir}")
    
    # Read and concatenate Parquet files
    logger.info("Reading and concatenating Parquet files...")
    snp_data = pd.concat([
        pd.read_parquet(file)[['phenotype_id', 'variant_id', 'slope', 'slope_se', 'pval_nominal', 'af']]
        for file in parquet_files
    ], ignore_index=True)
    
    # Write concatenated TSV
    logger.info(f"Writing concatenated TSV to: {concat_out}")
    snp_data.to_csv(concat_out, sep='\t', index=False)
    
    # Summarize SNPs per gene
    logger.info("Summarizing SNPs per gene...")
    snp_summary = snp_data.groupby('phenotype_id').size().reset_index(name='n_snps')
    
    # Write detailed SNPs per gene
    detailed_out = summary_out.replace('.tsv', '_detailed.tsv')
    snp_summary.to_csv(detailed_out, sep='\t', index=False)
    logger.info(f"Detailed SNPs per gene written to: {detailed_out}")
    
    # Compute summary statistics
    total_snps = len(snp_data)
    total_genes = snp_data['phenotype_id'].nunique()
    mean_snps = snp_summary['n_snps'].mean()
    median_snps = snp_summary['n_snps'].median()
    min_snps = snp_summary['n_snps'].min()
    max_snps = snp_summary['n_snps'].max()
    
    # Create summary table
    summary_table = pd.DataFrame({
        'Cell_Type': [cell_type],
        'Total_SNPs': [total_snps],
        'Total_Genes': [total_genes],
        'SNPs_Per_Gene_Mean': [mean_snps],
        'SNPs_Per_Gene_Median': [median_snps],
        'SNPs_Per_Gene_Min': [min_snps],
        'SNPs_Per_Gene_Max': [max_snps]
    })
    
    # Write summary to TSV
    logger.info(f"Writing summary to: {summary_out}")
    summary_table.to_csv(summary_out, sep='\t', index=False)
    
    return summary_table

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate and summarize TensorQTL nominal output")
    parser.add_argument("--nom_dir", required=True, help="Directory with TensorQTL Parquet files")
    parser.add_argument("--cell_type", required=True, help="Cell type name")
    parser.add_argument("--concat_out", required=True, help="Output concatenated TSV file")
    parser.add_argument("--summary_out", required=True, help="Output summary TSV file")
    args = parser.parse_args()
    
    concat_tensorqtl_nom_snps(args.nom_dir, args.cell_type, args.concat_out, args.summary_out)