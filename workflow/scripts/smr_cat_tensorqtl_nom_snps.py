#!/usr/bin/env python3
#--------------------------------------------------------------------------------------
#
#    Concatenate tensorQTL nominal output for SMR
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  10SMR | Rule: cat_tensorqtl_nom_snps (one job per cell type)
#            Upstream:   05TENSORQTL tensorqtl_nom (per-chromosome parquet files)
#            Downstream: create_query
#
# Purpose:   For one cell type:
#            - Read the tensorQTL nominal parquet files (chr 1-22) and keep
#              phenotype_id, variant_id, slope, slope_se, pval_nominal, af
#            - Write them as one TSV
#            - Summarise SNPs per gene (SMR's HEIDI test needs several per gene)
#
# Inputs:    --nom_dir      tensorQTL nominal output directory for the cell type
#            --cell_type    cell type name (used in the parquet file names)
#            --norm_method  tensorQTL normalisation method (used in the parquet file
#                           names; config tensorQTL: norm_methods[0])
#
# Outputs:   --concat_out   {cell_type}_nom.cis_qtl_pairs.tsv
#            --summary_out  {cell_type}_nom.cis_qtl_pairs_snp_summary.tsv: total SNPs,
#                           genes, and mean / median / min / max SNPs per gene
#            Also written:  <summary_out>_detailed.tsv: SNPs per gene
#
#--------------------------------------------------------------------------------------

import argparse
import glob
import logging
import os
import sys

import pandas as pd

## Set up logging  --------------------------------------------------------------------
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)


def concat_tensorqtl_nom_snps(nom_dir, cell_type, norm_method, concat_out, summary_out):

    ## Read and concatenate parquet files  --------------------------------------------
    # List Parquet files for chromosomes 1-22
    parquet_files = sorted(glob.glob(os.path.join(nom_dir, f"{cell_type}_{norm_method}_nom.cis_qtl_pairs.[1-9]*.parquet")))
    log.info(f"Found {len(parquet_files)} Parquet files for {cell_type}")

    if not parquet_files:
        raise FileNotFoundError(f"No Parquet files found in {nom_dir}")

    log.info("Reading and concatenating Parquet files...")
    snp_data = pd.concat([
        pd.read_parquet(file)[['phenotype_id', 'variant_id', 'slope', 'slope_se', 'pval_nominal', 'af']]
        for file in parquet_files
    ], ignore_index=True)

    log.info(f"Writing concatenated TSV to: {concat_out}")
    snp_data.to_csv(concat_out, sep='\t', index=False)

    ## Summarise SNPs per gene  -------------------------------------------------------
    log.info("Summarizing SNPs per gene...")
    snp_summary = snp_data.groupby('phenotype_id').size().reset_index(name='n_snps')

    # Write detailed SNPs per gene
    detailed_out = summary_out.replace('.tsv', '_detailed.tsv')
    snp_summary.to_csv(detailed_out, sep='\t', index=False)
    log.info(f"Detailed SNPs per gene written to: {detailed_out}")

    # Compute summary statistics
    total_snps = len(snp_data)
    total_genes = snp_data['phenotype_id'].nunique()
    mean_snps = snp_summary['n_snps'].mean()
    median_snps = snp_summary['n_snps'].median()
    min_snps = snp_summary['n_snps'].min()
    max_snps = snp_summary['n_snps'].max()

    summary_table = pd.DataFrame({
        'Cell_Type': [cell_type],
        'Total_SNPs': [total_snps],
        'Total_Genes': [total_genes],
        'SNPs_Per_Gene_Mean': [mean_snps],
        'SNPs_Per_Gene_Median': [median_snps],
        'SNPs_Per_Gene_Min': [min_snps],
        'SNPs_Per_Gene_Max': [max_snps]
    })

    log.info(f"Writing summary to: {summary_out}")
    summary_table.to_csv(summary_out, sep='\t', index=False)
    log.info("Done.")

    return summary_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate and summarize TensorQTL nominal output")
    parser.add_argument("--nom_dir", required=True, help="Directory with TensorQTL Parquet files")
    parser.add_argument("--cell_type", required=True, help="Cell type name")
    parser.add_argument("--norm_method", required=True, help="tensorQTL normalisation method")
    parser.add_argument("--concat_out", required=True, help="Output concatenated TSV file")
    parser.add_argument("--summary_out", required=True, help="Output summary TSV file")
    args = parser.parse_args()

    log.info("Concatenating tensorQTL nominal output for SMR ...")

    # Print what each variable is set to
    variables = vars(args)
    width = max(len(k) for k in variables)
    log.info("Variables")
    log.info("============================")
    for name, value in variables.items():
        log.info(f"{name:<{width}}  {value}")
    log.info("============================")

    concat_tensorqtl_nom_snps(args.nom_dir, args.cell_type, args.norm_method,
                              args.concat_out, args.summary_out)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
