#--------------------------------------------------------------------------------------
#
#    Extract genotype and expression data for the Fig 3 eQTL boxplots
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: eqtl_boxplots_py (single job)
#            Upstream:   04 exclude_SNPs (genotype VCF); 05 tensorqtl_prep_input
#                        (normalised expression BED per cell type)
#            Downstream: replication_fig (Fig 3E and 3G)
#
# Purpose:   For each SNP-gene pair in config (one per figure panel):
#            - Extract genotypes (GT, DS) for the SNP from the VCF with bcftools
#            - Extract the gene's normalised expression for the cell type
#            - Merge on sample and write one CSV per pair
#            - Optionally draw a seaborn check plot next to the CSV (not
#              publication quality; off by default)
#
# Inputs:    vcf          Genotype VCF
#            exp_<panel>  Expression BED for the pair's cell type
#
# Outputs:   <panel>      eqtl_data_{cell_type}_{snp}_{gene}.csv
#
# Notes:     Cell types use pipeline labels. Check plots are written as
#            eqtl_plot_{cell_type}_{snp}_{gene}_categorical.png and are not
#            tracked by Snakemake.
#
#--------------------------------------------------------------------------------------

import io
import logging
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Extracting eQTL boxplot data for Fig 3 ...")

##  Load functions and variables  -----------------------------------------------------

def extract_genotypes(genotype_file, snp_id):
    """Extract genotype data for a specific SNP from a VCF file using bcftools."""
    command = [
        "bcftools", "query",
        "-f", "[%SAMPLE\t%CHROM\t%ID\t%REF\t%ALT\t%GT\t%DS\n]",
        "-i", f'ID=="{snp_id}"',
        genotype_file
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"Error running bcftools for {snp_id}: {result.stderr}")
        return None
    genotypes_df = pd.read_csv(io.StringIO(result.stdout), sep="\t", header=None,
                               names=["Sample", "Chr", "rsID", "REF", "ALT", "GT", "DS"])
    return genotypes_df


def extract_expression(expression_file, gene_ensembl):
    """Extract expression data for a specific gene from a BED file."""
    if not os.path.exists(expression_file):
        log.error(f"Expression file not found: {expression_file}")
        return None
    expression_df = pd.read_csv(expression_file, sep="\t")
    gene_row = expression_df[expression_df['TargetID'] == gene_ensembl]
    if gene_row.empty:
        log.error(f"Gene {gene_ensembl} not found in {expression_file}.")
        return None
    sample_columns = expression_df.columns[4:]
    expression_values = gene_row.iloc[0, 4:].values
    expression_data = pd.DataFrame({
        'Sample': sample_columns,
        'Expression': expression_values.astype(float)
    })
    return expression_data


def merge_data(genotypes_df, expression_data):
    """Merge genotype and expression DataFrames on Sample column."""
    genotypes_df['Sample'] = genotypes_df['Sample'].astype(str)
    expression_data['Sample'] = expression_data['Sample'].astype(str)
    merged_df = pd.merge(genotypes_df, expression_data, on="Sample", how="inner")
    if merged_df.empty:
        log.error("No overlapping samples between genotype and expression data.")
        return None
    return merged_df


def plot_eqtl(merged_df, cell_type, snp_id, gene_symbol, output_file):
    """Generate a publication-quality eQTL plot."""
    sns.set(style="white")
    plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size': 14})

    def genotype_category(gt):
        if gt == '0|0':
            return 'Hom_ref'
        elif gt in ['0|1', '1|0']:
            return 'Het'
        elif gt == '1|1':
            return 'Hom_alt'
        else:
            return 'Unknown'
    
    merged_df['genotype_cat'] = merged_df['GT'].apply(genotype_category)
    merged_df = merged_df[merged_df['genotype_cat'] != 'Unknown']
    
    ax = sns.boxplot(x='genotype_cat', y='Expression', data=merged_df,
                     order=['Hom_ref', 'Het', 'Hom_alt'],
                     boxprops={'facecolor': 'lightblue', 'edgecolor': 'black'},
                     whiskerprops={'color': 'black'},
                     capprops={'color': 'black'},
                     medianprops={'color': 'black'})
    sns.stripplot(x='genotype_cat', y='Expression', data=merged_df, color='black',
                  size=3, jitter=True, order=['Hom_ref', 'Het', 'Hom_alt'])
    
    counts = merged_df['genotype_cat'].value_counts()
    y_max = merged_df['Expression'].max()
    for i, cat in enumerate(['Hom_ref', 'Het', 'Hom_alt']):
        if cat in counts:
            ax.text(i, y_max * 0.95, f'n={counts[cat]}',
                    horizontalalignment='center', fontsize=12)
    
    plt.xlabel("Genotype", labelpad=15)
    plt.ylabel(f"{gene_symbol} Expression", labelpad=15)
    plt.title(f"{cell_type}: {snp_id} vs. {gene_symbol}", fontsize=18, pad=20)
    
    for spine in plt.gca().spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


# Pairs, inputs and outputs
pairs       = snakemake.params.pairs
plot_checks = bool(snakemake.params.plot_checks)
vcf_file    = snakemake.input.vcf

variables = {"pairs": ", ".join(pairs), "plot_checks": plot_checks, "vcf_file": vcf_file}
for panel in pairs:
    variables[f"exp_{panel}"] = snakemake.input[f"exp_{panel}"]
    variables[f"out_{panel}"] = snakemake.output[panel]
width = max(len(k) for k in variables)
log.info("Variables")
log.info("============================")
for k, v in variables.items():
    log.info(f"{k:<{width}}  {v}")
log.info("============================")

##  Extract, merge and write each pair  -----------------------------------------------
for panel, pair in pairs.items():
    cell_type    = pair["cell_type"]
    snp_id       = pair["snp"]
    gene_ensembl = pair["gene_id"]
    gene_symbol  = pair["gene"]
    log.info(f"{panel}: {cell_type} {snp_id} {gene_symbol} ({gene_ensembl})")

    genotypes_df = extract_genotypes(vcf_file, snp_id)
    if genotypes_df is None:
        sys.exit(f"No genotypes extracted for {snp_id}")

    expression_data = extract_expression(snakemake.input[f"exp_{panel}"], gene_ensembl)
    if expression_data is None:
        sys.exit(f"No expression extracted for {gene_ensembl} in {cell_type}")

    merged_df = merge_data(genotypes_df, expression_data)
    if merged_df is None:
        sys.exit(f"No overlapping samples for {panel}")

    merged_data_file = snakemake.output[panel]
    os.makedirs(os.path.dirname(merged_data_file), exist_ok=True)
    merged_df.to_csv(merged_data_file, index=False)
    log.info(f"Saved merged data ({len(merged_df)} samples): {merged_data_file}")

    if plot_checks:
        output_file = os.path.join(os.path.dirname(merged_data_file),
                                   f"eqtl_plot_{cell_type}_{snp_id}_{gene_symbol}_categorical.png")
        plot_eqtl(merged_df, cell_type, snp_id, gene_symbol, output_file)
        log.info(f"Generated check plot: {output_file}")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
