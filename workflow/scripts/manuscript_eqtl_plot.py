import argparse
import subprocess
import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt
import os


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
        print(f"Error running bcftools for {snp_id}: {result.stderr}")
        return None
    genotypes_df = pd.read_csv(io.StringIO(result.stdout), sep="\t", header=None,
                               names=["Sample", "Chr", "rsID", "REF", "ALT", "GT", "DS"])
    return genotypes_df


def extract_expression(expression_file, gene_ensembl):
    """Extract expression data for a specific gene from a BED file."""
    if not os.path.exists(expression_file):
        print(f"Expression file not found: {expression_file}")
        return None
    expression_df = pd.read_csv(expression_file, sep="\t")
    gene_row = expression_df[expression_df['TargetID'] == gene_ensembl]
    if gene_row.empty:
        print(f"Gene {gene_ensembl} not found in {expression_file}.")
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
        print("No overlapping samples between genotype and expression data.")
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


def main():
    """Main function to generate eQTL plots from a pairs CSV file."""
    parser = argparse.ArgumentParser(description="Generate eQTL plots from a pairs CSV file.")
    parser.add_argument("--pairs_file", required=True, help="Path to CSV file with cell_type, variant, ensembl_id, symbol")
    parser.add_argument("--genotype_file", required=True, help="Path to VCF genotype file")
    parser.add_argument("--expression_dir", required=True, help="Directory containing {cell_type}_tmm.bed files")
    parser.add_argument("--output_dir", required=True, help="Directory to save the plots")
    
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Read pairs CSV
    pairs_df = pd.read_csv(args.pairs_file)
    required_columns = {'cell_type', 'variant', 'ensembl_id'}  # symbol is now optional
    if not required_columns.issubset(pairs_df.columns):
        raise ValueError("CSV must contain columns: cell_type, variant, ensembl_id (symbol is optional)")

    # Generate plots for each row
    for _, row in pairs_df.iterrows():
        cell_type = row['cell_type']
        snp_id = row['variant']
        gene_ensembl = row['ensembl_id']
        # Use symbol if present and non-empty, otherwise fall back to ensembl_id
        gene_symbol = row.get('symbol', gene_ensembl) if pd.notna(row.get('symbol')) else gene_ensembl
        
        # Construct expression file path using cell_type and expression_dir
        expression_file = os.path.join(args.expression_dir, f"{cell_type}_quantile.bed")
        
        # Extract and process data
        genotypes_df = extract_genotypes(args.genotype_file, snp_id)
        if genotypes_df is None:
            continue
        
        expression_data = extract_expression(expression_file, gene_ensembl)
        if expression_data is None:
            continue
        
        merged_df = merge_data(genotypes_df, expression_data)
        if merged_df is None:
            continue
        
        # Save merged data to file in output directory
        merged_data_file = os.path.join(args.output_dir, f"eqtl_data_{cell_type}_{snp_id}_{gene_symbol}.csv")
        merged_df.to_csv(merged_data_file, index=False)
        print(f"Saved merged data: {merged_data_file}")

        # Generate plot with cell_type in filename
        output_file = os.path.join(args.output_dir, f"eqtl_plot_{cell_type}_{snp_id}_{gene_symbol}_categorical.png")
        plot_eqtl(merged_df, cell_type, snp_id, gene_symbol, output_file)
        print(f"Generated plot: {output_file}")

    # Create empty done file for Snakemake
    done_file = os.path.join(args.output_dir, "plot_qtl.done")
    with open(done_file, 'w') as f:
        pass
    print(f"Created done file: {done_file}")


if __name__ == "__main__":
    main()

