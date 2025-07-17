#--------------------------------------------------------------------------------------
#
#    Munge expression and gene coord files for FUSION 
#
#--------------------------------------------------------------------------------------

# Option to run with subset of genes here for testing (L35-37)

## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nCreating plink ready expression data and gene coordinate file for FUSION weight calculation ...')

# Load required libraries
library(tidyverse)
library(biomaRt)

# Input and output paths
expr_in <- snakemake@input[[1]]
expr_out <- snakemake@output[['exp']]
coord_out <- snakemake@output[['coord']]
run_test <- TRUE # Add this to snakerule later to set dynamically

# Read expression file
message('Reading expression file ...')
expr_tbl <- read_tsv(expr_in, col_types=cols())

if (run_test) {
  message('Running with only a subset of genes for testing ...')
  message('Genes before: ', nrow(expr_tbl) - 1)
  expr_tbl <- head(expr_tbl)
  message('Genes after: ', nrow(expr_tbl) - 1)
}

# Extract sample IDs and expression data
coord_cols <- c("#Chr", "start", "end", "TargetID")
sample_cols <- setdiff(colnames(expr_tbl), coord_cols)

# Convert expression data to PLINK format
message('Converting expression data to PLINK format ...')
expr_plink_tbl <- expr_tbl %>%
  dplyr::select(TargetID, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "IID", values_to = "value") %>%
  pivot_wider(names_from = TargetID, values_from = value) %>%
  mutate(FID = 0) %>%
  dplyr::select(FID, IID, everything())

# Save PLINK phenotype file
message('Saving PLINK phenotype file ...')
write_tsv(expr_plink_tbl, expr_out)

# Connect to Ensembl BioMart
message('Setting BiomaRT mart ...')
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")

# Query gene coordinates
message('Get coords for genes ...')
gene_ids <- expr_tbl$TargetID
gene_coord <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
) %>%
  as_tibble() %>%
  rename(
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    gene_id = ensembl_gene_id
  ) %>%
  mutate(chr = str_remove(chr, "^chr")) %>%
  dplyr::select(chr, start, end, gene_id)

# Check for unmatched genes
missing_genes <- setdiff(gene_ids, gene_coord$gene_id)
if (length(missing_genes) > 0) {
  warning(str_glue("{length(missing_genes)} genes not found in BioMart: {str_c(head(missing_genes, 5), collapse = ', ')}..."))
}

# Save gene coordinate file
message('Saving gene coordinate file ...')
write_tsv(gene_coord, coord_out)

message(str_glue("Generated {expr_out} with {nrow(expression)} samples and {length(sample_cols)} genes\n"))
message(str_glue("Generated {coord_out} with {nrow(gene_coord)} genes\n"))