#--------------------------------------------------------------------------------------
#
#    Prep eQTL for SMR
#
#--------------------------------------------------------------------------------------

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
message('\n\nMunge eQTL, snp and gene level data for SMR input ...')

# Load packages
library(tidyverse)

# Input and output paths
eqtl_in <- snakemake@input[['eqtl']]
snps_in <- snakemake@input[['snps']]
genes_in <- snakemake@input[['genes']]
query_out <- snakemake@output[['query']]
genes_out <- snakemake@output[['gene_lst']]
fdr_thresh <- snakemake@params[['qval_thresh']]

# test code
# eqtl_tbl <- read_tsv('~/Desktop/test/smr_test/ExN-1_perm.cis_qtl.txt.gz')
# genes_tbl <- read_tsv('~/Desktop/test/smr_test/ExN-1_gene_meta.tsv') |>
#   rename(Probe_Chr = chromosome,
#          Orientation = strand,
#          Probe_bp = phenotype_pos, 
#          Gene = gene_id) |>
#   select(Probe_Chr, Probe_bp, Gene, Orientation)
# snps_tbl <- read_tsv('~/Desktop/test/smr_test/snp_positions.txt', 
#                       col_names = c('Chr', 'BP', 'SNP', 'A1', 'A2'))

# Read data
message('Loading data ...')
eqtl_tbl <- read_tsv(eqtl_in)
snps_tbl <- read_tsv(snps_in, col_names = c('Chr', 'BP', 'SNP', 'A1', 'A2'))
genes_tbl <- read_tsv(genes_in) |>
  rename(Probe_Chr = chromosome,
         Orientation = strand,
         Probe_bp = phenotype_pos, 
         Gene = gene_id) |>
  select(Probe_Chr, Probe_bp, Gene, Orientation)

# Filter for significant eQTLs (qval < threshold)
message('Filtering eQTL. FDR thresh set to: ', fdr_thresh)
eqtl_filt_tbl <- eqtl_tbl  %>%
  filter(qval < fdr_thresh) %>%
  select(
    SNP = variant_id,
    Gene = phenotype_id,
    b = slope,
    se = slope_se,
    p = pval_nominal, # nominal or qval?
    Freq = af
  )

# Merge with SNP position data
message('Merging SNP data ...')
snp_mrg_tbl <- eqtl_filt_tbl |>
  left_join(snps_tbl, by = "SNP") |>
  select(SNP, Chr, BP, A1, A2, Freq, Gene, b, se, p) |>
  mutate(Chr = str_replace(Chr, '^chr', ''))

# Merge with gene meta data
message('Merging Gene data ...')
smr_tbl <- snp_mrg_tbl |>
  left_join(genes_tbl, by = "Gene") |>
  mutate(Probe = Gene) |>
  relocate(SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr,
          Probe_bp, Gene, Orientation, b, se, p)

message('Running checks ...')
message(nrow(smr_tbl), ' eQTL after filtering and joining ...')
message('NAs in final dataframe? ', anyNA(smr_tbl))
message('Do SNP and Gene Chr annotations match? ', 
        sum(smr_tbl$Chr == smr_tbl$Probe_Chr) == nrow(smr_tbl))

# Write smr query file and gene list
message('Writing files ...')
write_tsv(smr_tbl, query_out)
genes_tbl <- smr_tbl |>
  select(Gene) |>
  distinct()
write_tsv(genes_tbl, genes_out)
message('Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
