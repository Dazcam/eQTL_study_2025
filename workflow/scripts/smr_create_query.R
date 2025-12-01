#--------------------------------------------------------------------------------------
#
#    Prep eQTL for SMR
#
#--------------------------------------------------------------------------------------

# Inputs
# - TensorQTL nominal eQTL summary statistics
# - SNP positions and frequencies from the hg38 build reference 
# - Gene metadata (generated in susie pipeline; prep_susie_gene_meta)

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
frq_in <- snakemake@input[['frq']]
query_out <- snakemake@output[['query']]
genes_out <- snakemake@output[['gene_lst']]
# fdr_thresh <- snakemake@params[['qval_thresh']] # Too stringent for SMR / HEIDI

# Read in data
message("\nLoading data ...\n")
message("eQTL loaded from: ", eqtl_in)
message("SNP ref info loaded from: ", snps_in)
message("Gene metadata loaded from: ", genes_in)
message("Allele frq info loaded from: ", frq_in)

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
frq_tbl <- read_table(frq_in, col_types = cols(.default = "c", MAF = "d")) |>
  select(CHR, SNP, A1, A2, MAF) |>
  rename(A1.ref = A1, A2.ref = A2)

# Filter for significant eQTLs (qval < threshold)
message('Extracting eQTL cols ...')
eqtl_filt_tbl <- eqtl_tbl  |>
  #  filter(qval < fdr_thresh) %>% # Too stringent; HEIDI fails for all SNPs
  #  filter(pval_nominal < 0.05) %>%  # Could use nominal p-value to include more SNPs
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

# Align eQTL alleles and frequencies to 1000G reference
message('Adjusting allele frequencies so Freq reflects frequency of current A1 (effect allele) ...')
snp_mrg_tbl <- snp_mrg_tbl %>%
  left_join(frq_tbl, by = c("Chr" = "CHR", "SNP")) %>%   # brings A1.ref, A2.ref, MAF
  mutate(
    allele_status = case_when(
      A1 == A1.ref & A2 == A2.ref ~ "aligned",
      A1 == A2.ref & A2 == A1.ref ~ "flipped",
      is.na(A1.ref)               ~ "missing_ref",
      TRUE                        ~ "mismatched"
    ),
    
    Freq = case_when(
      allele_status == "aligned"       ~ MAF,
      allele_status == "flipped"       ~ 1 - MAF,
      allele_status == "missing_ref"   ~ Freq,   # fall back to TensorQTL's own af
      TRUE                             ~ NA_real_
    )
  ) %>%
  select(-A1.ref, -A2.ref, -MAF, -allele_status) %>%   # drop helper columns
  filter(!is.na(Freq))

# Merge with gene meta data
message('Merging Gene data ...')
smr_tbl <- snp_mrg_tbl |>
  left_join(genes_tbl, by = "Gene") |>
  mutate(Probe = Gene) |>
  relocate(SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr,
           Probe_bp, Gene, Orientation, b, se, p)

message('Running checks ...')
if (anyNA(smr_tbl)) {
  message(sum(!complete.cases(smr_tbl)), " NAs in final df — rm rows with NA ...")
  smr_tbl <- smr_tbl |>
    drop_na()
} else {
  message("No NAs in final df — no action taken.")
}
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
