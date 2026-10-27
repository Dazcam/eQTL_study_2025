#--------------------------------------------------------------------------------------
#
#    Prep eQTL for SMR
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  10SMR | Rule: create_query (one job per cell type)
#            Upstream:   cat_tensorqtl_nom_snps (tensorQTL nominal pairs); cat_refs
#                        (1000G EUR hg38 .bim / .frq); 08SUSIE prep_susie_gene_meta
#            Downstream: create_besd
#
# Purpose:   Build the SMR query file for one cell type:
#            - Take SNP, gene, b, se, p and af from the tensorQTL nominal pairs
#            - Optionally drop extended MHC genes (remove_mhc in config)
#            - Add SNP position and alleles (A1, A2) from the 1000G .bim
#            - Set Freq to the 1000G frequency of A1 (align_freq_to_ref())
#            - Add gene position and strand; drop rows with any missing value
#
# Inputs:    eqtl         {cell_type}_nom.cis_qtl_pairs.tsv (concatenated nominal pairs)
#            snps         1000G.EUR.hg38.bim
#            genes        SuSiE gene metadata (chromosome, phenotype_pos, strand)
#            frq          1000G.EUR.hg38.frq
#            gene_lookup  hg38 gene lookup (positions for the xMHC filter)
#            remove_mhc   TRUE/FALSE: drop xMHC genes; region from config mhc_region
#
# Outputs:   query        {cell_type}_query.tsv (SMR --qfile format)
#            gene_lst     {cell_type}_gene_list.tsv (genes in the query file)
#
#--------------------------------------------------------------------------------------

## Set up logging for Snakemake
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") && length(snakemake@log) != 0) {
      log <- file(snakemake@log[[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
  log_smk()
}

message('\n\nMunge eQTL, snp and gene level data for SMR input ...')

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "functions.R"))
source(file.path(snakemake@scriptdir, "smr_functions.R"))

# Input and output paths
eqtl_in <- snakemake@input[['eqtl']]
snps_in <- snakemake@input[['snps']]
genes_in <- snakemake@input[['genes']]
frq_in <- snakemake@input[['frq']]
lookup_in <- snakemake@input[['gene_lookup']]
query_out <- snakemake@output[['query']]
genes_out <- snakemake@output[['gene_lst']]
remove_mhc <- isTRUE(snakemake@params[['remove_mhc']])
mhc_chr <- snakemake@params[['mhc_chr']]
mhc_start <- snakemake@params[['mhc_start']]
mhc_end <- snakemake@params[['mhc_end']]
# fdr_thresh <- snakemake@params[['qval_thresh']] # Too stringent for SMR / HEIDI

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_in", "snps_in", "genes_in", "frq_in", "lookup_in", "query_out",
               "genes_out", "remove_mhc", "mhc_region"),
  value    = c(eqtl_in, snps_in, genes_in, frq_in, lookup_in, query_out,
               genes_out, remove_mhc, paste0("chr", mhc_chr, ":", mhc_start, "-", mhc_end))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

## Load data  -------------------------------------------------------------------------
message('Loading data ...')
eqtl_tbl <- read_tsv(eqtl_in)
snps_tbl <- read_tsv(snps_in, col_names = c('Chr', 'SNP', 'CM', 'BP', 'A1', 'A2'))
genes_tbl <- read_tsv(genes_in) |>
  rename(Probe_Chr = chromosome,
         Orientation = strand,
         Probe_bp = phenotype_pos, 
         Gene = gene_id) |>
  select(Probe_Chr, Probe_bp, Gene, Orientation)
frq_tbl <- read_table(frq_in, col_types = cols(.default = "c", MAF = "d", NCHROBS = "i"))

## Build query table  -----------------------------------------------------------------
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

# Optionally drop extended MHC genes (config: smr: create_query: remove_mhc)
if (remove_mhc) {
  message('Removing xMHC genes (chr', mhc_chr, ':', mhc_start, '-', mhc_end, ') ...')
  lookup_tbl <- read_tsv(lookup_in, show_col_types = FALSE)
  mhc_genes <- get_mhc_genes(lookup_tbl, chr = mhc_chr, start = mhc_start, end = mhc_end)
  removed <- intersect(unique(eqtl_filt_tbl$Gene), mhc_genes)
  eqtl_filt_tbl <- eqtl_filt_tbl |> dplyr::filter(!Gene %in% mhc_genes)
  message('Removed ', length(removed), ' xMHC genes: ',
          paste(sort(removed), collapse = ', '))
} else {
  message('xMHC genes kept (remove_mhc: FALSE)')
}

# Merge with SNP position data from 100G hg38 .bim file
message('Merging SNP data ...')
snp_mrg_tbl <- eqtl_filt_tbl |>
  left_join(snps_tbl, by = "SNP") |>
  select(SNP, Chr, BP, A1, A2, Freq, Gene, b, se, p) |>
  mutate(Chr = str_replace(Chr, '^chr', ''))

# Align eQTL alleles and frequencies to 1000G reference
message('Adjusting allele frequencies so Freq reflects frequency of current A1 (effect allele) ...')
snp_mrg_tbl <- align_freq_to_ref(snp_mrg_tbl, frq_tbl, by = "SNP", freq_col = "Freq",
                                 keep_unmatched = TRUE, label = "eQTL")

# Merge with gene meta data
message('Merging Gene data ...')
smr_tbl <- snp_mrg_tbl |>
  left_join(genes_tbl, by = "Gene") |>
  mutate(Probe = Gene) |>
  select(SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr,
           Probe_bp, Gene, Orientation, b, se, p)

## Checks  ----------------------------------------------------------------------------
message('Running NA checks ...')
if (anyNA(smr_tbl)) {
  incomplete_n <- smr_tbl %>% filter(if_any(everything(), is.na)) %>% nrow()
  
  message("\n=== Rows with missing values (n = ", incomplete_n, ") ===")
  
  smr_tbl %>% 
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "column", values_to = "n_NA") %>%
    arrange(desc(n_NA)) %>%
    mutate(n_NA = scales::comma(n_NA)) %>%
    print(n = Inf)                     # ← this is the magic line
  
  message("======================================================\n")
  
  message("Dropping rows with any NA ...")
  smr_tbl <- smr_tbl %>% drop_na()
} else {
  message("No missing values in the final table.\n")
}
message('Do SNP and Gene Chr annotations match? ', 
        sum(smr_tbl$Chr == smr_tbl$Probe_Chr) == nrow(smr_tbl))

## Write outputs  ---------------------------------------------------------------------
message('Writing files ...')
write_tsv(smr_tbl, query_out)
genes_tbl <- smr_tbl |>
  select(Gene) |>
  distinct()
write_tsv(genes_tbl, genes_out)
message('Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
