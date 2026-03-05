#--------------------------------------------------------------------------------------
#
#    Generate eQTL nominal file for data sharing
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

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

message("\n\neQTL supplementry table for manuscript ...")

# -------------------------------------------------------------------------------------

library(tidyverse)

allele_file <- snakemake@params[["allele_file"]]
in_dir <- snakemake@params[["in_dir"]]
out_file <- snakemake@output[[1]]
out_dir <- paste0(dirname(out_file), '/')

norm_method <- 'quantile'
genPC <- 4

cell_types <- c("Glu-UL", "Glu-DL", "GABA", "NPC",
                "MG", "OPC", "Endo-Peri",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2",
                "NPC-0", "NPC-1", "NPC-2")

expPC_map <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30,
  "Glu-UL-0"     = 30,
  "Glu-UL-1"     = 50,
  "Glu-UL-2"     = 40,
  "Glu-DL-0"     = 40,
  "Glu-DL-1"     = 40,
  "Glu-DL-2"     = 30,
  "GABA-0"       = 30,
  "GABA-1"       = 40,
  "GABA-2"       = 40,
  "NPC-0"        = 50,
  "NPC-1"        = 30,
  "NPC-2"        = 20
)

gene_lookup_file <- '../resources/sheets/gene_lookup_hg38.tsv'
gene_lookup_tbl <- suppressMessages(read_tsv(gene_lookup_file)) |>
  select(ensembl_gene_id, external_gene_name)

# ----- 1. Load Genotype Metadata (pvar) for alleles -----
message("Loading pvar file...")
pvar <- read_tsv(allele_file, comment = "#", 
                 col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")) %>%
  dplyr::select(CHROM, POS, ID, REF, ALT)

# ----- 2. Iterate through Cell Types and build tbl for each cell type -----
for (cell_type in names(expPC_map)) {
  
  message('Generating eQTL nominal file for: ', cell_type)
  
  expPC <- expPC_map[[cell_type]]
  
  # Path construction
  log_file <- paste0(in_dir, cell_type, '/', cell_type, '_nom.cis_qtl_pairs.tsv')
  
  if (!file.exists(log_file)) {
    warning("File not found for: ", cell_type)
    next
  }
  
  # Load all nominal eQTL 
  message('\nLoading eQTL nominal file for:', cell_type)
  eqtl_tbl <- read_tsv(log_file, show_col_types = FALSE) %>%
    dplyr::rename(ensembl_id = phenotype_id, SNP = variant_id)
  
  if (nrow(eqtl_tbl) == 0) next
  
  # Cross-ref with pvar for REF/ALT and Coordinates
  message('Adding alleles ...')
  eqtl_enriched <- eqtl_tbl %>%
    inner_join(pvar, by = c("SNP" = "ID"))
  
  message('Adding chr prefix to eQTL tbl ...')
  eqtl_enriched <- eqtl_enriched %>%
    mutate(CHROM = ifelse(!str_detect(CHROM, "chr"), paste0("chr", CHROM), CHROM)) %>% 
    drop_na(CHROM, POS)
  
  # NA check
  message('Checking for NAs ...')
  na_count <- sum(is.na(eqtl_enriched))
  if (na_count > 0) {
    message(paste0("Warning: ", na_count, " NA values detected in ", cell_type, 
                   ". This usually means SNP IDs don't match between TensorQTL and PVAR."))}
  
  # Add cell type label, gene symbol and reorder columns
  message('Munginging tbl ...')
  eqtl_enriched <- eqtl_enriched %>%
    inner_join(gene_lookup_tbl, by = join_by(ensembl_id == ensembl_gene_id)) |>
    mutate(CHROM = str_remove(CHROM, "^chr")) |>
    dplyr::select(ensembl_id, symbol = external_gene_name, CHROM, 
                  SNP, POS, REF, ALT,  AF = af, slope, slope_se, pval_nominal)
  
  message('Any NAs in final tbl?', anyNA(eqtl_enriched))
  message("Processed ", cell_type, ": ", nrow(eqtl_enriched), " eQTLs")
  
  message('Writing file ...')
  write_tsv(eqtl_enriched, paste0(out_dir, cell_type, '_cis_eQTL_nominal.tsv.gz'))
  
}

message("Export Complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------