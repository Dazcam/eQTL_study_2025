#--------------------------------------------------------------------------------------
#
#    Generate eQTL supplementary table for manuscript
#
#--------------------------------------------------------------------------------------

# For each cell type extract sig. eQTL
# Add chr, position and allele info from genotype VCF
# Check for OCR overlap of each eQTL in Ziffra union of OCRs file

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
library(readxl)
library(GenomicRanges)
library(openxlsx)

peak_dir <- snakemake@params[["peak_dir"]]
allele_file <- snakemake@params[["allele_file"]]
in_dir <- snakemake@params[["in_dir"]]
out_file <- snakemake@output[[1]]

norm_method <- 'quantile'
genPC <- 4

expPC_map <- unlist(snakemake@params[["exp_pc_map"]])
cell_types <- names(expPC_map)

gene_lookup_file <- '../resources/sheets/gene_lookup_hg38.tsv'
gene_lookup_tbl <- suppressMessages(read_tsv(gene_lookup_file)) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ----- 1. Load ATAC-seq Peaks (Ziffra 2021) -----
message("Loading Ziffra peak data...")
peak_file <- file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx")

# Load the primary peak coordinates
peak_coords <- read_excel(peak_file, sheet = "ST2 AllPrimaryPeaks") %>%
  dplyr::select(seqnames, start, end, peak_name)

# Create a master GRanges object for fast overlap checking
master_peaks_gr <- makeGRangesFromDataFrame(peak_coords, keep.extra.columns = TRUE)

# ----- 2. Load Genotype Metadata (pvar) for alleles -----
message("Loading pvar file...")
pvar <- read_tsv(allele_file, comment = "#", 
                 col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")) %>%
  dplyr::select(CHROM, POS, ID, REF, ALT) %>%
  distinct(ID, .keep_all = TRUE)

# ----- 3. Iterate through Cell Types and Build Table -----
final_excel_list <- list()

for (cell_type in names(expPC_map)) {
  expPC <- expPC_map[[cell_type]]
  
  # Path construction
  log_dir <- paste0(in_dir, cell_type, '_', norm_method, '_genPC_', genPC, '_expPC_', expPC, '/')
  log_file <- paste0(log_dir, cell_type, '_', norm_method, '_perm.cis_qtl.txt.gz')
  
  if (!file.exists(log_file)) {
    warning("File not found for: ", cell_type)
    next
  }
  
  # Load significant eQTL
  message('\nLoading eQTL perm file for:', cell_type)
  eqtl_tbl <- read_tsv(log_file, show_col_types = FALSE) %>%
    dplyr::filter(qval < 0.05) %>%
    dplyr::select(ensembl_id = phenotype_id, SNP = variant_id, af, 
                  pval_nominal, pval_beta, qval, slope, slope_se)
  
  if (nrow(eqtl_tbl) == 0) next
  
  # Cross-reference with pvar for REF/ALT and Coordinates
  message('Adding alleles ...')
  n_before <- nrow(eqtl_tbl)
  eqtl_enriched <- eqtl_tbl %>%
    inner_join(pvar, by = c("SNP" = "ID"))
  
  if (nrow(eqtl_enriched) != n_before) {
    warning(cell_type, ": row count changed after pvar join (", n_before, " -> ", nrow(eqtl_enriched), "), check for duplicate SNP IDs")
  }  


  message('Adding chr prefix to eQTL tbl ...')
  eqtl_enriched <- eqtl_enriched %>%
    mutate(CHROM = ifelse(!str_detect(CHROM, "chr"), paste0("chr", CHROM), CHROM)) %>% 
    drop_na(CHROM, POS)
  
  # Define SNP positions as GRanges for overlap check
  message('Converting to GRanges ...')
  snp_gr <- GRanges(
    seqnames = eqtl_enriched$CHROM,
    ranges = IRanges(start = eqtl_enriched$POS, end = eqtl_enriched$POS)
  )
  
  # Identify SNPs overlapping any peak
  # countOverlaps > 0 returns TRUE if hit, FALSE otherwise
  message('Checking for SNP overlaps in snATACseq peaks ...')
  eqtl_enriched$In_Peak <- countOverlaps(snp_gr, master_peaks_gr) > 0
  
  # NA check
  message('Checking for NAs ...')
  na_count <- sum(is.na(eqtl_enriched))
  if (na_count > 0) {
    message(paste0("Warning: ", na_count, " NA values detected in ", cell_type, 
                   ". This usually means SNP IDs don't match between TensorQTL and PVAR."))
    # Optional: eqtl_enriched <- drop_na(eqtl_enriched) 
  }
  
  # Add cell type label, gene symbol and reorder columns
  message('Mungingfi tbl ...')
  n_before <- nrow(eqtl_tbl)
  eqtl_enriched <- eqtl_enriched %>%
    mutate(cell_type = !!cell_type) %>%
    inner_join(gene_lookup_tbl, by = join_by(ensembl_id == ensembl_gene_id)) |>
    mutate(CHROM = str_remove(CHROM, "^chr")) |>
    dplyr::select(cell_type, ensembl_id, symbol = external_gene_name, CHROM, 
                  SNP, POS, REF, ALT,  AF = af, slope, slope_se, pval_nominal, 
                  pval_beta, qval, In_Peak)
   
  if (nrow(eqtl_enriched) != n_before) {
    warning(cell_type, ": row count changed after gene lookup join (", n_before, " -> ", nrow(eqtl_enriched), "), check for duplicate Gene IDs")
  }
 
  message('Any NAs in final tbl?', anyNA(eqtl_enriched))
  
  final_excel_list[[cell_type]] <- eqtl_enriched
  message("Processed ", cell_type, ": ", nrow(eqtl_enriched), " eQTLs")
}

# Combine and Export
message("Writing to Excel file...")
write.xlsx(final_excel_list, 
           file = out_file, 
           overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))

message("Export Complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
