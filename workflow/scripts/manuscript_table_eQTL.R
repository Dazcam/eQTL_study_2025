#--------------------------------------------------------------------------------------
#
#    eQTL supplementary table (Supp Table 2)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: eqtl_tbl
#            Upstream:   05 tensorqtl_perm (quantile, genPC 4, expPC from exp_pc_map);
#                        05 prep_tensorQTL_input (pvar); Ziffra 2021 snATAC-seq peaks
#            Downstream: Supp Table 2
#
# Purpose:   For each cell type in exp_pc_map:
#            - Read significant eGenes (qval < 0.05) from the permutation results
#            - Add chromosome, position and alleles from the genotype pvar
#            - Flag lead SNPs inside a Ziffra 2021 primary peak (In_Peak)
#            - Add gene symbols; one sheet per cell type, manuscript labels
#
# Inputs:    in_dir            tensorQTL permutation results directory
#            allele_file       Genotype pvar (CHROM, POS, ID, REF, ALT)
#            peak_dir          Ziffra 2021 supplementary tables directory
#            exp_pc_map        Expression PCs per cell type
#            gene_lookup_file  Ensembl ID to gene symbol lookup
#
# Outputs:   out_file          xlsx, one sheet per cell type
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

message("\n\nBuilding eQTL supplementary table ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(readxl)
library(GenomicRanges)
library(openxlsx)

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
in_dir            <- snakemake@params[["in_dir"]]
allele_file       <- snakemake@params[["allele_file"]]
peak_dir          <- snakemake@params[["peak_dir"]]
gene_lookup_file  <- snakemake@params[["gene_lookup_file"]]
manuscript_labels <- snakemake@params[["labels"]]
expPC_map         <- unlist(snakemake@params[["exp_pc_map"]])
out_file          <- snakemake@output[[1]]

norm_method <- 'quantile'
genPC <- 4

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_dir", "allele_file", "peak_dir", "gene_lookup_file", "labels",
               "expPC_map", "norm_method", "genPC", "out_file"),
  value    = c(in_dir, allele_file, peak_dir, gene_lookup_file,
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               paste(names(expPC_map), expPC_map, sep = ":", collapse = ", "),
               norm_method, genPC, out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

gene_lookup_tbl <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

##  Ziffra 2021 snATAC-seq peaks  -----------------------------------------------------
message("Loading Ziffra peak data ...")
peak_file <- file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx")

peak_coords <- read_excel(peak_file, sheet = "ST2 AllPrimaryPeaks") |>
  dplyr::select(seqnames, start, end, peak_name)

master_peaks_gr <- makeGRangesFromDataFrame(peak_coords, keep.extra.columns = TRUE)

##  Genotype pvar (alleles and positions)  --------------------------------------------
message("Loading pvar file ...")
pvar <- read_tsv(allele_file, comment = "#", show_col_types = FALSE,
                 col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")) |>
  dplyr::select(CHROM, POS, ID, REF, ALT) |>
  distinct(ID, .keep_all = TRUE)

##  Build one table per cell type  ----------------------------------------------------
final_excel_list <- list()

for (cell_type in names(expPC_map)) {
  expPC <- expPC_map[[cell_type]]

  log_dir <- paste0(in_dir, cell_type, '_', norm_method, '_genPC_', genPC, '_expPC_', expPC, '/')
  log_file <- paste0(log_dir, cell_type, '_', norm_method, '_perm.cis_qtl.txt.gz')

  if (!file.exists(log_file)) {
    message("File not found for ", cell_type, ": ", log_file, " (skipped)")
    next
  }

  message('\nLoading eQTL perm file for: ', cell_type)
  eqtl_tbl <- read_tsv(log_file, show_col_types = FALSE) |>
    dplyr::filter(qval < 0.05) |>
    dplyr::select(ensembl_id = phenotype_id, SNP = variant_id, af,
                  pval_nominal, pval_beta, qval, slope, slope_se)

  if (nrow(eqtl_tbl) == 0) next

  # REF/ALT and coordinates from the pvar
  message('Adding alleles ...')
  n_before <- nrow(eqtl_tbl)
  eqtl_enriched <- eqtl_tbl |>
    inner_join(pvar, by = c("SNP" = "ID"))

  if (nrow(eqtl_enriched) != n_before) {
    message(cell_type, ": row count changed after pvar join (", n_before, " -> ",
            nrow(eqtl_enriched), "), check for duplicate SNP IDs")
  }

  message('Adding chr prefix to eQTL tbl ...')
  eqtl_enriched <- eqtl_enriched |>
    mutate(CHROM = ifelse(!str_detect(CHROM, "chr"), paste0("chr", CHROM), CHROM)) |>
    drop_na(CHROM, POS)

  # Lead SNP inside any Ziffra peak
  message('Checking for SNP overlaps in snATAC-seq peaks ...')
  snp_gr <- GRanges(
    seqnames = eqtl_enriched$CHROM,
    ranges = IRanges(start = eqtl_enriched$POS, end = eqtl_enriched$POS)
  )
  eqtl_enriched$In_Peak <- countOverlaps(snp_gr, master_peaks_gr) > 0

  message('Checking for NAs ...')
  na_count <- sum(is.na(eqtl_enriched))
  if (na_count > 0) {
    message(paste0("Warning: ", na_count, " NA values detected in ", cell_type,
                   ". This usually means SNP IDs don't match between TensorQTL and PVAR."))
  }

  # Cell type label, gene symbol and column order
  message('Formatting table ...')
  n_before <- nrow(eqtl_tbl)
  eqtl_enriched <- eqtl_enriched |>
    mutate(cell_type = relabel_cell_types(!!cell_type, manuscript_labels)) |>
    inner_join(gene_lookup_tbl, by = join_by(ensembl_id == ensembl_gene_id)) |>
    mutate(CHROM = str_remove(CHROM, "^chr")) |>
    dplyr::select(cell_type, ensembl_id, symbol = external_gene_name, CHROM,
                  SNP, POS, REF, ALT,  AF = af, slope, slope_se, pval_nominal,
                  pval_beta, qval, In_Peak)

  if (nrow(eqtl_enriched) != n_before) {
    message(cell_type, ": row count changed after gene lookup join (", n_before, " -> ",
            nrow(eqtl_enriched), "), check for duplicate Gene IDs")
  }

  message('Any NAs in final tbl? ', anyNA(eqtl_enriched))

  final_excel_list[[relabel_cell_types(cell_type, manuscript_labels)]] <- eqtl_enriched
  message("Processed ", cell_type, ": ", nrow(eqtl_enriched), " eQTLs")
}

##  Write xlsx  -----------------------------------------------------------------------
message("Writing to Excel file ...")
write.xlsx(final_excel_list,
           file = out_file,
           overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
