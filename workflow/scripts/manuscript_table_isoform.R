#--------------------------------------------------------------------------------------
#
#    eGenes with transcript-specific bulk eQTL (Supp Table 7)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: isoform_tbl
#            Upstream:   08 isoform_eqtl_comparison (per eGene-eSNP summary)
#            Downstream: Supp Table 7
#
# Purpose:   Prenatal eGene-eSNP pairs where, in O'Brien et al. 2018 bulk
#            transcript-level eQTL, the SNP affects only a subset of the gene's
#            transcripts (pattern isoform_specific, from 08), in the same
#            direction as the snRNA-seq effect (direction_concordant), and the
#            most significant transcript has P < bulk_p_thresh:
#            - Prenatal cell types listed in config cell_types order, manuscript
#              labels (config manuscript_labels)
#            - Sorted by gene symbol (case ignored, missing symbols last), then
#              by transcript P
#
# Inputs:    in_file            08 isoform_eqtl_comparison summary TSV
#            cell_types         config['cell_types'] (pipeline labels)
#            labels             config['manuscript_labels']
#            nominal_threshold  08 bulk P threshold for a significant transcript
#            bulk_p_thresh      P threshold for the most significant transcript
#
# Outputs:   out_file           xlsx, one sheet
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

message("\n\nBuilding the transcript-specific eGene table ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types, cell_type_levels)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
in_file           <- snakemake@input[[1]]
cell_types        <- unlist(snakemake@params[["cell_types"]])
manuscript_labels <- snakemake@params[["labels"]]
nominal_threshold <- as.numeric(snakemake@params[["nominal_threshold"]])
bulk_p_thresh     <- as.numeric(snakemake@params[["bulk_p_thresh"]])
out_file          <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_file", "cell_types", "labels", "nominal_threshold", "bulk_p_thresh", "out_file"),
  value    = c(in_file, paste(cell_types, collapse = ", "),
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               nominal_threshold, bulk_p_thresh, out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Filter  ---------------------------------------------------------------------------
summary_tbl <- read_tsv(in_file, show_col_types = FALSE) |>
  mutate(direction_concordant = as.logical(direction_concordant))
message("eGene-eSNP pairs in the 08 summary: ", nrow(summary_tbl))

filtered <- summary_tbl |>
  filter(pattern == "isoform_specific",
         direction_concordant %in% TRUE,
         most_sig_bulk_pval < bulk_p_thresh)
message("Isoform-specific, concordant, transcript P < ", bulk_p_thresh, ": ",
        nrow(filtered), " pairs, ", n_distinct(filtered$gene_id), " eGenes")

# Pipe-separated prenatal cell types, relabelled and in config order
order_cell_type_list <- function(x) {
  cts <- relabel_cell_types(str_split_1(x, fixed("|")), manuscript_labels)
  paste(cell_type_levels(cts, cell_types, manuscript_labels), collapse = "|")
}

##  Table  ----------------------------------------------------------------------------
n_sig_col <- paste0("N bulk transcripts P < ", nominal_threshold)

final_tbl <- filtered |>
  mutate(cell_types = map_chr(cell_types, order_cell_type_list)) |>
  arrange(is.na(gene_symbol), str_to_upper(gene_symbol), most_sig_bulk_pval) |>
  transmute(`Ensembl ID`                       = gene_id,
            `Gene symbol`                      = gene_symbol,
            `Variant ID`                       = variant_id,
            `Prenatal cell types`              = cell_types,
            `N prenatal cell types`            = n_cell_types,
            `Prenatal P`                       = snrna_pval_nominal,
            `Prenatal beta`                    = snrna_slope,
            `Prenatal beta SE`                 = snrna_slope_se,
            `N bulk transcripts tested`        = n_bulk_transcripts_tested,
            !!n_sig_col                        := n_bulk_transcripts_nominal,
            `Most significant bulk transcript` = most_sig_transcript,
            `Bulk transcript P`                = most_sig_bulk_pval,
            `Bulk transcript beta`             = most_sig_bulk_slope)

##  Write xlsx  -----------------------------------------------------------------------
message("\nWriting: ", out_file)
write.xlsx(list("Implicated eTranscripts" = final_tbl), file = out_file, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
