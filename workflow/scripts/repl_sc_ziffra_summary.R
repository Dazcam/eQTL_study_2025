#--------------------------------------------------------------------------------------
#
#    Combine Ziffra 2021 snATAC-seq peak enrichments across cell types
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  06QTL-REPLICATION | Rule: ziffra_summary (single job)
#            Upstream:   atac_enrich
#            Downstream: 06 report (Ziffra tab); 14 eqtl_qc_plt (Fig 2E)
#
# Purpose:   - Read the per-cell-type Ziffra enrichment tables
#            - Bind them into one table, one row per cell type x peak set
#
# Inputs:    enrich       atac_enrich outputs, one per cell type
#            cell_types   cell types in the same order as enrich (params)
#
# Outputs:   out_file     TSV: eqtl_cell_type, peak_cell_type, p_value,
#                         observed_overlap, mean_permuted_overlap,
#                         fold_enrichment, source_cell_type
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

message("\n\nCombining Ziffra snATAC-seq peak enrichments across cell types ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Input and output paths
enrich_files <- unlist(snakemake@input[["enrich"]])
cell_types   <- unlist(snakemake@params[["cell_types"]])
out_file     <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("n_enrich_files", "n_cell_types", "out_file"),
  value    = c(length(enrich_files), length(cell_types), out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

stopifnot(length(enrich_files) == length(cell_types))

##  Combine  ---------------------------------------------------------------------------
ziffra_tbl <- map2_dfr(enrich_files, cell_types, function(f, ct) {
  read_rds(f) |>
    mutate(source_cell_type = ct)
})

message("Rows: ", nrow(ziffra_tbl), " (", n_distinct(ziffra_tbl$eqtl_cell_type), " cell types)")

##  Write output  ----------------------------------------------------------------------
write_tsv(ziffra_tbl, out_file)
message("Written to: ", out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
