#--------------------------------------------------------------------------------------
#
#    SuSiE credible sets for data sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  18DATA-SHARING | Rule: data_susie (one job per cell type)
#            Upstream:   10SUSIE sort_susie (merged, sorted credible sets)
#            Downstream: data_susie_mk_tar
#
# Purpose:   For one cell type:
#            - Add a gene_symbol column directly after phenotype_id (Ensembl ID)
#            - Relabel any cell_type column to manuscript labels
#            - Keep all other rows and columns as they are, in the same order
#
# Inputs:    susie        Credible sets for the cell type (pipeline label)
#            gene_lookup  BioMart gene annotation (hg38)
#
# Outputs:   tsv          {cell_type}.susie.cred.hp.tsv.gz (manuscript label)
#
# Notes:     The {cell_type} wildcard is the manuscript label; cell_type_old is the
#            pipeline label used in the input path.
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

message("\n\nSuSiE credible sets for data sharing ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "functions.R"))
source(file.path(snakemake@scriptdir, "data_sharing_functions.R"))

# Input and output paths
susie_file       <- snakemake@input[["susie"]]
gene_lookup_file <- snakemake@input[["gene_lookup"]]
out_file         <- snakemake@output[["tsv"]]

# Settings
cell_type_old <- snakemake@params[["cell_type_old"]]
cell_type_new <- snakemake@wildcards[["cell_type"]]
labels        <- snakemake@params[["labels"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("susie_file", "gene_lookup_file", "out_file", "cell_type_old",
               "cell_type_new", "labels"),
  value    = c(susie_file, gene_lookup_file, out_file, cell_type_old, cell_type_new,
               paste(names(labels), unlist(labels), sep = " -> ", collapse = ", "))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

if (relabel_cell_types(cell_type_old, labels) != cell_type_new) {
  stop("Wildcard '", cell_type_new, "' is not the manuscript label of '", cell_type_old,
       "' (expected '", relabel_cell_types(cell_type_old, labels), "')")
}

## Gene lookup  ------------------------------------------------------------------------
gene_lookup <- read_gene_lookup(gene_lookup_file) |>
  dplyr::rename(gene_symbol = external_gene_name)

## Credible sets  ----------------------------------------------------------------------
susie_tbl <- read_tsv(susie_file, show_col_types = FALSE)
n_in <- nrow(susie_tbl)
message("Read ", n_in, " rows, ", ncol(susie_tbl), " columns: ", susie_file)

if (!"phenotype_id" %in% names(susie_tbl)) {
  stop("No phenotype_id column in ", susie_file)
}
if ("gene_symbol" %in% names(susie_tbl)) {
  stop("gene_symbol column already exists in ", susie_file)
}

final_tbl <- susie_tbl |>
  left_join(gene_lookup, by = join_by(phenotype_id == ensembl_gene_id)) |>
  relocate(gene_symbol, .after = phenotype_id)

if ("cell_type" %in% names(final_tbl)) {
  message("Relabelling values in the cell_type column")
  final_tbl <- final_tbl |> mutate(cell_type = relabel_cell_types(cell_type, labels))
}

## Checks: the join must not add, drop, reorder or change rows  ------------------------
if (nrow(final_tbl) != n_in) {
  stop("Row count changed after adding gene symbols: ", n_in, " -> ", nrow(final_tbl))
}
if (sum(duplicated(final_tbl)) > sum(duplicated(susie_tbl))) {
  stop("Duplicated rows introduced by adding gene symbols")
}
if (!identical(final_tbl$phenotype_id, susie_tbl$phenotype_id)) {
  stop("Row order changed after adding gene symbols")
}
if (!identical(setdiff(names(final_tbl), "gene_symbol"), names(susie_tbl))) {
  stop("Original columns changed after adding gene symbols")
}

n_no_symbol <- final_tbl |> filter(is.na(gene_symbol)) |> distinct(phenotype_id) |> nrow()
message(n_no_symbol, " of ", n_distinct(final_tbl$phenotype_id),
        " genes have no gene symbol (left as NA)")

## Write (readr compresses on the .gz extension)  --------------------------------------
write_tsv(final_tbl, out_file)
message("Wrote: ", out_file)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
