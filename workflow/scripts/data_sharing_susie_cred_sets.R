#--------------------------------------------------------------------------------------
#
#    Prepare SuSiE credible set output for data sharing (figshare)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  15DATA_SHARING | Rule: data_susie (one job per cell type)
#            Upstream:   08SUSIE sort_susie (merged, sorted credible sets)
#            Downstream: data_susie_mk_tar
#
# Purpose:   For one cell type:
#            - Add a gene_symbol column directly after phenotype_id (Ensembl ID)
#            - Name the output with the manuscript cell type label (Glu-UL -> Glu-A,
#              Glu-DL -> Glu-B, "-Q4-" stripped); relabel any cell_type column to match
#            - Keep all other rows and columns exactly as they are, in the same order
#
# Inputs:    susie        {cell_type_old}.mrgd.srtd.susie.cred.hp.txt.gz
#            gene_lookup  BioMart gene annotation (hg38)
#
# Outputs:   tsv          {cell_type}.susie.cred.hp.tsv.gz (manuscript label)
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

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

message("\n\nPreparing SuSiE credible sets for data sharing ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages(library(tidyverse))

# --- Set variables
susie_file       <- snakemake@input[["susie"]]
gene_lookup_file <- snakemake@input[["gene_lookup"]]
out_file         <- snakemake@output[["tsv"]]
cell_type_old    <- snakemake@params[["cell_type_old"]]
cell_type_new    <- snakemake@wildcards[["cell_type"]]


# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("susie_file", "gene_lookup_file", "out_file", "cell_type_old", "cell_type_new"),
  value    = c(susie_file, gene_lookup_file, out_file, cell_type_old, cell_type_new)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

message("=== ", cell_type_old, " -> ", cell_type_new, " ===")

# --- Relabel helper (display only): Glu-UL -> Glu-A, Glu-DL -> Glu-B, and strip
# "-Q4-" from pseudotime bin names. Same convention as manuscript_table_smr_sig.R.
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

if (relabel_cell_type(cell_type_old) != cell_type_new) {
  stop("Wildcard '", cell_type_new, "' is not the relabelled form of '", cell_type_old,
       "' (expected '", relabel_cell_type(cell_type_old), "')")
}
if (str_detect(cell_type_new, "UL|DL|Q4")) {
  stop("Old cell type naming remains in output label: ", cell_type_new)
}

# --- Gene lookup
gene_lookup_raw <- suppressMessages(read_tsv(gene_lookup_file, show_col_types = FALSE)) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct()

n_multi <- gene_lookup_raw |>
  count(ensembl_gene_id) |>
  filter(n > 1) |>
  nrow()
message("Gene lookup: ", n_multi, " Ensembl IDs with more than one distinct symbol; ",
        "keeping the first listed for each")

gene_lookup <- gene_lookup_raw |>
  distinct(ensembl_gene_id, .keep_all = TRUE) |>
  mutate(external_gene_name = na_if(external_gene_name, "")) |>
  rename(gene_symbol = external_gene_name)

# --- Credible sets
susie_tbl <- suppressMessages(read_tsv(susie_file, show_col_types = FALSE))
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
  final_tbl <- final_tbl |> mutate(cell_type = relabel_cell_type(cell_type))
}

# --- Checks: the join must not add or drop rows, and must not disturb the data
if (nrow(final_tbl) != n_in) {
  stop("Row count changed after adding gene symbols: ", n_in, " -> ", nrow(final_tbl))
}
if (sum(duplicated(final_tbl)) > sum(duplicated(susie_tbl))) {
  stop("Duplicated rows introduced by adding gene symbols")
}
if (!identical(final_tbl$phenotype_id, susie_tbl$phenotype_id)) {
  stop("Row order changed after adding gene symbols")
}
if (!identical(names(final_tbl)[-which(names(final_tbl) == "gene_symbol")], names(susie_tbl))) {
  stop("Original columns changed after adding gene symbols")
}

n_no_symbol <- final_tbl |> filter(is.na(gene_symbol)) |> distinct(phenotype_id) |> nrow()
message(n_no_symbol, " of ", n_distinct(final_tbl$phenotype_id),
        " genes have no gene symbol (left as NA)")
message("gene_symbol is column ", which(names(final_tbl) == "gene_symbol"),
        "; phenotype_id is column ", which(names(final_tbl) == "phenotype_id"))

# --- Write (readr compresses on the .gz extension)
write_tsv(final_tbl, out_file)
message("Wrote: ", out_file)

message("Done: ", cell_type_new)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
