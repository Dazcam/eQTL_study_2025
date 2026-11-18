#--------------------------------------------------------------------------------------
#
#    Prenatal eQTL tested in adult brain cell types (Supp Table 4)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: eqtl_specificity_tbl
#            Upstream:   15 classify_eqtl_specificity ($any_table)
#            Downstream: Supp Table 4
#
# Purpose:   One sheet per group of unique prenatal eGene-eSNP pairs, as
#            classified against the Jang et al. adult brain cell types by 15:
#            - eQTL not BF sig in adult cells: testable, not significant
#              (Bonferroni) in any adult cell type; sorted by adult P, highest first
#            - Shared eQTL: testable and significant in at least one adult cell type
#            - eGene not in adult dataset: gene not expressed in any adult cell type
#            Adult columns come from the adult cell type with the smallest P.
#            Shared and not-in-adult sheets are sorted by gene symbol (case
#            ignored, missing symbols last).
#
# Inputs:    in_file           15 classify_eqtl_specificity RDS
#            gene_lookup_file  Ensembl ID to gene symbol lookup
#
# Outputs:   out_file          xlsx, three sheets
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

message("\n\nBuilding the prenatal vs adult eQTL table ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Input and output paths
in_file          <- snakemake@input[["classification"]]
gene_lookup_file <- snakemake@input[["gene_lookup"]]
out_file         <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_file", "gene_lookup_file", "out_file"),
  value    = c(in_file, gene_lookup_file, out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read and annotate  ----------------------------------------------------------------
gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

any_table <- read_rds(in_file)$any_table |>
  left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) |>
  rename(gene_symbol = external_gene_name)

message("Pairs per category:")
any_table |> count(final_category) |> knitr::kable(format = "simple") |> cat(sep = "\n")

##  Sheets  ---------------------------------------------------------------------------
not_sig_adult <- any_table |>
  filter(final_category == "testable_specific") |>
  arrange(desc(random_p_adult)) |>
  transmute(`Ensembl ID`                       = phenotype_id,
            `Variant ID`                       = variant_id,
            `Gene symbol`                      = gene_symbol,
            `N adult cell types tested`        = n_jang_types_tested,
            `Most significant adult cell type` = best_jang_ct,
            `Adult P`                          = random_p_adult,
            `Adult beta`                       = beta_adult)

shared <- any_table |>
  filter(final_category == "testable_shared") |>
  arrange(is.na(gene_symbol), str_to_upper(gene_symbol), random_p_adult) |>
  transmute(`Ensembl ID`                       = phenotype_id,
            `Variant ID`                       = variant_id,
            `Gene symbol`                      = gene_symbol,
            `N adult cell types tested`        = n_jang_types_tested,
            `Most significant adult cell type` = best_jang_ct,
            `Adult P`                          = random_p_adult,
            `Adult beta`                       = beta_adult,
            `Adult lead SNP`                   = best_ct_is_jang_lead_snp)

not_in_adult <- any_table |>
  filter(!gene_expressed) |>
  arrange(is.na(gene_symbol), str_to_upper(gene_symbol), variant_id) |>
  transmute(`Ensembl ID`            = phenotype_id,
            `Variant ID`            = variant_id,
            `Gene symbol`           = gene_symbol,
            `Gene in adult dataset` = gene_expressed,
            `SNP genotyped`         = snp_genotyped)

sheets <- list("eQTL not BF sig in adult cells" = not_sig_adult,
               "Shared eQTL"                    = shared,
               "eGene not in adult dataset"     = not_in_adult)

message("\nRows per sheet:")
for (s in names(sheets)) {
  message("  ", s, ": ", nrow(sheets[[s]]), " pairs, ",
          n_distinct(sheets[[s]]$`Ensembl ID`), " eGenes")
}

##  Write xlsx  -----------------------------------------------------------------------
message("\nWriting: ", out_file)
write.xlsx(sheets, file = out_file, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
