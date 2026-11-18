#--------------------------------------------------------------------------------------
#
#    eQTL with larger effects in prenatal than adult brain cells (Supp Table 5)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: eqtl_heterogeneity_tbl
#            Upstream:   15 heterogeneity_test ($magnitude_consistent_pairs)
#            Downstream: Supp Table 5
#
# Purpose:   Prenatal eGene-eSNP pairs whose effect differs (Cochran's Q) from
#            every adult cell type they were tested in, with |prenatal beta| >
#            |adult beta| in all of them, as defined by 15:
#            - Prenatal cell type (lowest qval), number of adult comparisons and
#              the min / max of |prenatal beta| - |adult beta| (4 decimal places)
#            - Manuscript cell type labels (config manuscript_labels)
#            - Sorted by gene symbol (case ignored, missing symbols last)
#
# Inputs:    in_file   15 heterogeneity_test RDS
#            labels    config['manuscript_labels']
#
# Outputs:   out_file  xlsx, one sheet
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

message("\n\nBuilding the prenatal-larger eQTL table ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
in_file           <- snakemake@input[[1]]
manuscript_labels <- snakemake@params[["labels"]]
out_file          <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_file", "labels", "out_file"),
  value    = c(in_file,
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Table  ----------------------------------------------------------------------------
final_tbl <- read_rds(in_file)$magnitude_consistent_pairs |>
  arrange(is.na(gene_symbol), str_to_upper(gene_symbol), phenotype_id, variant_id) |>
  transmute(`Ensembl ID`                 = phenotype_id,
            `Gene symbol`                = gene_symbol,
            `Variant ID`                 = variant_id,
            `Prenatal cell type`         = relabel_cell_types(fetal_cell_type, manuscript_labels),
            `N adult comparisons`        = n_comparisons_tested,
            `Min effect size difference` = round(min_magnitude_diff, 4),
            `Max effect size difference` = round(max_magnitude_diff, 4))

message("Pairs: ", nrow(final_tbl), "; eGenes: ", n_distinct(final_tbl$`Ensembl ID`))

##  Write xlsx  -----------------------------------------------------------------------
message("\nWriting: ", out_file)
write.xlsx(list("Larger in prenatal" = final_tbl), file = out_file, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
