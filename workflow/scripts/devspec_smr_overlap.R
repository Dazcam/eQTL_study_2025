#--------------------------------------------------------------------------------------
#
#    Developmental specificity: SMR-significant SNP overlap with dev-specificity classification
#
#--------------------------------------------------------------------------------------

#  For every SMR result file (L1+L2 fetal cell types x 5 disorders), extracts the
#  SMR-and-HEIDI-significant (Gene, topSNP) pairs using the same per-file Bonferroni
#  + HEIDI threshold as smr_report.Rmd (p_SMR < p_smr/n_probes & p_HEIDI > p_heidi,
#  n_probes = unique probes tested WITHIN that cell_type/disorder file).
#
#  Joins these pairs against classify_eqtl_specificity's granular "any" tier
#  final_category (testable_shared / testable_specific / untestable_* /
#  testable_pair_untested) on (probeID, topSNP) == (phenotype_id, variant_id).
#
#  Retains every (cell_type, disorder, Gene, topSNP) row individually -- not
#  deduplicated -- since SMR selects its topSNP independently per cell type
#  (via each cell type's own besd file), so the same gene can have a different
#  topSNP, and therefore a different dev-specificity category, in different
#  cell types.
#
#  Pairs with no match in the dev-specificity universe are retained and flagged
#  explicitly ("not_in_dev_spec_universe") rather than dropped -- SMR's topSNP
#  selection is independent of TensorQTL's permutation-pass lead SNP, so a
#  non-trivial no-match rate is expected, not necessarily an error.
#
#  Step 3.5 of the developmental-specificity pipeline (standalone; reads
#  classify_eqtl_specificity's output, does not modify it).
#
#--------------------------------------------------------------------------------------

# Set up logging for Snakemake
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

message("\n\nExtracting SMR-significant SNPs and overlapping with dev-specificity classification ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
smr_files          <- snakemake@input[["smr_files"]]
classification_file <- snakemake@input[["classification"]]
gene_lookup_file   <- snakemake@input[["gene_lookup"]]
output_file        <- snakemake@output[[1]]
cell_types         <- snakemake@params[["cell_types"]]
gwas_list          <- snakemake@params[["gwas_list"]]
p_smr              <- as.numeric(snakemake@params[["p_smr"]])
p_heidi            <- as.numeric(snakemake@params[["p_heidi"]])

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("classification_file", "gene_lookup_file", "n_smr_files", "n_cell_types", "n_gwas",
               "p_smr", "p_heidi", "output_file"),
  value    = c(classification_file, gene_lookup_file, length(smr_files), length(cell_types), length(gwas_list),
               p_smr, p_heidi, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\nCell types (", length(cell_types), "): ", paste(cell_types, collapse = ", "))
message("Disorders (", length(gwas_list), "): ", paste(gwas_list, collapse = ", "))

# Functions ----------------------------------------------------------------------------

# Read one SMR file, apply the per-file Bonferroni + HEIDI threshold (identical
# logic to smr_report.Rmd's summarize_smr()/get_heidi() functions), return the
# significant rows only.
read_smr_significant <- function(cell_type, disorder, smr_dir) {
  file_path <- file.path(smr_dir, cell_type, paste0(cell_type, "_", disorder, ".smr"))

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  smr <- suppressMessages(read_tsv(file_path, show_col_types = FALSE))

  if (nrow(smr) == 0) {
    message("  ", cell_type, " / ", disorder, ": empty file (0 probes tested)")
    return(tibble())
  }

  n_probes <- n_distinct(smr$probeID)

  sig <- smr |>
    filter(p_SMR < (p_smr / n_probes), p_HEIDI > p_heidi & !is.na(p_HEIDI)) |>
    mutate(cell_type = cell_type, disorder = disorder, n_probes_tested = n_probes)

  message("  ", cell_type, " / ", disorder, ": ", n_probes, " probes tested, ",
          nrow(sig), " SMR+HEIDI significant")

  sig
}

# Main -----------------------------------------------------------------------------------

# Directory shared by all smr_files (same for every file, just take the common root)
smr_dir <- unique(dirname(dirname(smr_files)))
if (length(smr_dir) != 1) {
  stop("Expected a single common SMR output directory, got: ", paste(smr_dir, collapse = ", "))
}
message("\nSMR output directory: ", smr_dir)

message("\n--- Reading SMR files and applying significance threshold ---\n")

combinations <- expand.grid(cell_type = cell_types, disorder = gwas_list, stringsAsFactors = FALSE)

smr_sig_all <- map2_dfr(
  combinations$cell_type, combinations$disorder,
  ~ read_smr_significant(.x, .y, smr_dir)
)

message("\nTotal SMR+HEIDI significant (cell_type, disorder, Gene, topSNP) rows: ", nrow(smr_sig_all))

if (nrow(smr_sig_all) == 0) {
  stop("No SMR-significant rows found across any cell_type/disorder combination -- check p_smr/p_heidi and file paths.")
}

smr_sig_all <- smr_sig_all |>
  select(cell_type, disorder, probeID, Gene, topSNP,
         b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI, n_probes_tested)

# ========================================================================================
# Join against dev-specificity classification
# ========================================================================================

message("\n--- Loading dev-specificity classification ---\n")
classification <- read_rds(classification_file)
any_table <- classification$any_table
message("  any_table rows: ", nrow(any_table))

message("\n--- Joining SMR-significant pairs against any_table ---\n")

smr_overlap <- smr_sig_all |>
  left_join(
    any_table |> select(phenotype_id, variant_id, final_category),
    by = c("probeID" = "phenotype_id", "topSNP" = "variant_id")
  ) |>
  mutate(final_category = replace_na(final_category, "not_in_dev_spec_universe"))

n_no_match <- sum(smr_overlap$final_category == "not_in_dev_spec_universe")
message("SMR-significant pairs with NO match in the dev-specificity universe: ",
        n_no_match, " / ", nrow(smr_overlap),
        " (", round(100 * n_no_match / nrow(smr_overlap), 1), "%)")

# ========================================================================================
# Summary
# ========================================================================================

message("\n--- Overall final_category distribution among SMR-significant pairs ---\n")
overall_summary <- smr_overlap |>
  count(final_category, name = "n") |>
  mutate(prop = round(n / sum(n), 3)) |>
  arrange(desc(n))
print(overall_summary)

message("\n--- Breakdown by disorder ---\n")
smr_overlap |>
  count(disorder, final_category, name = "n") |>
  arrange(disorder, desc(n)) |>
  print(n = Inf)

message("\n--- Breakdown by cell type ---\n")
smr_overlap |>
  count(cell_type, final_category, name = "n") |>
  arrange(cell_type, desc(n)) |>
  print(n = Inf)

message("\n--- Unique genes / SNPs among SMR-significant pairs ---\n")
message("  Unique genes: ", n_distinct(smr_overlap$probeID))
message("  Unique topSNPs: ", n_distinct(smr_overlap$topSNP))
message("  Unique (gene, topSNP) pairs (collapsed across cell_type/disorder): ",
        nrow(distinct(smr_overlap, probeID, topSNP)))

# ========================================================================================
# Add gene symbols
# ========================================================================================

message("\nLoading gene symbol lookup: ", gene_lookup_file)
gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

smr_overlap <- smr_overlap |>
  left_join(gene_lookup, by = c("probeID" = "ensembl_gene_id")) |>
  rename(gene_symbol = external_gene_name) |>
  relocate(gene_symbol, .after = probeID)

n_missing_symbol <- sum(is.na(smr_overlap$gene_symbol))
if (n_missing_symbol > 0) {
  message("  WARNING: ", n_missing_symbol, " rows had no matching gene symbol in the lookup")
}

# ========================================================================================
# Write single tracked output
# ========================================================================================

write_tsv(smr_overlap, output_file)
message("\nSMR-significant / dev-specificity overlap table written to: ", output_file)
