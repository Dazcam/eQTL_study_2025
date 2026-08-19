#--------------------------------------------------------------------------------------
#
#    Developmental specificity: shared vs. fetal-specific eQTL classification
#
#--------------------------------------------------------------------------------------

#  For every fetal eGene-eSNP pair in the universe extracted by
#  devspec_extract_universe.R, classifies testability against Fujita 2024 adult
#  sn-eQTL data and, where testable, calls shared vs. fetal-specific using Fujita's
#  own published significance call (significant_by_2step_FDR).
#
#  "Any" tier (cell-type-agnostic, KEPT): operates on the 6,844 unique (gene, SNP)
#  pairs across all cell types. Tested against all 7 Fujita cell types pooled;
#  "shared" = significant in >=1 Fujita type. Produces the fetal-specific SNP list
#  for the OCR enrichment step (via LD clumping) and the fetal-specific input
#  population for the heterogeneity test (step C).
#
#  Step 2 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  ADDITION -- cell_specific tier retired:
#
#  The previous "cell_specific" tier (per fetal cell type, tested against its a
#  priori label-matched Fujita cell type only, with a sign-concordance check) has
#  been removed. It fed no downstream computation (its only consumer was the report's
#  cell_specific display) and its role -- a formal, cell-type-matched comparison
#  between fetal and adult effects -- is now served properly by the heterogeneity
#  test (step C), which tests effect-size difference directly rather than relying
#  on Fujita's own significance threshold and a bare sign check.
#
#  se added to the Fujita column selection, and pairs_all / gene_membership /
#  snp_membership -- already computed here to build the "any" tier, previously
#  discarded after use -- are now written to output. Step C reads these directly
#  (per-Fujita-cell-type beta/se and testability, keyed by the label-matched
#  cell type for each fetal row) rather than re-reading the Fujita files.
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

message("\n\nClassifying fetal eGene-eSNP pairs: testability + shared/specific vs. Fujita ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
fugita_files       <- snakemake@input[["fugita_files"]]
output_file        <- snakemake@output[[1]]
fugita_cell_types  <- snakemake@params[["fugita_cell_types"]]

names(fugita_files) <- fugita_cell_types

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "n_fugita_files", "output_file"),
  value    = c(eqtl_universe_file, length(fugita_files), output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# se added to col_select/col_types (see ADDITION note above) -- required by the
# heterogeneity test (C) for the adult side of Cochran's Q.
read_fugita_celltype <- function(file_path, fugita_ct, universe_genes, universe_snps) {
  message("Reading Fujita data for: ", fugita_ct, " (", file_path, ")")

  dat <- read_tsv(
    file_path,
    col_select = c(gene_id, snps, pvalue, beta, se, significant_by_2step_FDR, REF, ALT, ALT_AF),
    col_types = cols(gene_id = "c", snps = "c", pvalue = "d", beta = "d", se = "d",
                     significant_by_2step_FDR = "c", REF = "c", ALT = "c", ALT_AF = "d"),
    show_col_types = FALSE
  )

  gene_set <- unique(dat$gene_id)
  snp_set  <- unique(dat$snps)

  pair_lookup <- dat %>%
    filter(gene_id %in% universe_genes, snps %in% universe_snps) %>%
    distinct(gene_id, snps, .keep_all = TRUE) %>%
    mutate(fugita_cell_type = fugita_ct)

  message("  Genes tested: ", length(gene_set),
          " | SNPs tested: ", length(snp_set),
          " | Relevant pairs retained: ", nrow(pair_lookup))

  list(gene_set = gene_set, snp_set = snp_set, pair_lookup = pair_lookup)
}

# Main -----------------------------------------------------------------------------------

message("\nLoading fetal eGene/eSNP universe ...")
eqtl_universe <- read_rds(eqtl_universe_file)
message("  Rows: ", nrow(eqtl_universe), " | Unique eGenes: ", n_distinct(eqtl_universe$phenotype_id))

universe_genes <- unique(eqtl_universe$phenotype_id)
universe_snps  <- unique(eqtl_universe$variant_id)

# Deduplicate to unique (gene, SNP) pairs for the "any" tier
unique_pairs <- eqtl_universe %>%
  distinct(phenotype_id, variant_id)
message("  Unique (gene, SNP) pairs: ", nrow(unique_pairs))

# Loop over Fujita files ----------------------------------------------------------------
gene_sets    <- list()
snp_sets     <- list()
pair_lookups <- list()

for (ct in fugita_cell_types) {
  res <- read_fugita_celltype(fugita_files[[ct]], ct, universe_genes, universe_snps)
  gene_sets[[ct]]    <- res$gene_set
  snp_sets[[ct]]     <- res$snp_set
  pair_lookups[[ct]] <- res$pair_lookup

  rm(res)
  gc()
}

gene_set_any <- Reduce(union, gene_sets)
snp_set_any  <- Reduce(union, snp_sets)

pairs_all <- bind_rows(pair_lookups)

# Per-gene/SNP membership lookups. Used to build the "any" tier's testability check
# below, and (see ADDITION note) now also written to output for the heterogeneity
# test (C) to re-check testability in each row's specific label-matched Fujita
# cell type.
gene_membership <- map_dfr(names(gene_sets), ~ tibble(fugita_cell_type = .x, gene_id = gene_sets[[.x]]))
snp_membership  <- map_dfr(names(snp_sets),  ~ tibble(fugita_cell_type = .x, snps = snp_sets[[.x]]))

# =======================================================================================
# "ANY" TIER: 6,844 unique pairs, checked across all 7 Fujita cell types
# =======================================================================================
message("\n--- Building 'any' tier classification ---")

# Significance call: significant in >= 1 Fujita cell type
pairs_any_sig <- pairs_all %>%
  group_by(gene_id, snps) %>%
  summarise(
    significant_any = any(significant_by_2step_FDR == "Yes"),
    n_fugita_types_tested = n(),
    .groups = "drop"
  )

# Provenance: which Fujita cell type gave the smallest p-value (for QC/reporting
# only -- NOT used by the heterogeneity test, which looks up the label-matched
# cell type specifically via pairs_all instead; see ADDITION note above).
pairs_any_provenance <- pairs_all %>%
  group_by(gene_id, snps) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene_id, snps, best_fugita_ct = fugita_cell_type,
         pvalue_adult = pvalue, beta_adult = beta)

pairs_any <- left_join(pairs_any_sig, pairs_any_provenance, by = c("gene_id", "snps"))

# Join back to unique pairs for testability + classification
any_table <- unique_pairs %>%
  mutate(
    gene_expressed = phenotype_id %in% gene_set_any,
    snp_genotyped  = variant_id %in% snp_set_any
  ) %>%
  left_join(
    pairs_any,
    by = c("phenotype_id" = "gene_id", "variant_id" = "snps")
  ) %>%
  mutate(
    final_category = case_when(
      !gene_expressed & !snp_genotyped                        ~ "untestable_not_genotyped_not_expressed",
      snp_genotyped & !gene_expressed                         ~ "untestable_genotyped_not_expressed",
      !snp_genotyped & gene_expressed                         ~ "untestable_not_genotyped_expressed",
      gene_expressed & snp_genotyped & is.na(significant_any) ~ "testable_pair_untested",
      significant_any                                          ~ "testable_shared",
      !significant_any                                         ~ "testable_specific",
      TRUE                                                     ~ "uncategorised_check_logic"
    )
  )

message("\n'Any' tier summary:")
any_summary <- any_table %>%
  count(final_category, name = "n") %>%
  mutate(prop = round(n / sum(n), 3))
print(any_summary)

n_fetal_specific <- sum(any_table$final_category == "testable_specific")
message("Fetal-specific unique (gene, SNP) pairs: ", n_fetal_specific,
        " out of ", nrow(unique_pairs), " total unique pairs")

# Sanity check
if (any(any_table$final_category == "uncategorised_check_logic")) {
  message("WARNING: ", sum(any_table$final_category == "uncategorised_check_logic"),
          " 'any' tier rows fell through to uncategorised.")
}

# Write tracked output -------------------------------------------------------------------
write_rds(
  list(
    any_table       = any_table,
    any_summary     = any_summary,
    pairs_all       = pairs_all,
    gene_membership = gene_membership,
    snp_membership  = snp_membership
  ),
  output_file
)
message("\nSpecificity classification written to: ", output_file)

# Untracked TSVs for readability outside R -----------------------------------------------
output_dir <- dirname(output_file)
output_stem <- tools::file_path_sans_ext(basename(output_file))

write_tsv(any_table,       file.path(output_dir, paste0(output_stem, "_any_table.tsv")))
write_tsv(any_summary,     file.path(output_dir, paste0(output_stem, "_any_summary.tsv")))
write_tsv(pairs_all,       file.path(output_dir, paste0(output_stem, "_pairs_all.tsv")))
write_tsv(gene_membership, file.path(output_dir, paste0(output_stem, "_gene_membership.tsv")))
write_tsv(snp_membership,  file.path(output_dir, paste0(output_stem, "_snp_membership.tsv")))
message("\nAny-tier tables (TSV) written to: ", output_dir)
