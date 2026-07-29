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
#  Two analyses with different units:
#
#    - "any":           operates on the 6,844 unique (gene, SNP) pairs across all
#                        cell types. Tested against all 7 Fujita cell types pooled;
#                        "shared" = significant in >=1 Fujita type. Produces the
#                        fetal-specific SNP list for the OCR enrichment step.
#
#    - "cell_specific": operates per fetal cell type. Each cell type's own pairs
#                        tested against its a priori label-matched Fujita cell type
#                        only (e.g. GABA -> Inh). NPC excluded (no adult equivalent).
#                        Informative even where overlap is small.
#
#  Step 2 of the developmental-specificity pipeline (independent for now; designed
#  to fold into replication_pi1_enrichment.R later if retained for the paper).

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

# L1 cell types
cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

# A priori biological correspondence between fetal L1 cell types and Fujita adult
# cell types. NPC has no defensible adult equivalent (no progenitor population in
# the adult cortex), so it is excluded from cell_specific analysis entirely.
fugita_label_match <- c(
  "Glu-UL"    = "Exc",
  "Glu-DL"    = "Exc",
  "GABA"      = "Inh",
  "NPC"       = NA_character_,
  "MG"        = "Mic",
  "OPC"       = "OPC",
  "Endo-Peri" = "End"
)

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

# Strip L2 sub-cluster suffix to get parent L1 cell type
get_parent_cell_type <- function(cell_type, cell_types_L1) {
  if (cell_type %in% cell_types_L1) return(cell_type)
  sub("-[^-]+$", "", cell_type)
}

read_fugita_celltype <- function(file_path, fugita_ct, universe_genes, universe_snps) {
  message("Reading Fujita data for: ", fugita_ct, " (", file_path, ")")

  dat <- read_tsv(
    file_path,
    col_select = c(gene_id, snps, pvalue, beta, significant_by_2step_FDR),
    col_types = cols(gene_id = "c", snps = "c", pvalue = "d", beta = "d",
                      significant_by_2step_FDR = "c"),
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

# Per-gene/SNP membership lookups for cell_specific testability checks
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

# Provenance: which Fujita cell type gave the smallest p-value (for QC/reporting)
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
      significant_any & sign(pvalue_adult) >= 0               ~ ifelse(
        !is.na(beta_adult), "testable_shared", "testable_shared"
      ),
      significant_any                                          ~ "testable_shared",
      !significant_any                                         ~ "testable_specific",
      TRUE                                                     ~ "uncategorised_check_logic"
    )
  )

# No concordance check for the "any" tier -- there's no single fetal slope to compare
# against, since the same pair may appear in multiple fetal cell types with different
# effect sizes. Concordance is assessed in the cell_specific tier instead.

message("\n'Any' tier summary:")
any_summary <- any_table %>%
  count(final_category, name = "n") %>%
  mutate(prop = round(n / sum(n), 3))
print(any_summary)

n_fetal_specific <- sum(any_table$final_category == "testable_specific")
message("Fetal-specific unique (gene, SNP) pairs: ", n_fetal_specific,
        " out of ", nrow(unique_pairs), " total unique pairs")

# =======================================================================================
# "CELL_SPECIFIC" TIER: per fetal cell type, matched to label-matched Fujita type
# =======================================================================================
message("\n--- Building 'cell_specific' tier classification ---")

# Add parent cell type and label match to each row of the full universe
cs_table <- eqtl_universe %>%
  mutate(
    parent_cell_type = map_chr(cell_type, get_parent_cell_type, cell_types_L1 = cell_types_L1),
    ref_cell_type = fugita_label_match[parent_cell_type]
  )

# Exclude NPC (and NPC subtypes) -- no adult equivalent
cs_table_testable <- cs_table %>% filter(!is.na(ref_cell_type))
cs_table_npc      <- cs_table %>% filter(is.na(ref_cell_type))
message("  Rows with label match: ", nrow(cs_table_testable),
        " | NPC rows excluded: ", nrow(cs_table_npc))

# Check testability: gene expressed and SNP genotyped in the SPECIFIC matched Fujita type
cs_table_testable <- cs_table_testable %>%
  left_join(
    gene_membership %>% mutate(gene_expressed = TRUE),
    by = c("phenotype_id" = "gene_id", "ref_cell_type" = "fugita_cell_type")
  ) %>%
  left_join(
    snp_membership %>% mutate(snp_genotyped = TRUE),
    by = c("variant_id" = "snps", "ref_cell_type" = "fugita_cell_type")
  ) %>%
  mutate(
    gene_expressed = replace_na(gene_expressed, FALSE),
    snp_genotyped  = replace_na(snp_genotyped, FALSE)
  )

# Look up the exact pair in the matched Fujita cell type
cs_table_testable <- cs_table_testable %>%
  left_join(
    pairs_all %>%
      transmute(phenotype_id = gene_id, variant_id = snps,
                ref_cell_type = fugita_cell_type,
                significant_ref = significant_by_2step_FDR == "Yes",
                pvalue_adult = pvalue, beta_adult = beta),
    by = c("phenotype_id", "variant_id", "ref_cell_type")
  )

# Classify
cs_table_testable <- cs_table_testable %>%
  mutate(
    final_category = case_when(
      !gene_expressed & !snp_genotyped                          ~ "untestable_not_genotyped_not_expressed",
      snp_genotyped & !gene_expressed                           ~ "untestable_genotyped_not_expressed",
      !snp_genotyped & gene_expressed                           ~ "untestable_not_genotyped_expressed",
      gene_expressed & snp_genotyped & is.na(significant_ref)   ~ "testable_pair_untested",
      significant_ref & sign(slope) == sign(beta_adult)         ~ "testable_shared_concordant",
      significant_ref & sign(slope) != sign(beta_adult)         ~ "testable_shared_discordant",
      !significant_ref                                          ~ "testable_specific",
      TRUE                                                       ~ "uncategorised_check_logic"
    )
  )

# Add NPC rows back with their own category
cs_table_npc <- cs_table_npc %>%
  mutate(
    ref_cell_type = NA_character_,
    gene_expressed = NA,
    snp_genotyped = NA,
    significant_ref = NA,
    pvalue_adult = NA_real_,
    beta_adult = NA_real_,
    final_category = "no_adult_equivalent"
  )

cs_table <- bind_rows(cs_table_testable, cs_table_npc)

# Summary per cell type
message("\n'Cell-specific' tier summary (per cell type):")
cs_summary <- cs_table %>%
  count(cell_type, ref_cell_type, final_category, name = "n") %>%
  group_by(cell_type) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup()
print(cs_summary, n = Inf)

# Overall cell_specific summary
message("\n'Cell-specific' tier summary (all cell types combined):")
print(
  cs_table %>%
    count(final_category, name = "n") %>%
    mutate(prop = round(n / sum(n), 3))
)

# Sanity checks
if (any(any_table$final_category == "uncategorised_check_logic")) {
  message("WARNING: ", sum(any_table$final_category == "uncategorised_check_logic"),
          " 'any' tier rows fell through to uncategorised.")
}
if (any(cs_table$final_category == "uncategorised_check_logic")) {
  message("WARNING: ", sum(cs_table$final_category == "uncategorised_check_logic"),
          " 'cell_specific' tier rows fell through to uncategorised.")
}

# Write single tracked output ----------------------------------------------------------
write_rds(
  list(
    any_table = any_table,
    any_summary = any_summary,
    cell_specific_table = cs_table,
    cell_specific_summary = cs_summary
  ),
  output_file
)
message("\nSpecificity classification written to: ", output_file)
