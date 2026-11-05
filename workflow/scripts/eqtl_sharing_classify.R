#--------------------------------------------------------------------------------------
#
#    Classify fetal eQTL as shared with adult brain or fetal-specific
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: classify_eqtl_specificity (single job)
#            Upstream:   extract_eqtl_universe; Jang et al. 2026 SingleBrain full and
#                        top association files (all L1 + L2 cell types)
#            Downstream: heterogeneity_test; ocr_enrich; glm_enrich; report
#                        (Supp Table 4)
#
# Purpose:   For each unique fetal eGene-eSNP pair:
#            - Testable if the gene is expressed and the SNP genotyped in any Jang
#              cell type
#            - Shared if the pair's own Random_P is below a Bonferroni threshold
#              (0.05 / n Jang cell types / n testable pairs) in any Jang cell type;
#              otherwise fetal-specific
#            - Record the adult cell type with the smallest P and whether the
#              pair is Jang's own lead SNP there (informational only)
#
# Inputs:    eqtl_universe_file  fetal eGene-eSNP universe (extract_eqtl_universe)
#            jang_full_files     Jang full-association files, one per cell type
#            jang_top_files      Jang top-association files (lead-SNP identity only)
#            jang_cell_types     Jang cell types, same order as the files (params)
#
# Outputs:   output_file  rds list: any_table, any_summary, pairs_all,
#                         gene_membership, snp_membership, bonferroni_threshold,
#                         n_testable_pairs (plus untracked .tsv copies)
#
# Notes:     - Jang cell types exclude MiGA3 (separate microglia meta-analysis).
#            - Jang betas are sign-normalised to alt (align_jang_beta()) and
#              placeholder-cohort rows removed (filter_jang_artifacts()).
#            - n testable pairs is computed each run (3,479 for the paper run,
#              P < 4.106e-7); the Methods text gives 2,989 (P < 4.8e-7).
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

message("\n\nClassifying fetal eGene-eSNP pairs: testability + shared/specific vs. Jang et al. 2026 ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
jang_full_files    <- snakemake@input[["jang_full_files"]]
jang_top_files     <- snakemake@input[["jang_top_files"]]
output_file        <- snakemake@output[[1]]
jang_cell_types    <- snakemake@params[["jang_cell_types"]]

names(jang_full_files) <- jang_cell_types
names(jang_top_files)  <- jang_cell_types

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_universe_file", "n_jang_full_files", "n_jang_top_files", "output_file"),
  value    = c(eqtl_universe_file, length(jang_full_files), length(jang_top_files), output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  ------------------------------------------------------------------------

# Full-association file: every tested pair for one Jang cell type, filtered to the
# universe genes/SNPs; Random_P is the significance check for the pair
read_jang_full <- function(file_path, jang_ct, universe_genes, universe_snps) {
  message("Reading Jang full-association data for: ", jang_ct, " (", file_path, ")")

  dat <- read_tsv(
    file_path,
    col_select = c(feature, variant_id, ref, alt, Allele, fixed_beta, fixed_sd,
                    Fixed_P, Random_P),
    col_types = cols(feature = "c", variant_id = "c", ref = "c", alt = "c",
                      Allele = "c", fixed_beta = "d", fixed_sd = "d",
                      Fixed_P = "d", Random_P = "d"),
    show_col_types = FALSE
  )

  dat <- filter_jang_artifacts(dat, "Jang full-association")
  dat <- dat %>% select(-Fixed_P)

  # Strip GENCODE version suffix (Jang's feature is e.g. "ENSG00000131591.18";
  # our universe genes are bare ENSG IDs)
  dat <- dat %>%
    mutate(gene_id = sub("\\..*", "", feature)) %>%
    filter(str_detect(variant_id, "^rs"), str_detect(gene_id, "^ENSG"))

  gene_set <- unique(dat$gene_id)
  snp_set  <- unique(dat$variant_id)

  pair_lookup <- dat %>%
    filter(gene_id %in% universe_genes, variant_id %in% universe_snps) %>%
    distinct(gene_id, variant_id, .keep_all = TRUE) %>%
    mutate(
      # Sign-normalize: Jang's beta is relative to Allele, not always alt
      beta = align_jang_beta(fixed_beta, Allele, alt),
      se   = fixed_sd,
      random_p = Random_P,
      jang_cell_type = jang_ct
    ) %>%
    select(gene_id, variant_id, beta, se, random_p, ref, alt, jang_cell_type)

  message("  Genes tested: ", length(gene_set),
          " | SNPs tested: ", length(snp_set),
          " | Relevant pairs retained: ", nrow(pair_lookup))

  list(gene_set = gene_set, snp_set = snp_set, pair_lookup = pair_lookup)
}

# Top-association file: Jang's own lead SNP per gene, for the informational
# is_jang_lead_snp flag only (not used for significance)
read_jang_top <- function(file_path, jang_ct) {
  message("Reading Jang top-association data for: ", jang_ct, " (", file_path, ") -- lead-SNP identity only")

  dat <- read_tsv(
    file_path,
    col_select = c(feature, variant_id, Fixed_P, Random_P),
    col_types = cols(feature = "c", variant_id = "c", Fixed_P = "d", Random_P = "d"),
    show_col_types = FALSE
  )

  dat <- filter_jang_artifacts(dat, "Jang top-association")

  dat %>%
    mutate(gene_id = sub("\\..*", "", feature)) %>%
    filter(str_detect(variant_id, "^rs"), str_detect(gene_id, "^ENSG")) %>%
    transmute(gene_id, lead_variant_id = variant_id, jang_cell_type = jang_ct)
}

##  Load fetal universe  --------------------------------------------------------------

message("\nLoading fetal eGene/eSNP universe ...")
eqtl_universe <- read_rds(eqtl_universe_file)
message("  Rows: ", nrow(eqtl_universe), " | Unique eGenes: ", n_distinct(eqtl_universe$phenotype_id))

universe_genes <- unique(eqtl_universe$phenotype_id)
universe_snps  <- unique(eqtl_universe$variant_id)

# Deduplicate to unique (gene, SNP) pairs for the "any" tier
unique_pairs <- eqtl_universe %>%
  distinct(phenotype_id, variant_id)
message("  Unique (gene, SNP) pairs: ", nrow(unique_pairs))

##  Read Jang full-association files  -------------------------------------------------
message("\n--- Reading Jang full-association files (", length(jang_cell_types), " cell types) ---\n")

gene_sets    <- list()
snp_sets     <- list()
pair_lookups <- list()

for (ct in jang_cell_types) {
  res <- read_jang_full(jang_full_files[[ct]], ct, universe_genes, universe_snps)
  gene_sets[[ct]]    <- res$gene_set
  snp_sets[[ct]]     <- res$snp_set
  pair_lookups[[ct]] <- res$pair_lookup

  rm(res)
  gc()
}

gene_set_any <- Reduce(union, gene_sets)
snp_set_any  <- Reduce(union, snp_sets)

pairs_all <- bind_rows(pair_lookups)
message("\nTotal pairs_all rows (across all cell types): ", nrow(pairs_all))

# Per-gene/SNP membership lookups, used for the "any" tier's testability check,
# and written to output for the heterogeneity test.
gene_membership <- map_dfr(names(gene_sets), ~ tibble(jang_cell_type = .x, gene_id = gene_sets[[.x]]))
snp_membership  <- map_dfr(names(snp_sets),  ~ tibble(jang_cell_type = .x, variant_id = snp_sets[[.x]]))

##  Bonferroni threshold  -------------------------------------------------------------
# Testable = gene expressed AND SNP genotyped in any Jang cell type
message("\n--- Computing dynamic Bonferroni threshold ---\n")

n_testable_pairs <- unique_pairs %>%
  filter(phenotype_id %in% gene_set_any, variant_id %in% snp_set_any) %>%
  nrow()

n_jang_cell_types <- length(jang_cell_types)

bonferroni_threshold <- 0.05 / n_jang_cell_types / n_testable_pairs

message("  Testable pairs (gene_expressed & snp_genotyped both TRUE): ", n_testable_pairs)
message("  Jang cell types: ", n_jang_cell_types)
message("  Bonferroni threshold = 0.05 / ", n_jang_cell_types, " / ", n_testable_pairs,
        " = ", format(bonferroni_threshold, scientific = TRUE, digits = 4))

##  Jang lead-SNP identity (informational only)  --------------------------------------
message("\n--- Reading Jang top-association files (lead-SNP identity, informational only) ---\n")

top_lookups <- map(jang_cell_types, ~ read_jang_top(jang_top_files[[.x]], .x))
top_lookup_all <- bind_rows(top_lookups)
message("Total lead-SNP identity rows across all cell types: ", nrow(top_lookup_all))

pairs_all <- pairs_all %>%
  left_join(top_lookup_all, by = c("gene_id", "jang_cell_type")) %>%
  mutate(is_jang_lead_snp = !is.na(lead_variant_id) & variant_id == lead_variant_id) %>%
  select(-lead_variant_id)

n_lead_match <- sum(pairs_all$is_jang_lead_snp)
message("Rows in pairs_all that also match Jang's own lead SNP for that gene/cell type: ",
        n_lead_match, " / ", nrow(pairs_all),
        " (", round(100 * n_lead_match / nrow(pairs_all), 2), "%) -- informational, does not affect classification")

##  Classify unique pairs against all Jang cell types  --------------------------------
message("\n--- Building 'any' tier classification ---")

# Significant = the pair's own Random_P below the Bonferroni threshold
pairs_all <- pairs_all %>%
  mutate(significant_pair = !is.na(random_p) & random_p < bonferroni_threshold)

message("Rows significant (Random_P < ",
        format(bonferroni_threshold, scientific = TRUE, digits = 4), "): ",
        sum(pairs_all$significant_pair))
message("  Of which also Jang's own lead SNP: ",
        sum(pairs_all$significant_pair & pairs_all$is_jang_lead_snp))
message("  Of which NOT Jang's own lead SNP (significant via a non-lead SNP -- ",
        "exactly the case the old lead-SNP gate would have missed): ",
        sum(pairs_all$significant_pair & !pairs_all$is_jang_lead_snp))

# Significance call: significant_pair in >= 1 Jang cell type
pairs_any_sig <- pairs_all %>%
  group_by(gene_id, variant_id) %>%
  summarise(
    significant_any     = any(significant_pair),
    n_jang_types_tested  = n(),
    any_jang_lead_snp_match = any(is_jang_lead_snp),
    .groups = "drop"
  )

# Provenance: the Jang cell type with the smallest random_p for this pair
pairs_any_provenance <- pairs_all %>%
  group_by(gene_id, variant_id) %>%
  slice_min(random_p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene_id, variant_id, best_jang_ct = jang_cell_type,
         random_p_adult = random_p, beta_adult = beta,
         best_ct_is_jang_lead_snp = is_jang_lead_snp)

pairs_any <- left_join(pairs_any_sig, pairs_any_provenance, by = c("gene_id", "variant_id"))

# Join back to unique pairs for testability + classification
any_table <- unique_pairs %>%
  mutate(
    gene_expressed = phenotype_id %in% gene_set_any,
    snp_genotyped  = variant_id %in% snp_set_any
  ) %>%
  left_join(
    pairs_any,
    by = c("phenotype_id" = "gene_id", "variant_id" = "variant_id")
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

message("\n--- Lead-SNP match rate among testable_shared pairs (informational) ---")
shared_lead_match <- any_table %>%
  filter(final_category == "testable_shared") %>%
  count(best_ct_is_jang_lead_snp, name = "n")
print(shared_lead_match)

# Sanity check
if (any(any_table$final_category == "uncategorised_check_logic")) {
  message("WARNING: ", sum(any_table$final_category == "uncategorised_check_logic"),
          " 'any' tier rows fell through to uncategorised.")
}

##  Write output  ---------------------------------------------------------------------
write_rds(
  list(
    any_table            = any_table,
    any_summary          = any_summary,
    pairs_all            = pairs_all,
    gene_membership      = gene_membership,
    snp_membership       = snp_membership,
    bonferroni_threshold = bonferroni_threshold,
    n_testable_pairs     = n_testable_pairs
  ),
  output_file
)
message("\nSpecificity classification written to: ", output_file)

# Untracked TSVs for readability outside R
output_dir <- dirname(output_file)
output_stem <- tools::file_path_sans_ext(basename(output_file))

write_tsv(any_table,       file.path(output_dir, paste0(output_stem, "_any_table.tsv")))
write_tsv(any_summary,     file.path(output_dir, paste0(output_stem, "_any_summary.tsv")))
write_tsv(pairs_all,       file.path(output_dir, paste0(output_stem, "_pairs_all.tsv")))
write_tsv(gene_membership, file.path(output_dir, paste0(output_stem, "_gene_membership.tsv")))
write_tsv(snp_membership,  file.path(output_dir, paste0(output_stem, "_snp_membership.tsv")))
message("\nAny-tier tables (TSV) written to: ", output_dir)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
