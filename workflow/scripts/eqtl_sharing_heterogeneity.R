#--------------------------------------------------------------------------------------
#
#    Heterogeneity test (Cochran's Q) on fetal-specific eQTL
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: heterogeneity_test (single job)
#            Upstream:   extract_eqtl_universe; classify_eqtl_specificity;
#                        05TENSORQTL convert_genotypes (pvar)
#            Downstream: ocr_enrich; report (Supp Table 5)
#
# Purpose:   For each fetal-specific (gene, SNP) pair:
#            - Keep one fetal row: the cell type (L1 or L2) with the lowest qval
#            - Join it to every Jang et al. 2026 cell type the pair was tested in
#            - Harmonise alleles to our REF/ALT (palindromic pairs dropped)
#            - Cochran's Q on fetal vs adult beta; BH FDR pooled across all rows
#            - Call a pair heterogeneous if FDR < fdr_thresh in ALL its tested
#              adult comparisons, and magnitude-consistent if |fetal beta| >
#              |adult beta| in all of them
#
# Inputs:    eqtl_universe_file   fetal eGene-eSNP universe (extract_eqtl_universe)
#            classification_file  classify_eqtl_specificity output
#            pvar_file            our genotype pvar (REF/ALT per rsID)
#            gene_lookup_file     hg38 gene lookup (gene symbols)
#            fdr_thresh           FDR threshold for the Q test (params)
#
# Outputs:   output_file  rds list: usable, heterogeneous_pairs, heterogeneous_snps,
#                         magnitude_consistent_pairs, magnitude_inconsistent_pairs,
#                         stringent_magnitude, allele_status_summary, funnel
#                         (plus untracked .tsv copies)
#
# Notes:     - heterogeneous_snps (fed to ocr_enrich) are pairs heterogeneous in
#              all comparisons AND magnitude-consistent (218 in the paper run).
#            - Jang has no allele frequencies, so palindromic (A/T, C/G) pairs
#              cannot be resolved and are dropped.
#            - NPC is included; biologically_matched is a diagnostic flag only.
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

message("\n\nRunning heterogeneity test on fetal-specific eQTL vs Jang et al. 2026 (Step C, stringent 'all' rule + magnitude consistency) ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
eqtl_universe_file  <- snakemake@input[["eqtl_universe"]]
classification_file <- snakemake@input[["classification"]]
pvar_file           <- snakemake@input[["pvar_file"]]
gene_lookup_file    <- snakemake@input[["gene_lookup"]]
output_file         <- snakemake@output[[1]]
fdr_thresh          <- snakemake@params[["fdr_thresh"]]

if (is.null(fdr_thresh)) fdr_thresh <- 0.05

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_universe_file", "classification_file", "pvar_file",
               "gene_lookup_file", "fdr_thresh", "output_file"),
  value    = c(eqtl_universe_file, classification_file, pvar_file,
               gene_lookup_file, fdr_thresh, output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load universe and classification  -------------------------------------------------

message("\nLoading eQTL universe (A) and classification (B) ...")
eqtl_universe  <- read_rds(eqtl_universe_file)
classification <- read_rds(classification_file)
any_table <- classification$any_table
pairs_all <- classification$pairs_all

message("Loading gene symbol lookup ...")
gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

n_pairs_specific_B <- any_table %>% filter(final_category == "testable_specific") %>% nrow()
message("  B: unique (gene, SNP) pairs classified testable_specific: ", n_pairs_specific_B)

##  One fetal row per pair (lowest qval)  ---------------------------------------------

message("\n--- Collapsing to single most-significant fetal row per pair ---\n")

specific_pairs <- any_table %>%
  filter(final_category == "testable_specific") %>%
  distinct(phenotype_id, variant_id)

rows_specific_all <- eqtl_universe %>%
  semi_join(specific_pairs, by = c("phenotype_id", "variant_id"))

n_fetal_celltypes_significant <- rows_specific_all %>%
  distinct(phenotype_id, variant_id, cell_type) %>%
  count(phenotype_id, variant_id, name = "n_fetal_celltypes_significant")

fetal_ref_row <- rows_specific_all %>%
  group_by(phenotype_id, variant_id) %>%
  slice_min(qval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(n_fetal_celltypes_significant, by = c("phenotype_id", "variant_id"))

n_unique_pairs_specific <- nrow(fetal_ref_row)
message("  Unique (gene, SNP) pairs after collapsing to single fetal reference row: ",
        n_unique_pairs_specific)

if (n_unique_pairs_specific != n_pairs_specific_B) {
  message("  WARNING: pairs after collapse (", n_unique_pairs_specific,
          ") does not match B's testable_specific count (", n_pairs_specific_B,
          "). Check A/B are from the same run.")
}

message("\nDistribution of fetal cell types a pair was significant in (before collapse):")
n_fetal_celltypes_significant %>%
  mutate(n_bucket = if_else(n_fetal_celltypes_significant >= 3, "3+", as.character(n_fetal_celltypes_significant))) %>%
  count(n_bucket, name = "n_pairs") %>%
  arrange(n_bucket) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

message("\nFetal reference cell type chosen, by frequency:")
fetal_ref_row %>% count(cell_type, sort = TRUE) %>% print(n = Inf)

##  Join to every tested Jang cell type  ----------------------------------------------

message("\n--- Joining against all available Jang comparisons ---\n")

pairs_all_lookup <- pairs_all %>%
  select(phenotype_id = gene_id, variant_id, jang_cell_type,
         beta_jang = beta, se_jang = se,
         ref_jang = ref, alt_jang = alt)

joined <- fetal_ref_row %>%
  inner_join(pairs_all_lookup, by = c("phenotype_id", "variant_id"))

n_joined <- nrow(joined)
message("  C: joined rows (1 fetal row per pair x all tested Jang comparisons): ", n_joined)

joined <- joined %>%
  mutate(
    biologically_matched = case_when(
      is.na(jang_label_match) ~ FALSE,
      TRUE ~ jang_cell_type == jang_label_match
    )
  )

message("\nDistribution of adult (Jang) comparisons tested per pair:")
joined %>%
  count(phenotype_id, variant_id, name = "n_adult_tested") %>%
  mutate(n_bucket = if_else(n_adult_tested >= 4, "4+", as.character(n_adult_tested))) %>%
  count(n_bucket, name = "n_pairs") %>%
  arrange(n_bucket) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

##  Allele harmonisation  -------------------------------------------------------------

message("\n--- Allele harmonization (REF/ALT only, no AF check) ---\n")

pvar_lookup <- fread(
  pvar_file,
  select = c("ID", "REF", "ALT"),
  col.names = c("variant_id", "ref_my", "alt_my"),
  showProgress = FALSE
)

joined <- joined %>%
  rename(beta_my = slope, se_my = slope_se, af_my = af) %>%
  left_join(as_tibble(pvar_lookup), by = "variant_id")

joined <- joined %>%
  mutate(
    is_palindromic =
      (ref_my == "A" & alt_my == "T") | (ref_my == "T" & alt_my == "A") |
      (ref_my == "C" & alt_my == "G") | (ref_my == "G" & alt_my == "C"),
    allele_status        = get_allele_status(ref_my, alt_my, ref_jang, alt_jang),
    beta_jang_harmonized = harmonise_beta(beta_jang, allele_status)
  )

status_tbl <- joined %>% count(allele_status, name = "N") %>% mutate(prop = round(N / sum(N), 4))
message("Allele status breakdown:")
print(status_tbl)

n_allele_dropped <- sum(is.na(joined$beta_jang_harmonized))
usable <- joined %>% filter(!is.na(beta_jang_harmonized))
n_usable <- nrow(usable)
n_unique_pairs_usable <- usable %>% distinct(phenotype_id, variant_id) %>% nrow()

message("\n  C: dropped -- unresolved allele harmonization (mismatched/palindromic/missing): ",
        n_allele_dropped)
message("  C: usable rows entering Q-test: ", n_usable)
message("  C: unique (gene, SNP) pairs among usable rows: ", n_unique_pairs_usable)

n_pairs_lost_some_comparisons <- joined %>%
  distinct(phenotype_id, variant_id) %>%
  anti_join(
    usable %>%
      group_by(phenotype_id, variant_id) %>%
      summarise(n_usable = n(), .groups = "drop") %>%
      inner_join(
        joined %>% count(phenotype_id, variant_id, name = "n_total"),
        by = c("phenotype_id", "variant_id")
      ) %>%
      filter(n_usable == n_total),
    by = c("phenotype_id", "variant_id")
  ) %>%
  nrow()
message("  Pairs where >=1 (but not all) adult comparisons were dropped for allele ambiguity: ",
        n_pairs_lost_some_comparisons,
        " -- the 'all' rule below applies only across each pair's SURVIVING comparisons.")

##  Cochran's Q and FDR  --------------------------------------------------------------

message("\n--- Running Cochran's Q ---\n")

usable <- usable %>%
  mutate(
    delta_beta = beta_my - beta_jang_harmonized,
    se_delta   = sqrt(se_my^2 + se_jang^2),
    Q          = (delta_beta / se_delta)^2,
    p          = pchisq(Q, df = 1, lower.tail = FALSE),
    q          = p.adjust(p, method = "BH"),
    heterogeneous = q < fdr_thresh
  )

n_het_rows <- sum(usable$heterogeneous)
message("  Heterogeneous rows (FDR<", fdr_thresh, "): ", n_het_rows,
        " (", round(100 * n_het_rows / n_usable, 2), "%)")

##  Pair-level call: heterogeneous in all comparisons  --------------------------------

message("\n--- Aggregating to final (gene, SNP)-level heterogeneous call (STRINGENT: all()) ---\n")

heterogeneous_pairs <- usable %>%
  group_by(phenotype_id, variant_id) %>%
  summarise(
    fetal_cell_type              = first(cell_type),
    n_comparisons_tested          = n(),
    n_comparisons_heterogeneous    = sum(heterogeneous),
    heterogeneous_all               = all(heterogeneous),
    confirmed_by_matched_type        = any(biologically_matched & heterogeneous),
    matched_type_tested               = any(biologically_matched),
    .groups = "drop"
  )

message("\nPass/fail by number of comparisons tested (shows the stringency effect --",
        " pairs tested against more adult cell types are harder to pass):")
heterogeneous_pairs %>%
  mutate(n_bucket = if_else(n_comparisons_tested >= 4, "4+", as.character(n_comparisons_tested))) %>%
  group_by(n_bucket) %>%
  summarise(n_pairs = n(), n_passed_all = sum(heterogeneous_all), .groups = "drop") %>%
  arrange(n_bucket) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

n_het_pairs <- sum(heterogeneous_pairs$heterogeneous_all)
message("\n  Unique (gene, SNP) pairs heterogeneous in ALL tested comparisons: ", n_het_pairs,
        " / ", nrow(heterogeneous_pairs))
message("  Of which the biologically matched comparison was itself tested and confirmed: ",
        sum(heterogeneous_pairs$heterogeneous_all & heterogeneous_pairs$confirmed_by_matched_type))
message("  Of which NO biologically matched comparison was available at all (incl. all NPC pairs): ",
        sum(heterogeneous_pairs$heterogeneous_all & !heterogeneous_pairs$matched_type_tested))

##  Magnitude consistency  ------------------------------------------------------------

message("\n--- Magnitude consistency (|fetal beta| > |adult beta|) among stringent pairs ---\n")

stringent_magnitude <- usable %>%
  semi_join(heterogeneous_pairs %>% filter(heterogeneous_all),
            by = c("phenotype_id", "variant_id")) %>%
  mutate(
    magnitude_diff = abs(beta_my) - abs(beta_jang_harmonized),
    fetal_bigger    = magnitude_diff > 0
  )

pair_magnitude <- stringent_magnitude %>%
  group_by(phenotype_id, variant_id) %>%
  summarise(
    magnitude_consistent = all(fetal_bigger),
    min_magnitude_diff    = min(magnitude_diff),
    max_magnitude_diff    = max(magnitude_diff),
    .groups = "drop"
  )

heterogeneous_pairs <- heterogeneous_pairs %>%
  left_join(pair_magnitude, by = c("phenotype_id", "variant_id"))
# magnitude_consistent/min/max_magnitude_diff are only meaningful for
# heterogeneous_all == TRUE pairs; NA for pairs that never reached this stage.

n_pass_magnitude <- sum(pair_magnitude$magnitude_consistent)
n_stringent       <- nrow(pair_magnitude)
message("  ", n_pass_magnitude, " of ", n_stringent, " stringent pairs (",
        round(100 * n_pass_magnitude / n_stringent, 1),
        "%) also have |fetal beta| > |adult beta| in EVERY tested comparison.")

# Gene-symbol-annotated tables for the two groups
magnitude_consistent_pairs <- heterogeneous_pairs %>%
  filter(heterogeneous_all, magnitude_consistent) %>%
  left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) %>%
  rename(gene_symbol = external_gene_name) %>%
  relocate(gene_symbol, .after = phenotype_id)

magnitude_inconsistent_pairs <- heterogeneous_pairs %>%
  filter(heterogeneous_all, !magnitude_consistent) %>%
  left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) %>%
  rename(gene_symbol = external_gene_name) %>%
  relocate(gene_symbol, .after = phenotype_id)

message("  magnitude_consistent_pairs: ", nrow(magnitude_consistent_pairs), " rows")
message("  magnitude_inconsistent_pairs: ", nrow(magnitude_inconsistent_pairs), " rows")

n_missing_symbol_consistent <- sum(is.na(magnitude_consistent_pairs$gene_symbol))
n_missing_symbol_inconsistent <- sum(is.na(magnitude_inconsistent_pairs$gene_symbol))
if (n_missing_symbol_consistent > 0 | n_missing_symbol_inconsistent > 0) {
  message("  NOTE: ", n_missing_symbol_consistent, " / ", nrow(magnitude_consistent_pairs),
          " magnitude-consistent and ", n_missing_symbol_inconsistent, " / ",
          nrow(magnitude_inconsistent_pairs),
          " magnitude-inconsistent pairs had no matching gene symbol in the lookup.")
} else {
  message("  All rows in both tables matched a gene symbol in the lookup.")
}

# SNP set fed to ocr_enrich: heterogeneous in all comparisons AND magnitude-consistent
heterogeneous_snps <- heterogeneous_pairs %>%
  filter(heterogeneous_all, magnitude_consistent) %>%
  distinct(phenotype_id, variant_id)

message("\n  Final heterogeneous_snps fed to OCR enrichment (heterogeneous_all AND magnitude_consistent): ",
        nrow(heterogeneous_snps))

##  Funnel summary  -------------------------------------------------------------------

message("\n--- Population funnel ---\n")

funnel <- tibble(
  stage = c(
    "B: unique (gene, SNP) pairs classified testable_specific",
    "C: unique pairs after collapsing to single fetal reference row",
    "C: joined rows (1 fetal row/pair x all tested adult comparisons)",
    "C: dropped -- unresolved allele harmonization",
    "C: usable rows entering Q-test",
    "C: unique (gene, SNP) pairs among usable rows",
    "C: heterogeneous rows (FDR<thresh)",
    "C: unique (gene, SNP) pairs heterogeneous in ALL tested comparisons (stringent)",
    "C: of which also magnitude-consistent (final, fed to OCR)"
  ),
  N = c(n_pairs_specific_B, n_unique_pairs_specific, n_joined,
        n_allele_dropped, n_usable, n_unique_pairs_usable, n_het_rows, n_het_pairs,
        nrow(heterogeneous_snps))
)
print(funnel)

##  Write output  ---------------------------------------------------------------------

usable_out <- usable %>%
  select(cell_type, level, parent_cell_type, phenotype_id, variant_id,
         n_fetal_celltypes_significant,
         jang_cell_type, jang_label_match, biologically_matched,
         beta_my, se_my, af_my,
         beta_jang, beta_jang_harmonized, se_jang,
         ref_my, alt_my, ref_jang, alt_jang, allele_status,
         delta_beta, se_delta, Q, p, q, heterogeneous)

write_rds(
  list(
    usable                       = usable_out,
    heterogeneous_pairs           = heterogeneous_pairs,
    heterogeneous_snps             = heterogeneous_snps,
    magnitude_consistent_pairs      = magnitude_consistent_pairs,
    magnitude_inconsistent_pairs     = magnitude_inconsistent_pairs,
    stringent_magnitude                = stringent_magnitude,
    allele_status_summary                 = status_tbl,
    funnel                                   = funnel
  ),
  output_file
)
message("\nHeterogeneity test results written to: ", output_file)

# Untracked TSVs for readability outside R
output_dir  <- dirname(output_file)
output_stem <- tools::file_path_sans_ext(basename(output_file))

write_tsv(usable_out,                    file.path(output_dir, paste0(output_stem, "_usable.tsv")))
write_tsv(heterogeneous_pairs,           file.path(output_dir, paste0(output_stem, "_heterogeneous_pairs.tsv")))
write_tsv(magnitude_consistent_pairs,    file.path(output_dir, paste0(output_stem, "_magnitude_consistent_pairs.tsv")))
write_tsv(magnitude_inconsistent_pairs,  file.path(output_dir, paste0(output_stem, "_magnitude_inconsistent_pairs.tsv")))
write_tsv(stringent_magnitude,           file.path(output_dir, paste0(output_stem, "_stringent_magnitude.tsv")))
write_tsv(status_tbl,                    file.path(output_dir, paste0(output_stem, "_allele_status_summary.tsv")))
write_tsv(funnel,                        file.path(output_dir, paste0(output_stem, "_funnel.tsv")))
message("Tables (TSV) written to: ", output_dir)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
