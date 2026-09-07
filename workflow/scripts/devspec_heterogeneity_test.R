#--------------------------------------------------------------------------------------
#
#    Developmental specificity: heterogeneity test on fetal-specific eQTL (Step C)
#    (REVISED -- single most-significant fetal row per pair, ALL tested adult
#    comparisons must be heterogeneous)
#
#--------------------------------------------------------------------------------------

#  For every (gene, SNP) pair classified testable_specific by B, collapses to a
#  SINGLE fetal reference row per pair -- whichever fetal cell type (L1 or L2)
#  gives that pair's lowest permutation qval, i.e. where it is most significant.
#  That one fetal row is then tested against EVERY Jang et al. 2026 adult cell
#  type the pair was actually tested in (not just the a priori biologically-
#  matched type).
#
#  STRINGENT rule: a pair is called fetal-specific-and-heterogeneous only if it
#  is heterogeneous in ALL of its successfully-tested adult comparisons (not
#  just >=1). "Successfully tested" = survived allele harmonization; a
#  comparison dropped for unresolved allele ambiguity cannot be evaluated and is
#  excluded from both the count and the all() requirement -- logged separately
#  per pair so it's clear how many comparisons the "all" rule was actually
#  applied across.
#
#  This is a genuinely stricter bar than the "any" design: a pair tested against
#  more adult cell types has proportionally more chances to fail on just one.
#  The per-pair n_comparisons_tested x pass/fail breakdown is logged explicitly
#  so this effect is visible rather than hidden in a single headline count.
#
#  NPC included (no biological-match requirement for testing; biologically_matched
#  is retained as a per-row diagnostic flag only, always FALSE for NPC).
#
#  FDR (BH) pooled once across all usable rows (one row per pair per tested
#  adult comparison -- no fetal-side duplication now, since each pair
#  contributes exactly one fetal context).
#
#  Step 3 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  MIGRATION NOTE -- Fujita -> Jang et al. 2026 (SingleBrain):
#
#  All fugita_* naming renamed to jang_* throughout (jang_cell_type,
#  jang_label_match, beta_jang, se_jang, ref_jang, alt_jang,
#  beta_jang_harmonized), matching the rename already made in
#  devspec_extract_universe.R and devspec_classify_specificity.R.
#
#  Jang's data has no allele-frequency column (unlike Fujita's ALT_AF), so the
#  previous AF-based palindromic-SNP resolution (aligned/flipped disambiguation
#  via allele frequency for A/T, C/G pairs) has been REMOVED entirely -- not
#  substituted with an external reference (e.g. 1000G), since there is no
#  guarantee Jang's own allele/strand conventions were harmonized against 1000G,
#  and introducing that assumption to paper over a missing check would trade one
#  unvalidated assumption for another. Every palindromic pair is now
#  unconditionally categorised "palindromic_unresolved" and dropped from the
#  Q-test (beta_jang_harmonized = NA), regardless of allele frequency. This is
#  more conservative than the previous af_tol-based partial resolution -- the
#  af_tol parameter is no longer used by this script.
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

message("\n\nRunning heterogeneity test on fetal-specific eQTL vs Jang et al. 2026 (Step C, stringent 'all' rule) ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file  <- snakemake@input[["eqtl_universe"]]
classification_file <- snakemake@input[["classification"]]
pvar_file           <- snakemake@input[["pvar_file"]]
output_file         <- snakemake@output[[1]]
fdr_thresh          <- snakemake@params[["fdr_thresh"]]

if (is.null(fdr_thresh)) fdr_thresh <- 0.05

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "classification_file", "pvar_file",
               "fdr_thresh", "output_file"),
  value    = c(eqtl_universe_file, classification_file, pvar_file,
               fdr_thresh, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# ========================================================================================
# 1. Load A (universe) and B (classification)
# ========================================================================================

message("\nLoading eQTL universe (A) and classification (B) ...")
eqtl_universe  <- read_rds(eqtl_universe_file)
classification <- read_rds(classification_file)
any_table <- classification$any_table
pairs_all <- classification$pairs_all

n_pairs_specific_B <- any_table %>% filter(final_category == "testable_specific") %>% nrow()
message("  B: unique (gene, SNP) pairs classified testable_specific: ", n_pairs_specific_B)

# ========================================================================================
# 2. Collapse to a SINGLE fetal reference row per pair -- lowest qval (most
#    significant) across all fetal cell types (L1+L2, incl. NPC) it appears in
# ========================================================================================

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

# ========================================================================================
# 3. Join the single fetal reference row against EVERY Jang cell type the
#    pair was tested in
# ========================================================================================

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

# ========================================================================================
# 4. Allele harmonization (REF/ALT string comparison only -- no AF check; see
#    MIGRATION NOTE above)
# ========================================================================================

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
    allele_status = case_when(
      is.na(ref_my) | is.na(alt_my) | is.na(ref_jang) | is.na(alt_jang) ~ "missing_allele_info",
      is_palindromic                                                    ~ "palindromic_unresolved",
      ref_my == ref_jang & alt_my == alt_jang                           ~ "aligned",
      ref_my == alt_jang & alt_my == ref_jang                           ~ "flipped",
      TRUE                                                              ~ "mismatched"
    ),
    beta_jang_harmonized = case_when(
      allele_status == "flipped" ~ -beta_jang,
      allele_status == "aligned" ~  beta_jang,
      TRUE                        ~ NA_real_
    )
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

# ========================================================================================
# 5. Cochran's Q, FDR (pooled across all usable rows)
# ========================================================================================

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

# ========================================================================================
# 6. Aggregate to final pair-level call: heterogeneous in ALL tested comparisons
# ========================================================================================

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

heterogeneous_snps <- heterogeneous_pairs %>%
  filter(heterogeneous_all) %>%
  distinct(phenotype_id, variant_id)

# ========================================================================================
# 7. Funnel summary
# ========================================================================================

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
    "C: unique (gene, SNP) pairs heterogeneous in ALL tested comparisons (final, stringent)"
  ),
  N = c(n_pairs_specific_B, n_unique_pairs_specific, n_joined,
        n_allele_dropped, n_usable, n_unique_pairs_usable, n_het_rows, n_het_pairs)
)
print(funnel)

# ========================================================================================
# 8. Write output
# ========================================================================================

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
    usable                  = usable_out,
    heterogeneous_pairs      = heterogeneous_pairs,
    heterogeneous_snps        = heterogeneous_snps,
    allele_status_summary       = status_tbl,
    funnel                         = funnel
  ),
  output_file
)
message("\nHeterogeneity test results written to: ", output_file)

# Untracked TSVs for readability outside R -----------------------------------------------
output_dir  <- dirname(output_file)
output_stem <- tools::file_path_sans_ext(basename(output_file))

write_tsv(usable_out,          file.path(output_dir, paste0(output_stem, "_usable.tsv")))
write_tsv(heterogeneous_pairs, file.path(output_dir, paste0(output_stem, "_heterogeneous_pairs.tsv")))
write_tsv(status_tbl,          file.path(output_dir, paste0(output_stem, "_allele_status_summary.tsv")))
write_tsv(funnel,              file.path(output_dir, paste0(output_stem, "_funnel.tsv")))
message("Tables (TSV) written to: ", output_dir)
