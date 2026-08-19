#--------------------------------------------------------------------------------------
#
#    Developmental specificity: heterogeneity test on fetal-specific eQTL (Step C)
#
#--------------------------------------------------------------------------------------

#  For every (fetal cell type, gene, SNP) row in the universe (A) whose (gene, SNP)
#  pair was classified testable_specific by the "any" tier (B, i.e. significant in
#  fetal, present but not significant in any of the 7 Fujita cell types pooled),
#  looks up the JOINT (gene, SNP, matched Fujita cell type) row directly in B's
#  pairs_all table, harmonizes alleles, and tests whether the fetal and adult beta
#  differ by more than combined sampling error predicts (Cochran's Q):
#
#    Q = (beta_my - beta_fugita_harmonized)^2 / (se_my^2 + se_fugita^2) ~ chisq(df=1)
#
#  Matched Fujita cell type comes from A's fugita_label_match tag, already carried
#  per-row (L1 -> its own match; L2 -> its L1 parent's match). NPC rows (no adult
#  equivalent) are set aside, not tested, reported separately.
#
#  A pair failing the joint lookup for a given row (present in pairs_all under a
#  DIFFERENT Fujita cell type, or only under Ast/Oli, which no fetal type matches)
#  is dropped for that row -- not a marginal gene/SNP membership check, a direct
#  joint-pair lookup, since B's testable_specific already guarantees the pair was
#  jointly tested somewhere among the 7, just not necessarily in this row's
#  specific match.
#
#  Since A is L1+L2, a single (gene, SNP) pair can contribute more than one row
#  here (e.g. significant in both an L1 parent and its L2 child) -- tracked via
#  n_fetal_celltypes_tested, not collapsed.
#
#  FDR (BH) pooled across all usable rows. Classification is binary (heterogeneous
#  / not_heterogeneous) -- every row entering this test is already fetal-sig/
#  not-adult-sig by construction (that's what testable_specific means), so there's
#  no adult_sig axis left to cross against.
#
#  Standalone: does not feed D (LD clumping) or E (OCR enrichment), both of which
#  continue to read B's any_table directly, unchanged.
#
#  Step 3 (NEW) of the developmental-specificity pipeline.

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

message("\n\nRunning heterogeneity test on fetal-specific eQTL (Step C) ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file  <- snakemake@input[["eqtl_universe"]]
classification_file <- snakemake@input[["classification"]]
pvar_file           <- snakemake@input[["pvar_file"]]
output_file         <- snakemake@output[[1]]
af_tol              <- snakemake@params[["af_tol"]]
fdr_thresh          <- snakemake@params[["fdr_thresh"]]

if (is.null(af_tol))     af_tol <- 0.2
if (is.null(fdr_thresh)) fdr_thresh <- 0.05

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "classification_file", "pvar_file",
               "af_tol", "fdr_thresh", "output_file"),
  value    = c(eqtl_universe_file, classification_file, pvar_file,
               af_tol, fdr_thresh, output_file)
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
# 2. Build C's starting population: A's rows, restricted to B's testable_specific pairs
# ========================================================================================

specific_pairs <- any_table %>%
  filter(final_category == "testable_specific") %>%
  distinct(phenotype_id, variant_id)

rows_specific <- eqtl_universe %>%
  semi_join(specific_pairs, by = c("phenotype_id", "variant_id"))

n_rows_specific <- nrow(rows_specific)
n_unique_pairs_specific <- rows_specific %>% distinct(phenotype_id, variant_id) %>% nrow()

message("  C: rows entering (all fetal cell types, L1+L2): ", n_rows_specific)
message("  C: unique (gene, SNP) pairs represented: ", n_unique_pairs_specific)

if (n_unique_pairs_specific != n_pairs_specific_B) {
  message("  WARNING: unique pairs in C's starting population (", n_unique_pairs_specific,
          ") does not match B's testable_specific count (", n_pairs_specific_B,
          "). Check A/B are from the same run.")
}

# ========================================================================================
# 3. Set aside NPC (no adult equivalent)
# ========================================================================================

rows_npc      <- rows_specific %>% filter(is.na(fugita_label_match))
rows_testable <- rows_specific %>% filter(!is.na(fugita_label_match))

n_npc_excluded <- nrow(rows_npc)
message("  C: set aside -- no adult equivalent (NPC): ", n_npc_excluded)

# ========================================================================================
# 4. Joint (gene, SNP, matched Fujita cell type) lookup directly in pairs_all
# ========================================================================================

message("\n--- Joint pair lookup in matched Fujita cell type ---\n")

pairs_all_lookup <- pairs_all %>%
  select(phenotype_id = gene_id, variant_id = snps, fugita_cell_type,
         beta_fugita = beta, se_fugita = se,
         ref_fugita = REF, alt_fugita = ALT, af_fugita = ALT_AF)

joined <- rows_testable %>%
  left_join(
    pairs_all_lookup,
    by = c("phenotype_id", "variant_id", "fugita_label_match" = "fugita_cell_type")
  )

dropped   <- joined %>% filter(is.na(beta_fugita))
matched   <- joined %>% filter(!is.na(beta_fugita))

n_dropped <- nrow(dropped)
message("  C: dropped -- pair not tested in this row's matched Fujita cell type: ", n_dropped)

# Diagnostic: which Fujita cell type(s) the dropped pair WAS found in instead
if (n_dropped > 0) {
  dropped_found_elsewhere <- dropped %>%
    distinct(phenotype_id, variant_id) %>%
    left_join(
      pairs_all %>% select(phenotype_id = gene_id, variant_id = snps, fugita_cell_type),
      by = c("phenotype_id", "variant_id")
    ) %>%
    group_by(phenotype_id, variant_id) %>%
    summarise(found_in = paste(sort(unique(fugita_cell_type)), collapse = ", "),
              .groups = "drop")

  message("\n  Where dropped pairs were found instead (unique pairs):")
  dropped_found_elsewhere %>%
    count(found_in, name = "n_pairs") %>%
    arrange(desc(n_pairs)) %>%
    knitr::kable(format = "simple", align = "l") %>%
    print()
} else {
  dropped_found_elsewhere <- tibble()
}

# ========================================================================================
# 4b. Diagnostic: confirm drops are genuine coverage gaps, not a join key mismatch
# ========================================================================================

message("\n--- Diagnostic: verifying dropped rows are genuine coverage gaps ---\n")

# Check 1: every fugita_label_match value used in A should appear somewhere in
# B's pairs_all$fugita_cell_type -- a gross vocabulary mismatch would show up here.
labels_in_A <- sort(unique(na.omit(eqtl_universe$fugita_label_match)))
types_in_B  <- sort(unique(pairs_all$fugita_cell_type))

missing_labels <- setdiff(labels_in_A, types_in_B)

if (length(missing_labels) == 0) {
  message("  Check 1 (vocabulary): passed -- all fugita_label_match values in A",
          " (", paste(labels_in_A, collapse = ", "), ") are present in B's pairs_all.")
} else {
  message("  Check 1 (vocabulary): WARNING -- fugita_label_match value(s) not found",
          " anywhere in pairs_all: ", paste(missing_labels, collapse = ", "),
          ". These would never match regardless of coverage.")
}

# Check 2: for dropped rows only, redo the lookup with whitespace/case normalized
# on both sides. If this finds a match the original strict join (step 4) missed,
# that points to a join key mismatch (e.g. trailing whitespace, case difference)
# rather than a genuine absence of adult data.
normalize_label <- function(x) str_trim(toupper(x))

pairs_tested_in <- pairs_all %>%
  distinct(phenotype_id = gene_id, variant_id = snps, fugita_cell_type) %>%
  mutate(fugita_cell_type_norm = normalize_label(fugita_cell_type))

join_check <- dropped %>%
  distinct(phenotype_id, variant_id, fugita_label_match) %>%
  mutate(fugita_label_match_norm = normalize_label(fugita_label_match)) %>%
  left_join(
    pairs_tested_in,
    by = c("phenotype_id", "variant_id",
           "fugita_label_match_norm" = "fugita_cell_type_norm")
  )

suspicious <- join_check %>% filter(!is.na(fugita_cell_type))
n_suspicious <- nrow(suspicious)

if (n_suspicious == 0) {
  message("  Check 2 (join integrity): passed -- no dropped (gene, SNP) pair has a",
          " normalized-string match to its required Fujita cell type. Drops reflect",
          " genuine absence of testing in the matched adult cell type.")
} else {
  message("  Check 2 (join integrity): WARNING -- ", n_suspicious,
          " dropped (gene, SNP) pair(s) appear to have real data in their required",
          " Fujita cell type once whitespace/case is normalized -- possible join key",
          " mismatch. Inspect:")
  print(suspicious)
}

# ========================================================================================
# 5. Allele harmonization on matched rows
# ========================================================================================

message("\n--- Allele harmonization ---\n")

pvar_lookup <- fread(
  pvar_file,
  select = c("ID", "REF", "ALT"),
  col.names = c("variant_id", "ref_my", "alt_my"),
  showProgress = FALSE
)

matched <- matched %>%
  rename(beta_my = slope, se_my = slope_se, af_my = af) %>%
  left_join(as_tibble(pvar_lookup), by = "variant_id")

matched <- matched %>%
  mutate(
    is_palindromic =
      (ref_my == "A" & alt_my == "T") | (ref_my == "T" & alt_my == "A") |
      (ref_my == "C" & alt_my == "G") | (ref_my == "G" & alt_my == "C"),
    af_diff_aligned = abs(af_my - af_fugita),
    af_diff_flipped = abs(af_my - (1 - af_fugita)),
    allele_status = case_when(
      is.na(ref_my) | is.na(alt_my) | is.na(ref_fugita) | is.na(alt_fugita) ~ "missing_allele_info",
      !is_palindromic & ref_my == ref_fugita & alt_my == alt_fugita         ~ "aligned",
      !is_palindromic & ref_my == alt_fugita & alt_my == ref_fugita         ~ "flipped",
      !is_palindromic                                                      ~ "mismatched",
      is_palindromic & af_diff_aligned < af_tol & af_diff_flipped >= af_tol ~ "palindromic_af_aligned",
      is_palindromic & af_diff_flipped < af_tol & af_diff_aligned >= af_tol ~ "palindromic_af_flipped",
      TRUE                                                                 ~ "palindromic_ambiguous"
    ),
    beta_fugita_harmonized = case_when(
      allele_status %in% c("flipped", "palindromic_af_flipped") ~ -beta_fugita,
      allele_status %in% c("aligned", "palindromic_af_aligned")  ~  beta_fugita,
      TRUE                                                       ~ NA_real_
    )
  )

status_tbl <- matched %>% count(allele_status, name = "N") %>% mutate(prop = round(N / sum(N), 4))
message("Allele status breakdown:")
print(status_tbl)

n_allele_dropped <- sum(is.na(matched$beta_fugita_harmonized))
usable <- matched %>% filter(!is.na(beta_fugita_harmonized))
n_usable <- nrow(usable)
n_unique_pairs_usable <- usable %>% distinct(phenotype_id, variant_id) %>% nrow()

message("\n  C: dropped -- unresolved allele harmonization: ", n_allele_dropped)
message("  C: usable rows entering Q-test: ", n_usable)
message("  C: unique (gene, SNP) pairs among usable rows: ", n_unique_pairs_usable)

# ========================================================================================
# 6. Cochran's Q, FDR (pooled across all usable rows)
# ========================================================================================

message("\n--- Running Cochran's Q ---\n")

usable <- usable %>%
  mutate(
    delta_beta = beta_my - beta_fugita_harmonized,
    se_delta   = sqrt(se_my^2 + se_fugita^2),
    Q          = (delta_beta / se_delta)^2,
    p          = pchisq(Q, df = 1, lower.tail = FALSE),
    q          = p.adjust(p, method = "BH"),
    heterogeneous = q < fdr_thresh
  )

message("  Heterogeneous rows (FDR<", fdr_thresh, "): ", sum(usable$heterogeneous),
        " (", round(100 * mean(usable$heterogeneous), 2), "%)")

# ========================================================================================
# 7. Cross-cell-type duplication tracking
# ========================================================================================

celltypes_per_pair <- usable %>%
  distinct(phenotype_id, variant_id, cell_type) %>%
  count(phenotype_id, variant_id, name = "n_fetal_celltypes_tested")

usable <- usable %>% left_join(celltypes_per_pair, by = c("phenotype_id", "variant_id"))

message("\nDistribution of fetal cell types tested per (gene, SNP) pair, among usable rows:")
celltypes_per_pair %>%
  mutate(n_bucket = if_else(n_fetal_celltypes_tested >= 3, "3+", as.character(n_fetal_celltypes_tested))) %>%
  count(n_bucket, name = "n_pairs") %>%
  arrange(n_bucket) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

# ========================================================================================
# 8. Funnel summary
# ========================================================================================

message("\n--- Population funnel ---\n")

funnel <- tibble(
  stage = c(
    "B: unique (gene, SNP) pairs classified testable_specific",
    "C: rows entering (all fetal cell types, L1+L2)",
    "C: unique (gene, SNP) pairs represented",
    "C: set aside -- no adult equivalent (NPC)",
    "C: dropped -- pair not tested in matched Fujita cell type",
    "C: dropped -- unresolved allele harmonization",
    "C: usable rows entering Q-test",
    "C: unique (gene, SNP) pairs among usable rows"
  ),
  N = c(n_pairs_specific_B, n_rows_specific, n_unique_pairs_specific, n_npc_excluded,
        n_dropped, n_allele_dropped, n_usable, n_unique_pairs_usable)
)
print(funnel)

# ========================================================================================
# 9. Write output
# ========================================================================================

usable_out <- usable %>%
  select(cell_type, level, parent_cell_type, phenotype_id, variant_id,
         fugita_label_match, beta_my, se_my, af_my,
         beta_fugita, beta_fugita_harmonized, se_fugita, af_fugita,
         ref_my, alt_my, ref_fugita, alt_fugita, allele_status,
         delta_beta, se_delta, Q, p, q, heterogeneous, n_fetal_celltypes_tested)

write_rds(
  list(
    usable                  = usable_out,
    npc_excluded             = rows_npc,
    dropped_not_tested        = dropped,
    dropped_found_elsewhere    = dropped_found_elsewhere,
    allele_status_summary        = status_tbl,
    celltypes_per_pair             = celltypes_per_pair,
    funnel                           = funnel
  ),
  output_file
)
message("\nHeterogeneity test results written to: ", output_file)

# Untracked TSVs for readability outside R -----------------------------------------------
output_dir  <- dirname(output_file)
output_stem <- tools::file_path_sans_ext(basename(output_file))

write_tsv(usable_out,             file.path(output_dir, paste0(output_stem, "_usable.tsv")))
write_tsv(rows_npc,               file.path(output_dir, paste0(output_stem, "_npc_excluded.tsv")))
write_tsv(dropped_found_elsewhere, file.path(output_dir, paste0(output_stem, "_dropped_found_elsewhere.tsv")))
write_tsv(status_tbl,             file.path(output_dir, paste0(output_stem, "_allele_status_summary.tsv")))
write_tsv(funnel,                 file.path(output_dir, paste0(output_stem, "_funnel.tsv")))
message("Tables (TSV) written to: ", output_dir)
