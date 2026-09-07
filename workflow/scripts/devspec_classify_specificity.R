#--------------------------------------------------------------------------------------
#
#    Developmental specificity: shared vs. fetal-specific eQTL classification
#
#--------------------------------------------------------------------------------------

#  For every fetal eGene-eSNP pair in the universe extracted by
#  devspec_extract_universe.R, classifies testability against Jang et al. 2026
#  (SingleBrain) adult sn-eQTL data and, where testable, calls shared vs.
#  fetal-specific.
#
#  "Any" tier: operates on the unique (gene, SNP) pairs across all fetal cell
#  types. Tested against ALL Jang L1 + L2 cell types pooled (MiGA3 excluded --
#  a separate microglia-only meta-analysis, not part of the core L1+L2
#  hierarchy, and would double-count microglia signal alongside MG/MG1-4).
#  "Shared" = significant in >=1 Jang cell type. This is a broader reference set
#  than the previous Fujita-only (7 L1 types) design, by intent -- the point is
#  a more stringent bar for calling something fetal-specific.
#
#  Produces the fetal-specific SNP list for the OCR enrichment step (via LD
#  clumping) and the fetal-specific input population for the heterogeneity test
#  (step C).
#
#  Step 2 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  MIGRATION NOTE -- Fujita -> Jang et al. 2026 (SingleBrain):
#
#  Jang's files have no per-pair Yes/No significance flag equivalent to Fujita's
#  significant_by_2step_FDR. Significance is derived via a two-check design:
#
#    1. TESTABLE: the (gene, SNP) pair is present as a tested pair in that cell
#       type's full-association file (every tested pair, not just each gene's
#       lead SNP).
#    2. SIGNIFICANT: the pair is ALSO that gene's lead SNP in that cell type's
#       TOP-association file (one row per gene, its single most significant
#       SNP), AND that row's qval < 0.05 (qval = Random-effect P, FDR-adjusted
#       across variants within the gene, then Storey's q-value adjusted across
#       genes -- Jang's own two-stage correction, functionally the closest
#       equivalent to Fujita's significant_by_2step_FDR).
#
#  This is a MORE PERMISSIVE definition of "not significant" than Fujita's:
#  Fujita's flag was available for every tested pair independently, whereas here
#  a pair that IS tested (testable) but is NOT its gene's lead SNP has no qval
#  of its own to check at all -- it falls to testable_specific by default, not
#  because it was formally tested against the two-stage threshold and failed.
#  This will inflate testable_specific relative to the Fujita-based classification.
#  No new category is introduced for this middle case (decided deliberately, not
#  an oversight) -- the effect should instead be documented in the report's
#  methods text.
#
#  fixed_beta is used as the effect-size column (not Random_Z, which has no
#  corresponding beta/se on the natural scale -- see report methods text for the
#  fixed-vs-random rationale). Sign-normalized against the Allele column at read
#  time: Jang's beta is relative to the "Allele" column, which is NOT always ALT
#  (unlike Fujita, where beta was implicitly always ALT-relative) -- so beta is
#  flipped wherever Allele != alt, making pairs_all's beta always ALT-relative
#  downstream, exactly matching how Fujita's beta was already treated.
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

message("\n\nClassifying fetal eGene-eSNP pairs: testability + shared/specific vs. Jang et al. 2026 ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
jang_full_files    <- snakemake@input[["jang_full_files"]]
jang_top_files     <- snakemake@input[["jang_top_files"]]
output_file        <- snakemake@output[[1]]
jang_cell_types    <- snakemake@params[["jang_cell_types"]]

names(jang_full_files) <- jang_cell_types
names(jang_top_files)  <- jang_cell_types

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "n_jang_full_files", "n_jang_top_files", "output_file"),
  value    = c(eqtl_universe_file, length(jang_full_files), length(jang_top_files), output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# Full-association file: every tested pair for this cell type, filtered to our
# universe genes/SNPs. beta sign-normalized against Allele (see MIGRATION NOTE).
#
# ADDITION -- placeholder-cohort artifact filter (same fix as
# replication_pi1_enrichment.R's Jang branch): some cohorts encode "not
# tested for this pair" as beta=0/sd=0 rather than NA, giving that cohort
# infinite weight (1/SE^2) in the meta-analysis and forcing Fixed_P ==
# Random_P == 1 exactly. Removed entirely (not just beta/se nulled out) --
# per decision, every script reading Jang applies the identical filter, so
# there's one consistent set of rows omitted everywhere rather than a
# per-script judgement call about which use of the row is "affected".
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

  n_before_artifact <- nrow(dat)
  dat <- dat %>% filter(!(Fixed_P == 1 & Random_P == 1))
  n_removed_artifact <- n_before_artifact - nrow(dat)
  message("  Removed ", n_removed_artifact,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before_artifact, " (", round(100 * n_removed_artifact / n_before_artifact, 3), "%)")

  dat <- dat %>% select(-Fixed_P, -Random_P)

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
      beta = if_else(Allele == alt, fixed_beta, -fixed_beta),
      se   = fixed_sd,
      jang_cell_type = jang_ct
    ) %>%
    select(gene_id, variant_id, beta, se, ref, alt, jang_cell_type)

  message("  Genes tested: ", length(gene_set),
          " | SNPs tested: ", length(snp_set),
          " | Relevant pairs retained: ", nrow(pair_lookup))

  list(gene_set = gene_set, snp_set = snp_set, pair_lookup = pair_lookup)
}

# Top-association file: one row per gene (its lead SNP), with qval. This is the
# only place a two-stage-corrected significance value exists in Jang's data.
# Same artifact filter applied here too -- a gene where EVERY tested SNP is an
# artifact row would otherwise still surface as that gene's "lead SNP" with
# qval~1, which should correctly fall out as untested rather than be reported
# as a (trivially failing) significance check.
read_jang_top <- function(file_path, jang_ct) {
  message("Reading Jang top-association data for: ", jang_ct, " (", file_path, ")")

  dat <- read_tsv(
    file_path,
    col_select = c(feature, variant_id, qval, Fixed_P, Random_P),
    col_types = cols(feature = "c", variant_id = "c", qval = "d",
                      Fixed_P = "d", Random_P = "d"),
    show_col_types = FALSE
  )

  n_before_artifact <- nrow(dat)
  dat <- dat %>% filter(!(Fixed_P == 1 & Random_P == 1))
  n_removed_artifact <- n_before_artifact - nrow(dat)
  message("  Removed ", n_removed_artifact,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before_artifact)

  dat %>%
    mutate(gene_id = sub("\\..*", "", feature)) %>%
    filter(str_detect(variant_id, "^rs"), str_detect(gene_id, "^ENSG")) %>%
    transmute(gene_id, lead_variant_id = variant_id, qval, jang_cell_type = jang_ct)
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

# Loop over Jang full-association files -------------------------------------------------
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

# Per-gene/SNP membership lookups, same role as the Fujita version -- used for
# the "any" tier's testability check, and written to output for Step C.
gene_membership <- map_dfr(names(gene_sets), ~ tibble(jang_cell_type = .x, gene_id = gene_sets[[.x]]))
snp_membership  <- map_dfr(names(snp_sets),  ~ tibble(jang_cell_type = .x, variant_id = snp_sets[[.x]]))

# Loop over Jang top-association files ---------------------------------------------------
message("\n--- Reading Jang top-association files (", length(jang_cell_types), " cell types) ---\n")

top_lookups <- map(jang_cell_types, ~ read_jang_top(jang_top_files[[.x]], .x))
top_lookup_all <- bind_rows(top_lookups)
message("Total lead-SNP rows across all cell types: ", nrow(top_lookup_all))

# =======================================================================================
# "ANY" TIER: unique pairs, checked against ALL Jang cell types pooled
# =======================================================================================
message("\n--- Building 'any' tier classification ---")

# Is this pairs_all row (gene, snp, tested in cell type X) also cell type X's
# lead SNP for that gene, with qval < 0.05? (see MIGRATION NOTE for why this is
# a two-check design rather than a single per-pair flag)
pairs_all <- pairs_all %>%
  left_join(
    top_lookup_all,
    by = c("gene_id", "jang_cell_type")
  ) %>%
  mutate(
    is_lead_snp = !is.na(lead_variant_id) & variant_id == lead_variant_id,
    significant_lead = is_lead_snp & !is.na(qval) & qval < 0.05
  )

message("Rows in pairs_all that are their gene's lead SNP in that cell type: ",
        sum(pairs_all$is_lead_snp))
message("Rows additionally significant (qval<0.05): ", sum(pairs_all$significant_lead))

# Significance call: significant_lead in >= 1 Jang cell type
pairs_any_sig <- pairs_all %>%
  group_by(gene_id, variant_id) %>%
  summarise(
    significant_any    = any(significant_lead),
    n_jang_types_tested = n(),
    .groups = "drop"
  )

# Provenance: which Jang cell type gave the smallest qval among this pair's
# lead-SNP rows (for QC/reporting only). NA where the pair was never a lead SNP
# anywhere.
pairs_any_provenance <- pairs_all %>%
  filter(is_lead_snp) %>%
  group_by(gene_id, variant_id) %>%
  slice_min(qval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene_id, variant_id, best_jang_ct = jang_cell_type,
         qval_adult = qval, beta_adult = beta)

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

message("\nReminder: testable_specific here includes pairs that were tested but were\n",
        "never any cell type's lead SNP (no qval of their own to check) alongside pairs\n",
        "that WERE a lead SNP somewhere but failed qval<0.05 -- these two cases are not\n",
        "distinguished as separate categories (deliberate choice, see MIGRATION NOTE).\n",
        "This makes testable_specific more permissive than the prior Fujita-based\n",
        "classification; document this in the report.")

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
