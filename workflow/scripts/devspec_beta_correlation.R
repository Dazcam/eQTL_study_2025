#--------------------------------------------------------------------------------------
#
#    Developmental specificity: fetal-adult eQTL beta correlation
#
#--------------------------------------------------------------------------------------

#  Fetal significant eGene-eSNP pairs (qval<0.05) joined to matched Jang et al.
#  2026 (SingleBrain) adult betas, across 8 fetal-adult cell type comparisons
#  (NPC vs both Ext and IN per reviewer request). Identifies opposite-effect
#  pairs (sign(beta_my) != sign(beta_jang)) restricted to pairs ALSO significant
#  in Jang (two-check design, see MIGRATION NOTE), and annotates each with
#  allele-coding comparison against Jang's ref/alt.
#
#  Step 5 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  MIGRATION NOTE -- Fujita -> Jang et al. 2026 (SingleBrain):
#
#  1. Significance: Jang has no per-pair Yes/No flag equivalent to Fujita's
#     significant_by_2step_FDR. Same two-check design as classify_specificity.R
#     / devspec_heterogeneity_test.R: a pair is "significant in Jang" only if it
#     IS that gene's lead SNP in the top-association file AND qval < 0.05.
#     "Present but not the lead SNP" has no qval of its own and is treated as
#     not significant (same permissiveness caveat as elsewhere in this
#     migration -- document in report).
#
#  2. Effect size: fixed_beta/fixed_sd used (not Random_Z, which has no
#     corresponding beta/se on the natural scale -- see report methods text).
#     Sign-normalized against Allele at read time (Jang's beta is relative to
#     Allele, not always alt), so beta_jang is always alt-relative downstream,
#     matching how Fujita's beta was already implicitly treated.
#
#  3. Allele-frequency checks REMOVED ENTIRELY (not substituted with 1000G):
#     Jang has no AF column, and per the same decision made for Step C, 1000G
#     is not used as a substitute reference (no guarantee Jang's allele/strand
#     conventions were harmonized against it). For full consistency with Step
#     C, the previously-separate 1000G-vs-our-own-SNP cross-check
#     (allele_status_1000G / af_concordant_1000G) has ALSO been removed here,
#     even though it never depended on Jang's AF specifically -- kept it out
#     to avoid one dataset's harmonization being AF-free while this one still
#     used an external AF reference. Every palindromic pair is now
#     unconditionally "palindromic_unresolved" and dropped from opposite_effects
#     (beta_jang_harmonized = NA). AF_TOLERANCE is no longer used by this script.
#
#  4. Preloading: each of the 5 unique Jang adult cell types (Ext, IN, MG, OPC,
#     End) is read ONCE (full + top files), then joined per comparison -- the
#     previous Fujita version re-read the same file up to 3x per shared adult
#     cell type (e.g. Exc for Glu-UL, Glu-DL, and NPC); Jang's larger files
#     make that redundancy worth removing.
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

message("\n\nComputing fetal-adult (Jang et al. 2026) eQTL beta correlations and opposite-effects allele checks ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
jang_l1_files    <- snakemake@input[["jang_l1_files"]]    # full-assoc, 5 unique adult types
jang_l1_top_files <- snakemake@input[["jang_l1_top_files"]] # top-assoc, same 5 types
pvar_file        <- snakemake@input[["pvar_file"]]
output_file      <- snakemake@output[[1]]
nom_dir          <- snakemake@params[["nom_dir"]]

geno_pc     <- 4
norm_method <- "quantile"

exp_pc_map <- c(
  "Glu-UL"    = 50, "Glu-DL"    = 40, "GABA"      = 30,
  "NPC"       = 30, "MG"        = 30, "OPC"       = 30,
  "Endo-Peri" = 30
)

# 8 fetal-adult cell type comparisons (NPC vs both Ext and IN, no Ast).
# Jang label equivalents of the old Fujita mapping: Exc->Ext, Inh->IN,
# Mic->MG (NB: fetal "MG" and Jang "MG" now share the same string -- distinct
# columns (my_cell/jang_cell) keep this unambiguous throughout).
celltype_map <- tibble(
  my_cell   = c("Glu-UL", "Glu-DL", "GABA", "NPC",  "NPC",  "MG",  "OPC", "Endo-Peri"),
  jang_cell = c("Ext",    "Ext",    "IN",   "Ext",  "IN",   "MG",  "OPC", "End")
)

unique_jang_cts <- unique(celltype_map$jang_cell)
names(jang_l1_files)     <- unique_jang_cts
names(jang_l1_top_files) <- unique_jang_cts

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_jang_l1_files", "n_jang_l1_top_files", "pvar_file",
               "n_comparisons", "n_unique_adult_types", "output_file"),
  value    = c(length(jang_l1_files), length(jang_l1_top_files), pvar_file,
               nrow(celltype_map), length(unique_jang_cts), output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\nCell type comparisons:")
print(celltype_map)

# ========================================================================================
# 1. Load our REF/ALT from the pvar file (rsID-keyed lookup)
# ========================================================================================

message("\nLoading REF/ALT lookup from pvar file ...")
pvar_dt <- fread(
  pvar_file,
  select = c("ID", "REF", "ALT"),
  col.names = c("snp", "ref_my", "alt_my"),
  showProgress = FALSE
)
message("  Variants in pvar: ", nrow(pvar_dt))
setkey(pvar_dt, snp)

# ========================================================================================
# 2. Preload each unique Jang adult cell type ONCE (full + top files)
# ========================================================================================

message("\n--- Preloading Jang adult cell types (", length(unique_jang_cts), " unique) ---\n")

jang_full_tables <- list()
jang_top_tables  <- list()

for (jct in unique_jang_cts) {
  message("\nReading Jang full-association file for: ", jct, " (", jang_l1_files[[jct]], ")")

  full_dt <- fread(
    jang_l1_files[[jct]],
    select = c("feature", "variant_id", "ref", "alt", "Allele", "fixed_beta", "fixed_sd",
               "Fixed_P", "Random_P"),
    showProgress = FALSE
  )

  # ADDITION -- placeholder-cohort artifact filter (same fix applied
  # consistently across every script that reads Jang -- beta=0/sd=0
  # placeholder cohorts give infinite meta-analysis weight, forcing
  # Fixed_P==Random_P==1 exactly). Removed here before any other use.
  n_before_artifact <- nrow(full_dt)
  full_dt <- full_dt[!(Fixed_P == 1 & Random_P == 1)]
  n_removed_artifact <- n_before_artifact - nrow(full_dt)
  message("  Removed ", n_removed_artifact,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before_artifact, " (", round(100 * n_removed_artifact / n_before_artifact, 3), "%)")
  full_dt[, c("Fixed_P", "Random_P") := NULL]

  full_dt[, gene := sub("\\..*", "", feature)]
  full_dt <- full_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene)]

  # Sign-normalize against Allele (see MIGRATION NOTE)
  full_dt[, beta_jang := ifelse(Allele == alt, fixed_beta, -fixed_beta)]
  full_dt[, se_jang := fixed_sd]
  full_dt <- full_dt[, .(snp = variant_id, gene, beta_jang, se_jang,
                          ref_jang = ref, alt_jang = alt)]
  full_dt <- unique(full_dt, by = c("snp", "gene"))
  setkey(full_dt, snp, gene)
  message("  Full-association rows retained: ", nrow(full_dt))
  jang_full_tables[[jct]] <- full_dt

  message("Reading Jang top-association file for: ", jct, " (", jang_l1_top_files[[jct]], ")")
  top_dt <- fread(
    jang_l1_top_files[[jct]],
    select = c("feature", "variant_id", "qval", "Fixed_P", "Random_P"),
    showProgress = FALSE
  )

  n_before_artifact_top <- nrow(top_dt)
  top_dt <- top_dt[!(Fixed_P == 1 & Random_P == 1)]
  n_removed_artifact_top <- n_before_artifact_top - nrow(top_dt)
  message("  Removed ", n_removed_artifact_top,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before_artifact_top)
  top_dt[, c("Fixed_P", "Random_P") := NULL]

  top_dt[, gene := sub("\\..*", "", feature)]
  top_dt <- top_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene)]
  top_dt <- top_dt[, .(gene, lead_snp = variant_id, qval)]
  setkey(top_dt, gene)
  message("  Top-association (lead SNP) rows: ", nrow(top_dt))
  jang_top_tables[[jct]] <- top_dt

  gc()
}

# ========================================================================================
# 3. Run each fetal-adult comparison
# ========================================================================================

message("\n--- Running fetal-adult beta correlations ---\n")

paired_results   <- list()
opposite_results <- list()

for (i in seq_len(nrow(celltype_map))) {

  my_ct   <- celltype_map$my_cell[i]
  jang_ct <- celltype_map$jang_cell[i]
  expPC   <- exp_pc_map[[my_ct]]

  message("\n=== ", my_ct, " vs ", jang_ct, " ===")

  # --- Load fetal significant eGene-eSNP pairs (with af) ---
  file_dir  <- paste0(nom_dir, my_ct, "_", norm_method, "_genPC_", geno_pc, "_expPC_", expPC, "/")
  file_path <- paste0(file_dir, my_ct, "_", norm_method, "_perm.cis_qtl.txt.gz")

  message("Reading fetal perm file: ", file_path)
  pooled_my <- read_delim(file_path, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(snp = variant_id, gene = phenotype_id, beta_my = slope, af_my = af) %>%
    filter(str_detect(snp, "^rs"), str_detect(gene, "^ENSG"))

  message("  Sig. fetal eQTL pairs: ", nrow(pooled_my))

  # --- Join against the preloaded Jang full-association table ---
  pooled_my_dt <- as.data.table(pooled_my)
  setkey(pooled_my_dt, snp, gene)

  paired_betas <- jang_full_tables[[jang_ct]][pooled_my_dt, nomatch = 0L]
  message("  Overlapping SNP-gene pairs: ", nrow(paired_betas))

  # --- Two-check significance: is this pair the gene's Jang lead SNP, with
  # qval < 0.05? (see MIGRATION NOTE) ---
  top_lookup <- jang_top_tables[[jang_ct]]
  paired_betas <- top_lookup[paired_betas, on = "gene"]
  paired_betas[, is_lead_snp := !is.na(lead_snp) & snp == lead_snp]
  paired_betas[, significant_jang := is_lead_snp & !is.na(qval) & qval < 0.05]

  paired_betas <- as_tibble(paired_betas) %>%
    mutate(my_cell_type = my_ct, jang_cell_type = jang_ct) %>%
    select(snp, gene, beta_my, beta_jang, af_my, se_jang, ref_jang, alt_jang,
           significant_jang, my_cell_type, jang_cell_type)

  message("  Pairs significant in Jang (lead SNP, qval<0.05): ", sum(paired_betas$significant_jang))

  paired_results[[paste(my_ct, jang_ct, sep = "_vs_")]] <- paired_betas

  # --- Flag opposite effects, restricted to pairs significant in BOTH datasets ---
  opp <- paired_betas %>%
    filter(sign(beta_my) != sign(beta_jang), significant_jang)

  message("  Opposite-effect pairs (fetal AND Jang significant): ", nrow(opp))

  if (nrow(opp) > 0) {
    opp_dt <- as.data.table(opp)
    setkey(opp_dt, snp)
    opp_dt <- pvar_dt[opp_dt, on = "snp"]

    opp_annotated <- as_tibble(opp_dt) %>%
      mutate(
        is_palindromic =
          (ref_my == "A" & alt_my == "T") | (ref_my == "T" & alt_my == "A") |
          (ref_my == "C" & alt_my == "G") | (ref_my == "G" & alt_my == "C"),
        allele_status = case_when(
          is.na(ref_my) | is.na(alt_my) | is.na(ref_jang) | is.na(alt_jang) ~ "missing_allele_info",
          is_palindromic                                                     ~ "palindromic_unresolved",
          ref_my == ref_jang & alt_my == alt_jang                            ~ "aligned",
          ref_my == alt_jang & alt_my == ref_jang                            ~ "flipped",
          TRUE                                                               ~ "mismatched"
        ),
        beta_jang_harmonized = case_when(
          allele_status == "flipped" ~ -beta_jang,
          allele_status == "aligned" ~  beta_jang,
          TRUE                        ~ NA_real_
        )
      ) %>%
      select(my_cell_type, jang_cell_type, gene, snp,
             beta_my, beta_jang, beta_jang_harmonized,
             af_my, ref_my, alt_my, ref_jang, alt_jang,
             allele_status)

    aligned_or_flipped_n <- sum(opp_annotated$allele_status %in% c("aligned", "flipped"))
    message("  Allele status: ",
            paste(names(table(opp_annotated$allele_status)), table(opp_annotated$allele_status),
                  sep = "=", collapse = ", "))
    message("  Usable (harmonized) opposite-effect pairs: ", aligned_or_flipped_n, " / ", nrow(opp_annotated))

    opposite_results[[paste(my_ct, jang_ct, sep = "_vs_")]] <- opp_annotated
  }
}

# ========================================================================================
# 4. Combine and summarise
# ========================================================================================

message("\n--- Summary across all comparisons ---\n")

paired_combined <- bind_rows(paired_results)
opposite_combined <- bind_rows(opposite_results)

message("Total paired pairs across all comparisons: ", nrow(paired_combined))
message("Total opposite-effect pairs (fetal AND Jang significant): ", nrow(opposite_combined))

message("\nCorrelation per comparison (all paired pairs, regardless of significance):")
paired_combined %>%
  group_by(my_cell_type, jang_cell_type) %>%
  summarise(
    n_pairs = n(),
    n_also_sig_jang = sum(significant_jang),
    n_opposite_both_sig = sum(sign(beta_my) != sign(beta_jang) & significant_jang),
    pct_opposite_both_sig = round(100 * n_opposite_both_sig / n_pairs, 1),
    pearson_r = round(cor(beta_my, beta_jang, use = "complete.obs"), 3),
    .groups = "drop"
  ) %>%
  print(n = Inf)

if (nrow(opposite_combined) > 0) {
  message("\nAllele status summary, opposite-effect (both-significant) pairs only:")
  print(opposite_combined %>% count(my_cell_type, jang_cell_type, allele_status))
}

# ========================================================================================
# 5. Write output
# ========================================================================================

write_rds(
  list(
    paired_betas      = paired_combined,
    opposite_effects  = opposite_combined,
    celltype_map      = celltype_map
  ),
  output_file
)
message("\nBeta correlation and opposite-effects results written to: ", output_file)
