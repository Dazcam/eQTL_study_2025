#--------------------------------------------------------------------------------------
#
#    Fetal vs adult (Jang et al. 2026) eQTL effect-size concordance
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  06QTL-REPLICATION | Rule: beta_correlation (single job)
#            Upstream:   05TENSORQTL tensorqtl_perm (final expPC per cell type);
#                        Jang et al. 2026 SingleBrain L1 full + top association files
#            Downstream: 06 report (effect-size concordance tab; Fig 3F)
#
# Purpose:   - Join each fetal cell type's sig. eGene-eSNP pairs (qval < 0.05) to the
#              matched Jang adult betas, across 8 fetal-adult comparisons
#              (NPC tested against both Ext and IN)
#            - Flag a pair as sig. in Jang only if it is the gene's lead SNP in the
#              top-association file AND qval < 0.05 (two-check design)
#            - Report opposite-effect pairs (sign(beta_my) != sign(beta_jang)) among
#              pairs sig. in both datasets, annotated with an allele-coding check
#              against Jang's ref/alt
#
# Inputs:    jang_l1_files       Jang full-association files, 5 unique adult types
#            jang_l1_top_files   Jang top-association files, same 5 types
#            pvar_file           our genotype pvar (REF/ALT per rsID)
#            nom_dir             tensorQTL output dir (params)
#
# Outputs:   output_file   rds list: paired_betas, opposite_effects, celltype_map
#
# Notes:     - Jang placeholder-cohort rows (Fixed_P == Random_P == 1) are removed
#              before use (filter_jang_artifacts()).
#            - Jang beta is relative to `Allele`, not always alt; it is
#              sign-normalised to alt at read time (align_jang_beta()).
#            - Effect sizes use fixed_beta/fixed_sd; Random_Z has no matching
#              beta/se on the natural scale.
#            - No allele-frequency checks: Jang has no AF column and 1000G is not
#              used as a substitute, so palindromic pairs are unresolved and get
#              beta_jang_harmonized = NA.
#            - Each Jang adult cell type is read once and reused across comparisons.
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

message("\n\nComputing fetal-adult (Jang et al. 2026) eQTL beta correlations and opposite-effects allele checks ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
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

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("n_jang_l1_files", "n_jang_l1_top_files", "pvar_file", "nom_dir",
               "n_comparisons", "n_unique_adult_types", "output_file"),
  value    = c(length(jang_l1_files), length(jang_l1_top_files), pvar_file, nom_dir,
               nrow(celltype_map), length(unique_jang_cts), output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

message("\nCell type comparisons:")
print(celltype_map)

##  Load our REF/ALT from the pvar file (rsID-keyed lookup)  --------------------------
message("\nLoading REF/ALT lookup from pvar file ...")
pvar_dt <- fread(
  pvar_file,
  select = c("ID", "REF", "ALT"),
  col.names = c("snp", "ref_my", "alt_my"),
  showProgress = FALSE
)
message("  Variants in pvar: ", nrow(pvar_dt))
setkey(pvar_dt, snp)

##  Preload each Jang adult cell type once (full + top files)  ------------------------
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
  full_dt <- filter_jang_artifacts(full_dt, "Jang full-association")
  full_dt[, c("Fixed_P", "Random_P") := NULL]

  full_dt[, gene := sub("\\..*", "", feature)]
  full_dt <- full_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene)]

  # Sign-normalize against Allele (see Notes in header)
  full_dt[, beta_jang := align_jang_beta(fixed_beta, Allele, alt)]
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

  top_dt <- filter_jang_artifacts(top_dt, "Jang top-association")
  top_dt[, c("Fixed_P", "Random_P") := NULL]

  top_dt[, gene := sub("\\..*", "", feature)]
  top_dt <- top_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene)]
  top_dt <- top_dt[, .(gene, lead_snp = variant_id, qval)]
  setkey(top_dt, gene)
  message("  Top-association (lead SNP) rows: ", nrow(top_dt))
  jang_top_tables[[jct]] <- top_dt

  gc()
}

##  Run each fetal-adult comparison  --------------------------------------------------
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
  # qval < 0.05? (see Notes in header) ---
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
        allele_status        = get_allele_status(ref_my, alt_my, ref_jang, alt_jang),
        beta_jang_harmonized = harmonise_beta(beta_jang, allele_status)
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

##  Combine and summarise  ------------------------------------------------------------
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

##  Write output  ---------------------------------------------------------------------
write_rds(
  list(
    paired_betas      = paired_combined,
    opposite_effects  = opposite_combined,
    celltype_map      = celltype_map
  ),
  output_file
)
message("\nBeta correlation and opposite-effects results written to: ", output_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
