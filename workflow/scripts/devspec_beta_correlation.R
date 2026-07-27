#--------------------------------------------------------------------------------------
#
#    Developmental specificity: fetal-adult eQTL beta correlation
#
#--------------------------------------------------------------------------------------

#  Extends replication_beta_correlation_single.R's logic (significant fetal eGene-eSNP
#  pairs joined to matched Fujita adult betas) to 8 fetal-adult cell type pairs,
#  including NPC vs both Exc and Inh per reviewer request. Identifies opposite-effect
#  pairs (sign(beta_my) != sign(beta_fugita)) and annotates each with allele-coding
#  comparisons against Fujita's REF/ALT and 1000G EUR's A1/A2, following the
#  aligned/flipped/mismatched logic in smr_format_gwas.R. No correction is applied --
#  this is purely descriptive reporting of whether an apparent direction discordance
#  is consistent with an allele-coding artifact or not.
#
#  Step 5 of the developmental-specificity pipeline.

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

message("\n\nComputing fetal-adult eQTL beta correlations and opposite-effects allele checks ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
fugita_files    <- snakemake@input[["fugita_files"]]
thousandG_file  <- snakemake@input[["thousandG_frq"]]
pvar_file       <- snakemake@input[["pvar_file"]]
output_file     <- snakemake@output[[1]]
nom_dir         <- snakemake@params[["nom_dir"]]
fugita_dir      <- snakemake@params[["fugita_dir"]]

geno_pc     <- 4
norm_method <- "quantile"

exp_pc_map <- c(
  "Glu-UL"    = 50, "Glu-DL"    = 40, "GABA"      = 30,
  "NPC"       = 30, "MG"        = 30, "OPC"       = 30,
  "Endo-Peri" = 30
)

# 8 fetal-adult cell type comparisons (NPC vs both Exc and Inh, no Ast)
celltype_map <- tibble(
  my_cell     = c("Glu-UL", "Glu-DL", "GABA", "NPC",  "NPC",  "MG",  "OPC", "Endo-Peri"),
  fugita_cell = c("Exc",    "Exc",    "Inh",  "Exc",  "Inh",  "Mic", "OPC", "End")
)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_fugita_files", "thousandG_file", "pvar_file", "n_comparisons", "output_file"),
  value    = c(length(fugita_files), thousandG_file, pvar_file, nrow(celltype_map), output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\nCell type comparisons:")
print(celltype_map)

# ========================================================================================
# 1. Load 1000G EUR allele frequencies
# ========================================================================================

message("\nLoading 1000G EUR hg38 allele frequencies ...")
thousandG <- fread(thousandG_file, showProgress = FALSE) %>%
  rename(snp = SNP, a1_1000G = A1, a2_1000G = A2, maf_1000G = MAF) %>%
  select(snp, a1_1000G, a2_1000G, maf_1000G) %>%
  distinct(snp, .keep_all = TRUE)
message("  1000G variants loaded: ", nrow(thousandG))

# ========================================================================================
# 2. Load our REF/ALT from the pvar file (rsID-keyed lookup)
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
# 3. Run each fetal-adult comparison
# ========================================================================================

message("\n--- Running fetal-adult beta correlations ---\n")

paired_results    <- list()
opposite_results  <- list()

for (i in seq_len(nrow(celltype_map))) {

  my_ct <- celltype_map$my_cell[i]
  fu_ct <- celltype_map$fugita_cell[i]
  expPC <- exp_pc_map[[my_ct]]

  message("\n=== ", my_ct, " vs ", fu_ct, " ===")

  # --- Load fetal significant eGene-eSNP pairs (with af) ---
  log_dir   <- paste0(nom_dir, my_ct, "_", norm_method, "_genPC_", geno_pc, "_expPC_", expPC, "/")
  file_path <- paste0(log_dir, my_ct, "_", norm_method, "_perm.cis_qtl.txt.gz")

  message("Reading fetal perm file: ", file_path)
  pooled_my <- read_delim(file_path, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(snp = variant_id, gene = phenotype_id, beta_my = slope, af_my = af) %>%
    filter(str_detect(snp, "^rs"), str_detect(gene, "^ENSG")) %>%
    mutate(key = paste(snp, gene, sep = "_"))

  message("  Sig. fetal eQTL pairs: ", nrow(pooled_my))

  # --- Load matched Fujita file, dedup as in original script ---
  fugita_path <- paste0(fugita_dir, "celltype-eqtl-sumstats.", fu_ct, ".tsv")
  message("Reading Fujita file: ", fugita_path)

  fugita_df <- read_delim(
    fugita_path, delim = "\t", col_names = TRUE,
    na = c("", "NA", "na", "-"), show_col_types = FALSE
  ) %>%
    select(snp = snps, gene = gene_id, beta_fugita = beta,
           ref_fugita = REF, alt_fugita = ALT, af_fugita = ALT_AF) %>%
    mutate(beta_fugita = as.numeric(beta_fugita)) %>%
    filter(str_detect(snp, "^rs"), str_detect(gene, "^ENSG"), !is.na(beta_fugita))

  message("  Total SNP-gene pairs in Fujita ", fu_ct, ": ", nrow(fugita_df))

  setDT(fugita_df)

  # Dedup within this Fujita file: largest |beta| per (snp, gene)
  pooled_fugita <- fugita_df[
    , .SD[which.max(abs(beta_fugita))],
    by = .(snp, gene)
  ][, .(key = paste(snp, gene, sep = "_"), snp, gene,
        beta_fugita, ref_fugita, alt_fugita, af_fugita)] %>%
    as_tibble()

  rm(fugita_df)
  gc()

  message("  Fujita pairs after dedup: ", nrow(pooled_fugita))

  # --- Join fetal and Fujita on key ---
  paired_betas <- pooled_my %>%
    inner_join(
      pooled_fugita %>% select(key, beta_fugita, ref_fugita, alt_fugita, af_fugita),
      by = "key"
    ) %>%
    select(snp, gene, beta_my, beta_fugita, af_my,
           ref_fugita, alt_fugita, af_fugita) %>%
    mutate(my_cell_type = my_ct, fugita_cell_type = fu_ct)

  message("  Overlapping SNP-gene pairs: ", nrow(paired_betas))

  paired_results[[paste(my_ct, fu_ct, sep = "_vs_")]] <- paired_betas

  # --- Flag opposite effects ---
  opp <- paired_betas %>%
    filter(sign(beta_my) != sign(beta_fugita))

  message("  Opposite-effect pairs: ", nrow(opp))

  if (nrow(opp) > 0) {
    # Join in our REF/ALT from pvar
    opp_dt <- as.data.table(opp)
    setkey(opp_dt, snp)
    opp_dt <- pvar_dt[opp_dt, on = "snp"]

    opp_annotated <- as_tibble(opp_dt) %>%
      left_join(thousandG, by = "snp") %>%
      mutate(
        allele_status_fugita = case_when(
          ref_my == ref_fugita & alt_my == alt_fugita ~ "aligned",
          ref_my == alt_fugita & alt_my == ref_fugita  ~ "flipped",
          TRUE                                          ~ "mismatched"
        ),
        allele_status_1000G = case_when(
          is.na(a1_1000G) | is.na(a2_1000G)            ~ NA_character_,
          ref_my == a1_1000G & alt_my == a2_1000G       ~ "aligned",
          ref_my == a2_1000G & alt_my == a1_1000G       ~ "flipped",
          TRUE                                          ~ "mismatched"
        )
      ) %>%
      select(my_cell_type, fugita_cell_type, gene, snp,
             beta_my, beta_fugita,
             af_my, ref_my, alt_my,
             ref_fugita, alt_fugita, af_fugita,
             a1_1000G, a2_1000G, maf_1000G,
             allele_status_fugita, allele_status_1000G)

    opposite_results[[paste(my_ct, fu_ct, sep = "_vs_")]] <- opp_annotated

    message("  Allele status (vs Fujita): ",
            paste(names(table(opp_annotated$allele_status_fugita)),
                  table(opp_annotated$allele_status_fugita), sep = "=", collapse = ", "))
    message("  Allele status (vs 1000G): ",
            paste(names(table(opp_annotated$allele_status_1000G, useNA = "ifany")),
                  table(opp_annotated$allele_status_1000G, useNA = "ifany"), sep = "=", collapse = ", "))
  }
}

# ========================================================================================
# 4. Combine and summarise
# ========================================================================================

message("\n--- Summary across all comparisons ---\n")

paired_combined <- bind_rows(paired_results)
opposite_combined <- bind_rows(opposite_results)

message("Total paired pairs across all comparisons: ", nrow(paired_combined))
message("Total opposite-effect pairs across all comparisons: ", nrow(opposite_combined))

message("\nCorrelation per comparison:")
paired_combined %>%
  group_by(my_cell_type, fugita_cell_type) %>%
  summarise(
    n_pairs = n(),
    n_opposite = sum(sign(beta_my) != sign(beta_fugita)),
    pct_opposite = round(100 * n_opposite / n_pairs, 1),
    pearson_r = round(cor(beta_my, beta_fugita, use = "complete.obs"), 3),
    .groups = "drop"
  ) %>%
  print(n = Inf)

if (nrow(opposite_combined) > 0) {
  message("\nAllele status summary (vs Fujita), opposite-effect pairs only:")
  print(opposite_combined %>% count(my_cell_type, fugita_cell_type, allele_status_fugita))

  message("\nAllele status summary (vs 1000G), opposite-effect pairs only:")
  print(opposite_combined %>% count(my_cell_type, fugita_cell_type, allele_status_1000G))
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
