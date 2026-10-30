#--------------------------------------------------------------------------------------
#
#    Developmental specificity: GLM comparison of two eQTL categories' relative
#    concentration in a peak set (Li et al. 2023 adult cortical union OCRs).
#
#--------------------------------------------------------------------------------------

#  Unlike the bootRanges-based enrichment tests (ocr_enrich / chip_enrich), this
#  test has NO genomic background model -- it is a direct comparison of two
#  proportions (fraction of fetal_specific SNPs overlapping peaks vs. fraction
#  of shared SNPs overlapping peaks), via logistic regression rather than a
#  plain Fisher's/chi-square test, to allow future covariate adjustment.
#
#  This answers a DIFFERENT question to the permutation tests: not "is this
#  set enriched over a genomically realistic null" but "do these two sets
#  differ from each other in overlap rate". Treat as complementary to, not a
#  replacement for, the background-aware enrichment tests.
#
#  Current scope: fetal_specific vs. shared only, no covariates. No LD-score
#  or other covariate included -- SNPs are already LD-pruned (independent_snps),
#  which addresses row-independence for the regression; it does not address
#  potential confounding of overlap rate itself (e.g. distance-to-TSS, local
#  peak density), which is a separate question left open for now.
#
#  Peak set: Li et al. 2023 (cortical union), reusing the same Table S6/S7
#  loading + filtering logic as the retired ocr_enrich script.

#--------------------------------------------------------------------------------------

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

message("\n\nGLM comparison: fetal_specific vs. shared eQTL concentration in Li et al. 2023 adult cortical union OCRs ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

# Input and output paths ---------------------------------------------------------------
classification_file    <- snakemake@input[["classification"]]
independent_snps_file  <- snakemake@input[["independent_snps"]]
snp_lookup_files        <- snakemake@input[["snp_lookups"]]
li_s6_file             <- snakemake@input[["li_s6"]]
li_s7_file             <- snakemake@input[["li_s7"]]
output_file            <- snakemake@output[[1]]

message("\nVariables")
cat("============================")
tibble(
  variable = c("classification_file", "independent_snps_file", "n_snp_lookups",
               "li_s6_file", "li_s7_file", "output_file"),
  value    = c(classification_file, independent_snps_file, length(snp_lookup_files),
               li_s6_file, li_s7_file, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

all_cortical_types <- c(
  "ITL23_1", "ITL23_2", "ITL23_3", "ITL23_4", "ITL23_5", "ITL23_6",
  "ITL34", "ITL4_1", "ITL4_2", "ITL45_1", "ITL45_2",
  "ITV1C_1", "ITV1C_2", "ITV1C_3",
  "ITL5_1", "ITL5_2", "ITL5_3", "ITL5_4",
  "ITL6_1_1", "ITL6_1_2", "ITL6_2_1", "ITL6_2_2",
  "CT_1", "CT_2", "L6B_1", "L6B_2", "NP_1", "NP_2", "NP_3", "ET",
  "SST_1", "SST_2", "SST_3", "SST_4", "SST_5", "SST_CHODL",
  "PVALB_1", "PVALB_2", "PVALB_3", "PVALB_4", "PV_ChCs",
  "VIP_1", "VIP_2", "VIP_3", "VIP_4", "VIP_5", "VIP_6", "VIP_7",
  "SNCG_1", "SNCG_2", "SNCG_3", "SNCG_4", "SNCG_5",
  "LAMP5", "LAMP5_LHX6",
  "MGC_1", "MGC_2",
  "OPC", "COP",
  "EC", "SMC",
  "ASCT_1", "ASCT_2", "ASCT_3",
  "OGC_1", "OGC_2", "OGC_3"
)

# Main -----------------------------------------------------------------------------------

# 1. Load fetal-specific + shared SNPs from classification, restricted to the
#    LD-independent universe. ------------------------------------------------------------
message("\nLoading classification output (B) ...")
classification <- read_rds(classification_file)
any_table <- classification$any_table

message("\nRestricting to LD-independent eGene-eSNP pairs ...")
independent_universe <- read_rds(independent_snps_file)

n_before_prune <- nrow(any_table)
any_table <- any_table %>%
  semi_join(independent_universe, by = c("phenotype_id", "variant_id"))
message("  any_table rows before LD pruning: ", n_before_prune)
message("  any_table rows after LD pruning: ", nrow(any_table))

fetal_specific_snps <- any_table %>%
  filter(final_category == "testable_specific") %>%
  pull(variant_id) %>%
  unique()
message("  Fetal-specific unique SNPs: ", length(fetal_specific_snps))

shared_snps <- any_table %>%
  filter(final_category == "testable_shared") %>%
  pull(variant_id) %>%
  unique()
message("  Shared unique SNPs: ", length(shared_snps))

# 2. Build SNP coordinate lookup from existing per-cell-type files ----------------------
message("\nBuilding SNP coordinate lookup from ", length(snp_lookup_files), " existing files ...")

snp_coords <- map_dfr(snp_lookup_files, read_rds) %>%
  distinct(refsnp_id, .keep_all = TRUE)
message("  Total unique SNPs with coordinates: ", nrow(snp_coords))

make_snp_gr <- function(snp_ids, coords) {
  coords %>%
    filter(refsnp_id %in% snp_ids) %>%
    distinct(refsnp_id, .keep_all = TRUE) %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chr_name",
      start.field    = "chrom_start",
      end.field      = "chrom_end",
      keep.extra.columns = TRUE
    )
}

specific_gr <- make_snp_gr(fetal_specific_snps, snp_coords)
shared_gr   <- make_snp_gr(shared_snps, snp_coords)

message("  Fetal-specific SNPs with coordinates: ", length(specific_gr))
message("  Shared SNPs with coordinates: ", length(shared_gr))

# 3. Load Li et al. 2023 adult peaks, cortical union only -------------------------------
message("\nLoading Li et al. 2023 adult brain cCREs (cortical union only) ...")

options(scipen = 999)

li_s6 <- read_tsv(li_s6_file, col_names = c("chr", "start", "end", "ccre_id"),
                   col_types = "ciic", show_col_types = FALSE)
li_s7 <- read_tsv(li_s7_file, col_names = c("ccre_id", "cell_types"),
                   col_types = "cc", show_col_types = FALSE)

ccre_ids_union <- li_s7 %>%
  filter(map_lgl(cell_types, ~ any(str_split_1(.x, ",") %in% all_cortical_types))) %>%
  pull(ccre_id)

adult_cortical_union_gr <- li_s6 %>%
  filter(ccre_id %in% ccre_ids_union) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )

message("  adult_cortical_union: ", length(adult_cortical_union_gr), " peaks")

# 4. Build the per-SNP GLM input table ---------------------------------------------------
#    Full SNP sets (not restricted to overlappers) -- in_category is the outcome,
#    is_set_A the predictor (1 = fetal_specific, 0 = shared). No covariates for now.
message("\nBuilding per-SNP overlap table ...")

glm_data <- bind_rows(
  tibble(
    variant_id = specific_gr$refsnp_id,
    is_set_A   = 1L,
    in_category = as.integer(overlapsAny(specific_gr, adult_cortical_union_gr))
  ),
  tibble(
    variant_id = shared_gr$refsnp_id,
    is_set_A   = 0L,
    in_category = as.integer(overlapsAny(shared_gr, adult_cortical_union_gr))
  )
)

message("  Total SNPs in GLM input: ", nrow(glm_data))
message("  Set A (fetal_specific): ", sum(glm_data$is_set_A == 1),
        " (", sum(glm_data$is_set_A == 1 & glm_data$in_category == 1), " overlap peaks)")
message("  Set B (shared): ", sum(glm_data$is_set_A == 0),
        " (", sum(glm_data$is_set_A == 0 & glm_data$in_category == 1), " overlap peaks)")

# 5. Fit the model -------------------------------------------------------------------------
message("\n--- Fitting logistic regression (no covariates) ---\n")

model <- glm(in_category ~ is_set_A, data = glm_data, family = binomial)
print(summary(model))

# Manual extraction of OR + 95% CI (avoids depending on broom, which errored
# here -- likely not installed in this container).
coefs <- coef(summary(model))
est   <- coefs["is_set_A", "Estimate"]
se    <- coefs["is_set_A", "Std. Error"]
pval  <- coefs["is_set_A", "Pr(>|z|)"]

results <- tibble(
  comparison      = "fetal_specific_vs_shared",
  peak_set        = "li_adult_cortical_union",
  n_set_A         = sum(glm_data$is_set_A == 1),
  n_set_B         = sum(glm_data$is_set_A == 0),
  n_overlap_set_A = sum(glm_data$is_set_A == 1 & glm_data$in_category == 1),
  n_overlap_set_B = sum(glm_data$is_set_A == 0 & glm_data$in_category == 1),
  odds_ratio      = exp(est),
  ci_low          = exp(est - 1.96 * se),
  ci_high         = exp(est + 1.96 * se),
  p.value         = pval
)

message("\n--- Results summary ---\n")
print(results)

# Write tracked output ----------------------------------------------------------------
write_rds(list(model = model, results = results, glm_data = glm_data), output_file)
message("\nGLM enrichment comparison written to: ", output_file)

# Untracked TSV for readability outside R -----------------------------------------------
output_tsv <- sub("\\.rds$", ".tsv", output_file)
write_tsv(results, output_tsv)
message("GLM enrichment comparison (TSV) written to: ", output_tsv)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
