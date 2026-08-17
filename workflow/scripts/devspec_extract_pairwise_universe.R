#--------------------------------------------------------------------------------------
#
#    Developmental specificity: pairwise fetal-adult tested-pair background
#
#--------------------------------------------------------------------------------------

#  Builds the background of all tested (significant AND non-significant) SNP-gene
#  pairs for 9 pairwise fetal (L1) - Fujita adult cell type comparisons, by
#  inner-joining full nominal cis-eQTL summary statistics on each side.
#
#  Also attaches, per pair:
#    - fetal_sig / adult_sig  (each dataset's own significance call)
#    - allele_status / beta_fugita_harmonized  (REF/ALT harmonization; palindromic
#      SNPs resolved via allele frequency where string matching is uninformative)
#
#  Nothing is dropped here -- unusable pairs (allele_status = mismatched /
#  palindromic_ambiguous / missing_allele_info) are retained with
#  beta_fugita_harmonized = NA, so the full allele-status breakdown can be reviewed
#  before the heterogeneity test (step 2) filters them out.
#
#  Comparisons (9): Glu-UL-Exc, Glu-DL-Exc, GABA-Inh, MG-Mic, OPC-OPC,
#  Endo-Peri-End, NPC-Exc, NPC-Inh, NPC-OPC. L1 only for this first pass.
#
#  Memory strategy: the 5 unique Fujita adult files (+ REF/ALT/AF/sig columns) and
#  the fetal permutation lead-SNP files (small) are preloaded once. The 7 unique
#  fetal nominal files are large (~2-3GB each) and read one at a time, joined
#  against whichever adult table(s) apply, then discarded before the next fetal
#  cell type is read. The fetal pvar REF/ALT lookup is also small and read once.
#
#  Step 1 (revised) of the developmental-specificity pipeline. Standalone --
#  writes to its own output path, does not touch extract_universe.R or
#  classify_specificity.R's existing outputs.

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

message("\n\nBuilding pairwise fetal-adult tested-pair background (all tested pairs) ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
nom_files    <- snakemake@input[["nom_files"]]     # large, per fetal L1 cell type
perm_files   <- snakemake@input[["perm_files"]]    # small, per fetal L1 cell type
fugita_files <- snakemake@input[["fugita_files"]]  # per Fujita cell type
pvar_file    <- snakemake@input[["pvar_file"]]     # single file, fetal REF/ALT lookup
output_file  <- snakemake@output[[1]]

fetal_cell_types  <- snakemake@params[["fetal_cell_types"]]
fugita_cell_types <- snakemake@params[["fugita_cell_types"]]
af_tol             <- snakemake@params[["af_tol"]]

if (is.null(af_tol)) af_tol <- 0.2

names(nom_files)    <- fetal_cell_types
names(perm_files)   <- fetal_cell_types
names(fugita_files) <- fugita_cell_types

# Fixed biological pairwise comparison map -- not a tunable pipeline parameter,
# hence hardcoded here rather than in config, consistent with fugita_label_match /
# celltype_map elsewhere in this pipeline. NPC compared against all adult types
# except those from a different embryonic lineage (MG = yolk-sac mesoderm-derived
# microglia; Endo-Peri = mesoderm-derived vasculature).
comparison_map <- tibble(
  my_cell     = c("Glu-UL", "Glu-DL", "GABA", "MG",  "OPC", "Endo-Peri",
                  "NPC",  "NPC",  "NPC"),
  fugita_cell = c("Exc",    "Exc",    "Inh",  "Mic", "OPC", "End",
                  "Exc",  "Inh",  "OPC")
)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_fetal_files", "n_fugita_files", "n_comparisons", "af_tol", "output_file"),
  value    = c(length(nom_files), length(fugita_files), nrow(comparison_map),
               af_tol, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\nPairwise comparisons:")
print(comparison_map)

# ========================================================================================
# 1. Preload fetal REF/ALT lookup (single small file, read once)
# ========================================================================================

message("\n--- Preloading fetal REF/ALT lookup from pvar file ---\n")

pvar_lookup <- fread(
  pvar_file,
  select = c("ID", "REF", "ALT"),
  col.names = c("variant_id", "ref_my", "alt_my"),
  showProgress = FALSE
)
setkey(pvar_lookup, variant_id)
message("  Variants in pvar lookup: ", nrow(pvar_lookup))

# ========================================================================================
# 2. Preload fetal permutation lead-SNP calls (for fetal_sig flag)
# ========================================================================================

message("\n--- Preloading fetal permutation lead-SNP calls ---\n")

perm_tables <- list()
for (fct in fetal_cell_types) {
  message("Reading perm file for: ", fct, " (", perm_files[[fct]], ")")
  perm_dt <- fread(perm_files[[fct]], select = c("phenotype_id", "variant_id", "qval"),
                    showProgress = FALSE)
  perm_dt[, fetal_sig := !is.na(qval) & qval < 0.05]
  perm_tables[[fct]] <- perm_dt[, .(phenotype_id, variant_id, fetal_sig)]
  setkey(perm_tables[[fct]], phenotype_id, variant_id)
  message("  Lead-SNP rows: ", nrow(perm_tables[[fct]]),
          " | fetal_sig TRUE: ", sum(perm_tables[[fct]]$fetal_sig))
  rm(perm_dt); gc()
}

# ========================================================================================
# 3. Preload all unique Fujita adult nominal files (REF/ALT/AF/significance included)
# ========================================================================================

message("\n--- Preloading Fujita adult nominal files ---\n")

fugita_target_types <- unique(comparison_map$fugita_cell)
fugita_tables <- list()

for (fct in fugita_target_types) {

  message("Reading Fujita file for: ", fct, " (", fugita_files[[fct]], ")")

  fug_dt <- fread(
    fugita_files[[fct]],
    select = c("gene_id", "snps", "beta", "se", "significant_by_2step_FDR",
               "REF", "ALT", "ALT_AF"),
    showProgress = FALSE
  )

  fug_dt[, gene_id := sub("\\..*", "", gene_id)]  # strip ENSG version suffix if present

  fug_dt <- fug_dt[
    !is.na(beta) & !is.na(se) & se > 0 &
    grepl("^rs", snps) & grepl("^ENSG", gene_id)
  ]

  setnames(fug_dt,
           c("gene_id", "snps", "beta", "se", "significant_by_2step_FDR",
             "REF", "ALT", "ALT_AF"),
           c("phenotype_id", "variant_id", "beta_fugita", "se_fugita", "adult_sig_raw",
             "ref_fugita", "alt_fugita", "af_fugita"))

  fug_dt[, adult_sig := adult_sig_raw == "Yes"]
  fug_dt[, adult_sig_raw := NULL]

  setkey(fug_dt, phenotype_id, variant_id)
  message("  Rows after filtering: ", nrow(fug_dt))

  fugita_tables[[fct]] <- fug_dt
  rm(fug_dt); gc()
}

# ========================================================================================
# 4. Loop over unique fetal cell types, join against applicable adult table(s),
#    attach fetal_sig, harmonize alleles
# ========================================================================================

message("\n--- Processing fetal cell types (one at a time) ---\n")

background_list <- list()
allele_summary_list <- list()

for (fct in fetal_cell_types) {

  message("\n=== Reading fetal nominal file for: ", fct, " (", nom_files[[fct]], ") ===")

  fetal_dt <- fread(
    nom_files[[fct]],
    select = c("phenotype_id", "variant_id", "slope", "slope_se", "af"),
    showProgress = FALSE
  )

  fetal_dt <- fetal_dt[
    !is.na(slope) & !is.na(slope_se) & slope_se > 0 &
    grepl("^rs", variant_id) & grepl("^ENSG", phenotype_id)
  ]

  setnames(fetal_dt, c("slope", "slope_se", "af"), c("beta_my", "se_my", "af_my"))
  setkey(fetal_dt, phenotype_id, variant_id)

  message("  Rows after filtering: ", nrow(fetal_dt))

  applicable_adult <- comparison_map %>% filter(my_cell == fct) %>% pull(fugita_cell)
  message("  Comparisons for ", fct, ": ", paste(applicable_adult, collapse = ", "))

  for (act in applicable_adult) {

    comparison_label <- paste(fct, act, sep = "_vs_")
    message("  Joining: ", comparison_label)

    joined <- fetal_dt[fugita_tables[[act]], nomatch = 0L]
    joined <- perm_tables[[fct]][joined, on = c("phenotype_id", "variant_id")]
    joined[is.na(fetal_sig), fetal_sig := FALSE]
    joined <- pvar_lookup[joined, on = "variant_id"]

    message("    Tested-in-both pairs: ", nrow(joined))

    # --- Allele harmonization ---
    joined[, is_palindromic :=
             (ref_my == "A" & alt_my == "T") | (ref_my == "T" & alt_my == "A") |
             (ref_my == "C" & alt_my == "G") | (ref_my == "G" & alt_my == "C")]

    joined[, af_diff_aligned := abs(af_my - af_fugita)]
    joined[, af_diff_flipped := abs(af_my - (1 - af_fugita))]

    joined[, allele_status := fcase(
      is.na(ref_my) | is.na(alt_my) | is.na(ref_fugita) | is.na(alt_fugita),
        "missing_allele_info",
      !is_palindromic & ref_my == ref_fugita & alt_my == alt_fugita,
        "aligned",
      !is_palindromic & ref_my == alt_fugita & alt_my == ref_fugita,
        "flipped",
      !is_palindromic,
        "mismatched",
      is_palindromic & af_diff_aligned < af_tol & af_diff_flipped >= af_tol,
        "palindromic_af_aligned",
      is_palindromic & af_diff_flipped < af_tol & af_diff_aligned >= af_tol,
        "palindromic_af_flipped",
      default = "palindromic_ambiguous"
    )]

    joined[, beta_fugita_harmonized := fcase(
      allele_status %in% c("flipped", "palindromic_af_flipped"), -beta_fugita,
      allele_status %in% c("aligned", "palindromic_af_aligned"),  beta_fugita,
      default = NA_real_
    )]

    # QC-only flag on string-resolved pairs -- not used to drop anything
    joined[allele_status == "aligned", af_consistent := af_diff_aligned < af_tol]
    joined[allele_status == "flipped", af_consistent := af_diff_flipped < af_tol]

    status_tbl <- joined[, .N, by = allele_status][order(-N)]
    status_tbl[, prop := round(N / nrow(joined), 4)]
    message("    Allele status breakdown:")
    print(status_tbl)

    n_af_inconsistent <- joined[allele_status %in% c("aligned", "flipped") &
                                 af_consistent == FALSE, .N]
    message("    'aligned'/'flipped' pairs with AF cross-check disagreement: ",
            n_af_inconsistent)

    allele_summary_list[[comparison_label]] <- as_tibble(status_tbl) %>%
      mutate(comparison = comparison_label)

    joined_tbl <- as_tibble(joined) %>%
      mutate(my_cell_type = fct, fugita_cell_type = act) %>%
      select(my_cell_type, fugita_cell_type, phenotype_id, variant_id,
             beta_my, se_my, af_my, fetal_sig,
             beta_fugita, beta_fugita_harmonized, se_fugita, af_fugita, adult_sig,
             ref_my, alt_my, ref_fugita, alt_fugita,
             allele_status, af_consistent)

    background_list[[comparison_label]] <- joined_tbl
  }

  rm(fetal_dt); gc()
}

# ========================================================================================
# 5. Summary and write output
# ========================================================================================

message("\n--- Summary: tested-in-both pairs and allele status per comparison ---\n")

comparison_summary <- tibble(
  comparison = names(background_list),
  n_tested   = map_int(background_list, nrow)
)
print(comparison_summary)

allele_summary <- bind_rows(allele_summary_list) %>%
  select(comparison, allele_status, N, prop)

message("\nAllele status breakdown, all comparisons:")
allele_summary %>%
  knitr::kable(format = "simple") %>%
  print()

n_usable <- sum(map_int(background_list, ~ sum(!is.na(.x$beta_fugita_harmonized))))
n_total  <- sum(comparison_summary$n_tested)
message("\nPairs with a usable harmonized beta (aligned/flipped, resolved): ",
        n_usable, " / ", n_total, " (", round(100 * n_usable / n_total, 2), "%)")

write_rds(
  list(
    background          = background_list,
    comparison_map       = comparison_map,
    comparison_summary    = comparison_summary,
    allele_summary        = allele_summary
  ),
  output_file
)
message("\nPairwise tested-pair background written to: ", output_file)
