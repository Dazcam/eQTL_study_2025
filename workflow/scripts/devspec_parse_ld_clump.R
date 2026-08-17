#--------------------------------------------------------------------------------------
#
#    Developmental specificity: LD-independence pre-processing (step 3 of 3)
#
#--------------------------------------------------------------------------------------

#  Final step of the LD-independence sub-pipeline (see devspec_prep_ld_clump_input.R
#  for steps 1-2). Reads every per-gene plink .clumped file produced by run_ld_clump,
#  keeps each clump's index SNP (the representative of one independent signal) and
#  drops everything absorbed into it (SP2), combines these with the singleton-gene
#  SNPs (genes with only 1 mapped eSNP, which never needed clumping), and writes the
#  final pruned, LD-independent (eGene, eSNP) universe consumed by ocr_enrich.
#
#  A gene can contribute more than one row to the final universe if plink split its
#  SNPs into more than one clump -- i.e. more than one independent signal at that
#  gene. This is intentional: independent signals at the same gene are kept as
#  separate rows; only SNPs absorbed into the same clump as a more significant SNP
#  are treated as redundant and dropped.

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

message("\n\nParsing LD clump results and building final pruned eGene-eSNP universe ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
clump_dir           <- snakemake@input[["clump_dir"]]
run_ld_clump_sentinel <- snakemake@input[["run_ld_clump_sentinel"]]
singleton_snps_file <- snakemake@input[["singleton_snps"]]
output_file         <- snakemake@output[[1]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("clump_dir", "singleton_snps_file", "output_file"),
  value    = c(clump_dir, singleton_snps_file, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# Count entries in a clump's SP2 field ("NONE", or comma-separated "rsID(1),rsID(1),...")
count_sp2 <- function(sp2_str) {
  if (is.na(sp2_str) || sp2_str == "NONE") return(0L)
  length(str_split(sp2_str, ",")[[1]])
}

# Read a single .clumped file. Whitespace-delimited with variable padding -- base
# read.table handles this correctly regardless of run-length, unlike a fixed-delimiter
# reader.
read_clump_file <- function(f) {
  gene_id <- tools::file_path_sans_ext(basename(f))

  df <- tryCatch(
    read.table(f, header = TRUE, stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) {
      message("  WARNING: failed to parse ", f, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(df) || nrow(df) == 0) {
    message("  WARNING: no clump rows found in ", f)
    return(tibble())
  }

  df %>%
    as_tibble() %>%
    mutate(phenotype_id = gene_id,
           n_sp2 = map_int(SP2, count_sp2)) %>%
    select(phenotype_id, variant_id = SNP, n_sp2)
}

# Main -----------------------------------------------------------------------------------

# --- 1. Load singleton-gene SNPs (no clumping needed) ----------------------------------
message("\nLoading singleton-gene SNPs ...")
singleton_snps <- read_rds(singleton_snps_file)
message("  Singleton-gene (eGene, eSNP) pairs: ", nrow(singleton_snps))

# --- 2. Parse all .clumped files --------------------------------------------------------
clump_files <- list.files(clump_dir, pattern = "\\.clumped$", full.names = TRUE)
message("\nFound ", length(clump_files), " .clumped files in: ", clump_dir)

clump_results <- map_dfr(clump_files, read_clump_file)

n_genes_parsed <- n_distinct(clump_results$phenotype_id)
message("  Genes successfully parsed: ", n_genes_parsed, " (expected: ", length(clump_files), ")")

if (n_genes_parsed != length(clump_files)) {
  warning("Number of genes parsed does not match number of .clumped files found -- ",
          "check WARNING messages above for files that failed to parse or were empty.")
}

# --- 3. Independent-signal diagnostics --------------------------------------------------
message("\n--- Independent-signal diagnostics (multi-SNP genes) ---\n")

signals_per_gene <- clump_results %>% count(phenotype_id, name = "n_signals")

message("Distribution of independent signals per multi-SNP gene:")
signals_per_gene %>%
  count(n_signals, name = "n_genes") %>%
  arrange(n_signals) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

n_genes_1_signal  <- sum(signals_per_gene$n_signals == 1)
n_genes_2plus_signals <- sum(signals_per_gene$n_signals > 1)

message("\nMulti-SNP genes collapsing to exactly 1 independent signal: ", n_genes_1_signal)
message("Multi-SNP genes retaining 2+ independent signals: ", n_genes_2plus_signals)

n_index_snps    <- nrow(clump_results)
n_absorbed_snps <- sum(clump_results$n_sp2)
n_multi_gene_input_snps <- n_index_snps + n_absorbed_snps

message("\nIndex (kept) SNPs from multi-SNP genes: ", n_index_snps)
message("Absorbed (dropped, LD-redundant) SNPs from multi-SNP genes: ", n_absorbed_snps)
message("Total multi-SNP-gene SNPs accounted for (index + absorbed): ", n_multi_gene_input_snps)

# --- 4. Combine into final pruned universe ----------------------------------------------
message("\n--- Building final pruned, LD-independent eGene-eSNP universe ---\n")

independent_universe <- bind_rows(
  singleton_snps,
  clump_results %>% select(phenotype_id, variant_id)
)

n_final <- nrow(independent_universe)
n_final_genes <- n_distinct(independent_universe$phenotype_id)

message("Singleton-gene pairs: ", nrow(singleton_snps))
message("Multi-SNP-gene index (independent) SNPs: ", n_index_snps)
message("---------------------------------------------------------")
message("FINAL PRUNED UNIVERSE -- total independent (eGene, eSNP) pairs: ", n_final)
message("FINAL PRUNED UNIVERSE -- total unique eGenes represented: ", n_final_genes)
message("---------------------------------------------------------")

n_dropped_total <- n_absorbed_snps
message("\nTotal SNPs dropped as LD-redundant: ", n_dropped_total,
        " (", round(100 * n_dropped_total / (nrow(singleton_snps) + n_multi_gene_input_snps), 1),
        "% of the original universe)")

# Write single tracked output ----------------------------------------------------------
write_rds(independent_universe, output_file)
message("\nFinal pruned, LD-independent eGene-eSNP universe written to: ", output_file)
