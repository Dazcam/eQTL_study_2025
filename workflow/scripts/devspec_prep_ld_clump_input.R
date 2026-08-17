#--------------------------------------------------------------------------------------
#
#    Developmental specificity: LD-independence pre-processing (step 1 of 3)
#
#--------------------------------------------------------------------------------------

#  Reviewer request: before treating fetal (eGene, eSNP) pairs as distinct signals in
#  the OCR enrichment analysis, confirm that genes with more than one mapped eSNP
#  (across the 19 fetal cell types) don't actually reflect the same underlying signal
#  via LD between those eSNPs. This is step 1 of a 3-step LD-independence sub-pipeline:
#
#    1. devspec_prep_ld_clump_input.R (this script) -- read eqtl_universe, compute a
#       single ranking p-value per (eGene, eSNP) pair (min qval across cell types),
#       split genes into singleton (1 SNP, no clumping needed) vs multi-SNP (>1 SNP,
#       candidates for LD clumping), and write one plink association file per
#       multi-SNP gene.
#    2. run_ld_clump (plink, separate envmodules rule -- kept out of the R container
#       to avoid mixing singularity/envmodules dependency management in one rule) --
#       runs plink --clump per gene against the cohort's own genotypes.
#    3. devspec_parse_ld_clump.R -- combines clump index SNPs + singleton-gene SNPs
#       into the final pruned, LD-independent (eGene, eSNP) universe used by
#       ocr_enrich.
#
#  Same diagnostic figures as previously eyeballed in devspec_extract_universe.R's log
#  (unique pairs, unique genes, SNPs-per-gene distribution, singleton vs multi-SNP gene
#  counts) are reproduced here so the two logs can be cross-checked against each other.

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

message("\n\nPreparing LD-clumping input (min-qval ranking, per-gene association files) ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
assoc_dir          <- snakemake@output[["assoc_dir"]]
singleton_snps_file <- snakemake@output[["singleton_snps"]]
snp_list_file       <- snakemake@output[["snp_list"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "assoc_dir", "singleton_snps_file", "snp_list_file"),
  value    = c(eqtl_universe_file, assoc_dir, singleton_snps_file, snp_list_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Main -----------------------------------------------------------------------------------

message("\nLoading eQTL universe ...")
eqtl_universe <- read_rds(eqtl_universe_file)

# --- Unique (eGene, eSNP) pairs and SNPs-per-gene diagnostics (matches the log output
# from the earlier devspec_extract_universe.R diagnostic addition) --------------------
message("\n--- SNPs-per-gene diagnostic ---\n")

unique_gene_snp_pairs <- eqtl_universe %>%
  distinct(phenotype_id, variant_id)

message("Unique (eGene, eSNP) pairs across all cell types: ", nrow(unique_gene_snp_pairs))
message("Unique eGenes represented: ", n_distinct(unique_gene_snp_pairs$phenotype_id))

snps_per_gene <- unique_gene_snp_pairs %>%
  count(phenotype_id, name = "n_snps")

message("\nDistribution of SNPs per gene:")
snps_per_gene %>%
  count(n_snps, name = "n_genes") %>%
  arrange(n_snps) %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

n_single_snp_genes <- sum(snps_per_gene$n_snps == 1)
n_multi_snp_genes  <- sum(snps_per_gene$n_snps > 1)

message("\nGenes with exactly 1 SNP: ", n_single_snp_genes)
message("Genes with >1 SNP (candidates for LD clumping): ", n_multi_snp_genes)

# --- Min qval across cell types, per (eGene, eSNP) pair --------------------------------
# A given (phenotype_id, variant_id) pair can appear under multiple cell types with
# different qvals; take the minimum (strongest evidence) as the pair's ranking p-value
# for plink --clump's index-SNP selection.
message("\nComputing min qval across cell types per (eGene, eSNP) pair ...")

gene_snp_pvals <- eqtl_universe %>%
  group_by(phenotype_id, variant_id) %>%
  summarise(min_qval = min(qval, na.rm = TRUE), .groups = "drop")

message("  Ranked (eGene, eSNP) pairs: ", nrow(gene_snp_pvals))

# --- Split singleton vs multi-SNP genes -------------------------------------------------
singleton_genes <- snps_per_gene %>% filter(n_snps == 1) %>% pull(phenotype_id)
multi_genes     <- snps_per_gene %>% filter(n_snps > 1)  %>% pull(phenotype_id)

singleton_snps <- gene_snp_pvals %>%
  filter(phenotype_id %in% singleton_genes) %>%
  select(phenotype_id, variant_id)

message("\nSingleton-gene (eGene, eSNP) pairs (pass through, no clumping needed): ",
        nrow(singleton_snps))

# --- Write one plink association file per multi-SNP gene ------------------------------
message("\nWriting per-gene plink association files to: ", assoc_dir)
dir.create(assoc_dir, recursive = TRUE, showWarnings = FALSE)

n_written <- 0
for (g in multi_genes) {
  gene_assoc <- gene_snp_pvals %>%
    filter(phenotype_id == g) %>%
    select(SNP = variant_id, P = min_qval)

  out_path <- file.path(assoc_dir, paste0(g, ".assoc"))
  write_tsv(gene_assoc, out_path)
  n_written <- n_written + 1
}

message("  Association files written: ", n_written, " (expected: ", length(multi_genes), ")")

if (n_written != length(multi_genes)) {
  warning("Number of association files written does not match number of multi-SNP genes.")
}

# --- Write singleton-gene SNPs and full candidate SNP list -----------------------------
write_rds(singleton_snps, singleton_snps_file)
message("\nSingleton-gene (eGene, eSNP) pairs written to: ", singleton_snps_file)

# Full candidate SNP list (unique variant_id across the whole universe), for the
# plink --extract subsetting step in run_ld_clump -- avoids re-reading the full
# genome-wide bfile once per gene.
all_candidate_snps <- unique_gene_snp_pairs %>% distinct(variant_id) %>% pull(variant_id)
writeLines(all_candidate_snps, snp_list_file)
message("Full candidate SNP list (", length(all_candidate_snps), " SNPs) written to: ", snp_list_file)

message("\nLD-clumping input preparation complete.")
