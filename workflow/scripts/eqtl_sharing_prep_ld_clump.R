#--------------------------------------------------------------------------------------
#
#    Prepare LD-clumping input for fetal eGene-eSNP pairs
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: prep_ld_clump_input (single job)
#            Upstream:   extract_eqtl_universe
#            Downstream: run_ld_clump (plink --clump per gene); parse_ld_clump
#
# Purpose:   Check that genes with more than one eSNP across cell types are not
#            counting the same signal twice (LD between eSNPs):
#            - Rank each unique (eGene, eSNP) pair by its minimum qval across cell
#              types
#            - Split genes into singletons (one eSNP; passed straight through) and
#              multi-SNP genes (clumped)
#            - Write one plink association file per multi-SNP gene and the full
#              candidate SNP list for plink --extract
#
# Inputs:    eqtl_universe_file  fetal eGene-eSNP universe (extract_eqtl_universe)
#
# Outputs:   assoc_dir            {gene}.assoc (SNP, P = min qval), multi-SNP genes
#            singleton_snps_file  rds: singleton-gene (phenotype_id, variant_id) pairs
#            snp_list_file        all candidate SNPs, one per line
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

message("\n\nPreparing LD-clumping input (min-qval ranking, per-gene association files) ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
assoc_dir          <- snakemake@output[["assoc_dir"]]
singleton_snps_file <- snakemake@output[["singleton_snps"]]
snp_list_file       <- snakemake@output[["snp_list"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_universe_file", "assoc_dir", "singleton_snps_file", "snp_list_file"),
  value    = c(eqtl_universe_file, assoc_dir, singleton_snps_file, snp_list_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load universe  --------------------------------------------------------------------

message("\nLoading eQTL universe ...")
eqtl_universe <- read_rds(eqtl_universe_file)

##  SNPs-per-gene diagnostic (as in the extract_eqtl_universe log)  -------------------
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

##  Rank pairs by minimum qval  -------------------------------------------------------
# A pair can appear in several cell types; its minimum qval ranks it for
# plink --clump's index-SNP selection
message("\nComputing min qval across cell types per (eGene, eSNP) pair ...")

gene_snp_pvals <- eqtl_universe %>%
  group_by(phenotype_id, variant_id) %>%
  summarise(min_qval = min(qval, na.rm = TRUE), .groups = "drop")

message("  Ranked (eGene, eSNP) pairs: ", nrow(gene_snp_pvals))

##  Split singleton and multi-SNP genes  ----------------------------------------------
singleton_genes <- snps_per_gene %>% filter(n_snps == 1) %>% pull(phenotype_id)
multi_genes     <- snps_per_gene %>% filter(n_snps > 1)  %>% pull(phenotype_id)

singleton_snps <- gene_snp_pvals %>%
  filter(phenotype_id %in% singleton_genes) %>%
  select(phenotype_id, variant_id)

message("\nSingleton-gene (eGene, eSNP) pairs (pass through, no clumping needed): ",
        nrow(singleton_snps))

##  Write per-gene association files  -------------------------------------------------
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

##  Write output  ---------------------------------------------------------------------
write_rds(singleton_snps, singleton_snps_file)
message("\nSingleton-gene (eGene, eSNP) pairs written to: ", singleton_snps_file)

# Candidate SNP list for plink --extract in run_ld_clump (subsets the genome-wide
# bfile once rather than per gene)
all_candidate_snps <- unique_gene_snp_pairs %>% distinct(variant_id) %>% pull(variant_id)
writeLines(all_candidate_snps, snp_list_file)
message("Full candidate SNP list (", length(all_candidate_snps), " SNPs) written to: ", snp_list_file)

message("\nLD-clumping input preparation complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
