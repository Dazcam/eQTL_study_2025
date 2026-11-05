#--------------------------------------------------------------------------------------
#
#    Extract the fetal eGene-eSNP universe (L1 + L2)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: extract_eqtl_universe (single job)
#            Upstream:   05TENSORQTL tensorqtl_perm (final expPC per cell type)
#            Downstream: classify_eqtl_specificity; heterogeneity_test;
#                        prep_ld_clump_input; eqtl_effect_sizes
#
# Purpose:   - Take the lead eSNP of every sig. eGene (qval < 0.05) in each L1 and L2
#              cell type
#            - Tag each row with its level (L1/L2), L1 parent cell type and the
#              matched Jang et al. 2026 adult cell type (via the parent; NA for NPC)
#            - Log diagnostics: row counts and SNPs per gene (sizes the LD check)
#
# Inputs:    qtl_perm_files    tensorQTL permutation output, one per cell type
#            cell_types        cell types in the same order as qtl_perm_files (params)
#            cell_types_L1     fetal L1 cell types (params)
#            jang_label_match  fetal L1 -> matched Jang L1 type; NPC null (params)
#            eqtl_fdr          qval threshold for sig. eGenes (params)
#
# Outputs:   output_universe  rds: variant_id, phenotype_id, slope, slope_se, af,
#                             qval, cell_type, level, parent_cell_type,
#                             jang_label_match (plus an untracked .tsv copy)
#
# Notes:     NPC has no adult equivalent, so jang_label_match is NA.
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

message("\n\nExtracting fetal eGene/eSNP universe (L1 + L2) ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load functions
source(file.path(snakemake@scriptdir, "eqtl_sharing_functions.R"))

# Input and output paths
qtl_perm_files  <- snakemake@input[["qtl_perm_files"]]
output_universe <- snakemake@output[[1]]
cell_types      <- snakemake@params[["cell_types"]]
geno_pc         <- snakemake@params[["geno_pc"]]
norm_method     <- snakemake@params[["norm_method"]]
eqtl_fdr        <- snakemake@params[["eqtl_fdr"]]
cell_types_L1   <- unlist(snakemake@params[["cell_types_L1"]])

names(qtl_perm_files) <- cell_types  # preserve cell_type <-> file correspondence

# A priori matched Jang et al. 2026 adult cell type for each fetal L1 cell type
# (config null -> NA: NPC has no adult equivalent)
jang_label_match <- snakemake@params[["jang_label_match"]]
jang_label_match <- setNames(
  vapply(jang_label_match, function(x) if (is.null(x)) NA_character_ else x, character(1)),
  names(jang_label_match)
)

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("n_cell_types", "cell_types_L1", "geno_pc", "norm_method", "eqtl_fdr",
               "output_universe"),
  value    = c(length(cell_types), paste(cell_types_L1, collapse = ", "), geno_pc,
               norm_method, eqtl_fdr, output_universe)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  ------------------------------------------------------------------------
# Read FDR<0.05 eGene-eSNP pairs for a single cell type
# (tensorQTL perm output already reports one lead eSNP per eGene)
# slope_se is kept for the heterogeneity test (Cochran's Q).
read_sig_eqtl <- function(cell_type, file_path) {
  message("Reading sig. eQTL for: ", cell_type, " (", file_path, ")")

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  read_tsv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(qval) & qval < eqtl_fdr) %>%
    select(variant_id, phenotype_id, slope, slope_se, af, qval) %>%
    mutate(cell_type = cell_type)
}

##  Read sig. eQTL  -------------------------------------------------------------------
eqtl_universe <- map2_dfr(names(qtl_perm_files), qtl_perm_files, read_sig_eqtl)

# Dedup safeguard: one row per (cell_type, phenotype_id)
n_before <- nrow(eqtl_universe)
eqtl_universe <- eqtl_universe %>% distinct(cell_type, phenotype_id, .keep_all = TRUE)
n_after <- nrow(eqtl_universe)

if (n_before != n_after) {
  message("WARNING: removed ", n_before - n_after,
          " duplicate (cell_type, phenotype_id) rows.")
} else {
  message("No duplicate (cell_type, phenotype_id) rows found.")
}

# Tag level, parent cell type, and a priori Jang label-match
eqtl_universe <- eqtl_universe %>%
  mutate(
    level = ifelse(cell_type %in% cell_types_L1, "L1", "L2"),
    parent_cell_type = map_chr(cell_type, get_parent_cell_type, cell_types_L1 = cell_types_L1),
    jang_label_match = jang_label_match[parent_cell_type]
  )

# Sanity check against manuscript's reported ~3,121 unique eGene figure
n_unique_egenes <- n_distinct(eqtl_universe$phenotype_id)
message("\nTotal unique eGenes across all ", length(cell_types),
        " cell types (L1+L2): ", n_unique_egenes,
        " (manuscript reports ~3,121 -- confirm these match)")

message("\nRows per cell type:")
print(eqtl_universe %>% count(cell_type, level))

##  Diagnostics  ----------------------------------------------------------------------
message("\n--- Row-count relationship diagnostic ---\n")

n_rows          <- nrow(eqtl_universe)
n_unique_genes  <- n_distinct(eqtl_universe$phenotype_id)
n_unique_pairs  <- eqtl_universe %>% distinct(phenotype_id, variant_id) %>% nrow()

tibble(
  quantity = c("Total rows (cell_type x gene)",
               "Unique eGenes",
               "Unique (eGene, eSNP) pairs"),
  value    = c(n_rows, n_unique_genes, n_unique_pairs)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()

if (!(n_unique_genes <= n_unique_pairs && n_unique_pairs <= n_rows)) {
  message("WARNING: expected ordering n_unique_genes <= n_unique_pairs <= n_rows",
          " does not hold -- check for unexpected duplication or a join/dedup bug.")
} else {
  message("Row-count relationship holds as expected (n_unique_genes <= n_unique_pairs <= n_rows).")
}

message("\n--- SNPs-per-gene diagnostic (for OCR LD-independence check) ---\n")

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
message("Genes with >1 SNP (candidates for LD check): ", n_multi_snp_genes)

# Genes with >1 SNP are the candidates for LD clumping (log only)
multi_snp_gene_detail <- unique_gene_snp_pairs %>%
  semi_join(snps_per_gene %>% filter(n_snps > 1), by = "phenotype_id") %>%
  arrange(phenotype_id, variant_id)

message("\nSNPs mapped to genes with >1 SNP (first 20 rows shown):")
print(head(multi_snp_gene_detail, 20))

##  Write output  ---------------------------------------------------------------------
write_rds(eqtl_universe, output_universe)
message("\nEgene/eSNP universe written to: ", output_universe)

# Untracked TSV for readability outside R
output_universe_tsv <- sub("\\.rds$", ".tsv", output_universe)
write_tsv(eqtl_universe, output_universe_tsv)
message("Egene/eSNP universe (TSV) written to: ", output_universe_tsv)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
