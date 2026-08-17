#--------------------------------------------------------------------------------------
#
#    Developmental specificity: fetal eGene/eSNP universe extraction
#
#--------------------------------------------------------------------------------------

#  Builds the master long-format table of FDR<0.05 eGene-eSNP pairs across all
#  L1 + L2 fetal cell types, tagging each row with:
#    - level (L1 / L2)
#    - parent_cell_type (L2 -> parent L1 label; L1 -> itself)
#    - fugita_label_match (a priori adult Fujita cell type, via parent for L2;
#      NA where no defensible adult equivalent exists, e.g. NPC)
#
#  Step 1 of the developmental-specificity pipeline (independent for now;
#  designed to fold into replication_pi1_enrichment.R later if retained for the paper).
#
#--------------------------------------------------------------------------------------
#
#  ADDITION -- SNPs-per-gene diagnostic (LD-independence check):
#
#  Before pruning eSNPs for the OCR enrichment analysis, we need to know how many
#  of the ~6,844 unique (eGene, eSNP) pairs across all 19 cell types actually
#  involve more than one distinct eSNP mapped to the same eGene -- since those are
#  the only cases where LD between the two eSNPs is even a possibility, and the
#  only cases that would need pruning before being treated as "independent" signals.
#  This step reports the distribution of SNPs-per-gene (1, 2, 3, 4+) and lists which
#  genes/SNPs fall into the >1 bucket, purely to size the problem before deciding on
#  an LD-testing approach (pairwise plink --r2 vs per-gene --clump). It does not
#  change eqtl_universe's contents or its written output -- diagnostic only, printed
#  to the log.
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

message("\n\nExtracting fetal eGene/eSNP universe (L1 + L2) ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
qtl_perm_files  <- snakemake@input[["qtl_perm_files"]]
output_universe <- snakemake@output[[1]]
cell_types      <- snakemake@params[["cell_types"]]
geno_pc         <- snakemake@params[["geno_pc"]]
norm_method     <- snakemake@params[["norm_method"]]

names(qtl_perm_files) <- cell_types  # preserve cell_type <-> file correspondence

# L1 cell types -- fixed by cell-type annotation, not a tunable pipeline parameter.
cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

# A priori biological correspondence between fetal L1 cell types and Fujita adult
# cell types. Fixed mapping -- not a tunable pipeline parameter, hence hardcoded
# here rather than in config. NPC has no defensible adult equivalent (no progenitor
# population in the adult cortex), so it is excluded from label/best-match
# comparisons and only included in the "any" matching tier downstream.
fugita_label_match <- c(
  "Glu-UL"    = "Exc",
  "Glu-DL"    = "Exc",
  "GABA"      = "Inh",
  "NPC"       = NA_character_,
  "MG"        = "Mic",
  "OPC"       = "OPC",
  "Endo-Peri" = "End"
)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_cell_types", "geno_pc", "norm_method", "output_universe"),
  value    = c(length(cell_types), geno_pc, norm_method, output_universe)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# Strip L2 sub-cluster suffix to get parent L1 cell type
# (same convention used in tensorqtl_report.Rmd's egenes_L2 %>% mutate(parent = ...))
get_parent_cell_type <- function(cell_type, cell_types_L1) {
  if (cell_type %in% cell_types_L1) return(cell_type)
  sub("-[^-]+$", "", cell_type)
}

# Read FDR<0.05 eGene-eSNP pairs for a single cell type
# (tensorQTL perm output already reports one lead eSNP per eGene)
read_sig_eqtl <- function(cell_type, file_path) {
  message("Reading sig. eQTL for: ", cell_type, " (", file_path, ")")

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  read_tsv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(qval) & qval < 0.05) %>%
    select(variant_id, phenotype_id, slope, qval) %>%
    mutate(cell_type = cell_type)
}

# Main -----------------------------------------------------------------------------------

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

# Tag level, parent cell type, and a priori Fujita label-match
eqtl_universe <- eqtl_universe %>%
  mutate(
    level = ifelse(cell_type %in% cell_types_L1, "L1", "L2"),
    parent_cell_type = map_chr(cell_type, get_parent_cell_type, cell_types_L1 = cell_types_L1),
    fugita_label_match = fugita_label_match[parent_cell_type]
  )

# Sanity check against manuscript's reported ~3,121 unique eGene figure
n_unique_egenes <- n_distinct(eqtl_universe$phenotype_id)
message("\nTotal unique eGenes across all ", length(cell_types),
        " cell types (L1+L2): ", n_unique_egenes,
        " (manuscript reports ~3,121 -- confirm these match)")

message("\nRows per cell type:")
print(eqtl_universe %>% count(cell_type, level))

# --- SNPs-per-gene diagnostic (see ADDITION note at top of file) ----------------------
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

# Full per-gene SNP listing for genes with >1 SNP -- this is the direct input
# list for whichever LD-testing approach is chosen next (pairwise plink --r2
# vs per-gene --clump). Diagnostic only, printed to the log -- not written to
# any output file at this stage.
multi_snp_gene_detail <- unique_gene_snp_pairs %>%
  semi_join(snps_per_gene %>% filter(n_snps > 1), by = "phenotype_id") %>%
  arrange(phenotype_id, variant_id)

message("\nSNPs mapped to genes with >1 SNP (first 20 rows shown):")
print(head(multi_snp_gene_detail, 20))

# Write single tracked output ----------------------------------------------------------
write_rds(eqtl_universe, output_universe)
message("\nEgene/eSNP universe written to: ", output_universe)
