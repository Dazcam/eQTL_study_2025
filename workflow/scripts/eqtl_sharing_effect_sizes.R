#--------------------------------------------------------------------------------------
#
#    eQTL effect sizes across fetal and adult cell types; subcluster-specific eQTL
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: eqtl_effect_sizes (single job)
#            Upstream:   extract_eqtl_universe; 12SMR nominal export ({ct}_nom.cis_qtl_pairs);
#                        Jang et al. 2026 L1 full-association files
#            Downstream: report (effect-size and subcluster tabs);
#                        16 subcluster_specific_tbl (Supp Table 3);
#                        17 family_heatmaps_fig (Supp Figs 3-6)
#
# Purpose:   For every unique fetal eGene-eSNP pair:
#            - Get the nominal slope and P in each fetal cell type (L1 + L2) and in
#              each Jang L1 cell type
#            - Build a wide effect-size matrix (fetal L1 x Jang L1)
#            - Find subcluster-specific eGenes (significant in an L2 subcluster but
#              not its L1 parent) and supercluster-specific eGenes (significant in
#              an L1 type but none of its subclusters)
#
# Inputs:    eqtl_universe_file  fetal eGene-eSNP universe (extract_eqtl_universe)
#            nom_files           fetal nominal pairs, one per cell type (L1 + L2)
#            jang_l1_files       Jang L1 full-association files
#            gene_lookup_file    hg38 gene lookup (gene symbols)
#            cell_types          fetal cell types, same order as nom_files (params)
#            cell_types_L1       fetal L1 cell types (params)
#            jang_l1_types       Jang L1 types, same order as jang_l1_files (params)
#
# Outputs:   output_file  rds list: effect_size_matrix, effect_size_long,
#                         subcluster_specific, supercluster_specific
#
# Notes:     Jang effect sizes are fixed_beta, sign-normalised to alt
#            (align_jang_beta()); Jang nominal P is Random_P.
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

message("\n\nExtracting eQTL effect sizes across fetal (L1) and Jang adult (L1) cell types ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))
source(file.path(snakemake@scriptdir, "eqtl_sharing_functions.R"))

# Input and output paths
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
nom_files          <- snakemake@input[["nom_files"]]        # all 19 fetal cell types
jang_l1_files      <- snakemake@input[["jang_l1_files"]]    # Jang's 7 L1-pooled types only
output_file        <- snakemake@output[[1]]
cell_types         <- snakemake@params[["cell_types"]]       # all 19 fetal (L1+L2)
jang_l1_types      <- snakemake@params[["jang_l1_types"]]    # 7: Ast, End, Ext, IN, MG, OD, OPC

names(nom_files)     <- cell_types
names(jang_l1_files) <- jang_l1_types

gene_lookup_file   <- snakemake@input[["gene_lookup"]]  # symbols for display; matching is on ENSG
cell_types_L1      <- unlist(snakemake@params[["cell_types_L1"]])

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_universe_file", "n_nom_files", "n_jang_l1_files",
               "gene_lookup_file", "output_file"),
  value    = c(eqtl_universe_file, length(nom_files), length(jang_l1_files),
               gene_lookup_file, output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load gene lookup and universe  ----------------------------------------------------

message("Loading gene symbol lookup ...")
gene_lookup <- fread(gene_lookup_file, select = c("ensembl_gene_id", "external_gene_name"),
                     showProgress = FALSE) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
message("  Gene symbols loaded: ", nrow(gene_lookup))

message("Loading eGene/eSNP universe ...")
eqtl_universe <- read_rds(eqtl_universe_file)
message("  Rows: ", nrow(eqtl_universe), " | Unique eGenes: ", n_distinct(eqtl_universe$phenotype_id))

# Unique pairs for the effect size matrix
unique_pairs <- eqtl_universe %>%
  distinct(phenotype_id, variant_id) %>%
  left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) %>%
  rename(gene_symbol = external_gene_name)

message("  Unique (gene, SNP) pairs: ", nrow(unique_pairs))
message("  Pairs with gene symbol: ", sum(!is.na(unique_pairs$gene_symbol)))

# Convert to data.table for fast keyed joins
pairs_dt <- as.data.table(unique_pairs[, c("phenotype_id", "variant_id")])
setkey(pairs_dt, phenotype_id, variant_id)

##  Fetal effect sizes (all L1 + L2 cell types)  --------------------------------------

message("\n--- Extracting fetal nominal effect sizes (all 19 cell types) ---\n")

fetal_effects <- list()

for (ct in cell_types) {
  message("Reading nominal file for: ", ct, " (", nom_files[[ct]], ")")

  nom_dt <- fread(
    nom_files[[ct]],
    select = c("phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"),
    showProgress = FALSE
  )
  setkey(nom_dt, phenotype_id, variant_id)

  matched <- nom_dt[pairs_dt, nomatch = 0L]
  message("  Matched pairs: ", nrow(matched))

  fetal_effects[[ct]] <- matched[, .(phenotype_id, variant_id,
                                      slope, slope_se, pval_nominal,
                                      cell_type = ct)]

  rm(nom_dt)
  gc()
}

fetal_effects_long <- rbindlist(fetal_effects)
message("\nTotal fetal effect rows: ", nrow(fetal_effects_long))

##  Jang adult L1 effect sizes  -------------------------------------------------------

message("\n--- Extracting Jang adult (L1-pooled) nominal effect sizes ---\n")

jang_effects <- list()

for (jct in jang_l1_types) {
  message("Reading Jang L1 file for: ", jct, " (", jang_l1_files[[jct]], ")")

  jang_dt <- fread(
    jang_l1_files[[jct]],
    select = c("feature", "variant_id", "ref", "alt", "Allele", "fixed_beta", "fixed_sd",
               "Fixed_P", "Random_P"),
    showProgress = FALSE
  )

  # Drop placeholder-cohort rows (Fixed_P == Random_P == 1)
  jang_dt <- filter_jang_artifacts(jang_dt, "Jang full-association")
  jang_dt[, Fixed_P := NULL]

  jang_dt[, gene_id := sub("\\..*", "", feature)]
  jang_dt <- jang_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene_id)]

  # Jang beta is relative to Allele; align to alt
  jang_dt[, beta := align_jang_beta(fixed_beta, Allele, alt)]

  setkey(jang_dt, gene_id, variant_id)

  pairs_jang <- copy(pairs_dt)
  setnames(pairs_jang, c("phenotype_id", "variant_id"), c("gene_id", "variant_id"))
  setkey(pairs_jang, gene_id, variant_id)

  matched <- jang_dt[pairs_jang, nomatch = 0L]
  message("  Matched pairs: ", nrow(matched))

  jang_effects[[jct]] <- matched[, .(phenotype_id = gene_id, variant_id,
                                      slope = beta, slope_se = fixed_sd,
                                      pval_nominal = Random_P,
                                      cell_type = paste0("Jang_", jct))]

  rm(jang_dt, pairs_jang)
  gc()
}

jang_effects_long <- rbindlist(jang_effects)
message("\nTotal Jang effect rows: ", nrow(jang_effects_long))

##  Wide effect-size matrix (fetal L1 x Jang L1)  -------------------------------------

message("\n--- Building wide effect-size matrix (fetal L1 x Jang L1-pooled only) ---\n")

fetal_l1_effects_long <- fetal_effects_long[cell_type %in% cell_types_L1]
message("  Fetal L1 rows retained for matrix: ", nrow(fetal_l1_effects_long),
        " (of ", nrow(fetal_effects_long), " total fetal rows -- L2 excluded from matrix, still",
        " retained in fetal_effects_long for subcluster analysis below)")

all_effects_long_matrix <- rbind(fetal_l1_effects_long, jang_effects_long)

effect_size_matrix <- unique_pairs %>%
  left_join(
    as_tibble(all_effects_long_matrix) %>%
      select(phenotype_id, variant_id, cell_type, slope) %>%
      pivot_wider(names_from = cell_type, values_from = slope, names_prefix = "slope_"),
    by = c("phenotype_id", "variant_id")
  )

message("  Effect size matrix: ", nrow(effect_size_matrix), " rows x ",
        ncol(effect_size_matrix), " columns")

# Quick correlation check for logging
slope_cols <- grep("^slope_", names(effect_size_matrix), value = TRUE)
if (length(slope_cols) >= 2) {
  slope_mat <- as.matrix(effect_size_matrix[, slope_cols])
  cor_mat <- cor(slope_mat, use = "pairwise.complete.obs")
  message("\nEffect-size correlation matrix (pairwise complete):")
  print(round(cor_mat, 3))
}

##  Subcluster-specific eQTL (L2 sig., parent L1 not sig.)  ---------------------------

message("\n--- Identifying subcluster-specific eQTL ---\n")

cell_types_L2 <- setdiff(cell_types, cell_types_L1)

# Get the set of FDR<0.05 eGenes per cell type from the universe
sig_per_ct <- eqtl_universe %>%
  distinct(cell_type, phenotype_id, variant_id)

subcluster_specific <- list()

for (l2 in cell_types_L2) {
  parent <- get_parent_cell_type(l2, cell_types_L1)

  # Pairs significant in L2
  l2_pairs <- sig_per_ct %>% filter(cell_type == l2)
  # Pairs significant in parent L1
  parent_pairs <- sig_per_ct %>% filter(cell_type == parent)

  # L2-specific: gene significant in L2 but not in parent L1
  l2_specific_genes <- setdiff(l2_pairs$phenotype_id, parent_pairs$phenotype_id)

  if (length(l2_specific_genes) == 0) {
    message("  ", l2, " → ", parent, ": no subcluster-specific eGenes")
    next
  }

  message("  ", l2, " → ", parent, ": ", length(l2_specific_genes), " subcluster-specific eGenes")

  # Get the L2 perm stats for these genes (slope/qval from permutation lead SNP)
  l2_stats <- eqtl_universe %>%
    filter(cell_type == l2, phenotype_id %in% l2_specific_genes) %>%
    select(phenotype_id, variant_id, slope_L2 = slope, qval_L2 = qval)

  # L2's own nominal stats (SE, raw p-value) for the same lead SNP, sourced from
  # fetal_effects_long (already built in step 2 for every cell type) rather than
  # re-reading nominal files (same source as parent_nom below)
  l2_nom <- fetal_effects_long[
    cell_type == l2 & phenotype_id %in% l2_specific_genes
  ]

  # Nominal stats from parent L1 for the same (gene, SNP) pairs
  parent_nom <- fetal_effects_long[
    cell_type == parent & phenotype_id %in% l2_specific_genes
  ]

  l2_with_parent <- l2_stats %>%
    left_join(
      as_tibble(l2_nom) %>%
        select(phenotype_id, variant_id,
               se_L2 = slope_se, pval_nom_L2 = pval_nominal),
      by = c("phenotype_id", "variant_id")
    ) %>%
    left_join(
      as_tibble(parent_nom) %>%
        select(phenotype_id, variant_id,
               slope_parent = slope, se_parent = slope_se, pval_nom_parent = pval_nominal),
      by = c("phenotype_id", "variant_id")
    ) %>%
    left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) %>%
    rename(gene_symbol = external_gene_name) %>%
    mutate(
      cell_type_L2 = l2,
      parent_L1 = parent
    )

  subcluster_specific[[l2]] <- l2_with_parent
}

subcluster_tbl <- bind_rows(subcluster_specific)
message("\nTotal subcluster-specific eGene-eSNP rows: ", nrow(subcluster_tbl))

# Also identify L1-specific (in L1 but no L2 subtype has it)
message("\n--- Identifying supercluster-specific eQTL (L1 only, no L2) ---\n")

supercluster_specific <- list()

for (l1 in cell_types_L1) {
  l2_children <- cell_types_L2[sapply(cell_types_L2, get_parent_cell_type, cell_types_L1 = cell_types_L1) == l1]

  if (length(l2_children) == 0) next

  l1_genes <- sig_per_ct %>% filter(cell_type == l1) %>% pull(phenotype_id) %>% unique()
  l2_genes <- sig_per_ct %>% filter(cell_type %in% l2_children) %>% pull(phenotype_id) %>% unique()

  l1_only <- setdiff(l1_genes, l2_genes)

  if (length(l1_only) > 0) {
    message("  ", l1, ": ", length(l1_only), " eGenes significant in L1 but no L2 subtype")

    supercluster_specific[[l1]] <- eqtl_universe %>%
      filter(cell_type == l1, phenotype_id %in% l1_only) %>%
      select(phenotype_id, variant_id, slope_L1 = slope, qval_L1 = qval) %>%
      left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) %>%
      rename(gene_symbol = external_gene_name) %>%
      mutate(parent_L1 = l1)
  }
}

supercluster_tbl <- bind_rows(supercluster_specific)
message("Total supercluster-specific eGene-eSNP rows: ", nrow(supercluster_tbl))

##  Summary  --------------------------------------------------------------------------

message("\n--- Summary ---\n")

if (nrow(subcluster_tbl) > 0) {
  message("Subcluster-specific eGenes by L2 cell type:")
  print(subcluster_tbl %>% count(cell_type_L2, parent_L1, name = "n_specific_eGenes"))

  message("\nMedian nominal p-value of subcluster-specific eGenes in parent L1:")
  print(
    subcluster_tbl %>%
      group_by(cell_type_L2) %>%
      summarise(
        median_pval_parent = median(pval_nom_parent, na.rm = TRUE),
        prop_same_direction = mean(sign(slope_L2) == sign(slope_parent), na.rm = TRUE),
        prop_nom_sig_parent = mean(pval_nom_parent < 0.05, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

##  Write output  ---------------------------------------------------------------------

write_rds(
  list(
    effect_size_matrix      = effect_size_matrix,
    effect_size_long        = as_tibble(rbind(fetal_effects_long, jang_effects_long)),
    subcluster_specific     = subcluster_tbl,
    supercluster_specific   = supercluster_tbl
  ),
  output_file
)
message("\nEffect sizes and subcluster analysis written to: ", output_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
