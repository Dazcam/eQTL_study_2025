#--------------------------------------------------------------------------------------
#
#    Developmental specificity: eQTL effect sizes across fetal and adult cell types
#
#--------------------------------------------------------------------------------------

#  For every unique (eGene, eSNP) pair in the universe (6,844 pairs), extracts
#  the nominal effect size (slope) and p-value from each fetal L1 and L2 cell
#  type's nominal cis-eQTL output, and from Jang et al. 2026 adult cell types.
#
#  Produces:
#    1. effect_size_matrix — wide table for correlation heatmap. SCOPE CHANGE
#       (see MIGRATION NOTE): fetal L1 only (7 types) x Jang L1-pooled only (7
#       types), no L2 on either side, for now.
#    2. subcluster_specific_eqtl — L2 eGenes not significant in parent L1,
#       with nominal stats (slope, SE, nominal p-value) from BOTH the L2
#       subcluster and its parent L1, for the report's subcluster gene table
#    3. supercluster_specific — L1 eGenes not significant in any L2 subtype
#    4. de_crossref — whether subcluster-specific genes are DE markers
#
#  (2)-(4) are unaffected by the Fujita -> Jang migration: they compare fetal
#  L1 against fetal L2 directly and never touch adult data.
#
#  Step 6 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  MIGRATION NOTE -- Fujita -> Jang et al. 2026 (SingleBrain):
#
#  effect_size_matrix's SCOPE is deliberately narrowed from the previous
#  19-fetal (L1+L2) x 7-Fujita design to 7-fetal-L1 x 7-Jang-L1-pooled, per
#  request -- L2 excluded from both sides for now (two report heatmaps:
#  fetal-L1-vs-fetal-L1 population-ordered/no dendrogram, and fetal-L1-vs-
#  adult-L1 only). fetal_effects_long is still built across all 19 fetal cell
#  types unchanged, since subcluster_specific/supercluster_specific still need
#  L2 nominal stats -- only the wide matrix used for the heatmap is narrowed.
#
#  Adult side: Jang's 7 L1-pooled cell types (Ast, End, Ext, IN, MG, OD, OPC),
#  not the 28 L2 subtypes. fixed_beta used as the effect-size column (matches
#  beta_correlation.R/classify_specificity.R's fixed-vs-random rationale --
#  see report methods text), sign-normalized against Allele (Jang's beta is
#  relative to Allele, not always alt). Random_P used for pval_nominal (the
#  informational nominal p-value shown alongside slope), consistent with using
#  the random-effects side for evidence/significance framing throughout this
#  migration and the fixed-effects side for magnitude.
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

message("\n\nExtracting eQTL effect sizes across fetal (L1) and Jang adult (L1) cell types ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
})

# Input and output paths ---------------------------------------------------------------
eqtl_universe_file <- snakemake@input[["eqtl_universe"]]
nom_files          <- snakemake@input[["nom_files"]]        # all 19 fetal cell types
jang_l1_files      <- snakemake@input[["jang_l1_files"]]    # Jang's 7 L1-pooled types only
output_file        <- snakemake@output[[1]]
cell_types         <- snakemake@params[["cell_types"]]       # all 19 fetal (L1+L2)
jang_l1_types      <- snakemake@params[["jang_l1_types"]]    # 7: Ast, End, Ext, IN, MG, OD, OPC

names(nom_files)     <- cell_types
names(jang_l1_files) <- jang_l1_types

# Gene symbol lookup — all matching on ENSG, symbols for display only.
gene_lookup_file <- "../resources/sheets/gene_lookup_hg38.tsv"

# DE marker table (hardcoded path — fixed analysis artifact, not a tunable parameter)
de_file <- "../results/13MANUSCRIPT_PLOTS_TABLES/tables/diff_exp_genes_L1_and_L2_clusters.xlsx"

# L1 cell types and parent mapping
cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

get_parent_cell_type <- function(ct, L1s) {
  if (ct %in% L1s) return(ct)
  sub("-[^-]+$", "", ct)
}

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_universe_file", "n_nom_files", "n_jang_l1_files",
               "gene_lookup_file", "de_file", "output_file"),
  value    = c(eqtl_universe_file, length(nom_files), length(jang_l1_files),
               gene_lookup_file, de_file, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# ========================================================================================
# 1. Load gene lookup and universe
# ========================================================================================

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

# ========================================================================================
# 2. Extract slopes from fetal nominal files (ALL 19 cell types, L1+L2 -- still
#    needed in full for the subcluster/supercluster analysis in section 5)
# ========================================================================================

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

# ========================================================================================
# 3. Extract effect sizes from Jang adult L1-pooled files (7 types only)
# ========================================================================================

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

  # ADDITION -- placeholder-cohort artifact filter (same fix applied
  # consistently across every script that reads Jang -- see
  # replication_pi1_enrichment.R / devspec_classify_specificity.R for the
  # original diagnosis: beta=0/sd=0 placeholder cohorts give infinite meta-
  # analysis weight, forcing Fixed_P==Random_P==1 exactly). Removed here
  # before any other use of these rows.
  n_before_artifact <- nrow(jang_dt)
  jang_dt <- jang_dt[!(Fixed_P == 1 & Random_P == 1)]
  n_removed_artifact <- n_before_artifact - nrow(jang_dt)
  message("  Removed ", n_removed_artifact,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before_artifact, " (", round(100 * n_removed_artifact / n_before_artifact, 3), "%)")
  jang_dt[, Fixed_P := NULL]

  jang_dt[, gene_id := sub("\\..*", "", feature)]
  jang_dt <- jang_dt[grepl("^rs", variant_id) & grepl("^ENSG", gene_id)]

  # Sign-normalize against Allele (see MIGRATION NOTE)
  jang_dt[, beta := ifelse(Allele == alt, fixed_beta, -fixed_beta)]

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

# ========================================================================================
# 4. Build wide effect-size matrix for correlation heatmap: fetal L1 ONLY x
#    Jang L1-pooled ONLY (L2 excluded from both sides -- see MIGRATION NOTE)
# ========================================================================================

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

# ========================================================================================
# 5. Subcluster-specific eQTL (L2 significant, parent L1 not significant)
#    -- UNCHANGED, fetal-only, no adult data involved
# ========================================================================================

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
  # re-reading nominal files -- same source as parent_nom below, for consistency.
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
  l2_children <- cell_types_L2[sapply(cell_types_L2, get_parent_cell_type, L1s = cell_types_L1) == l1]

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

# ========================================================================================
# 6. Cross-reference with DE markers -- UNCHANGED
# ========================================================================================

message("\n--- Cross-referencing with DE marker genes ---\n")

de_markers <- tibble()

if (file.exists(de_file)) {
  de_sheets <- excel_sheets(de_file)
  message("  DE file sheets: ", paste(de_sheets, collapse = ", "))

  de_markers <- map_dfr(de_sheets, function(sheet) {
    tryCatch({
      read_excel(de_file, sheet = sheet) %>%
        mutate(de_cell_type = sheet)
    }, error = function(e) {
      message("  WARNING: failed to read sheet '", sheet, "': ", e$message)
      tibble()
    })
  })

  message("  Total DE rows loaded: ", nrow(de_markers))

  gene_col <- intersect(c("gene", "gene_id", "Gene", "gene_name",
                           "phenotype_id", "ensembl_gene_id"), names(de_markers))
  if (length(gene_col) > 0) {
    gene_col <- gene_col[1]
    message("  Using DE gene column: ", gene_col)

    if (nrow(subcluster_tbl) > 0) {
      subcluster_tbl <- subcluster_tbl %>%
        mutate(
          is_de_marker = map2_lgl(phenotype_id, cell_type_L2, function(gene, ct) {
            any(de_markers[[gene_col]] == gene & de_markers$de_cell_type == ct)
          })
        )
      message("  Subcluster-specific eGenes that are DE markers: ",
              sum(subcluster_tbl$is_de_marker, na.rm = TRUE), " / ", nrow(subcluster_tbl))
    }
  } else {
    message("  WARNING: could not identify gene column in DE table. Columns: ",
            paste(names(de_markers), collapse = ", "))
  }
} else {
  message("  WARNING: DE file not found at: ", de_file)
}

# ========================================================================================
# 7. Summary stats for logging
# ========================================================================================

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

# ========================================================================================
# 8. Write output
# ========================================================================================

write_rds(
  list(
    effect_size_matrix      = effect_size_matrix,
    effect_size_long        = as_tibble(rbind(fetal_effects_long, jang_effects_long)),
    subcluster_specific     = subcluster_tbl,
    supercluster_specific   = supercluster_tbl,
    de_markers_raw          = de_markers
  ),
  output_file
)
message("\nEffect sizes and subcluster analysis written to: ", output_file)
