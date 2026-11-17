#--------------------------------------------------------------------------------------
#
#    SuSiE fine-mapping diagnostics per cell type (Supp Table 11)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: susie_diag_tbl
#            Upstream:   10 sort_susie (merged, sorted cred.hp files);
#                        10 prep_susie_gene_meta (gene TSS)
#            Downstream: Supp Table 11
#
# Purpose:   One row per cell type summarising SuSiE fine-mapping, computed
#            from the cred.hp files as in the susie_report.Rmd diagnostics
#            section:
#              Genes with >=1 CS           Fine-mapped genes with at least one credible set
#              Mean CS / Gene              Mean credible sets per gene with >=1 CS
#              % Genes, 1 / 2 / >=3 CS     Share of genes with >=1 CS having 1, 2 or >=3 sets
#              Median CS Size              Median variants per credible set
#              Median Distance to TSS (bp) Median |index variant - TSS| across credible sets
#
# Inputs:    susie_files       cred.hp file per cell type (same order as cell_types)
#            gene_meta_files   Gene meta file per cell type (phenotype_pos = TSS)
#
# Outputs:   xlsx              One sheet; rows: L1 cell types, L2 subclusters by
#                              parent, then pseudotime bins; manuscript labels
#
# Notes:     The index variant is the highest-PIP variant in a credible set.
#            Percentages use genes with >=1 credible set as the denominator.
#            Stops if a cell type has no file or no credible sets.
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

message("\n\nBuilding SuSiE fine-mapping diagnostics table ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
cell_types        <- snakemake@params[["cell_types"]]
manuscript_labels <- snakemake@params[["labels"]]
susie_files       <- snakemake@input[["susie_files"]]
gene_meta_files   <- snakemake@input[["gene_meta_files"]]
xlsx_out          <- snakemake@output[["xlsx"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("cell_types", "labels", "susie_files", "gene_meta_files", "xlsx_out"),
  value    = c(paste(cell_types, collapse = ", "),
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               paste(susie_files, collapse = ", "), paste(gene_meta_files, collapse = ", "),
               xlsx_out)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

if (length(susie_files) != length(cell_types) || length(gene_meta_files) != length(cell_types)) {
  stop("Number of input files (", length(susie_files), " susie, ", length(gene_meta_files),
       " gene meta) doesn't match number of cell types (", length(cell_types),
       ") -- check the expand(..., cell_type = ...) lists match `cell_types` exactly.")
}

##  Index variants per credible set  --------------------------------------------------
# Highest-PIP variant per credible set, its distance to the TSS and the number
# of credible sets per gene. As the report's get_index_variants(), but stops if
# a file is missing so no cell type drops out of the table.
get_index_variants <- function(cell_type, susie_path, meta_path) {
  if (!file.exists(susie_path) || !file.exists(meta_path)) {
    stop("  [", cell_type, "] missing susie or gene meta file: ", susie_path, " / ", meta_path)
  }
  message("  [", cell_type, "] reading ", susie_path)

  susie_tbl <- read_tsv(susie_path, show_col_types = FALSE)
  gene_meta <- read_tsv(meta_path, show_col_types = FALSE) |>
    dplyr::select(phenotype_id, phenotype_pos)

  susie_tbl |>
    group_by(cs_id) |>
    slice_max(pip, n = 1, with_ties = FALSE) |>
    ungroup() |>
    dplyr::select(phenotype_id, cs_id, cs_size, pos) |>
    left_join(gene_meta, by = "phenotype_id") |>
    mutate(dist_to_tss = abs(pos - phenotype_pos)) |>
    group_by(phenotype_id) |>
    mutate(n_cs_for_gene = n()) |>
    ungroup() |>
    mutate(cell_type = cell_type)
}

cs_summary_all <- pmap_dfr(
  list(cell_types, susie_files, gene_meta_files),
  get_index_variants
)

# Every cell type has credible sets and every index variant has a TSS
missing_ct <- setdiff(cell_types, unique(cs_summary_all$cell_type))
if (length(missing_ct) > 0) {
  stop("No credible sets found for: ", paste(missing_ct, collapse = ", "))
}
if (anyNA(cs_summary_all$dist_to_tss)) {
  stop("Missing TSS distance for some index variants (phenotype_id not in gene meta file)")
}

##  Summaries per cell type (definitions as in the report)  ---------------------------
cs_per_gene <- cs_summary_all |>
  distinct(cell_type, phenotype_id, n_cs_for_gene) |>
  group_by(cell_type) |>
  summarise(
    n_genes_with_cs    = n_distinct(phenotype_id),
    mean_cs_per_gene   = mean(n_cs_for_gene),
    pct_genes_1_cs     = mean(n_cs_for_gene == 1) * 100,
    pct_genes_2_cs     = mean(n_cs_for_gene == 2) * 100,
    pct_genes_3plus_cs = mean(n_cs_for_gene >= 3) * 100,
    .groups = "drop"
  )

cs_size <- cs_summary_all |>
  group_by(cell_type) |>
  summarise(median_cs_size = median(cs_size), .groups = "drop")

dist_to_tss <- cs_summary_all |>
  group_by(cell_type) |>
  summarise(median_dist_to_tss = median(dist_to_tss), .groups = "drop")

final_tbl <- cs_per_gene |>
  left_join(cs_size, by = "cell_type") |>
  left_join(dist_to_tss, by = "cell_type")

##  Manuscript labels and row order  --------------------------------------------------
final_tbl <- final_tbl |>
  mutate(cell_type = relabel_cell_types(cell_type, manuscript_labels))

old_names <- str_detect(final_tbl$cell_type, "UL|DL|Q4")
if (any(old_names)) {
  stop("Old cell type names remain after relabelling: ",
       paste(final_tbl$cell_type[old_names], collapse = ", "))
}

# Row order: major types, then subclusters grouped by parent, then pseudotime bins
major <- c("Glu-A", "Glu-B", "GABA", "NPC", "OPC", "MG", "Endo-Peri")

final_tbl <- final_tbl |>
  mutate(
    tier   = case_when(cell_type %in% major ~ 1L,
                       str_detect(cell_type, "^NPC-to-") ~ 3L,
                       TRUE ~ 2L),
    parent = case_when(tier == 1L ~ cell_type,
                       tier == 2L ~ str_remove(cell_type, "-[0-9]+$"),
                       tier == 3L ~ str_extract(cell_type, "Glu-[AB]")),
    parent_ord = match(parent, major)
  )

if (anyNA(final_tbl$parent_ord)) {
  stop("Could not assign a parent cell type for: ",
       paste(final_tbl$cell_type[is.na(final_tbl$parent_ord)], collapse = ", "))
}

final_tbl <- final_tbl |>
  arrange(tier, parent_ord, cell_type) |>
  dplyr::select(cell_type, n_genes_with_cs, mean_cs_per_gene, pct_genes_1_cs, pct_genes_2_cs,
         pct_genes_3plus_cs, median_cs_size, median_dist_to_tss) |>
  mutate(across(c(mean_cs_per_gene, pct_genes_1_cs, pct_genes_2_cs, pct_genes_3plus_cs),
                ~ round(.x, 1)),
         median_dist_to_tss = round(median_dist_to_tss))

names(final_tbl) <- c("Cell Type", "Genes with \u22651 CS", "Mean CS / Gene", "% Genes, 1 CS",
                      "% Genes, 2 CS", "% Genes, \u22653 CS", "Median CS Size",
                      "Median Distance to TSS (bp)")

message("Summarised ", nrow(final_tbl), " cell types")

##  Write xlsx  -----------------------------------------------------------------------
message("Writing xlsx: ", xlsx_out)
write.xlsx(
  list(`SuSiE diagnostics` = final_tbl),
  file = xlsx_out,
  overwrite = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
