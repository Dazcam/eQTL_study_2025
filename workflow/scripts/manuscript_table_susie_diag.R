#--------------------------------------------------------------------------------------
#
#    Generate SuSiE fine-mapping diagnostics table for manuscript (Supplementary Table 9)
#
#--------------------------------------------------------------------------------------

# One row per cell type, summarising SuSiE fine-mapping performance (reviewer 2
# response). Recomputed here directly from the sorted, merged cred.hp files so
# the table is independent of the SuSiE report; the logic mirrors the
# "Fine-mapping diagnostics (reviewer response)" section of susie_report.Rmd.
#
# Columns:
#   Cell Type                   Cell type (Glu-UL/Glu-DL relabelled Glu-A/Glu-B, "-Q4-"
#                               stripped from pseudotime bin names)
#   Genes with >=1 CS           Number of fine-mapped genes with at least one credible set
#   Mean CS / Gene              Mean number of credible sets per gene with >=1 CS
#   % Genes, 1 CS               % of genes with >=1 CS that have exactly 1 credible set
#   % Genes, 2 CS               % of genes with >=1 CS that have exactly 2 credible sets
#   % Genes, >=3 CS             % of genes with >=1 CS that have 3 or more credible sets
#   Median CS Size              Median number of variants per credible set
#   Median Distance to TSS (bp) Median absolute distance between the credible set index
#                               variant and the gene TSS, across all credible sets
#
# Notes for the table legend:
#   - A credible set's index variant is the variant with the highest PIP in that set.
#   - Percentages use genes with >=1 credible set as the denominator, not genes tested.
#   - Distance to TSS is |index variant position - phenotype_pos| from the gene meta file.
#
# Row order: major cell types (Glu-A, Glu-B, GABA, NPC, OPC, MG, Endo-Peri), then
# subclusters grouped by parent in the same order, then pseudotime bins.

## Info  ------------------------------------------------------------------------------

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

message("\n\nGenerating SuSiE fine-mapping diagnostics table for the manuscript ...")

# -------------------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# --- Set variables
cell_types      <- snakemake@params[["cell_types"]]
susie_files     <- snakemake@input[["susie_files"]]
gene_meta_files <- snakemake@input[["gene_meta_files"]]
xlsx_out        <- snakemake@output[["xlsx"]]

# cell_types      <- c("Glu-UL", "Glu-DL", "NPC", ...)
# susie_files     <- "../results/08SUSIE/susie/{cell_type}/{cell_type}.mrgd.srtd.susie.cred.hp.txt.gz"
# gene_meta_files <- "../results/08SUSIE/susie_input/{cell_type}_gene_meta.tsv"

message("=== SuSiE diagnostics table: ", length(cell_types), " cell types ===")

if (length(susie_files) != length(cell_types) || length(gene_meta_files) != length(cell_types)) {
  stop("Number of input files (", length(susie_files), " susie, ", length(gene_meta_files),
       " gene meta) doesn't match number of cell types (", length(cell_types),
       ") -- check the expand(..., cell_type = ...) lists match `cell_types` exactly.")
}

# --- Index variant (highest-PIP variant per credible set) with distance to TSS
# and number of credible sets per gene, for one cell type. Same as the report's
# get_index_variants(), but errors rather than skipping if a file is missing so
# no cell type can silently drop out of the manuscript table.
get_index_variants <- function(cell_type, susie_path, meta_path) {
  if (!file.exists(susie_path) || !file.exists(meta_path)) {
    stop("  [", cell_type, "] missing susie or gene meta file: ", susie_path, " / ", meta_path)
  }
  message("  [", cell_type, "] reading ", susie_path)

  susie_tbl <- suppressMessages(read_tsv(susie_path, show_col_types = FALSE))
  gene_meta <- suppressMessages(read_tsv(meta_path, show_col_types = FALSE)) |>
    select(phenotype_id, phenotype_pos)

  susie_tbl |>
    group_by(cs_id) |>
    slice_max(pip, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(phenotype_id, cs_id, cs_size, pos) |>
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

# --- Checks: every cell type has credible sets, and every index variant has a TSS
missing_ct <- setdiff(cell_types, unique(cs_summary_all$cell_type))
if (length(missing_ct) > 0) {
  stop("No credible sets found for: ", paste(missing_ct, collapse = ", "))
}
if (anyNA(cs_summary_all$dist_to_tss)) {
  stop("Missing TSS distance for some index variants (phenotype_id not in gene meta file)")
}

# --- Summaries (same definitions as the report)
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

# --- Relabel helper (display only): Glu-UL -> Glu-A, Glu-DL -> Glu-B, and strip
# "-Q4-" from pseudotime bin names, e.g. "NPC-to-Glu-DL-Q4-Bin1" ->
# "NPC-to-Glu-B-Bin1". Same convention as manuscript_table_smr_sig.R.
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

final_tbl <- final_tbl |>
  mutate(cell_type = relabel_cell_type(cell_type))

old_names <- str_detect(final_tbl$cell_type, "UL|DL|Q4")
if (any(old_names)) {
  stop("Old cell type names remain after relabelling: ",
       paste(final_tbl$cell_type[old_names], collapse = ", "))
}

# --- Row order: major types, then subclusters grouped by parent, then pseudotime bins
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
  select(cell_type, n_genes_with_cs, mean_cs_per_gene, pct_genes_1_cs, pct_genes_2_cs,
         pct_genes_3plus_cs, median_cs_size, median_dist_to_tss) |>
  mutate(across(c(mean_cs_per_gene, pct_genes_1_cs, pct_genes_2_cs, pct_genes_3plus_cs),
                ~ round(.x, 1)),
         median_dist_to_tss = round(median_dist_to_tss))

names(final_tbl) <- c("Cell Type", "Genes with \u22651 CS", "Mean CS / Gene", "% Genes, 1 CS",
                      "% Genes, 2 CS", "% Genes, \u22653 CS", "Median CS Size",
                      "Median Distance to TSS (bp)")

message("Summarised ", nrow(final_tbl), " cell types")

message("Writing xlsx (Supplementary Table 9): ", xlsx_out)
write.xlsx(
  list(`SuSiE diagnostics` = final_tbl),
  file = xlsx_out,
  overwrite = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

message("Done: ", xlsx_out)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
