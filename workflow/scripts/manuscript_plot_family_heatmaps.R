#--------------------------------------------------------------------------------------
#
#    Generate subcluster-specific family heatmaps for manuscript
#
#--------------------------------------------------------------------------------------

# L1 eQTL beta hetmaps across the parent and its 3
# L2 children for every subcluster-specific (gene, SNP) pair in that family.
# Black border marks the cell where the pair was actually classified
# subcluster-specific. 

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

message("\n\nGenerating subcluster-specific family heatmaps for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Set variables
eqtl_effects_file <- snakemake@params[["eqtl_effects"]]
gene_lookup_file   <- snakemake@params[["gene_lookup"]]

out_files <- list(
  "Glu-A" = snakemake@output[["glu_a"]],
  "Glu-B" = snakemake@output[["glu_b"]],
  "GABA"  = snakemake@output[["gaba"]],
  "NPC"   = snakemake@output[["npc"]]
)

# eqtl_effects_file <- "../results/19DEV-SPECIFICITY/eqtl_effects/eqtl_effect_sizes.rds"
# gene_lookup_file  <- "../resources/sheets/gene_lookup_hg38.tsv"

# --- Read and unpack
eqtl_effects   <- read_rds(eqtl_effects_file)
subcluster_tbl <- eqtl_effects$subcluster_specific
effect_long    <- eqtl_effects$effect_size_long

gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

# --- Relabel cell type names for plotting: Glu-UL -> Glu-A, Glu-DL -> Glu-B
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

subcluster_tbl <- subcluster_tbl |>
  mutate(parent_L1    = relabel_cell_type(parent_L1),
         cell_type_L2 = relabel_cell_type(cell_type_L2))

effect_long <- effect_long |>
  mutate(cell_type = relabel_cell_type(cell_type))

# --- Family definitions (new Glu-A/Glu-B names)
family_order <- list(
  "Glu-A" = c("Glu-A", "Glu-A-0", "Glu-A-1", "Glu-A-2"),
  "Glu-B" = c("Glu-B", "Glu-B-0", "Glu-B-1", "Glu-B-2"),
  "GABA"  = c("GABA", "GABA-0", "GABA-1", "GABA-2"),
  "NPC"   = c("NPC", "NPC-0", "NPC-1", "NPC-2")
)

# --- Build the beta matrix + significance mask for one family (unchanged
# logic from the report's build_family_matrix())
build_family_matrix <- function(family, subcluster_tbl, effect_long, gene_lookup) {
  members <- family_order[[family]]

  fam_pairs <- subcluster_tbl |>
    filter(parent_L1 == family) |>
    distinct(phenotype_id, variant_id, cell_type_L2) |>
    left_join(gene_lookup, by = c("phenotype_id" = "ensembl_gene_id")) |>
    mutate(pair_label = paste0(coalesce(external_gene_name, phenotype_id), " (", variant_id, ")"))

  if (nrow(fam_pairs) == 0) return(NULL)

  slopes <- effect_long |>
    filter(cell_type %in% members,
           phenotype_id %in% fam_pairs$phenotype_id,
           variant_id %in% fam_pairs$variant_id) |>
    inner_join(fam_pairs |> select(phenotype_id, variant_id, pair_label),
               by = c("phenotype_id", "variant_id"))

  slope_wide <- slopes |>
    select(pair_label, cell_type, slope) |>
    distinct(pair_label, cell_type, .keep_all = TRUE) |>
    pivot_wider(names_from = cell_type, values_from = slope)

  for (m in members) if (!(m %in% names(slope_wide))) slope_wide[[m]] <- NA_real_

  mat <- as.matrix(slope_wide[, members])
  rownames(mat) <- slope_wide$pair_label

  sig_mat <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(fam_pairs))) {
    pl <- fam_pairs$pair_label[i]
    sig_ct <- fam_pairs$cell_type_L2[i]
    if (pl %in% rownames(sig_mat) && sig_ct %in% colnames(sig_mat)) {
      sig_mat[pl, sig_ct] <- TRUE
    }
  }

  list(mat = mat, sig_mat = sig_mat)
}

col_fun_family <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

build_heatmap_from_data <- function(fam_data) {
  Heatmap(
    fam_data$mat,
    name = "Slope",
    col = col_fun_family,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = "grey90",
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 10),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (fam_data$sig_mat[i, j]) {
        grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 2, fill = NA))
      }
    }
  )
}

# --- Render each family to its own full-page TIFF. Height scales with row
# count (the same 0.18in/row rule the report used for Glu-UL, applied
# consistently to all four here) so dense families don't come out squashed;
# width is fixed since every family always has the same 4 columns.
width_in <- 6

for (family in names(family_order)) {
  fam_data <- build_family_matrix(family, subcluster_tbl, effect_long, gene_lookup)

  if (is.null(fam_data)) {
    message("No subcluster-specific genes for ", family, " -- skipping.")
    next
  }

  height_in <- max(8, 0.18 * nrow(fam_data$mat))
  ht <- build_heatmap_from_data(fam_data)

  message("Writing ", family, " family heatmap (", nrow(fam_data$mat), " pairs, ",
          round(height_in, 1), "in tall) -> ", out_files[[family]])

  tiff(
    filename = out_files[[family]],
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  draw(ht)
  dev.off()
}

message("Export complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
