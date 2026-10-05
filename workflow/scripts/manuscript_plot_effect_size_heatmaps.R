#--------------------------------------------------------------------------------------
#
#    Generate effect-size correlation heatmaps for manuscript
#
#--------------------------------------------------------------------------------------

#   A. Fetal -> Fetal : discovery (FDR<0.05 in one fetal cell type) x replication
#                        (same pair's slope in every other fetal cell type, L1+L2,
#                        no significance cutoff on replication) -- asymmetric
#   B. Fetal L1 -> Adult (Jang et al. 2026) : fetal L1 discovery x the 7 Jang
#                        L1-pooled replication cell types -- asymmetric
#
#--------------------------------------------------------------------------------------

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

message("\n\nGenerating effect-size correlation heatmaps for the manuscript ...")

# -------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# --- Set variables
eqtl_effect_replication_file <- snakemake@params[["eqtl_effect_replication"]]

fetal_vs_fetal_file <- snakemake@output[["fetal_vs_fetal"]]
fetal_vs_jang_file  <- snakemake@output[["fetal_vs_jang"]]

# eqtl_effect_replication_file <- "../results/19DEV-SPECIFICITY/eqtl_effect_replication/eqtl_effect_replication.rds"

# --- Read and unpack
effect_repl <- read_rds(eqtl_effect_replication_file)

message("Discovery FDR threshold used upstream: ", effect_repl$settings$fdr_thresh)
message("Fetal cell types (L1+L2): ", length(effect_repl$settings$fetal_cell_types))
message("Fetal L1 cell types: ", paste(effect_repl$settings$fetal_l1, collapse = ", "))
message("Jang cell types: ", paste(effect_repl$settings$jang_cell_types, collapse = ", "))

# --- Relabel cell type names for plotting: Glu-UL -> Glu-A, Glu-DL -> Glu-B
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

# Diverging red-white-blue scale, r in [-1, 1] -- kept from the previous
# version of this script (ComplexHeatmap's colorRamp2 equivalent for ggplot)
corr_colours <- c("#2166AC", "white", "#B2182B")

# Long-format helper: matrix (rows = discovery, cols = replication) -> tibble
mat_to_long <- function(mat) {
  as_tibble(mat, rownames = "discovery") |>
    pivot_longer(-discovery, names_to = "replication", values_to = "r")
}

# Shared heatmap builder, styled after plot_jang_heatmap()/plot_int_heatmap()
# in manuscript_plot_replication.R -- same borders/labels/theme, diverging
# fill scale instead of pi1's sequential one.
plot_corr_heatmap <- function(df, x_lab, y_lab, cell_text_size, axis_text_size,
                              base_size = 13) {
  ggplot(df, aes(x = replication, y = discovery, fill = r)) +
    geom_tile(color = "black", lwd = 1.1, linetype = 1) +
    geom_text(aes(label = ifelse(is.na(r), "NA", sprintf("%.2f", r))),
              color = "black", size = cell_text_size) +
    scale_fill_gradientn(
      colours = corr_colours,
      limits = c(-1, 1),
      na.value = "grey80",
      name = "Pearson r"
    ) +
    coord_equal() +
    labs(x = x_lab, y = y_lab) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = axis_text_size),
      axis.text.y = element_text(face = "bold", size = axis_text_size),
      panel.grid = element_blank()
    )
}

# ========================================================================================
# A. Fetal -> Fetal (19 L1+L2 types, discovery x replication)
# ========================================================================================

message("\n--- Fetal -> fetal heatmap ---\n")

fetal_family_order <- c(
  "Glu-A", "Glu-A-0", "Glu-A-1", "Glu-A-2",
  "Glu-B", "Glu-B-0", "Glu-B-1", "Glu-B-2",
  "GABA", "GABA-0", "GABA-1", "GABA-2",
  "NPC", "NPC-0", "NPC-1", "NPC-2",
  "MG", "OPC", "Endo-Peri"
)

fetal_long <- mat_to_long(effect_repl$fetal_vs_fetal$pearson) |>
  mutate(discovery = relabel_cell_type(discovery),
         replication = relabel_cell_type(replication))

fetal_order_present <- fetal_family_order[fetal_family_order %in% unique(fetal_long$discovery)]
message("Fetal cell types present (expected 19): ", length(fetal_order_present))

fetal_long <- fetal_long |>
  mutate(discovery = factor(discovery, levels = rev(fetal_order_present)),
         replication = factor(replication, levels = fetal_order_present))

# 19x19 grid -- smaller text than the 7x7 fetal-adult grid so labels stay
# legible without shrinking to unreadable size; canvas widened to compensate
# (9x9in was sized for ComplexHeatmap's tighter grid.text layout)
ht_fetal <- plot_corr_heatmap(
  fetal_long,
  x_lab = "Prenatal cell type (replication)",
  y_lab = "Prenatal cell type (discovery)",
  cell_text_size = 2.6,
  axis_text_size = 9,
  base_size = 12
)

message("Writing fetal-vs-fetal correlation heatmap -> ", fetal_vs_fetal_file)
ggsave(
  filename = fetal_vs_fetal_file,
  plot = ht_fetal,
  width = 12, height = 12, units = "in",
  dpi = 600, device = "tiff", compression = "lzw"
)

# ========================================================================================
# B. Fetal L1 -> Adult (Jang et al. 2026, 7x7)
# ========================================================================================

message("\n--- Fetal L1 -> Jang heatmap ---\n")

fetal_l1_order <- c("Glu-A", "Glu-B", "GABA", "NPC", "MG", "OPC", "Endo-Peri")
jang_l1_order  <- effect_repl$settings$jang_cell_types  # already unprefixed: Ext, IN, MG, OPC, End, Ast, OD

adult_long <- mat_to_long(effect_repl$fetal_vs_adult$pearson) |>
  mutate(discovery = relabel_cell_type(discovery))

adult_order_present <- fetal_l1_order[fetal_l1_order %in% unique(adult_long$discovery)]
jang_order_present  <- jang_l1_order[jang_l1_order %in% unique(adult_long$replication)]
message("Fetal L1 cell types present (expected 7): ", length(adult_order_present))
message("Jang cell types present (expected 7): ", length(jang_order_present))

adult_long <- adult_long |>
  mutate(discovery = factor(discovery, levels = rev(adult_order_present)),
         replication = factor(replication, levels = jang_order_present))

ht_cross <- plot_corr_heatmap(
  adult_long,
  x_lab = "Adult cell type (replication)",
  y_lab = "Prenatal cell type (discovery)",
  cell_text_size = 4,
  axis_text_size = 13,
  base_size = 13
)

message("Writing fetal-vs-Jang correlation heatmap -> ", fetal_vs_jang_file)
ggsave(
  filename = fetal_vs_jang_file,
  plot = ht_cross,
  width = 7, height = 6, units = "in",
  dpi = 600, device = "tiff", compression = "lzw"
)

message("Export complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
