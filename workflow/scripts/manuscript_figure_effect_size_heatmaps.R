#--------------------------------------------------------------------------------------
#
#    eQTL effect-size correlation heatmaps (Supp Figs 7 and 8)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: effect_size_heatmaps_fig
#            Upstream:   15 eqtl_effect_replication (Pearson r matrices)
#            Downstream: Supp Figs 7 and 8
#
# Purpose:   Plot the Pearson r of eQTL slopes computed in 15 (rows: discovery
#            cell type, FDR < 0.05 eQTL; columns: replication cell type, same
#            SNP-gene pairs, no significance cutoff); asymmetric matrices
#            - Supp Fig 7: fetal -> fetal, 19 L1 and L2 cell types
#            - Supp Fig 8: fetal L1 -> adult L1 (Jang et al.), 7 x 7
#            Manuscript cell type labels (config manuscript_labels); fetal cell
#            types in config cell_types order, each L1 type followed by its subclusters
#
# Inputs:    eqtl_effect_replication  15 eqtl_effect_replication RDS
#            cell_types               config['cell_types'] (pipeline labels)
#
# Outputs:   fetal_vs_fetal           Supp Fig 7 (TIFF, 600 dpi, LZW)
#            fetal_vs_jang            Supp Fig 8 (TIFF, 600 dpi, LZW)
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

message("\n\nPlotting effect-size correlation heatmaps ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load functions (relabel_cell_types, cell_type_levels)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
eqtl_effect_replication_file <- snakemake@params[["eqtl_effect_replication"]]
manuscript_labels            <- snakemake@params[["labels"]]
cell_types                   <- unlist(snakemake@params[["cell_types"]])
fetal_vs_fetal_file          <- snakemake@output[["fetal_vs_fetal"]]
fetal_vs_jang_file           <- snakemake@output[["fetal_vs_jang"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_effect_replication_file", "labels", "cell_types", "fetal_vs_fetal_file",
               "fetal_vs_jang_file"),
  value    = c(eqtl_effect_replication_file, paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               paste(cell_types, collapse = ", "),
               fetal_vs_fetal_file, fetal_vs_jang_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read effect-size correlations  ----------------------------------------------------
effect_repl <- read_rds(eqtl_effect_replication_file)

message("Discovery FDR threshold used upstream: ", effect_repl$settings$fdr_thresh)
message("Fetal cell types (L1+L2): ", length(effect_repl$settings$fetal_cell_types))
message("Fetal L1 cell types: ", paste(effect_repl$settings$fetal_l1, collapse = ", "))
message("Jang cell types: ", paste(effect_repl$settings$jang_cell_types, collapse = ", "))

##  Functions  ------------------------------------------------------------------------
# Diverging red-white-blue scale, r in [-1, 1]
corr_colours <- c("#2166AC", "white", "#B2182B")

# Long-format helper: matrix (rows = discovery, cols = replication) -> tibble
mat_to_long <- function(mat) {
  as_tibble(mat, rownames = "discovery") |>
    pivot_longer(-discovery, names_to = "replication", values_to = "r")
}

# Heatmap builder, styled as plot_jang_heatmap() / plot_int_heatmap() in
# manuscript_figure_replication.R, with a diverging fill scale
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

##  Supp Fig 7: fetal -> fetal (19 L1 and L2 cell types)  -----------------------------
message("\n--- Fetal -> fetal heatmap ---\n")

fetal_long <- mat_to_long(effect_repl$fetal_vs_fetal$pearson) |>
  mutate(discovery = relabel_cell_types(discovery, manuscript_labels),
         replication = relabel_cell_types(replication, manuscript_labels))

# Each L1 type followed by its subclusters; pseudotime bins not shown
fetal_order_present <- cell_type_levels(fetal_long$discovery, cell_types, manuscript_labels)
fetal_order_present <- fetal_order_present[!str_detect(fetal_order_present, "-to-")]
message("Fetal cell types present (expected 19): ", length(fetal_order_present))

fetal_long <- fetal_long |>
  mutate(discovery = factor(discovery, levels = rev(fetal_order_present)),
         replication = factor(replication, levels = fetal_order_present))

# 19 x 19 grid: smaller text and a larger canvas than the 7 x 7 heatmap
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

##  Supp Fig 8: fetal L1 -> adult L1 (Jang et al., 7 x 7)  ----------------------------
message("\n--- Fetal L1 -> Jang heatmap ---\n")

jang_l1_order  <- effect_repl$settings$jang_cell_types  # already unprefixed: Ext, IN, MG, OPC, End, Ast, OD

adult_long <- mat_to_long(effect_repl$fetal_vs_adult$pearson) |>
  mutate(discovery = relabel_cell_types(discovery, manuscript_labels))

adult_order_present <- cell_type_levels(adult_long$discovery, cell_types, manuscript_labels)
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

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
