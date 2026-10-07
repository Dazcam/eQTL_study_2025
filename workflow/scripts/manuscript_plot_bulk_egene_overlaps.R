#--------------------------------------------------------------------------------------
#
#    Generate bulk eGene overlap supplementary figure for manuscript
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  14MANUSCRIPT_PLOTS | Rule: bulk_egene_overlaps_plt
#            Upstream:   extract_unique_egenes (per-cell-type eGenes and their
#                        overlap with O'Brien 2018 bulk eQTL)
#            Downstream: None (supplementary figure)
#
# Purpose:   Per cell type, total eGenes (dark bar) and the number also significant
#            in O'Brien 2018 bulk fetal brain eQTL (light bar), grouped in three
#            blocks: L1, L2, pseudotime bins.
#            A: gene-level bulk eQTL overlap
#            B: transcript-level bulk eQTL overlap
#            - Keeps only cell types in config['cell_types']
#            - Relabels cell types (Glu-UL -> Glu-A, Glu-DL -> Glu-B, drop "-Q4-")
#
# Inputs:    egenes_per_celltype  One row per (cell_type, eGene), with level
#            overlap              One row per unique eGene: in_gene_bulk,
#                                 in_transcript_bulk
#            cell_types           config['cell_types'] (pipeline labels)
#
# Outputs:   out_file             Two-panel figure (TIFF, 600 dpi, LZW)
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

message("\n\nGenerating bulk eGene overlap plot for the manuscript ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(cowplot)

# Input and output paths
egenes_per_celltype_file <- snakemake@params[["egenes_per_celltype"]]
overlap_file             <- snakemake@params[["overlap"]]
cell_types               <- unlist(snakemake@params[["cell_types"]])
out_file                 <- snakemake@output[[1]]

# egenes_per_celltype_file <- "../results/.../egenes_per_celltype.tsv"
# overlap_file             <- "../results/.../overlap.tsv"
# cell_types               <- c("Glu-UL", "Glu-DL", "NPC", ...)   # config['cell_types']
# out_file                 <- "../results/14MANUSCRIPT_PLOTS/bulk_egene_overlaps.tiff"

# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("egenes_per_celltype_file", "overlap_file", "cell_types", "out_file"),
  value    = c(egenes_per_celltype_file, overlap_file,
               paste(cell_types, collapse = ", "), out_file)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

# Custom colour palette: identical to the LDSR figure so cell types share colours
# across every manuscript figure
custom_palette <- c(
  'Glu-A'               = '#4363d8',
  'Glu-B'               = '#00B6EB',
  'NPC'                 = '#FF5959',
  'GABA'                = '#3CBB75FF',
  'Endo-Peri'           = '#B200ED',
  'MG'                  = '#F58231',
  'OPC'                 = '#FDE725FF',
  'Trajectory-to-Glu-A' = '#E91E8C',
  'Trajectory-to-Glu-B' = '#2F4F4F'
)

# Relabel helper: Glu-UL -> Glu-A, Glu-DL -> Glu-B, strip "-Q4-" from bin names,
# e.g. "NPC-to-Glu-UL-Q4-Bin1" -> "NPC-to-Glu-A-Bin1"
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

# Dependency-free lighten (base R), avoids needing colorspace for one plot
lighten_colour <- function(hex, amount = 0.55) {
  rgb_mat   <- col2rgb(hex)
  rgb_light <- rgb_mat + (255 - rgb_mat) * amount
  rgb(rgb_light[1, ], rgb_light[2, ], rgb_light[3, ], maxColorValue = 255)
}

## Read, filter and join  ---------------------------------------------------------------
# egenes_per_celltype: one row per (cell_type, eGene)
# overlap: one row per unique eGene across all cell types (no cell_type column)
egenes_per_celltype <- read_tsv(egenes_per_celltype_file, show_col_types = FALSE)
overlap             <- read_tsv(overlap_file, show_col_types = FALSE)

# Keep only config cell types (pipeline labels, so filter before relabelling)
file_cell_types <- unique(egenes_per_celltype$cell_type)
message("Cell types in eGene file: ", paste(sort(file_cell_types), collapse = ", "))

not_in_config <- setdiff(file_cell_types, cell_types)
if (length(not_in_config) > 0) {
  message("Dropping cell types not in config['cell_types']: ",
          paste(not_in_config, collapse = ", "))
}
not_in_file <- setdiff(cell_types, file_cell_types)
if (length(not_in_file) > 0) {
  message("WARNING: config cell types absent from the eGene file: ",
          paste(not_in_file, collapse = ", "))
}

egenes_per_celltype <- egenes_per_celltype |>
  filter(cell_type %in% cell_types)

if (nrow(egenes_per_celltype) == 0) {
  stop("No rows left after filtering to config['cell_types']. ",
       "Check the cell_types param is passed and that names match the eGene file.")
}

joined <- egenes_per_celltype |>
  left_join(overlap, by = "phenotype_id")

n_unmatched <- sum(is.na(joined$in_gene_bulk))
if (n_unmatched > 0) {
  message("WARNING: ", n_unmatched, " (cell_type, phenotype_id) rows had no matching overlap ",
          "record. Check both files came from the same extract_unique_egenes run.")
}

celltype_overlap <- joined |>
  group_by(cell_type, level) |>
  summarise(
    n_egenes                  = n_distinct(phenotype_id),
    n_gene_bulk_overlap       = sum(in_gene_bulk, na.rm = TRUE),
    n_transcript_bulk_overlap = sum(in_transcript_bulk, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(cell_type = relabel_cell_type(cell_type))

message("\neGene counts per cell type:")
celltype_overlap |> knitr::kable(format = "simple", align = "l") |> print()

## Colour groups and x-axis order  ------------------------------------------------------
# Trajectory check comes first, since "NPC-to-Glu-A-Bin1" contains "Glu-A" and "NPC"
celltype_overlap <- celltype_overlap |>
  mutate(
    is_trajectory = str_detect(cell_type, "-to-"),
    main_type = case_when(
      is_trajectory & str_detect(cell_type, "to-Glu-A") ~ "Trajectory-to-Glu-A",
      is_trajectory & str_detect(cell_type, "to-Glu-B") ~ "Trajectory-to-Glu-B",
      str_detect(cell_type, "Glu-A")     ~ "Glu-A",
      str_detect(cell_type, "Glu-B")     ~ "Glu-B",
      str_detect(cell_type, "GABA")      ~ "GABA",
      str_detect(cell_type, "NPC")       ~ "NPC",
      str_detect(cell_type, "OPC")       ~ "OPC",
      str_detect(cell_type, "MG")        ~ "MG",
      str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
      TRUE ~ cell_type
    )
  )

# Three level blocks (L1, L2, pseudotime), family order within each block.
# Entries not present after the config filter are dropped below.
l1_order <- c("Glu-A", "Glu-B", "GABA", "NPC", "Endo-Peri", "MG", "OPC")
l2_order <- c("Glu-A-0", "Glu-A-1", "Glu-A-2",
              "Glu-B-0", "Glu-B-1", "Glu-B-2",
              "GABA-0", "GABA-1", "GABA-2",
              "NPC-0", "NPC-1", "NPC-2")
pt_order <- c(sort(grep("^NPC-to-Glu-A", celltype_overlap$cell_type, value = TRUE)),
              sort(grep("^NPC-to-Glu-B", celltype_overlap$cell_type, value = TRUE)))

cell_order <- c(l1_order, l2_order, pt_order)
cell_order <- cell_order[cell_order %in% celltype_overlap$cell_type]

unordered <- setdiff(celltype_overlap$cell_type, cell_order)
if (length(unordered) > 0) {
  stop("Cell types not covered by the x-axis ordering: ", paste(unordered, collapse = ", "))
}

celltype_overlap <- celltype_overlap |>
  mutate(cell_type = factor(cell_type, levels = cell_order))

# x position: 1 unit between bars, plus block_gap units at the two level
# boundaries (L1 -> L2, L2 -> pseudotime) only
block_gap <- 1
block_order <- celltype_overlap |>
  distinct(cell_type, level) |>
  mutate(level = factor(level, levels = c("L1", "L2", "pseudotime"))) |>
  arrange(level, cell_type) |>
  mutate(
    block_change = level != lag(level, default = first(level)),
    x_pos        = row_number() + cumsum(block_change) * block_gap
  )

celltype_overlap <- celltype_overlap |>
  left_join(block_order |> dplyr::select(cell_type, x_pos), by = "cell_type")

## Plots  -------------------------------------------------------------------------------
# Broadly matches the LDSR figure's base_theme; no plot titles
base_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "none",
    plot.title         = element_blank(),
    plot.margin        = margin(t = 18, r = 6, b = 6, l = 6)
  )

bar_width <- 0.9

# Two solid bar layers: total eGenes (dark) drawn first, bulk overlap (light,
# shorter) drawn on top, so it reads as the lower portion of one bar
plot_overlap_bar <- function(dat, overlap_col) {
  dat <- dat |> rename(overlap_n = all_of(overlap_col))
  light_fill <- lighten_colour(custom_palette[as.character(dat$main_type)], amount = 0.55)

  ggplot(dat, aes(x = x_pos)) +
    geom_col(aes(y = n_egenes, fill = main_type), width = bar_width, colour = "black") +
    geom_col(aes(y = overlap_n), fill = light_fill, width = bar_width, colour = "black") +
    geom_text(aes(y = n_egenes, label = n_egenes), vjust = -0.6, size = 3) +
    geom_text(aes(y = overlap_n, label = overlap_n), vjust = -0.6, size = 2.8, colour = "grey20") +
    scale_x_continuous(breaks = block_order$x_pos, labels = block_order$cell_type) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = custom_palette, name = "Cell Type") +
    base_theme +
    labs(x = NULL, y = "eGenes")
}

plot_A <- plot_overlap_bar(celltype_overlap, "n_gene_bulk_overlap")         # gene-level
plot_B <- plot_overlap_bar(celltype_overlap, "n_transcript_bulk_overlap")   # transcript-level

# Single column, A above B; labels sit in the top margin reserved by plot.margin
final_plot <- plot_grid(
  plot_A, plot_B,
  labels     = c("A", "B"),
  label_size = 20,
  label_x    = 0,
  label_y    = 1,
  hjust      = -0.3,
  vjust      = 1.1,
  ncol       = 1,
  align      = "v"
) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

## Save  --------------------------------------------------------------------------------
message("\nWriting: ", out_file)
ggsave(
  filename    = out_file,
  plot        = final_plot,
  width       = 8,
  height      = 10,
  units       = "in",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
