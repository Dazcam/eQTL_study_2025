#--------------------------------------------------------------------------------------
#
#    Generate bulk eGene overlap supplementary figure for manuscript
#
#--------------------------------------------------------------------------------------

# A: eGenes overlapping O'Brien 2018 gene-level bulk eQTL
# B: eGenes overlapping O'Brien 2018 transcript-level bulk eQTL

## Info  ------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)

# --- Set variables
egenes_per_celltype_file <- snakemake@params[["egenes_per_celltype"]]
overlap_file             <- snakemake@params[["overlap"]]
out_file                 <- snakemake@output[[1]]

# egenes_per_celltype_file <- "../results/.../egenes_per_celltype.tsv"
# overlap_file             <- "../results/.../overlap.tsv"
# out_file                 <- "../results/13MANUSCRIPT_PLOTS_TABLES/bulk_egene_overlaps.tiff"

# --- Custom colour palette -- identical to the LDSR report's (new Glu-A/Glu-B
# names), so cell types read the same colour across every manuscript figure.
custom_palette <- c(
  'Glu-A' = '#4363d8',
  'Glu-B' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF',
  'Trajectory-to-Glu-A' = '#E91E8C',
  'Trajectory-to-Glu-B' = '#2F4F4F'
)

# --- Relabel cell type names for plotting: Glu-UL -> Glu-A, Glu-DL -> Glu-B,
# and drop "-Q4-" from pseudotime trajectory bin names, e.g.
# "NPC-to-Glu-UL-Q4-Bin1" -> "NPC-to-Glu-A-Bin1"
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

# --- Dependency-free lighten (base R col2rgb/rgb) -- avoids requiring the
# colorspace package just for this one plot.
lighten_colour <- function(hex, amount = 0.55) {
  rgb_mat <- col2rgb(hex)
  rgb_light <- rgb_mat + (255 - rgb_mat) * amount
  rgb(rgb_light[1, ], rgb_light[2, ], rgb_light[3, ], maxColorValue = 255)
}

# --- Read and join
# egenes_per_celltype: one row per (cell_type, eGene) -- which cell type(s)
# each eGene was significant in.
# overlap: one row per unique eGene ACROSS all cell types -- whether it's
# gene-/transcript-bulk significant, no cell_type column at all.
egenes_per_celltype <- read_tsv(egenes_per_celltype_file, show_col_types = FALSE)
overlap             <- read_tsv(overlap_file, show_col_types = FALSE)

joined <- egenes_per_celltype |>
  left_join(overlap, by = "phenotype_id")

n_unmatched <- sum(is.na(joined$in_gene_bulk))
if (n_unmatched > 0) {
  warning(n_unmatched, " (cell_type, phenotype_id) rows had no matching overlap record -- ",
          "check egenes_per_celltype and overlap came from the same extract_unique_egenes run.")
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

# --- Family / trajectory detection (new Glu-A/Glu-B names). Trajectory check
# BEFORE the family cascade, since "NPC-to-Glu-A-Q4-Bin1" contains "Glu-A" as
# a substring and would otherwise be mis-classified.
celltype_overlap <- celltype_overlap |>
  mutate(
    is_trajectory = str_detect(cell_type, "-to-"),
    main_type = case_when(
      is_trajectory & str_detect(cell_type, "to-Glu-A") ~ "Trajectory-to-Glu-A",
      is_trajectory & str_detect(cell_type, "to-Glu-B") ~ "Trajectory-to-Glu-B",
      str_detect(cell_type, "Glu-A") ~ "Glu-A",
      str_detect(cell_type, "Glu-B") ~ "Glu-B",
      str_detect(cell_type, "GABA") ~ "GABA",
      str_detect(cell_type, "NPC") ~ "NPC",
      str_detect(cell_type, "OPC") ~ "OPC",
      str_detect(cell_type, "MG") ~ "MG",
      str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
      TRUE ~ cell_type
    )
  )

# --- x-axis order: 3 LEVEL blocks (L1, L2, pseudotime), family order within
# each block (not family-interleaved across the whole axis).
l1_order <- c("Glu-A", "Glu-B", "GABA", "NPC", "Endo-Peri", "MG", "OPC")
l2_order <- c("Glu-A-0", "Glu-A-1", "Glu-A-2",
              "Glu-B-0", "Glu-B-1", "Glu-B-2",
              "GABA-0", "GABA-1", "GABA-2",
              "NPC-0", "NPC-1", "NPC-2")
pt_order <- c(sort(grep("^NPC-to-Glu-A", celltype_overlap$cell_type, value = TRUE)),
              sort(grep("^NPC-to-Glu-B", celltype_overlap$cell_type, value = TRUE)))

cell_order <- c(l1_order, l2_order, pt_order)
cell_order <- cell_order[cell_order %in% celltype_overlap$cell_type]

celltype_overlap <- celltype_overlap |>
  mutate(cell_type = factor(cell_type, levels = cell_order))

# x position: sequential within a block, with a gap ONLY at the 2 block
# boundaries (L1->L2, L2->pseudotime) -- not at every family transition.
# block_gap sets how many x-units wide that boundary gap is (bar spacing
# within a block is always 1 unit, regardless of this value).
block_gap <- 1
block_order <- celltype_overlap |>
  arrange(cell_type) |>
  distinct(cell_type, level) |>
  mutate(level = factor(level, levels = c("L1", "L2", "pseudotime"))) |>
  arrange(level, cell_type) |>
  mutate(
    block_change = level != lag(level, default = first(level)),
    x_pos = row_number() + cumsum(block_change) * block_gap
  )

celltype_overlap <- celltype_overlap |>
  left_join(block_order |> select(cell_type, x_pos), by = "cell_type")

# --- Base theme -- broadly matches the LDSR figure's base_theme (same
# base_size, text sizes, grid/legend/strip treatment), no plot titles.
base_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(t = 18, r = 6, b = 6, l = 6)
  )

# Two SOLID bar layers, dark drawn first, light drawn second on top (light
# bar is shorter, so it renders as the bottom portion of the same bar) --
# both with black outlines. No plot titles.
bar_width <- 0.9

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

# --- Plot A: gene-level overlap
plot_A <- plot_overlap_bar(celltype_overlap, "n_gene_bulk_overlap")

# --- Plot B: transcript-level overlap
plot_B <- plot_overlap_bar(celltype_overlap, "n_transcript_bulk_overlap")

# --- Combine into final figure -- single column, A above B. Labels sit in
# the top margin reserved by plot.margin above, so they no longer overlap
# the bars; a white border is added around the whole figure.
final_plot <- plot_grid(
  plot_A, plot_B,
  labels = c("A", "B"),
  label_size = 20,
  label_x = 0,
  label_y = 1,
  hjust = -0.3,
  vjust = 1.1,
  ncol = 1,
  align = "v"
) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# High-res TIFF
ggsave(
  filename = out_file,
  plot = final_plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

message("Export complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
