#--------------------------------------------------------------------------------------
#
#    Generate MaxCPP S-LDSR plot for manuscript (Figure 5)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  14MANUSCRIPT_PLOTS | Rule: ldsr_plt
#            Upstream:   09SLDSR strat_bl_v12 (MaxCPP summary, all cell types x GWAS)
#            Downstream: None. Same processing as manuscript_table_sldsr.R; keep in sync.
#
# Purpose:   S-LDSR enrichment bar charts, one facet per disorder:
#            A: L1 populations   B: L2 subclusters   C: pseudotime trajectory bins
#            - Keeps only cell types in config['cell_types']
#            - Relabels cell types (Glu-UL -> Glu-A, Glu-DL -> Glu-B, drop "-Q4-")
#            - Dashed line: per-level Bonferroni (0.05 / n cell types in that level,
#              counted from the data); dotted line: nominal P = 0.05
#
# Inputs:    in_dir      ldsr_strat_hg38_bl_v12.maxCPP.summary.tsv
#            cell_types  config['cell_types'] (pipeline labels)
#
# Outputs:   out_file    Combined three-panel figure (PDF)
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

message("\n\nGenerating MaxCPP S-LDSR plot for the manuscript ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(cowplot)

# Input and output paths
in_dir     <- snakemake@params[["in_dir"]]
cell_types <- unlist(snakemake@params[["cell_types"]])
out_file   <- snakemake@output[[1]]

# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_dir", "cell_types", "out_file"),
  value    = c(in_dir, paste(cell_types, collapse = ", "), out_file)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

# Custom colour palette (trajectory colours match the SLDSR report Rmd and are
# kept distinct from every main-cluster colour)
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

## Read and prepare table  -------------------------------------------------------------
ldsr_tbl <- read_tsv(paste0(in_dir, 'ldsr_strat_hg38_bl_v12.maxCPP.summary.tsv'),
                     show_col_types = FALSE) |>
  separate(Category, into = c("annot", "file_name"), sep = "/", remove = FALSE) |>
  separate(file_name, into = c("cell_type", "suffix"), sep = "_hg38", remove = TRUE) |>
  mutate(disorder = toupper(str_extract(suffix, "scz|bpd|mdd|adhd|ocd"))) |>
  dplyr::select(-Category, -suffix) |>
  relocate(cell_type, disorder) |>
  mutate(ldsr = if_else(`Coefficient_z-score` > 0,
                        -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) |>
  # Trajectory rows are flagged BEFORE the digit-based L1/L2 split, since bin
  # names contain digits and would otherwise be classed as L2 subclusters
  mutate(is_trajectory = str_detect(cell_type, "-to-")) |>
  mutate(level = case_when(
    is_trajectory                ~ 3L,
    str_detect(cell_type, "\\d") ~ 2L,
    TRUE                         ~ 1L
  )) |>
  filter(disorder != 'PTSD') |>
  filter(cell_type %in% cell_types) |>          # config cell types (pipeline labels)
  mutate(disorder = recode(disorder,
                           "SCZ" = "Schizophrenia",
                           "BPD" = "Bipolar Disorder")) |>
  mutate(disorder = factor(disorder,
                           levels = c("Schizophrenia", "Bipolar Disorder", "MDD", "ADHD", "OCD")))

not_in_file <- setdiff(cell_types, unique(ldsr_tbl$cell_type))
if (length(not_in_file) > 0) {
  message("WARNING: config cell types absent from the summary file: ",
          paste(not_in_file, collapse = ", "))
}

ldsr_tbl <- ldsr_tbl |>
  mutate(cell_type = relabel_cell_type(cell_type))

## Cell type order and colour groups  --------------------------------------------------
# Plots A/B: main clusters + subclusters. Patterns are anchored ("-[0-9]+$") so
# they can't catch trajectory names (e.g. "^NPC-" would match "NPC-to-Glu-A-Bin1")
cell_order_main <- c(
  "Glu-A", sort(grep("^Glu-A-[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  "Glu-B", sort(grep("^Glu-B-[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  "GABA",  sort(grep("^GABA-[0-9]+$",  ldsr_tbl$cell_type, value = TRUE)),
  "NPC",   sort(grep("^NPC-[0-9]+$",   ldsr_tbl$cell_type, value = TRUE)),
  "OPC", "MG", "Endo-Peri"
) |> unique()

# Plot C: trajectory bins, Glu-A before Glu-B, each in Bin1-4 order
# (Glu-B entries match nothing unless that trajectory is in config['cell_types'])
cell_order_traj <- c(
  sort(grep("^NPC-to-Glu-A-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  sort(grep("^NPC-to-Glu-B-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE))
) |> unique()

cell_order <- c(cell_order_main, cell_order_traj)

# Assign main cluster to subclusters / trajectory bins (trajectory check first,
# since e.g. "NPC-to-Glu-A-Bin1" would otherwise match the "NPC" or "Glu-A" pattern)
ldsr_tbl <- ldsr_tbl |>
  mutate(cell_type = factor(cell_type, levels = rev(cell_order))) |>
  mutate(main_type = case_when(
    is_trajectory & str_detect(cell_type, "to-Glu-A") ~ "Trajectory-to-Glu-A",
    is_trajectory & str_detect(cell_type, "to-Glu-B") ~ "Trajectory-to-Glu-B",
    str_detect(cell_type, "Glu-A")     ~ "Glu-A",
    str_detect(cell_type, "Glu-B")     ~ "Glu-B",
    str_detect(cell_type, "GABA")      ~ "GABA",
    str_detect(cell_type, "NPC")       ~ "NPC",
    str_detect(cell_type, "OPC")       ~ "OPC",
    str_detect(cell_type, "MG")        ~ "MG",
    str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ as.character(cell_type)
  ))

## Significance thresholds  -------------------------------------------------------------
# Per-level Bonferroni; denominator counted from the cell types actually plotted
n_tests <- ldsr_tbl |> distinct(level, cell_type) |> count(level)
message("\nBonferroni denominators per level:")
n_tests |> knitr::kable(format = "simple", align = "l") |> print()

thresh_level1  <- -log10(0.05 / n_tests$n[n_tests$level == 1])
thresh_level2  <- -log10(0.05 / n_tests$n[n_tests$level == 2])
thresh_level3  <- -log10(0.05 / n_tests$n[n_tests$level == 3])
thresh_nominal <- -log10(0.05)

## Plots  -------------------------------------------------------------------------------
base_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    panel.spacing.x    = unit(2, "lines")
  )

# One panel per level; only the Bonferroni line differs
plot_level <- function(lvl, thresh) {
  ldsr_tbl |>
    filter(level == lvl) |>
    ggplot(aes(x = ldsr, y = cell_type, fill = main_type)) +
    geom_col(width = 0.7, colour = 'black') +
    facet_grid(. ~ disorder, scales = "free_y", space = "free_y") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    geom_vline(xintercept = thresh, linetype = "dashed", color = "black") +
    geom_vline(xintercept = thresh_nominal, linetype = "dotted", color = "black") +
    scale_fill_manual(values = custom_palette) +
    base_theme +
    labs(x = expression(-log[10](P)), y = "Cell Type") +
    coord_cartesian(xlim = c(0, 4))
}

plot_A <- plot_level(1, thresh_level1)   # L1 populations
plot_B <- plot_level(2, thresh_level2)   # L2 subclusters
plot_C <- plot_level(3, thresh_level3)   # Pseudotime trajectory bins

# rel_heights are approximate, sized by row count per panel (7 / 12 / 4 rows);
# adjust if the balance looks off once rendered
final_plot <- plot_grid(
  plot_A, plot_B, plot_C,
  labels      = c("A", "B", "C"),
  label_size  = 20,
  ncol        = 1,
  align       = "v",
  rel_heights = c(0.3, 0.45, 0.2)
)

## Save  --------------------------------------------------------------------------------
message("\nWriting: ", out_file)
ggsave(
  filename = out_file,
  plot     = final_plot,
  width    = 10,     # inches, full-page width
  height   = 8.5,    # reduced for 4-row Plot C; adjust after rendering
  units    = "in",
  device   = "pdf",
  dpi      = 300
)

# High-res PNG
# ggsave(filename = out_file, plot = final_plot, width = 10, height = 7.5,
#        units = "in", dpi = 300, device = "png", bg = "white")

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
