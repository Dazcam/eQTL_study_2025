#--------------------------------------------------------------------------------------
#
#    Generate TSS-distance density figure for manuscript
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  14MANUSCRIPT_PLOTS | Rule: tss_density_plt
#            Upstream:   17PSEUDOTIME per-bin tensorQTL perm results (NPC-to-Glu-UL)
#            Downstream: None (manuscript figure)
#
# Purpose:   Distance-to-TSS density of significant eQTL (qval < 0.05) for the
#            NPC to Glu-A trajectory, one line per Q4 pseudotime bin.
#            Logic and styling as in pseudotime_report.Rmd Section 4
#            "TSS distance -- combined". The Glu-B (NPC-to-Glu-DL) trajectory is
#            no longer shown.
#
# Inputs:    pseudotime_dir  Per-bin tensorQTL perm results per trajectory
#            exp_pc_map      Final expression PC count per bin
#
# Outputs:   out_file        Single-panel figure (TIFF, 600 dpi, LZW)
#
# Notes:     geno_pc (4) is fixed to the final eQTL setting used throughout the
#            manuscript.
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

message("\n\nGenerating TSS-distance density figure for the manuscript ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Input and output paths
pseudotime_dir <- snakemake@params[["pseudotime_dir"]]
exp_pc_map     <- snakemake@params[["exp_pc_map"]]
out_file       <- snakemake@output[[1]]

# pseudotime_dir <- "../results/17PSEUDOTIME/"
# exp_pc_map     <- list("NPC-to-Glu-UL-Q4-Bin1" = 20, ...)   # config['exp_pc_map']
# out_file       <- "../results/14MANUSCRIPT_PLOTS/tss_density.tiff"

trajectory  <- "NPC-to-Glu-UL"   # pipeline label
panel_title <- "NPC to Glu-A"
geno_pc     <- 4
quantile_n  <- 4

# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("pseudotime_dir", "out_file", "trajectory", "panel_title",
               "geno_pc", "quantile_n"),
  value    = c(pseudotime_dir, out_file, trajectory, panel_title,
               geno_pc, quantile_n)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

bin_palette <- setNames(
  colorRampPalette(c("#FDE725", "#21908C", "#440154"))(6),
  as.character(1:6)
)

# Significant eQTL TSS distances for one bin (same logic as the report's
# read_perm_tss()); NULL if the perm file is missing
read_perm_tss <- function(bin) {
  pc <- exp_pc_map[[sprintf("%s-Q%d-Bin%d", trajectory, quantile_n, bin)]]
  perm_file <- file.path(
    pseudotime_dir, trajectory, "tensorqtl", "perm",
    sprintf("Q%d_bin%d_genPC%d_expPC%d", quantile_n, bin, geno_pc, pc),
    sprintf("Q%d_bin%d_perm.cis_qtl.txt.gz", quantile_n, bin)
  )
  if (!file.exists(perm_file)) {
    message("WARNING: not found: ", perm_file)
    return(NULL)
  }
  message("Reading: ", perm_file)
  read_tsv(perm_file, show_col_types = FALSE) |>
    filter(!is.na(qval) & qval < 0.05, !is.na(start_distance)) |>
    mutate(distance_kb = start_distance / 1000, bin = bin)
}

## Load data  ---------------------------------------------------------------------------
message("Loading Q", quantile_n, " TSS distance data for ", trajectory, " ...")
tss_tbl <- map_dfr(seq_len(quantile_n), read_perm_tss)

if (nrow(tss_tbl) == 0) {
  stop("No significant eQTL loaded for ", trajectory, ". Check pseudotime_dir and exp_pc_map.")
}

message("\nSignificant eQTL per bin:")
tss_tbl |> count(bin) |> knitr::kable(format = "simple", align = "l") |> print()

n_outside <- sum(abs(tss_tbl$distance_kb) > 600)
if (n_outside > 0) {
  message("Note: ", n_outside, " eQTL beyond +/-600 kb are outside the x-axis limits ",
          "and removed by ggplot")
}

## Plot  --------------------------------------------------------------------------------
final_plot <- ggplot(tss_tbl, aes(x = distance_kb, color = as.factor(bin))) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = bin_palette, name = "Bin") +
  scale_x_continuous(limits = c(-600, 600), breaks = seq(-600, 600, 200)) +
  labs(title = panel_title, x = "Distance to TSS (kb)", y = "Density") +
  theme_minimal(base_size = 12) +
  theme(
    plot.margin      = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(face = "bold")
  )

## Save  --------------------------------------------------------------------------------
message("\nWriting: ", out_file)
ggsave(
  filename    = out_file,
  plot        = final_plot,
  width       = 8,
  height      = 5,     # single panel (was 10 for two stacked panels)
  units       = "in",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
