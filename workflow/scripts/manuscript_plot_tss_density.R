#--------------------------------------------------------------------------------------
#
#    Generate combined TSS-distance density figure for manuscript
#
#--------------------------------------------------------------------------------------

# Two stacked panels, one per trajectory: eQTL distance-to-TSS density,
# coloured by Q4 pseudotime bin. Logic and styling carried over unchanged
# from pseudotime_report.Rmd's Section 4 "TSS distance -- combined" tab.
#
#   A: Glu-A to NPC (NPC-to-Glu-UL)
#   B: Glu-B to NPC (NPC-to-Glu-DL)

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

message("\n\nGenerating combined TSS-distance density figure for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)

# --- Set variables
pseudotime_dir <- snakemake@params[["pseudotime_dir"]]
exp_pc_map     <- snakemake@params[["exp_pc_map"]]
trajectories   <- snakemake@params[["trajectories"]]
out_file       <- snakemake@output[[1]]

# pseudotime_dir <- "../results/17PSEUDOTIME/"
# out_file       <- "../results/13MANUSCRIPT_PLOTS_TABLES/plots/tss_density_combined.tiff"

geno_pc <- 4

# --- Panel titles: user-specified labels (not a reversal of trajectory direction)
panel_titles <- c(
  "NPC-to-Glu-UL" = "NPC to Glu-A",
  "NPC-to-Glu-DL" = "NPC to Glu-B"
)

bin_palette <- setNames(
  colorRampPalette(c("#FDE725", "#21908C", "#440154"))(6),
  as.character(1:6)
)

# --- Read Q4 bin eQTL TSS distances for one trajectory/bin (same logic as
# the report's read_perm_tss())
read_perm_tss <- function(pseudotime_dir, trajectory, bin, exp_pc_map,
                           geno_pc = 4, quantile_n = 4) {
  pc_key <- sprintf("%s-Q%d-Bin%d", trajectory, quantile_n, bin)
  pc <- exp_pc_map[[pc_key]]
  perm_file <- file.path(
    pseudotime_dir, trajectory, "tensorqtl", "perm",
    sprintf("Q%d_bin%d_genPC%d_expPC%d", quantile_n, bin, geno_pc, pc),
    sprintf("Q%d_bin%d_perm.cis_qtl.txt.gz", quantile_n, bin)
  )
  if (!file.exists(perm_file)) { message("Not found: ", perm_file); return(NULL) }
  read_tsv(perm_file, show_col_types = FALSE) |>
    filter(!is.na(qval) & qval < 0.05, !is.na(start_distance)) |>
    mutate(distance_kb = start_distance / 1000, trajectory = trajectory, bin = bin)
}

message("Loading Q4 TSS distance data...")
tss_q4 <- map_dfr(trajectories, function(traj) {
  map_dfr(1:4, function(b) read_perm_tss(pseudotime_dir, traj, b, exp_pc_map, geno_pc = geno_pc))
})

# --- Combined density plot for one trajectory (unchanged styling from the report)
plot_tss_combined <- function(traj_data, title_label) {
  ggplot(traj_data, aes(x = distance_kb, color = as.factor(bin))) +
    geom_density() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = bin_palette, name = "Bin") +
    scale_x_continuous(limits = c(-600, 600), breaks = seq(-600, 600, 200)) +
    labs(title = title_label, x = "Distance to TSS (kb)", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )
}

plot_A <- plot_tss_combined(tss_q4 |> filter(trajectory == "NPC-to-Glu-UL"),
                             panel_titles[["NPC-to-Glu-UL"]])
plot_B <- plot_tss_combined(tss_q4 |> filter(trajectory == "NPC-to-Glu-DL"),
                             panel_titles[["NPC-to-Glu-DL"]])

final_plot <- plot_grid(plot_A, plot_B, labels = c("A", "B"), label_size = 20,
                         ncol = 1, align = "v")

message("Writing TSS density figure -> ", out_file)
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
