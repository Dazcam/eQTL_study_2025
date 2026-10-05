#--------------------------------------------------------------------------------------
#
#    Generate eQTL p1 replication plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: Pie chart of shared eGenes accross cell types
# B: Upset Plot
# C: Internal Pi1 heatmap
# D: Fetal vs. adult Pi1 heatmap (Jang et al. 2026; Glu-A, Glu-B, GABA, NPC only)
# E: NEUROD1 (rs115583150) eQTL boxplot -- Glu-B
# F: Glu-A vs Ext beta correlation
# G: PLBD2 (rs12825284) eQTL boxplot -- Glu-A

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

message("\n\nGenerating replication plot for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)
library(UpSetR)
library(grid)

# --- Set variables
in_dir <- snakemake@params[['in_dir']]
internal_dir <- snakemake@params[['internal_dir']]
jang_pi1_dir <- snakemake@params[['jang_pi1_dir']]
jang_beta_file <- snakemake@input[['beta_file']]
out_file <- snakemake@output[[1]]

# in_dir <- '../results/05TENSORQTL/tensorqtl_perm/'
# internal_dir <- "../results/06QTL-REPLICATION/internal/"
# jang_pi1_dir <- "../results/19DEV-SPECIFICITY/pi1_jang/"
# jang_beta_file <- "../results/19DEV-SPECIFICITY/<step5_output>.rds"
# out_dir <- "../results/14MANUSCRIPT_PLOTS/"

# NOTE: these old-format names (Glu-UL / Glu-DL) drive file paths for the eQTL
# input files (Panel A) and the internal Pi1 comparison (Panel C), which are
# still stored under the old naming on disk. Do not rename these - the
# relabel_cell_type() helper below handles renaming for display only.
cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG")

# Exp PC map (old names - keys file paths for Panels A and C)
expPC_map <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30
)

# --- Relabel map for cell type name changes, for display/plotting only
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

cell_types_new <- relabel_cell_type(cell_types)

# Define custom color palette (new Glu-A/Glu-B keys, since this is applied
# only to already-relabelled data downstream)
custom_palette <- c(
  'Glu-A' = '#4363d8',
  'Glu-B' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF'
)

eqtl_list <- list()
genPC <- 4
norm_method <- 'quantile'

gene_lookup <- read_tsv('../resources/sheets/gene_lookup_hg38.tsv') |>
  select(gene = ensembl_gene_id, symbol = external_gene_name)


# ----- Pie chart shared eQTL amoung L1 cell types (Panel A - unchanged) -----
for (cell_type in cell_types) {
  
  expPC <- expPC_map[[cell_type]]
  
  f <- file.path(
    in_dir,
    paste0(cell_type, "_", norm_method, "_genPC_", genPC, "_expPC_", expPC),
    paste0(cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz")  
  )
  
  if (file.exists(f)) {
    dat <- read_tsv(f, show_col_types = FALSE) %>%
      filter(!is.na(start_distance)) %>%
      mutate(distance_kb = start_distance / 1000) %>%
      filter(qval < 0.05)   # significance cutoff
    
    if (nrow(dat) > 0) {
      eqtl_list[[cell_type]] <- dat %>% mutate(cell_type = cell_type)
    }
  }
}

eqtls_all <- bind_rows(eqtl_list)

# Relabel here (feeds both Panel A's counts, which are cell-type-agnostic,
# and Panel B's set names, which do need the new labels)
gene_cell <- eqtls_all %>%
  distinct(phenotype_id, cell_type) %>%
  mutate(cell_type = relabel_cell_type(cell_type))

gene_counts <- gene_cell %>%
  count(phenotype_id) %>%  
  count(n)  

label_threshold <- 0.08

pie_dat <- gene_counts %>%
  mutate(label = paste0(n, " cell type", ifelse(n > 1, "s", "")),
         count = nn) %>%
  arrange(n) %>%
  mutate(fraction = count / sum(count),
         ymax = cumsum(fraction),
         ymin = lag(ymax, default = 0),
         midpoint = (ymax + ymin) / 2,
         label = factor(label, levels = label)) |>
  mutate(large_slice = fraction >= label_threshold)

mono_cols <- colorRampPalette(c("#08306B", "#DEEBF7"))(nrow(pie_dat))

# --- Build the pie as real Cartesian polygons instead of geom_rect() +
# coord_polar(). coord_polar applies its transform AFTER stat/geom
# computation, so geom_text_repel (which repels labels in PRE-transform
# space) can't account for it -- that mismatch is what was misaligning the
# small-slice labels and their connector lines. Building genuine (x, y)
# points up front means there is only one coordinate system, so repel's own
# leader lines (segment.color below) line up with the labels by construction.
# Calibrated to match the original orientation exactly: y = 0 sits at
# 12 o'clock, increasing y goes clockwise -- theta = pi/2 - 2*pi*fraction
# (ggplot2's actual coord_polar(theta = "y") default direction/start).
n_arc_pts <- 100

make_wedge <- function(row, r = 1) {
  theta <- pi / 2 - 2 * pi * c(row$ymin, seq(row$ymin, row$ymax, length.out = n_arc_pts), row$ymax)
  tibble(
    label = row$label,
    x = c(0, r * cos(theta), 0),
    y = c(0, r * sin(theta), 0)
  )
}

wedge_polygons <- pie_dat %>%
  split(seq_len(nrow(.))) %>%
  map_dfr(make_wedge)

# Inner labels (large slices) at r = 0.6, same radial position as the
# original's x = 0.5 xmax/xmin = 0/1 rect gave visually
inner_labels <- pie_dat %>%
  filter(large_slice) %>%
  mutate(theta = pi / 2 - 2 * pi * midpoint,
         x = 0.6 * cos(theta),
         y = 0.6 * sin(theta))

# Outer labels (small slices): fully deterministic fan-out, NOT geom_text_repel.
# ggrepel's placement is an adaptive physics simulation that depends on the
# panel's actual rendered aspect ratio -- it can look correct in one render
# and come out broken in another (e.g. once embedded in the full multi-panel
# figure at a different width/height than a standalone test). Spacing labels
# out by hand, left-to-right in the same clockwise order the slices appear,
# removes that instability entirely: the same input always produces the same
# layout, regardless of how the panel ends up sized.
outer_labels <- pie_dat %>%
  filter(!large_slice) %>%
  arrange(midpoint) %>%
  mutate(
    theta  = pi / 2 - 2 * pi * midpoint,
    x_edge = 1 * cos(theta),
    y_edge = 1 * sin(theta),
    x_label = seq(-0.4 * (n() - 1) / 2, 0.4 * (n() - 1) / 2, length.out = n()),
    y_label = 1.3
  )

pie_chart <- ggplot() +
  geom_polygon(data = wedge_polygons, aes(x = x, y = y, group = label, fill = label),
               color = "black", linewidth = 0.2) +
  geom_text(data = inner_labels, aes(x = x, y = y, label = count),
            color = "white", size = 5, fontface = "bold") +
  geom_segment(data = outer_labels,
               aes(x = x_edge, y = y_edge, xend = x_label, yend = y_label - 0.06),
               color = "grey60", linewidth = 0.4) +
  geom_text(data = outer_labels, aes(x = x_label, y = y_label, label = count),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = mono_cols) +
  coord_equal(clip = "off") +
  theme_void(base_size = 15) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin = margin(t = 10, r = 0, b = 10, l = 0),
        plot.margin = margin(t = 60, r = 30, b = 30, l = 30))

# --- Upset Plot (Panel B - relabelled set names) ----
gene_by_cell <- gene_cell %>%
  pivot_wider(names_from = cell_type, values_from = cell_type,
              values_fill = 0, values_fn = function(x) 1) %>%
  column_to_rownames("phenotype_id")

# Explicit set order (matching every other manuscript figure's L1 cell-type
# order) instead of UpSetR's default size-based sort. This also fixes
# sets.bar.color below, which was mismatched (NPC/GABA and OPC/Endo-Peri
# swapped) because the colour vector's order didn't match UpSetR's own
# (size-sorted) set order when no explicit `sets=` was given.
# Confirmed from the rendered output: UpSetR draws the first-listed set at the
# BOTTOM and the last-listed at the TOP, so the vector below is written
# bottom-to-top (Endo-Peri first -> bottom, Glu-A last -> top).
set_order <- c("Endo-Peri", "MG", "OPC", "NPC", "GABA", "Glu-B", "Glu-A")

# Render to PNG
tmp_upset <- tempfile(fileext = ".png")

# UpSetR's outer margins are fixed in absolute size, not proportional to the
# device -- so shrinking the canvas (2800 -> 2000) compressed the chart itself
# while the margins stayed the same size, making the whitespace proportionally
# WORSE, not better. Going the other way -- widening well past the original
# 2800 -- gives the chart more room relative to those fixed margins and
# spreads it out as intended. Height (top/bottom whitespace) still unchanged.
png(tmp_upset, width = 3400, height = 1800, res = 300)

upset(
  gene_by_cell,
  sets           = set_order,
  keep.order     = TRUE,
  nsets          = length(set_order),
  order.by       = "freq",
  nintersects    = 20,
  sets.bar.color = custom_palette[set_order],
  point.size     = 3.8,
  line.size      = 2,
  # text.scale: c(intersection size title, intersection size tick labels, 
  #               set size title, set size tick labels, 
  #               set names, numbers above bars)
  text.scale     = c(2.2, 2, 1.8, 1.8, 1.8, 1.7)
)

dev.off()

# Re-import as grob for cowplot
img <- png::readPNG(tmp_upset)

upset_grob <- grid::rasterGrob(
  img,
  interpolate = TRUE
)

upset_plt <- ggplotify::as.ggplot(upset_grob) +
  theme(
    plot.margin = margin(t = 40, r = 0, b = 40, l = 0, unit = "pt")
  )

# --- Internal pi1 heatmap (Panel C - relabelled + new axis titles) -----
read_pi1_results <- function(ct, ref_ct) {
  
  expPC <- expPC_map[[ct]]
  ref_expPC <- expPC_map[[ref_ct]]
  
  file_path <- paste0(internal_dir, ct, "_vs_", ref_ct, "_quantile_genPC_4_expPC_", 
                      expPC, "_expPCref_", ref_expPC,"_pi1_results_tbl.rds")
  if (file.exists(file_path)) {
    pi1_list <- read_rds(file_path)
    
    tibble(
      cell_type = ct,
      ref_cell_type = ref_ct,
      pi1_forward = pi1_list$forward$pi1,
      pi1_reverse = pi1_list$reverse$pi1,
      prop_replicating_forward = pi1_list$forward$prop_replicating,
      prop_replicating_reverse = pi1_list$reverse$prop_replicating,
      prop_same_direction_forward = pi1_list$forward$prop_same_direction,
      prop_same_direction_reverse = pi1_list$reverse$prop_same_direction
    )
  } else {
    message("File missing: ", file_path)
    NULL
  }
}

# Generate all combinations (old names - matches file paths on disk)
combinations <- expand.grid(cell_type = cell_types, ref_cell_type = cell_types, stringsAsFactors = FALSE)

# Read all pi1 results
pi1_result_tbl <- map2_dfr(combinations$cell_type, combinations$ref_cell_type, read_pi1_results)

# Prepare the full square data without manual mirroring
pi1_square_tbl <- pi1_result_tbl %>%
  mutate(
    # Use the specific Forward result for every unique combination provided by expand.grid
    pi1_final = ifelse(cell_type == ref_cell_type, 1.0, pi1_forward)
  ) %>%
  select(ref = cell_type, query = ref_cell_type, pi1 = pi1_final) %>%
  mutate(
    ref = relabel_cell_type(ref),
    query = relabel_cell_type(query),
    # rev() ensures the first cell type is at the top (standard matrix view)
    query = factor(query, levels = rev(cell_types_new)), 
    ref = factor(ref, levels = cell_types_new)
  )

# Updated heatmap function
plot_int_heatmap <- function(df) {
  ggplot(df, aes(x = query, y = ref, fill = pi1)) +
    geom_tile(color = "black", lwd = 1.1, linetype = 1) +
    geom_text(aes(label = ifelse(is.na(pi1), "NA", sprintf("%.2f", pi1))),
              color = "black", size = 3.5) +
    scale_fill_gradientn(
      colours = c("white", "yellow", "red"),
      limits = c(0.5, 1.0),
      na.value = "grey80",
      name = expression(pi[1])
    ) +
    coord_equal() +
    labs(x = "Prenatal cell type (replication)", y = "Prenatal cell type (discovery)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

# Generate only the forward heatmap (no title)
pi1_int_heatmap <- plot_int_heatmap(pi1_square_tbl)


# --- Jang pi1 heatmap (Panel D - replaces Fugita; restricted to Glu-A, Glu-B,
# GABA, NPC rows in that order, all 7 Jang adult reference types as columns) -----
jang_L1_cell_types <- c("Glu-UL", "Glu-DL", "GABA", "NPC")  # old names - match file names on disk
jang_adult_cell_types <- c("Ast", "End", "Ext", "IN", "MG", "OD", "OPC")

read_jang_pi1_results <- function(ct, ref_ct, jang_dir) {
  file_path <- paste0(jang_dir, ct, "_vs_", ref_ct, "_pi1_results_tbl.rds")
  if (file.exists(file_path)) {
    pi1_list <- read_rds(file_path)
    
    tibble(
      cell_type = ct,
      ref_cell_type = ref_ct,
      pi1_forward = pi1_list$forward$pi1,
      pi1_reverse = pi1_list$reverse$pi1
    )
  } else {
    message("File missing: ", file_path)
    NULL
  }
}

jang_combinations <- expand.grid(cell_type = jang_L1_cell_types, ref_cell_type = jang_adult_cell_types,
                                 stringsAsFactors = FALSE)

pi1_jang_tbl <- map2_dfr(jang_combinations$cell_type, jang_combinations$ref_cell_type,
                         read_jang_pi1_results, jang_dir = jang_pi1_dir)

pi1_jang_square_tbl <- pi1_jang_tbl %>%
  select(cell_type, ref_cell_type, pi1 = pi1_forward) %>%
  mutate(
    cell_type = relabel_cell_type(cell_type),
    # rev() so Glu-A appears at the top (standard matrix view)
    cell_type = factor(cell_type, levels = rev(relabel_cell_type(jang_L1_cell_types))),
    ref_cell_type = factor(ref_cell_type, levels = jang_adult_cell_types)
  )

# Function to generate the Jang heatmap
plot_jang_heatmap <- function(df) {
  ggplot(df, aes(x = ref_cell_type, y = cell_type, fill = pi1)) +
    geom_tile(color = "black", lwd = 1.1, linetype = 1) +
    geom_text(aes(label = ifelse(is.na(pi1), "NA", sprintf("%.2f", pi1))),
              color = "black", size = 3.5) +
    scale_fill_gradientn(
      colours = c("white", "yellow", "red"),
      limits = c(0, 1.0),
      na.value = "grey80",
      name = expression(pi[1])
    ) +
    coord_equal() +
    labs(x = "Adult cell type (replication)", y = "Prenatal cell type (discovery)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(50, 50, 50, 50, unit = "pt")
    ) 
}

pi1_jang_heatmap <- plot_jang_heatmap(pi1_jang_square_tbl)


### --- beta correlation plt (Panel E, right-hand side - Jang data) -----
jang_beta_list <- read_rds(jang_beta_file)
paired_betas_all <- jang_beta_list$paired_betas

# genes dropped from the table entirely (not just left unlabelled)
exclude_genes <- c("ABCC8", "CLHC1")

make_jang_beta_plot <- function(paired_betas_all, my_ct, jang_ct, gene_lookup,
                                exclude_genes = character(0), label_genes = NULL,
                                ct_labels) {
  
  paired_betas <- paired_betas_all %>%
    filter(my_cell_type == my_ct, jang_cell_type == jang_ct) %>%
    left_join(gene_lookup, by = "gene") %>%
    mutate(label_text = ifelse(is.na(symbol) | symbol == "NA", gene, symbol)) %>%
    filter(!(label_text %in% exclude_genes))
  
  ## Correlation
  cor_val <- cor(
    paired_betas$beta_my,
    paired_betas$beta_jang,
    use = "complete.obs",
    method = "pearson"
  )
  
  cor_label <- sprintf("r = %.2f", cor_val)
  
  # --- Point classification: two-check "significant_jang" convention (see
  # devspec_beta_correlation.R / devspec_report.Rmd), replacing the old
  # sign-mismatch-only "discordant" rule. A pair is only highlighted as
  # concordant/discordant if it is ALSO significant in Jang (i.e. it IS that
  # gene's lead SNP in Jang's top-association file AND qval < 0.05) -
  # everything else, including sign-mismatched pairs that aren't
  # Jang-significant, falls into a single grey "other" bucket.
  paired_betas <- paired_betas |>
    mutate(
      concordant = sign(beta_my) == sign(beta_jang),
      point_category = case_when(
        significant_jang & concordant  ~ "concordant_both_sig",
        significant_jang & !concordant ~ "discordant_both_sig",
        TRUE                            ~ "other"
      ),
      should_label = if (!is.null(label_genes)) label_text %in% label_genes else FALSE
    )
  
  ggplot(paired_betas, aes(x = beta_jang, y = beta_my)) +
    geom_point(
      data = filter(paired_betas, point_category == "other"),
      alpha = 0.2, color = "grey70", size = 1.8
    ) +
    geom_point(
      data = filter(paired_betas, point_category == "discordant_both_sig"),
      color = "#d32f2f", size = 3, alpha = 0.6
    ) +
    geom_point(
      data = filter(paired_betas, point_category == "concordant_both_sig"),
      color = "#1976d2", size = 3, alpha = 0.6
    ) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 0.9) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.5
    ) +
    geom_text_repel(
      data = filter(paired_betas, should_label),
      aes(label = label_text),
      size = 3.1,
      fontface = "bold",
      box.padding = 0.45,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = 25,
      force = 1.5,
      force_pull = 0.8,
      direction = "both",
      seed = 2025
    ) +
#    annotate(
#      "text", x = -1.25, y = Inf, label = cor_label,
#      hjust = 0, vjust = 1.8, size = 5
#    ) +
    labs(
      x = substitute("Adult" ~ x ~ beta, list(x = ct_labels[2])),
      y = substitute("Prenatal" ~ x ~ beta, list(x = ct_labels[1]))
    ) +
    # Centred, fixed y-range so this panel reads consistently alongside the
    # boxplot next to it, rather than auto-scaling to the data
    coord_cartesian(clip = "off", ylim = c(-2, 2)) +
    theme_minimal(base_size = 13) +
    theme(
      plot.margin = margin(25, 50, 22, 50, unit = "pt"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      title = element_blank(),
      subtitle = element_blank(),
    )
}

beta_gluA_plt <- make_jang_beta_plot(paired_betas_all, "Glu-UL", "Ext", gene_lookup,
                                     exclude_genes = exclude_genes,
                                     label_genes = "PLBD2",
                                     ct_labels = c("Glu-A", "Ext"))

### --- eQTL boxplots (Panel E left-hand side + Panel F) -----
# Shared builder for every per-gene boxplot in this figure (PLBD2 for Panel E,
# NEUROD1 + EMX1 for Panel F). Box colour matches the discovery cell type's
# colour used throughout every other manuscript figure (custom_palette).
# Style, per manuscript request:
#   - no panel border
#   - no "(log2)" in the y-axis label
#   - larger x-axis (genotype) tick text
#   - genotype and n= on two separate lines, not a plotmath subscript
make_eqtl_boxplot <- function(csv_file, gene, rsid, cell_type_label, fill_colour) {

  # Data file is named with the OLD cell-type naming (matching what's on
  # disk), same convention as every other file path in this pipeline - only
  # display labels (cell_type_label passed in) use the new Glu-A/Glu-B names.
  dat <- read_csv(csv_file, show_col_types = FALSE,
                  col_types = cols(
                    Sample = col_character(),
                    Chr    = col_character(),
                    rsID   = col_character(),
                    REF    = col_character(),  # explicit: readr misreads a single-letter
                    ALT    = col_character(),  # column like REF="T" as logical TRUE
                    GT     = col_character(),
                    .default = col_guess()
                  ))

  ref_allele <- unique(dat$REF)[1]
  alt_allele <- unique(dat$ALT)[1]

  dat <- dat |>
    mutate(
      Genotype = case_when(
        GT %in% c("0|0", "0/0") ~ paste0(ref_allele, ref_allele),
        GT %in% c("0|1", "1|0", "0/1", "1/0") ~ paste0(ref_allele, alt_allele),
        GT %in% c("1|1", "1/1") ~ paste0(alt_allele, alt_allele),
        TRUE ~ NA_character_
      ),
      Genotype = factor(
        Genotype,
        levels = c(
          paste0(ref_allele, ref_allele),
          paste0(ref_allele, alt_allele),
          paste0(alt_allele, alt_allele)
        )
      )
    )

  # Genotype counts feed the two-line x-axis labels ("TT" / "n=95") - only
  # genotypes actually present in the data get a label; scale_x_discrete's
  # default drop = TRUE already excludes any unused factor level from the axis
  geno_counts <- dat |>
    count(Genotype, name = "n") |>
    mutate(label = paste0(Genotype, "\nn=", n))

  geno_labels <- setNames(geno_counts$label, as.character(geno_counts$Genotype))

  boxplot_base <- ggplot(dat, aes(x = Genotype, y = Expression)) +
    geom_boxplot(width = 0.5, outlier.size = 2, fill = fill_colour, colour = "black") +
    scale_x_discrete(labels = geno_labels) +
    scale_y_continuous(breaks = scales::breaks_width(1)) +  # whole-number breaks only
    labs(title = cell_type_label, x = "Genotype", y = paste0(gene, " Expression")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      plot.margin = margin(25, 50, 22, 50, unit = "pt")  # matches make_jang_beta_plot's margin
      # panel.border intentionally omitted -- no black border around the plot
    )

  # rsID label centred beneath the panel, matching the original style
  ggdraw(boxplot_base) +
    draw_label(rsid, x = 0.5, y = 0.02, fontface = "bold", size = 11)
}

glu_a_colour <- custom_palette[["Glu-A"]]
glu_b_colour <- custom_palette[["Glu-B"]]

# File paths use the OLD Glu-DL naming on disk (see make_eqtl_boxplot note above)
plbd2_boxplot_plt <- make_eqtl_boxplot(
  csv_file = "reports/14MANUSCRIPT_PLOTS/eqtl_boxplots/eqtl_data_Glu-UL_rs12825284_PLBD2.csv",
  gene = "PLBD2", rsid = "rs12825284", cell_type_label = "Glu-A", fill_colour = glu_a_colour
)

neurod1_boxplot_plt <- make_eqtl_boxplot(
  csv_file = "reports/14MANUSCRIPT_PLOTS/eqtl_boxplots/eqtl_data_Glu-DL_rs115583150_NEUROD1.csv",
  gene = "NEUROD1", rsid = "rs115583150", cell_type_label = "Glu-B", fill_colour = glu_b_colour
)

### --- plot -----
# Final plot
top_row <- plot_grid(pie_chart, upset_plt, labels = c("A", "B"),
  label_size = 24, ncol = 2,rel_widths = c(1, 1.3))

# Stack heatmaps
heatmaps_stacked <- plot_grid(
  pi1_int_heatmap,
  pi1_jang_heatmap,
  ncol = 1,
  rel_heights = c(1, 0.8),   # C slightly taller
  labels = c("C", "D"),
  label_size = 24,
  align = "v"
)

# Row 1: E (NEUROD1) + F (Glu-A vs Ext beta correlation), individually labelled
ef_row <- plot_grid(neurod1_boxplot_plt, beta_gluA_plt,
                     labels = c("E", "F"), label_size = 24,
                     ncol = 2, rel_widths = c(0.9, 1.1))

# Row 2: G (PLBD2 boxplot) centred under E and F.
# G keeps the same width as E (0.9 / 2 = 0.45 of the row), with equal blank
# space (0.275) on each side, so it sits in the middle of the row.
g_row <- plot_grid(NULL, plbd2_boxplot_plt, NULL,
                   labels = c("", "G", ""), label_size = 24,
                   ncol = 3, rel_widths = c(0.275, 0.45, 0.275))

# NOTE: rel_heights below (1, 1) gives the single G panel the same row
# height as the paired E/F row above -- reasonable starting point, but
# since G alone ends up wider than either E or F individually, it may be
# worth revisiting once rendered.
betas_stacked <- plot_grid(
  ef_row, g_row,
  ncol = 1,
  rel_heights = c(1, 1)
)

# Bottom row: heatmaps + betas
# NOTE: betas_stacked now holds two 2-panel rows instead of two single
# panels, so it needs more horizontal room than before -- rel_widths below
# is a starting guess, check the render and adjust.
bottom_row <- plot_grid(
  heatmaps_stacked,
  betas_stacked,
  ncol = 2,
  rel_widths = c(0.8, 1.2), 
  align = "h"
)
final_plt <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1.4)     
)

ggsave(
  filename = out_file,
  plot = final_plt,
  width = 20,   # widened from 16 -- betas_stacked is now a 2x2 grid of panels, not 2 stacked
  height = 14,
  units = "in",
  device = "pdf",
  dpi = 300
)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
