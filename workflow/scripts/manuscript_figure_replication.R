#--------------------------------------------------------------------------------------
#
#    eQTL sharing between prenatal and adult cell populations (Fig 3)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: replication_fig
#            Upstream:   05 tensorqtl_perm (final norm method, genPC, expPCs); 07 internal pi1;
#                        15 pi1_jang and beta correlation (Jang et al. adult eQTL);
#                        eqtl_boxplots_py (17)
#            Downstream: Fig 3
#
# Purpose:   Panels for the 7 L1 cell types (manuscript labels):
#            A  Pie chart: eGenes (FDR < 0.05) by the number of L1 cell types
#               they are found in
#            B  UpSet plot of eGenes shared between L1 cell types (the 20 most
#               frequent intersections)
#            C  pi1 heatmap, prenatal discovery x prenatal replication (forward
#               estimates; diagonal set to 1)
#            D  pi1 heatmap, prenatal discovery (Glu-A, Glu-B, GABA, NPC) x the
#               7 Jang adult cell types
#            E  eQTL boxplot (config pair fig3e: NEUROD1, rs115583150, Glu-B)
#            F  Glu-A vs adult Ext eQTL slopes for Glu-A FDR < 0.05 pairs; blue /
#               red: also significant in Jang (lead SNP, qval < 0.05) with the
#               same / opposite sign; PLBD2 labelled
#            G  eQTL boxplot (config pair fig3g: PLBD2, rs12825284, Glu-A)
#
# Inputs:    in_dir         tensorQTL permutation results directory
#            internal_dir   07 internal pi1 results
#            jang_pi1_dir   15 Jang pi1 results
#            beta_file      15 beta correlation RDS ($paired_betas)
#            boxplot_e/g    Boxplot CSVs from eqtl_boxplots_py
#            gene_lookup    Ensembl ID to gene symbol
#
# Outputs:   out_file       Fig 3 (PDF)
#
# Notes:     Cell types, expPC values, genPC, norm method and the Jang adult cell
#            types come from config; the panel D discovery rows are set here.
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

message("\n\nPlotting eQTL sharing figure ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)
library(UpSetR)
library(grid)

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
in_dir            <- snakemake@params[['in_dir']]
internal_dir      <- snakemake@params[['internal_dir']]
jang_pi1_dir      <- snakemake@params[['jang_pi1_dir']]
gene_lookup_file  <- snakemake@params[['gene_lookup']]
manuscript_labels <- snakemake@params[["labels"]]
boxplot_pairs     <- snakemake@params[['boxplot_pairs']]
jang_beta_file    <- snakemake@input[['beta_file']]
boxplot_e_file    <- snakemake@input[['boxplot_e']]
boxplot_g_file    <- snakemake@input[['boxplot_g']]
expPC_map         <- unlist(snakemake@params[["exp_pc_map"]])
genPC             <- snakemake@params[["geno_pc"]]
norm_method       <- snakemake@params[["norm_method"]]
eqtl_fdr          <- as.numeric(snakemake@params[["eqtl_fdr"]])
out_file          <- snakemake@output[[1]]

# L1 cell types (pipeline labels, config order: these build the input file
# paths and set the heatmap order)
l1_types   <- unlist(snakemake@params[["l1_types"]])
cell_types <- unlist(snakemake@params[["cell_types"]])
cell_types <- cell_types[cell_types %in% l1_types]

# Panel D: prenatal discovery rows (NPC and developing neurons) and the Jang
# adult L1 replication columns (config jang_l1_types)
jang_L1_cell_types    <- c("Glu-UL", "Glu-DL", "GABA", "NPC")
jang_adult_cell_types <- unlist(snakemake@params[["jang_types"]])

cell_types_new <- relabel_cell_types(cell_types, manuscript_labels)

# Colour palette (manuscript labels)
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

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_dir", "internal_dir", "jang_pi1_dir", "gene_lookup_file", "labels",
               "jang_beta_file", "boxplot_e_file", "boxplot_g_file", "boxplot_pairs",
               "cell_types", "jang_L1_cell_types", "jang_adult_cell_types",
               "expPC_map", "norm_method", "genPC", "eqtl_fdr", "out_file"),
  value    = c(in_dir, internal_dir, jang_pi1_dir, gene_lookup_file, paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               jang_beta_file, boxplot_e_file, boxplot_g_file,
               paste(names(boxplot_pairs), map_chr(boxplot_pairs, ~ paste(unlist(.x), collapse = " ")),
                     sep = ": ", collapse = "; "),
               paste(cell_types, collapse = ", "), paste(jang_L1_cell_types, collapse = ", "),
               paste(jang_adult_cell_types, collapse = ", "),
               paste(names(expPC_map[cell_types]), expPC_map[cell_types], sep = ":", collapse = ", "),
               norm_method, genPC, eqtl_fdr, out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  dplyr::select(gene = ensembl_gene_id, symbol = external_gene_name)

##  Panel A: pie chart of eGenes by number of L1 cell types  -------------------------
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
      filter(qval < eqtl_fdr)
    
    if (nrow(dat) > 0) {
      eqtl_list[[cell_type]] <- dat %>% mutate(cell_type = cell_type)
    }
  }
}

eqtls_all <- bind_rows(eqtl_list)

# Manuscript labels (needed for the panel B set names)
gene_cell <- eqtls_all %>%
  distinct(phenotype_id, cell_type) %>%
  mutate(cell_type = relabel_cell_types(cell_type, manuscript_labels))

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

# Pie drawn as Cartesian polygons (not coord_polar) so label positions and
# connector lines share one coordinate system. Slices start at 12 o'clock
# and run clockwise: theta = pi/2 - 2*pi*fraction.
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

# Inner labels (large slices) at r = 0.6
inner_labels <- pie_dat %>%
  filter(large_slice) %>%
  mutate(theta = pi / 2 - 2 * pi * midpoint,
         x = 0.6 * cos(theta),
         y = 0.6 * sin(theta))

# Outer labels (small slices): fixed fan-out above the pie, left to right in
# slice order, so the layout does not depend on the rendered panel size
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

##  Panel B: UpSet plot  --------------------------------------------------------------
gene_by_cell <- gene_cell %>%
  pivot_wider(names_from = cell_type, values_from = cell_type,
              values_fill = 0, values_fn = function(x) 1) %>%
  column_to_rownames("phenotype_id")

# Fixed set order (keep.order = TRUE) so sets.bar.color matches the sets.
# UpSetR draws the first set at the bottom, so config order is reversed.
set_order <- rev(cell_types_new)

# UpSetR is base graphics: render to a temporary PNG and re-import as a grob.
# Its margins are fixed in absolute size, so a wide canvas leaves less
# whitespace.
tmp_upset <- tempfile(fileext = ".png")
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

# Re-import as a grob for cowplot
img <- png::readPNG(tmp_upset)

upset_grob <- grid::rasterGrob(
  img,
  interpolate = TRUE
)

upset_plt <- ggplotify::as.ggplot(upset_grob) +
  theme(
    plot.margin = margin(t = 40, r = 0, b = 40, l = 0, unit = "pt")
  )

##  Panel C: prenatal pi1 heatmap  ----------------------------------------------------
read_pi1_results <- function(ct, ref_ct) {
  
  expPC <- expPC_map[[ct]]
  ref_expPC <- expPC_map[[ref_ct]]
  
  file_path <- paste0(internal_dir, ct, "_vs_", ref_ct, "_", norm_method, "_genPC_", genPC, "_expPC_",
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

# All discovery x replication pairs (pipeline labels: these build the file paths)
combinations <- expand.grid(cell_type = cell_types, ref_cell_type = cell_types, stringsAsFactors = FALSE)

# Read all pi1 results
pi1_result_tbl <- map2_dfr(combinations$cell_type, combinations$ref_cell_type, read_pi1_results)

# Forward pi1 for every pair; 1 on the diagonal
pi1_square_tbl <- pi1_result_tbl %>%
  mutate(
    pi1_final = ifelse(cell_type == ref_cell_type, 1.0, pi1_forward)
  ) %>%
  dplyr::select(ref = cell_type, query = ref_cell_type, pi1 = pi1_final) %>%
  mutate(
    ref = relabel_cell_types(ref, manuscript_labels),
    query = relabel_cell_types(query, manuscript_labels),
    # rev() ensures the first cell type is at the top (standard matrix view)
    query = factor(query, levels = rev(cell_types_new)), 
    ref = factor(ref, levels = cell_types_new)
  )

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

pi1_int_heatmap <- plot_int_heatmap(pi1_square_tbl)


##  Panel D: prenatal vs adult (Jang) pi1 heatmap  ------------------------------------

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
  dplyr::select(cell_type, ref_cell_type, pi1 = pi1_forward) %>%
  mutate(
    cell_type = relabel_cell_types(cell_type, manuscript_labels),
    # rev() so Glu-A appears at the top (standard matrix view)
    cell_type = factor(cell_type, levels = rev(relabel_cell_types(jang_L1_cell_types, manuscript_labels))),
    ref_cell_type = factor(ref_cell_type, levels = jang_adult_cell_types)
  )

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


##  Panel F: Glu-A vs adult Ext eQTL slopes  ------------------------------------------
jang_beta_list <- read_rds(jang_beta_file)
paired_betas_all <- jang_beta_list$paired_betas

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
  
  # Point colour, as in repl_sc_beta_correlation.R / eqtl_sharing_report.Rmd:
  # pairs also significant in Jang (lead SNP, qval < 0.05) are blue (same
  # sign) or red (opposite sign); all other pairs are grey
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
    labs(
      x = substitute("Adult" ~ x ~ beta, list(x = ct_labels[2])),
      y = substitute("Prenatal" ~ x ~ beta, list(x = ct_labels[1]))
    ) +
    # Fixed y-range, centred on 0
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

##  Panels E and G: eQTL boxplots  ----------------------------------------------------
# Box colour is the discovery cell type's colour; x labels give the genotype
# and its n on two lines
make_eqtl_boxplot <- function(csv_file, gene, rsid, cell_type_label, fill_colour) {

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

  # Genotype counts for the two-line x-axis labels ("TT" / "n=95"); only
  # genotypes present in the data are shown
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
      plot.margin = margin(25, 50, 22, 50, unit = "pt")  # as make_jang_beta_plot
    )

  # rsID label centred beneath the panel
  ggdraw(boxplot_base) +
    draw_label(rsid, x = 0.5, y = 0.02, fontface = "bold", size = 11)
}

glu_a_colour <- custom_palette[["Glu-A"]]
glu_b_colour <- custom_palette[["Glu-B"]]

plbd2_boxplot_plt <- make_eqtl_boxplot(
  csv_file = boxplot_g_file,
  gene = boxplot_pairs$fig3g$gene, rsid = boxplot_pairs$fig3g$snp,
  cell_type_label = "Glu-A", fill_colour = glu_a_colour
)

neurod1_boxplot_plt <- make_eqtl_boxplot(
  csv_file = boxplot_e_file,
  gene = boxplot_pairs$fig3e$gene, rsid = boxplot_pairs$fig3e$snp,
  cell_type_label = "Glu-B", fill_colour = glu_b_colour
)

##  Assemble and save  ----------------------------------------------------------------
# A and B
top_row <- plot_grid(pie_chart, upset_plt, labels = c("A", "B"),
  label_size = 24, ncol = 2,rel_widths = c(1, 1.3))

# C above D
heatmaps_stacked <- plot_grid(
  pi1_int_heatmap,
  pi1_jang_heatmap,
  ncol = 1,
  rel_heights = c(1, 0.8),
  labels = c("C", "D"),
  label_size = 24,
  align = "v"
)

# E and F side by side
ef_row <- plot_grid(neurod1_boxplot_plt, beta_gluA_plt,
                     labels = c("E", "F"), label_size = 24,
                     ncol = 2, rel_widths = c(0.9, 1.1))

# G centred below, the same width as E (0.45 of the row)
g_row <- plot_grid(NULL, plbd2_boxplot_plt, NULL,
                   labels = c("", "G", ""), label_size = 24,
                   ncol = 3, rel_widths = c(0.275, 0.45, 0.275))

betas_stacked <- plot_grid(
  ef_row, g_row,
  ncol = 1,
  rel_heights = c(1, 1)
)

# C/D beside E/F/G
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

message("Writing: ", out_file)
ggsave(
  filename = out_file,
  plot = final_plt,
  width = 20,
  height = 14,
  units = "in",
  device = "pdf",
  dpi = 300
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
