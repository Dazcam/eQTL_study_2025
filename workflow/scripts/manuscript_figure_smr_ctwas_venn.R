#--------------------------------------------------------------------------------------
#
#    SMR vs cTWAS gene Venn diagrams per GWAS trait (Fig 6)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: smr_ctwas_venn_fig
#            Upstream:   16 smr_ctwas_venn_tbl (gene membership table per trait)
#            Downstream: Fig 6
#
# Purpose:   For each GWAS trait, a Venn diagram of SMR vs cTWAS genes that
#            lists every gene in its region (symbols A-Z, then Ensembl IDs;
#            one column for regions with <= single_col_threshold genes, two
#            otherwise):
#            - SCZ as a single panel
#            - BPD, MDD, ADHD and OCD as a 2 x 2 grid, in that order
#
# Inputs:    genes_tbls     venn_{gwas}_genes.tsv, one per trait
#            traits         config gwas keys (must include scz, bpd, mdd, adhd, ocd)
#
# Outputs:   plot_scz       SCZ Venn (PNG, 300 dpi)
#            plot_grid      BPD / MDD / ADHD / OCD grid (PNG, 300 dpi)
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

message("\n\nPlotting SMR vs cTWAS Venn diagrams ...")

##  Load packages and variables  ------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(ggvenn)
  library(cowplot)
})

# Inputs and outputs
traits        <- unlist(snakemake@params[["traits"]])
genes_tbls    <- unlist(snakemake@input[["genes_tbls"]])
out_plot_scz  <- snakemake@output[["plot_scz"]]
out_plot_grid <- snakemake@output[["plot_grid"]]

# Plot settings: SCZ single panel
plot_width       <- 16
plot_height      <- 14
text_size        <- 3.6
# 2 x 2 grid: every panel lists all its genes, so the canvas is larger than
# 2x a single panel
grid_plot_width  <- 32
grid_plot_height <- 28
grid_text_size   <- 4.2
# Regions with this many genes or fewer are drawn as one column, not two
single_col_threshold <- 8

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("traits", "genes_tbls", "out_plot_scz", "out_plot_grid",
               "plot_width x height", "text_size", "grid_plot_width x height",
               "grid_text_size", "single_col_threshold"),
  value    = c(paste(traits, collapse = ", "), paste(genes_tbls, collapse = ", "),
               out_plot_scz, out_plot_grid, paste(plot_width, plot_height, sep = " x "),
               text_size, paste(grid_plot_width, grid_plot_height, sep = " x "),
               grid_text_size, single_col_threshold)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  ------------------------------------------------------------------------
# Gene membership table for one trait (16 smr_ctwas_venn_tbl). Columns: gene, symbol, in_SMR, in_cTWAS, category (both / SMR_only / cTWAS_only)
read_membership <- function(gwas, genes_tbls) {
  f <- genes_tbls[str_detect(basename(genes_tbls), paste0("_", gwas, "_"))]
  if (length(f) != 1) stop("Expected one gene table for ", gwas, ", found ", length(f))
  message("Reading gene table: ", f)
  read_tsv(f, show_col_types = FALSE,
           col_types = cols(.default = col_character(),
                            in_SMR = col_logical(), in_cTWAS = col_logical()))
}

# Sort (symbols alphabetically, then Ensembl IDs) and pack into two-column
# row strings for the Venn diagram
sort_symbols_first <- function(display_names) {
  if (length(display_names) == 0) return(character(0))
  is_ensg <- str_detect(display_names, "^ENSG\\d+$")
  symbols <- sort(display_names[!is_ensg])
  # Ensembl IDs have a fixed digit count, so a string sort is numeric order
  ensg <- sort(display_names[is_ensg])
  c(symbols, ensg)
}

make_two_col_rows <- function(display_names) {
  ordered <- sort_symbols_first(display_names)
  n <- length(ordered)
  if (n == 0) return(character(0))
  if (n %% 2 == 1) ordered <- c(ordered, "")
  col1 <- ordered[seq(1, length(ordered), by = 2)]
  col2 <- ordered[seq(2, length(ordered), by = 2)]
  width1 <- max(nchar(col1))
  ifelse(col2 == "",
         col1,
         sprintf(paste0("%-", width1, "s   %s"), col1, col2))
}

# One column for regions with <= single_col_threshold genes, two otherwise
make_rows <- function(display_names, single_col_threshold) {
  ordered <- sort_symbols_first(display_names)
  if (length(ordered) <= single_col_threshold) {
    ordered
  } else {
    make_two_col_rows(display_names)
  }
}

# Venn plot for one trait
run_gwas_venn <- function(gwas, genes_tbls, text_size, single_col_threshold) {

  message("--- ", gwas, " ---")
  membership <- read_membership(gwas, genes_tbls)

  both_disp       <- membership$symbol[membership$category == "both"]
  smr_only_disp   <- membership$symbol[membership$category == "SMR_only"]
  ctwas_only_disp <- membership$symbol[membership$category == "cTWAS_only"]

  rows_both       <- make_rows(both_disp, single_col_threshold)
  rows_smr_only   <- make_rows(smr_only_disp, single_col_threshold)
  rows_ctwas_only <- make_rows(ctwas_only_disp, single_col_threshold)

  # Pre-built row strings are the set elements; rows_both is in both sets, so
  # ggvenn places it in the overlap
  gene_sets <- list(
    SMR   = c(rows_smr_only, rows_both),
    cTWAS = c(rows_ctwas_only, rows_both)
  )

  # White background, trait acronym as the title, every gene shown
  p <- ggvenn(
    gene_sets,
    show_elements   = TRUE,
    label_sep       = "\n",
    show_percentage = FALSE,
    fill_color      = c("#66c2a5", "#fc8d62"),
    text_size       = text_size,
    max_elements    = Inf,    # ggvenn defaults to 6 and truncates - show everything
    text_truncate   = FALSE
  ) +
    labs(title = toupper(gwas)) +
    theme(
      # Gene lists use a monospace font to keep the columns aligned; the
      # title stays sans
      plot.title       = element_text(hjust = 0.5, size = 16, face = "bold", family = "sans"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      text             = element_text(family = "mono")  # keeps column padding aligned
    )

  message(sprintf(
    "%s: both=%d, SMR_only=%d, cTWAS_only=%d",
    gwas, length(both_disp), length(smr_only_disp), length(ctwas_only_disp)
  ))

  list(plot = p)
}

##  Venn diagrams  --------------------------------------------------------------------
# SCZ uses text_size; the grid traits use grid_text_size
results <- map(traits, function(gwas) {
  run_gwas_venn(gwas, genes_tbls,
               text_size = if (gwas == "scz") text_size else grid_text_size,
               single_col_threshold = single_col_threshold)
})
names(results) <- traits

# SCZ: single panel
if (!("scz" %in% traits)) {
  stop("'scz' must be included in `traits` to produce the standalone SCZ plot.")
}
message("Saving SCZ plot: ", out_plot_scz)
ggsave(out_plot_scz, results[["scz"]]$plot,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")

# BPD, MDD, ADHD, OCD: 2 x 2 grid, in that order
grid_traits <- c("bpd", "mdd", "adhd", "ocd")
missing_traits <- setdiff(grid_traits, traits)
if (length(missing_traits) > 0) {
  stop("Traits required for the 2x2 grid are missing from `traits`: ",
       paste(missing_traits, collapse = ", "))
}

grid_plt <- plot_grid(
  plotlist = map(grid_traits, ~ results[[.x]]$plot),
  ncol = 2, nrow = 2
)

message("Saving combined BPD/MDD/ADHD/OCD grid plot: ", out_plot_grid)
ggsave(out_plot_grid, grid_plt,
       width = grid_plot_width, height = grid_plot_height, dpi = 300, bg = "white")

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
