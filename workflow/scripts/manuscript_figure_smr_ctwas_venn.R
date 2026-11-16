# manuscript_figure_smr_ctwas_venn.R
#
# Compare significant genes between SMR and cTWAS across all 5 GWAS traits in
# a single run. Produces:
#   - one Venn diagram plot for SCZ alone
#   - one combined 2x2 grid Venn diagram plot for BPD, MDD, ADHD, OCD (in
#     that order)
#   (gene sets read from the 16 smr_ctwas_venn_tbl membership tables)
#
# Runs as a Snakemake `script:`; paths and cell types come from the rule.
#
# CHANGE FROM PREVIOUS VERSION: no longer driven by a per-gwas Snakemake
# wildcard - all 5 traits are processed in this single script invocation,
# looping internally rather than being invoked once per trait.

# ---- Logging ---------------------------------------------------------------
# Snakemake does NOT auto-redirect console output to `log:` for `script:`
# rules (that only happens automatically for `shell:` rules). We sink()
# explicitly, before anything else runs, so every message/warning/error
# actually lands in the log file instead of vanishing.
if (exists("snakemake")) {
  log_conn <- file(snakemake@log[[1]], open = "wt")
  sink(log_conn, type = "output")
  sink(log_conn, type = "message")
}

suppressMessages({
  library(tidyverse)
  library(ggvenn)
  library(cowplot)
})

# ---- Parameters ------------------------------------------------------------
# Traits and the per-trait SMR vs cTWAS gene tables (16 smr_ctwas_venn_tbl)
# come from the rule; the plotting settings are fixed here.
traits           <- unlist(snakemake@params[["traits"]])
genes_tbls       <- unlist(snakemake@input[["genes_tbls"]])
out_plot_scz     <- snakemake@output[["plot_scz"]]
# single-panel sizing (SCZ) - same as the original single-gwas script
plot_width       <- 16
plot_height      <- 14
text_size        <- 3.6
# combined 2x2 grid sizing (BPD, MDD, ADHD, OCD): each panel still shows its
# full unabridged gene list (max_elements = Inf), so the combined canvas
# needs to be substantially larger than a single panel, not just 2x its
# dimensions.
grid_plot_width  <- 32
grid_plot_height <- 28
grid_text_size   <- 4.2
# Regions with this many genes or fewer render as a single column instead of
# two - without this, a sparse region (e.g. ADHD/OCD's 4-5 genes) still gets
# packed into 2 columns and ends up a tiny, hard-to-read block lost in the
# middle of its circle, regardless of text_size. Tuned against toy_venn.R.
single_col_threshold <- 8

# Only the final grid plot is tracked as a Snakemake output (it's written
# last, so its existence implies everything above it succeeded); everywhere
# else this script runs the same whether invoked via Snakemake or standalone.
out_plot_grid <- if (exists("snakemake")) {
  snakemake@output[["plot_grid"]]
} else {
  "venn_grid_bpd_mdd_adhd_ocd.png"
}

message("=== smr_ctwas_venn_fig: ", paste(traits, collapse = ", "), " ===")

# ---- 1. Gene membership tables (16 smr_ctwas_venn_tbl), one per trait ------
# Columns: gene, symbol, in_SMR, in_cTWAS, category (both / SMR_only / cTWAS_only)
read_membership <- function(gwas, genes_tbls) {
  f <- genes_tbls[str_detect(basename(genes_tbls), paste0("_", gwas, "_"))]
  if (length(f) != 1) stop("Expected one gene table for ", gwas, ", found ", length(f))
  message("Reading gene table: ", f)
  read_tsv(f, show_col_types = FALSE,
           col_types = cols(.default = col_character(),
                            in_SMR = col_logical(), in_cTWAS = col_logical()))
}

# ---- 4. Sort (symbols alphabetically, then Ensembl IDs) and pack into ------
#         two-column "row" strings for the Venn diagram
sort_symbols_first <- function(display_names) {
  if (length(display_names) == 0) return(character(0))
  is_ensg <- str_detect(display_names, "^ENSG\\d+$")
  symbols <- sort(display_names[!is_ensg])
  # Standard Ensembl IDs are zero-padded to a fixed digit count, so a plain
  # string sort already gives numeric order - no need to parse digits out.
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

# Sparse regions (<= single_col_threshold genes) render as a single column;
# denser regions keep the two-column packing above. See toy_venn.R for the
# comparison that motivated this - text_size alone can't fix a sparse region
# looking lost, since 2-column packing centres a handful of genes into a
# tiny block regardless of font size.
make_rows <- function(display_names, single_col_threshold) {
  ordered <- sort_symbols_first(display_names)
  if (length(ordered) <= single_col_threshold) {
    ordered
  } else {
    make_two_col_rows(display_names)
  }
}

# ---- 5. Per-trait Venn plot ------------------------------------------------
run_gwas_venn <- function(gwas, genes_tbls, text_size, single_col_threshold) {

  message("--- ", gwas, " ---")
  membership <- read_membership(gwas, genes_tbls)

  both_disp       <- membership$symbol[membership$category == "both"]
  smr_only_disp   <- membership$symbol[membership$category == "SMR_only"]
  ctwas_only_disp <- membership$symbol[membership$category == "cTWAS_only"]

  rows_both       <- make_rows(both_disp, single_col_threshold)
  rows_smr_only   <- make_rows(smr_only_disp, single_col_threshold)
  rows_ctwas_only <- make_rows(ctwas_only_disp, single_col_threshold)

  # Feed ggvenn pre-built row strings (single- or two-column, per region) as
  # the set elements. Since rows_both is identical in both vectors, ggvenn's
  # internal intersect() still correctly resolves it as the shared region.
  gene_sets <- list(
    SMR   = c(rows_smr_only, rows_both),
    cTWAS = c(rows_ctwas_only, rows_both)
  )

  # Venn diagram: white background, centred acronym-only title, every gene
  # name shown (no truncation), single- or two-column per region depending
  # on how many genes it holds
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
      # text = element_text(family = "mono") below is needed to keep the
      # two-column gene-list padding aligned, but plot.title must NOT
      # inherit it (that was the earlier bug making the title look
      # different from the SMR/cTWAS set-name labels, which ggvenn draws
      # separately and were never affected by this theme in the first
      # place) - explicitly reset plot.title back to a normal sans family.
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

# ---- 6. Run every trait --------------------------------------------------------
# SCZ keeps the single-panel text_size; the other 4 (destined for the 2x2
# grid) use grid_text_size instead.
results <- map(traits, function(gwas) {
  run_gwas_venn(gwas, genes_tbls,
               text_size = if (gwas == "scz") text_size else grid_text_size,
               single_col_threshold = single_col_threshold)
})
names(results) <- traits

# ---- 7. SCZ: standalone single-panel plot ---------------------------------
if (!("scz" %in% traits)) {
  stop("'scz' must be included in `traits` to produce the standalone SCZ plot.")
}
message("Saving SCZ plot: ", out_plot_scz)
ggsave(out_plot_scz, results[["scz"]]$plot,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")

# ---- 8. BPD, MDD, ADHD, OCD: combined 2x2 grid (in that order) -----------
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

message("Done: ", paste(traits, collapse = ", "))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
