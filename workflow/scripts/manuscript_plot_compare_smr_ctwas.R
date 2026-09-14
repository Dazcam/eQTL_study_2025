# manuscript_plot_compare_smr_ctwas.R
#
# Compare significant genes between SMR and cTWAS across all 5 GWAS traits in
# a single run. Produces:
#   - one Venn diagram plot for SCZ alone
#   - one combined 2x2 grid Venn diagram plot for BPD, MDD, ADHD, OCD (in
#     that order)
#   - one gene-membership TSV per trait (Ensembl ID + symbol + region),
#     unchanged from the per-gwas version
#
# Runs either as a Snakemake `script:` (uses the snakemake@ object) or
# standalone from the command line / RStudio (edit the ELSE block below).
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
# Fixed paths/settings for this script, hardcoded directly rather than routed
# through Snakemake params: this rule is a pipeline dead-end (final
# manuscript figures - nothing downstream reads these files) and none of
# these values vary between runs, so there's no need for config-driven
# flexibility here. Edit these directly if a path changes.
traits           <- c("scz", "bpd", "mdd", "adhd", "ocd")
smr_tbl_dir      <- "../results/13MANUSCRIPT_PLOTS_TABLES/tables/"
ctwas_in_dir     <- "../results/12CTWAS/"
weights_dir      <- "../results/12CTWAS/weights"
cell_types       <- c(
  "Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG",
  "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
  "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
  "GABA-0", "GABA-1", "GABA-2",
  "NPC-0", "NPC-1", "NPC-2",
  "NPC-to-Glu-DL-Q4-Bin1", "NPC-to-Glu-DL-Q4-Bin2", "NPC-to-Glu-DL-Q4-Bin3", "NPC-to-Glu-DL-Q4-Bin4",
  "NPC-to-Glu-UL-Q4-Bin1", "NPC-to-Glu-UL-Q4-Bin2", "NPC-to-Glu-UL-Q4-Bin3", "NPC-to-Glu-UL-Q4-Bin4"
)
gene_lookup_file <- "../resources/sheets/gene_lookup_hg38.tsv"
genes_out_tmpl   <- "../results/13MANUSCRIPT_PLOTS_TABLES/tables/venn_{gwas}_genes.tsv" # untracked by Snakemake, written silently
out_plot_scz     <- "../results/13MANUSCRIPT_PLOTS_TABLES/plots/venn_scz.png"           # untracked by Snakemake, written silently
ctwas_pip_thresh <- 0.8
ctwas_use_bf     <- TRUE    # BF-corrected significance is the real criterion; PIP>0.8 alone was the report's bug
ctwas_bf_thresh  <- 0.05
use_symbols      <- TRUE
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

message("=== compare_smr_ctwas: ", paste(traits, collapse = ", "), " ===")

# ---- 1. SMR significant genes (union across cell types) ------------------
# Reuses the tables smr_report.Rmd already writes: {tbl_dir}/{gwas}_smr.tsv
get_smr_genes <- function(gwas, smr_tbl_dir) {
  f <- file.path(smr_tbl_dir, paste0(gwas, "_smr.tsv"))
  if (!file.exists(f)) {
    message("No SMR sig table for ", gwas, " (", f, ") - treating as empty set")
    return(character(0))
  }
  message("Reading SMR table: ", f)
  tbl <- tryCatch(
    suppressMessages(read_tsv(f, show_col_types = FALSE)),
    error = function(e) stop("Failed to read SMR table (", f, "): ", conditionMessage(e))
  )
  if (!"Ensembl ID" %in% names(tbl)) {
    stop("SMR table ", f, " has no 'Ensembl ID' column - columns present: ",
         paste(names(tbl), collapse = ", "))
  }
  unique(tbl[["Ensembl ID"]])
}

# ---- 2. cTWAS significant genes (union across cell types) ----------------
# Same filter ctwas_report.Rmd uses to build its displayed table:
# group != 'SNP' & susie_pip > 0.8 & in a credible set (!is.na(cs))
`%||%` <- function(a, b) if (is.null(a)) b else a

# Correction factor for the optional BF-corrected significance path, matching
# the "Correction" tab in ctwas_report.Rmd: the number of non-empty
# .wgt.RDat FUSION weight files for a cell type is the count of genes that
# were actually testable there (hsq_p < 0.01) - i.e. the number of tests,
# NOT the number of genes that came out significant.
get_bf_factor <- function(ct, weights_dir) {
  ct_dir <- file.path(weights_dir, ct)
  if (!dir.exists(ct_dir)) {
    message("  [", ct, "] weights dir not found for BF factor: ", ct_dir)
    return(NA_integer_)
  }
  files <- list.files(ct_dir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
  files <- files[file.info(files)$size > 0]
  n <- length(files)
  if (n == 0) {
    message("  [", ct, "] zero non-empty .wgt.RDat files found in ", ct_dir)
    return(NA_integer_)
  }
  n
}

get_ctwas_genes <- function(gwas, ctwas_in_dir, weights_dir, cell_types,
                             pip_thresh, use_bf, bf_thresh) {
  bf_skipped_cts <- character(0)  # cell types skipped for missing/empty weights

  genes <- map(cell_types, function(ct) {
    f <- file.path(ctwas_in_dir, "output", paste0("ctwas_", ct, "_", gwas, "_ctwas.rds"))
    if (!file.exists(f)) {
      message("  [", ct, "] file not found: ", f)
      return(character(0))
    }
    message("  [", ct, "] reading ", f)
    res <- tryCatch(
      read_rds(f),
      error = function(e) stop("Failed to read RDS for cell type '", ct, "' (", f, "): ",
                                conditionMessage(e))
    )
    if (!is.null(res$status) && res$status == "skipped") {
      message("  [", ct, "] skipped: ", res$reason %||% "no reason given")
      return(character(0))
    }
    if (is.null(res$finemap_res) || nrow(res$finemap_res) == 0) {
      message("  [", ct, "] no finemap_res rows")
      return(character(0))
    }

    fm <- tryCatch({
      as_tibble(res$finemap_res) |>
        filter(group != "SNP", susie_pip > pip_thresh, !is.na(cs)) |>
        separate(id, into = c("gene", "id"), sep = "\\|") |>
        separate(id, into = c("cell_type", "id"), sep = "_")
    }, error = function(e) {
      stop("Failed to parse finemap_res for cell type '", ct, "' (", f, "): ",
           conditionMessage(e))
    })

    if (use_bf) {
      if (nrow(fm) == 0) {
        # Nothing survived the PIP/credible-set filter for this cell type -
        # no point computing a correction factor, this cell type just
        # contributes zero genes.
        return(character(0))
      }
      bf_factor <- get_bf_factor(ct, weights_dir)
      if (is.na(bf_factor)) {
        # A cell type can legitimately have no FUSION weights at all (too
        # few cells to model cis-heritability, etc.) - especially likely
        # among the finer subcluster/pseudotime-bin cell types. That's a
        # normal outcome, not a pipeline error, so we skip this cell type
        # (it contributes zero genes) rather than aborting the whole GWAS
        # run. The message() above already recorded why, for the log.
        message("  [", ct, "] skipping BF correction (see reason above) - contributes 0 genes")
        bf_skipped_cts <<- c(bf_skipped_cts, ct)
        return(character(0))
      }
      fm <- fm |>
        mutate(pval = 2 * pnorm(-abs(z)),
               pval_bf = pmin(pval * bf_factor, 1)) |>
        filter(pval_bf < bf_thresh)
    }
    message("  [", ct, "] ", length(unique(fm$gene)), " significant genes")
    unique(fm$gene)
  })

  if (use_bf && length(bf_skipped_cts) == length(cell_types) && length(cell_types) > 0) {
    message("WARNING: ALL ", length(cell_types), " cell types were skipped for BF ",
            "correction (missing/empty weights) for ", gwas, " - this looks like a ",
            "'weights_dir' path problem rather than normal per-cell-type variation. ",
            "Check the path: ", weights_dir)
  } else if (use_bf && length(bf_skipped_cts) > 0) {
    message("  BF correction skipped for ", length(bf_skipped_cts), "/", length(cell_types),
            " cell types (no weights found): ", paste(bf_skipped_cts, collapse = ", "))
  }

  unique(unlist(genes))
}

# ---- 3. Ensembl -> symbol lookup, falling back to the Ensembl ID ----------
to_display <- function(ids, use_symbols, gene_lookup_file) {
  if (length(ids) == 0) return(character(0))
  if (!use_symbols || is.null(gene_lookup_file) || !file.exists(gene_lookup_file)) {
    if (use_symbols) message("Gene lookup file not found (", gene_lookup_file,
                              ") - falling back to Ensembl IDs")
    return(ids)
  }
  lookup <- suppressMessages(read_tsv(gene_lookup_file, show_col_types = FALSE)) |>
    select(ensembl_gene_id, external_gene_name) |>
    distinct()
  tibble(ensembl_gene_id = ids) |>
    left_join(lookup, by = "ensembl_gene_id") |>
    mutate(external_gene_name = case_when(
      is.na(external_gene_name) | external_gene_name == "" ~ ensembl_gene_id,
      TRUE ~ external_gene_name
    )) |>
    pull(external_gene_name)
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

# ---- 5. Per-trait pipeline: gene sets -> Venn plot + membership table -----
# Encapsulates the original single-gwas script's full flow (steps 1-6 there)
# as a function, so it can be run once per trait and the resulting plot
# object either saved standalone (SCZ) or combined into a grid
# (BPD/MDD/ADHD/OCD) without duplicating the gene-retrieval/display logic.
run_gwas_venn <- function(gwas, smr_tbl_dir, ctwas_in_dir, weights_dir, cell_types,
                          ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh,
                          use_symbols, gene_lookup_file, text_size,
                          single_col_threshold) {

  message("--- ", gwas, " ---")
  message("Getting SMR genes...")
  smr_genes   <- get_smr_genes(gwas, smr_tbl_dir)
  message("Getting cTWAS genes...")
  ctwas_genes <- get_ctwas_genes(gwas, ctwas_in_dir, weights_dir, cell_types,
                                  ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh)
  message("SMR genes: ", length(smr_genes), " | cTWAS genes: ", length(ctwas_genes))

  both_disp       <- to_display(intersect(smr_genes, ctwas_genes), use_symbols, gene_lookup_file)
  smr_only_disp   <- to_display(setdiff(smr_genes, ctwas_genes), use_symbols, gene_lookup_file)
  ctwas_only_disp <- to_display(setdiff(ctwas_genes, smr_genes), use_symbols, gene_lookup_file)

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

  # Gene membership table (Ensembl ID + symbol, always unabridged)
  membership <- tibble(gene = union(smr_genes, ctwas_genes)) |>
    mutate(
      symbol   = to_display(gene, use_symbols, gene_lookup_file),
      in_SMR   = gene %in% smr_genes,
      in_cTWAS = gene %in% ctwas_genes,
      category = case_when(
        in_SMR  & in_cTWAS  ~ "both",
        in_SMR  & !in_cTWAS ~ "SMR_only",
        !in_SMR & in_cTWAS  ~ "cTWAS_only"
      )
    ) |>
    arrange(category, gene)

  message(sprintf(
    "%s: SMR=%d, cTWAS=%d, both=%d, SMR_only=%d, cTWAS_only=%d",
    gwas, length(smr_genes), length(ctwas_genes),
    sum(membership$category == "both"),
    sum(membership$category == "SMR_only"),
    sum(membership$category == "cTWAS_only")
  ))

  list(plot = p, membership = membership)
}

# ---- 6. Run every trait, write its membership table -----------------------
# SCZ keeps the single-panel text_size; the other 4 (destined for the 2x2
# grid) use grid_text_size instead.
results <- map(traits, function(gwas) {
  run_gwas_venn(gwas, smr_tbl_dir, ctwas_in_dir, weights_dir, cell_types,
               ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh,
               use_symbols, gene_lookup_file,
               text_size = if (gwas == "scz") text_size else grid_text_size,
               single_col_threshold = single_col_threshold)
})
names(results) <- traits

for (gwas in traits) {
  genes_out <- str_replace(genes_out_tmpl, fixed("{gwas}"), gwas)
  write_tsv(results[[gwas]]$membership, genes_out)
  message("Wrote membership table: ", genes_out)
}

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
