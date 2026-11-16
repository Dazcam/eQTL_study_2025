#--------------------------------------------------------------------------------------
#
#    Functions for the tensorQTL pipeline report
#
#--------------------------------------------------------------------------------------
#
# Functions: tqtl_dir()                  results directory for one tensorQTL run
#            parse_perm_log()            summary numbers from a cis permutation log
#            summarise_logs()            parse_perm_log() over cell types, expPCs, norms
#            read_sig_eqtls()            significant lead eQTL per cell type
#            read_egenes()               significant eGenes for one cell type
#            read_indep()                conditionally independent eQTL per cell type
#            gt_style()                  house style for the report's gt tables
#            plot_sensitivity()          eGenes by expPC, one panel per norm method
#            sensitivity_best_gt()       best expPC / norm method per cell type
#            sensitivity_all_gt()        eGenes for every expPC / norm method
#            plot_density()              TSS distance density, one cell type
#            plot_density_combined()     TSS distance density, all cell types
#            plot_egene_counts()         eGenes per cell type (bar chart)
#            indep_counts_gt()           independent eQTL counts per rank
#            indep_tss_gt()              median TSS distance per rank
#            inputs_gt()                 genes, samples, variants of the runs used
#
# Packages:  the caller must have dplyr, tidyr, purrr, readr, stringr, tibble,
#            ggplot2, cowplot and gt installed. The read and summarise functions
#            need only dplyr, purrr, readr, stringr and tibble.
#
# Used in:   tensorqtl_report.Rmd   all
#
#--------------------------------------------------------------------------------------

# Directory of one tensorQTL run, e.g. <base_dir>/NPC_quantile_genPC_4_expPC_30/.
# base_dir must end in "/".
tqtl_dir <- function(base_dir, cell_type, norm_method, geno_pc, exp_pc) {
  paste0(base_dir, cell_type, "_", norm_method, "_genPC_", geno_pc, "_expPC_", exp_pc, "/")
}

# Summary numbers from one tensorQTL cis (permutation) log; one-row tibble.
# The cis-window is read as the last number on its line (tensorQTL prints
# '\u00b1500,000'), so the result does not depend on the locale.
parse_perm_log <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)

  gene_line      <- stringr::str_subset(lines, "phenotypes tested:")
  sample_line    <- stringr::str_subset(lines, "samples")
  covariate_line <- utils::tail(stringr::str_subset(lines, "covariates"), 1)
  variants_line  <- stringr::str_subset(lines, "variants")
  qtl_line       <- stringr::str_subset(lines, "QTL phenotypes @ FDR 0.05:")
  window_line    <- stringr::str_subset(lines, "cis-window:")
  run_line       <- stringr::str_subset(lines, "Time elapsed:")
  beta_line      <- stringr::str_subset(lines, "Correlation between Beta-approximated and empirical p-values:")
  min_p_line     <- stringr::str_subset(lines, "min p-value threshold @ FDR 0.05:")
  prop_sig_line  <- stringr::str_subset(lines, "Proportion of significant phenotypes")

  min_p <- as.numeric(stringr::str_extract(min_p_line, "[0-9\\.eE-]+$"))

  tibble::tibble(
    genes_tested  = as.integer(stringr::str_extract(gene_line, "\\d+$")),
    samples       = as.integer(stringr::str_extract(sample_line, "\\d+")),
    covariates    = as.integer(stringr::str_extract(covariate_line, "\\d+(?= covariates)")),
    variants      = as.integer(stringr::str_extract(variants_line, "\\d+")),
    qtl_fdr_0.05  = as.integer(stringr::str_extract(qtl_line, "\\d+$")),
    cis_window_kb = as.integer(stringr::str_extract(gsub(",", "", window_line), "\\d+$")) %/% 1000L,
    beta_cor      = as.numeric(stringr::str_extract(beta_line, "[0-9\\.eE-]+$")),
    min_p         = formatC(min_p, format = "e", digits = 2),
    prop_sig      = as.integer(stringr::str_extract(prop_sig_line, "\\d+$")),
    run_time_min  = as.integer(stringr::str_extract(run_line, "\\d+"))
  )
}

# parse_perm_log() for every cell type x expPC x norm method combination; one row
# per log found (missing logs are skipped with a message). Rows are ordered with
# cell type varying fastest, then expPC, then norm method.
summarise_logs <- function(cell_types, in_dir, exp_pcs, geno_pc, norm_methods) {
  grid <- expand.grid(cell_type = cell_types, expPC = exp_pcs, norm_method = norm_methods,
                      stringsAsFactors = FALSE)

  purrr::pmap(grid, function(cell_type, expPC, norm_method) {
    log_file <- paste0(tqtl_dir(in_dir, cell_type, norm_method, geno_pc, expPC),
                       cell_type, "_", norm_method, "_perm.tensorQTL.cis.log")
    if (!file.exists(log_file)) {
      message("Log file not found: ", log_file)
      return(NULL)
    }
    dplyr::bind_cols(
      tibble::tibble(cell_type = cell_type, genPC = geno_pc, expPC = expPC,
                     norm_methods = norm_method),
      parse_perm_log(log_file))
  }) |>
    dplyr::bind_rows()
}

# Lead eQTL at qval < fdr with a TSS distance, one tensorQTL permutation file per
# cell type (expPC from exp_pc_map). Adds distance_kb and cell_type; missing files
# and cell types with no significant eQTL are skipped.
read_sig_eqtls <- function(cell_types, in_dir, exp_pc_map, norm_method, geno_pc, fdr) {
  purrr::map(cell_types, function(cell_type) {
    f <- paste0(tqtl_dir(in_dir, cell_type, norm_method, geno_pc, exp_pc_map[[cell_type]]),
                cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz")
    if (!file.exists(f)) {
      message("File not found: ", f)
      return(NULL)
    }
    dat <- readr::read_tsv(f, show_col_types = FALSE) |>
      dplyr::filter(!is.na(start_distance)) |>
      dplyr::mutate(distance_kb = start_distance / 1000) |>
      dplyr::filter(qval < fdr)
    if (nrow(dat) == 0) return(NULL)
    dplyr::mutate(dat, cell_type = cell_type)
  }) |>
    dplyr::bind_rows()
}

# Distinct eGenes (qval < fdr) for one cell type: phenotype_id, cell_type.
# A missing file gives an empty tibble.
read_egenes <- function(cell_type, in_dir, norm_method, geno_pc, exp_pc, fdr) {
  f <- paste0(tqtl_dir(in_dir, cell_type, norm_method, geno_pc, exp_pc),
              cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz")
  if (!file.exists(f)) return(tibble::tibble(cell_type = cell_type, eGene = character(0)))

  readr::read_tsv(f, show_col_types = FALSE) |>
    dplyr::filter(!is.na(qval) & qval < fdr) |>
    dplyr::distinct(phenotype_id) |>
    dplyr::mutate(cell_type = cell_type)
}

# Conditionally independent eQTL (tensorQTL cis_independent), all cell types.
read_indep <- function(cell_types, ind_dir, exp_pc_map, norm_method, geno_pc) {
  purrr::map(cell_types, function(cell_type) {
    f <- paste0(tqtl_dir(ind_dir, cell_type, norm_method, geno_pc, exp_pc_map[[cell_type]]),
                cell_type, "_", norm_method, "_independent.cis_independent_qtl.txt.gz")
    readr::read_tsv(f, show_col_types = FALSE) |>
      dplyr::mutate(cell_type = cell_type)
  }) |>
    dplyr::bind_rows()
}

# House style for the report's gt tables: bold column labels, black rules.
# subtitle = TRUE also sets the subtitle font size (sensitivity and settings tables).
gt_style <- function(tbl, subtitle = FALSE) {
  tbl <- gt::tab_style(tbl, style = gt::cell_text(weight = "bold"),
                       locations = gt::cells_column_labels())
  opts <- list(
    table.font.size = 14,
    heading.title.font.size = 16,
    heading.title.font.weight = "bold",
    heading.subtitle.font.size = if (subtitle) 14 else NULL,
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.top.color = "black",
    table.border.bottom.color = "transparent",
    table_body.hlines.color = "grey80",
    table_body.border.bottom.color = "black"
  )
  do.call(gt::tab_options, c(list(data = tbl), purrr::compact(opts)))
}

# eGenes (FDR < 0.05, from the logs) by expPC, one panel per norm method.
plot_sensitivity <- function(summary_tbl, palette, y_max) {
  make_plot <- function(norm) {
    ggplot2::ggplot(dplyr::filter(summary_tbl, norm_methods == norm),
                    ggplot2::aes(x = expPC, y = qtl_fdr_0.05, color = cell_type, group = cell_type)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::labs(
        x = "Expression PCs",
        y = "eGene count (FDR < 0.05)",
        color = "Cell Type",
        title = paste("No. of Sig. QTL by Expr PCs\nNorm method:", norm)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black"),
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      ) + ggplot2::ylim(0, y_max)
  }
  plots <- lapply(unique(summary_tbl$norm_methods), make_plot)
  cowplot::plot_grid(plotlist = plots, ncol = 2)
}

# Best expPC / norm method (most eGenes) per cell type.
sensitivity_best_gt <- function(summary_tbl) {
  summary_tbl |>
    dplyr::select(cell_type, expPC, norm_methods, qtl_fdr_0.05) |>
    dplyr::group_by(cell_type) |>
    dplyr::slice_max(qtl_fdr_0.05, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    gt::gt() |>
    gt::tab_header(title = "Table 1. Best param configuration per cell type") |>
    gt::cols_label(
      cell_type = "Cell Type",
      expPC = "Expression PCs",
      norm_methods = "Normalisation method",
      qtl_fdr_0.05 = "No. of eGenes"
    ) |>
    gt_style(subtitle = TRUE)
}

# eGenes for every expPC / norm method; one column per cell type.
sensitivity_all_gt <- function(summary_tbl) {
  summary_tbl |>
    dplyr::select(cell_type, expPC, norm_methods, qtl_fdr_0.05) |>
    tidyr::pivot_wider(
      id_cols     = c(expPC, norm_methods),
      names_from  = cell_type,
      values_from = qtl_fdr_0.05,
      names_sort  = TRUE) |>
    dplyr::arrange(dplyr::desc(expPC)) |>
    gt::gt() |>
    gt::tab_header(title = "Table 2. All param configurations ") |>
    gt::cols_label(
      expPC = "Expression PCs",
      norm_methods = "Normalisation method"
    ) |>
    gt_style(subtitle = TRUE)
}

# Shared theme for the density and count plots.
report_theme <- function() {
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(face = "bold")
  )
}

# TSS distance density for one cell type's significant eQTL (distance_kb).
plot_density <- function(eqtl_tbl, title) {
  ggplot2::ggplot(eqtl_tbl, ggplot2::aes(x = distance_kb)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::annotate("text", x = 100, y = Inf,
                      label = paste("eQTLs:", nrow(eqtl_tbl)),
                      hjust = 0, vjust = 1.5, size = 4, fontface = "bold") +
    ggplot2::scale_x_continuous(limits = c(-600, 600),
                                breaks = seq(-600, 600, by = 200),
                                labels = seq(-600, 600, by = 200)) +
    ggplot2::labs(x = "Distance to TSS (kb)", y = "Density", title = title) +
    ggplot2::theme_minimal() +
    report_theme()
}

# TSS distance density, one line per cell type.
plot_density_combined <- function(eqtls_all, palette) {
  ggplot2::ggplot(eqtls_all, ggplot2::aes(x = distance_kb, color = cell_type)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_x_continuous(limits = c(-600, 600), breaks = seq(-600, 600, 200)) +
    ggplot2::labs(x = "Distance to TSS (kb)", y = "Density") +
    ggplot2::theme_minimal() +
    report_theme()
}

# eGenes per cell type; gene_cell has one row per phenotype_id x cell_type.
plot_egene_counts <- function(gene_cell, palette) {
  dplyr::count(gene_cell, cell_type) |>
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(cell_type, n), y = n, fill = cell_type)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Cell type", y = "Number of eGenes") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    report_theme() +
    ggplot2::theme(legend.position = "none")
}

# Independent eQTL counts per cell type and rank (primary, secondary, tertiary).
indep_counts_gt <- function(all_eqtls, title) {
  all_eqtls |>
    dplyr::count(cell_type, rank) |>
    tidyr::pivot_wider(names_from = rank, values_from = n, values_fill = 0) |>
    dplyr::rename(primary = `1`, secondary = `2`, tertiary = `3`) |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::cols_label(
      cell_type = "Cell Type",
      primary   = "Primary eQTL",
      secondary = "Secondary eQTL",
      tertiary  = "Tertiary eQTL"
    ) |>
    gt::fmt_number(columns = c(primary, secondary, tertiary), decimals = 0, use_seps = TRUE) |>
    gt_style()
}

# Median absolute TSS distance (kb) per cell type and rank.
indep_tss_gt <- function(all_eqtls, title) {
  all_eqtls |>
    dplyr::mutate(tss_distance = abs(start_distance)) |>
    dplyr::group_by(cell_type, rank) |>
    dplyr::summarise(median_tss_distance = stats::median(tss_distance, na.rm = TRUE) / 1000,
                     .groups = "drop") |>
    tidyr::pivot_wider(names_from = rank, values_from = median_tss_distance, values_fill = 0) |>
    dplyr::rename(primary = `1`, secondary = `2`, tertiary = `3`) |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::cols_label(
      cell_type = "Cell Type",
      primary   = "Primary eQTL",
      secondary = "Secondary eQTL",
      tertiary  = "Tertiary eQTL"
    ) |>
    gt::fmt_number(columns = c(primary, secondary, tertiary), decimals = 0, use_seps = TRUE) |>
    gt_style()
}

# Genes tested, samples and variants per cell type, from the logs of the runs used
# for results (final norm method, expPC from exp_pc_map).
inputs_gt <- function(summary_tbl, exp_pc_map, final_norm, title) {
  exp_pc_tbl <- tibble::tibble(cell_type = names(exp_pc_map), expPC = as.numeric(exp_pc_map))

  summary_tbl |>
    dplyr::filter(norm_methods == final_norm) |>
    dplyr::inner_join(exp_pc_tbl, by = c("cell_type", "expPC")) |>
    dplyr::select(cell_type, genes_tested, samples, variants) |>
    dplyr::distinct(cell_type, .keep_all = TRUE) |>
    dplyr::arrange(cell_type) |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::cols_label(
      cell_type = "Cell Type",
      genes_tested = "Genes Tested",
      samples = "Samples",
      variants = "Variants"
    ) |>
    gt::fmt_number(columns = c(genes_tested, samples, variants), decimals = 0, use_seps = TRUE) |>
    gt_style(subtitle = TRUE)
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
