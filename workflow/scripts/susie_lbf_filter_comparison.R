#--------------------------------------------------------------------------------------
#
#    Compare SuSiE credible set retention with vs without an LBF filter
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08SUSIE | Rule: susie_lbf_comparison
#            Upstream:   sort_susie
#            Downstream: none (reviewer response)
#
# Purpose:   - Load the credible set variants for every cell type. The cred.hp files
#              are the baseline: susie_run.R has already removed low-purity CSs
#              (min.abs.corr < 0.5) and any CS sharing a variant with an earlier CS
#              for the same gene
#            - Collapse to one row per credible set
#            - Count CSs and eGenes retained with an extra filter of
#              cs_log10bf > lbf_thresh, and of ln(BF) > lbf_thresh
#            - Plot the cs_log10bf distribution per cell type
#
# Inputs:    input_files  *.mrgd.srtd.susie.cred.hp.txt.gz, one per cell type
#                         (variant-level: one row per variant per credible set)
#
# Outputs:   out_overall  TSV of CS / eGene counts per filter, all cell types combined
#            out_by_ct    TSV of CS / eGene counts per filter and cell type
#            out_full     gzipped TSV of all input rows with filter flags
#            out_plot     PNG histogram of cs_log10bf per cell type
#
# Notes:     cs_log10bf is log10(BF): susie_run.R converts susieR's natural-log LBF
#            with log10(exp(lbf)). The reviewer's "LBF > 2" gave no base, so both
#            interpretations are reported.
#            out_by_ct, out_full and out_plot are derived from the out_overall path
#            and are not declared in the rule.
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

message("\n\nComparing SuSiE credible sets with vs without LBF filter ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages(library(tidyverse))

# Input and output paths
input_files <- snakemake@input[["susie_files"]]
out_overall <- snakemake@output[["overall"]]
out_by_ct   <- str_replace(out_overall, "_overall\\.tsv$", "_by_celltype.tsv")
out_full    <- str_replace(out_overall, "_overall\\.tsv$", "_full_annotated.tsv.gz")
out_plot    <- str_replace(out_overall, "_overall\\.tsv$", "_distribution.png")
lbf_thresh  <- as.numeric(snakemake@params[["lbf_thresh"]])

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("input_files", "out_overall", "out_by_ct", "out_full", "out_plot",
               "lbf_thresh"),
  value    = c(paste(input_files, collapse = ", "), out_overall, out_by_ct, out_full,
               out_plot, lbf_thresh)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load credible sets  -----------------------------------------------------------------
cs_all <- map_dfr(input_files, function(path) {
  ct <- basename(path) |> str_extract("^[^.]+")
  message("Loading ", ct, ": ", path)
  read_tsv(path, show_col_types = FALSE) |>
    mutate(cell_type = ct)
})

required_cols <- c("phenotype_id", "cs_id", "cs_log10bf", "cs_size")
missing_cols  <- setdiff(required_cols, colnames(cs_all))
if (length(missing_cols) > 0) {
  stop("Input file(s) missing expected column(s): ", paste(missing_cols, collapse = ", "),
       "\nFound columns: ", paste(colnames(cs_all), collapse = ", "))
}

cs_all <- cs_all |>
  mutate(
    lbf_natural      = cs_log10bf * log(10),
    passes_lbf_log10 = cs_log10bf  > lbf_thresh,
    passes_lbf_ln    = lbf_natural > lbf_thresh)

# One row per credible set, so counts are not inflated by CS size
cs_level <- cs_all |>
  distinct(cell_type, cs_id, .keep_all = TRUE) |>
  dplyr::select(cell_type, phenotype_id, cs_id, cs_index, cs_size, cs_log10bf,
                lbf_natural, passes_lbf_log10, passes_lbf_ln)

baseline       <- cs_level
with_lbf_log10 <- baseline |> dplyr::filter(passes_lbf_log10)
with_lbf_ln    <- baseline |> dplyr::filter(passes_lbf_ln)

##  Per cell type summary  --------------------------------------------------------------
summarise_set <- function(df, label) {
  df |>
    group_by(cell_type) |>
    summarise(
      n_credible_sets = n(),
      n_eGenes        = n_distinct(phenotype_id),
      median_cs_size  = median(cs_size),
      .groups = "drop") |>
    mutate(filter = label)
}

summary_df <- bind_rows(
  summarise_set(baseline,       "baseline_no_lbf_filter"),
  summarise_set(with_lbf_log10, "baseline_plus_LBFlog10_gt_thresh"),
  summarise_set(with_lbf_ln,    "baseline_plus_LBFln_gt_thresh"))

message("\nCS / eGene retention by filter:")
print(summary_df)

##  Overall summary  --------------------------------------------------------------------
overall_df <- bind_rows(
  baseline       |> mutate(filter = "baseline_no_lbf_filter"),
  with_lbf_log10 |> mutate(filter = "baseline_plus_LBFlog10_gt_thresh"),
  with_lbf_ln    |> mutate(filter = "baseline_plus_LBFln_gt_thresh")) |>
  group_by(filter) |>
  summarise(
    n_credible_sets = n(),
    n_eGenes        = n_distinct(phenotype_id),
    pct_of_baseline = round(100 * n() / nrow(baseline), 1),
    .groups = "drop")

message("\nOverall totals (all cell types combined):")
print(overall_df)

message("\ncs_log10bf distribution among purity-passed CSs:")
print(summary(baseline$cs_log10bf))

##  Sanity checks  ----------------------------------------------------------------------
message("\nSanity checks:")
message("Total credible sets (baseline): ", nrow(baseline))
message("NA cs_log10bf: ", sum(is.na(baseline$cs_log10bf)))
message("Inf cs_log10bf: ", sum(is.infinite(baseline$cs_log10bf)))
message("Min cs_log10bf: ", signif(min(baseline$cs_log10bf, na.rm = TRUE), 4),
        "  (ln(BF) > ", lbf_thresh, " equals log10(BF) > ",
        signif(lbf_thresh / log(10), 4), ")")
message("Max finite cs_log10bf: ",
        signif(max(baseline$cs_log10bf[is.finite(baseline$cs_log10bf)]), 4))
message("CSs with ln(BF) <= ", lbf_thresh, ": ",
        sum(baseline$lbf_natural <= lbf_thresh, na.rm = TRUE))
message("CSs with log10(BF) <= ", lbf_thresh, ": ",
        sum(baseline$cs_log10bf <= lbf_thresh, na.rm = TRUE))

##  Write tables  -----------------------------------------------------------------------
write_tsv(summary_df, out_by_ct)
write_tsv(overall_df, out_overall)
write_tsv(cs_all, out_full)

##  Plot  -------------------------------------------------------------------------------
p <- ggplot(baseline, aes(x = cs_log10bf)) +
  geom_histogram(bins = 60, fill = "#1f77b4", alpha = 0.8) +
  geom_vline(xintercept = lbf_thresh, linetype = "dashed", color = "#d62728") +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(
    title = "Distribution of credible-set log10(BF), purity-passed CSs",
    subtitle = paste0("Dashed line: cs_log10bf = ", lbf_thresh,
                      " (confirm log10 vs natural-log base against reviewer's intent)"),
    x = "cs_log10bf", y = "Count of credible sets") +
  theme_minimal(base_size = 11)

ggsave(out_plot, p, width = 10, height = 6, dpi = 300)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
