#--------------------------------------------------------------------------------------
#
#    Generate MaxCPP S-LDSR enrichment tables for manuscript
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  13MANUSCRIPT_TABLES | Rule: sldsr_tbl
#            Upstream:   09SLDSR strat_bl_v12 (MaxCPP summary, all cell types x GWAS)
#
# Purpose:   Reproduce the table behind the MaxCPP S-LDSR figure (panels A-C):
#            - Parse cell type and disorder from the annotation Category
#            - Relabel cell types (Glu-UL -> Glu-A, Glu-DL -> Glu-B, drop "-Q4-")
#            - Keep only cell types listed in config['cell_types']
#            - Add the plotted -log10(P), level (L1 / L2 / Pseudotime), main cell
#              type, per-level Bonferroni threshold and significance flags
#
# Inputs:    in_file     ldsr_strat_hg38_bl_v12.maxCPP.summary.tsv
#            cell_types  config['cell_types'] (pipeline labels)
#
# Outputs:   xlsx        One sheet per disorder (sheet name = full disorder name)
#
# Notes:     P is one-sided on the coefficient z-score (enrichment > 0), as in the
#            plot. The plotted -log10(P) is set to 0 where z <= 0. Bonferroni is
#            per level: 0.05 / number of cell types in that level, counted from
#            the data after filtering to config['cell_types'].
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

message("\n\nGenerating MaxCPP S-LDSR tables for the manuscript ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Input and output paths
in_file    <- snakemake@input[[1]]
xlsx_out   <- snakemake@output[[1]]
cell_types <- unlist(snakemake@params[["cell_types"]])

# in_file    <- "../results/09SLDSR/strat_bl_v12/ldsr_strat_hg38_bl_v12.maxCPP.summary.tsv"
# xlsx_out   <- "../results/13MANUSCRIPT_TABLES/supplementary_table_X.xlsx"
# cell_types <- c("Glu-UL", "Glu-DL", "NPC", ...)   # config['cell_types']

# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_file", "xlsx_out", "cell_types"),
  value    = c(in_file, xlsx_out, paste(cell_types, collapse = ", "))) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

alpha <- 0.05

# Relabel helper: Glu-UL -> Glu-A, Glu-DL -> Glu-B, strip "-Q4-" from bin names
# (same convention as the rest of the manuscript pipeline)
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

# Full disorder names for xlsx sheets (same mapping as manuscript_table_smr_sig.R)
disorder_labels <- c(
  SCZ  = "Schizophrenia",
  BPD  = "Bipolar disorder",
  MDD  = "Major depressive disorder",
  ADHD = "ADHD",
  OCD  = "Obsessive-compulsive disorder"
)

## Read and prepare table  -------------------------------------------------------------
message("Reading: ", in_file)
ldsr_tbl <- read_tsv(in_file, show_col_types = FALSE) |>
  separate(Category, into = c("annot", "file_name"), sep = "/", remove = FALSE) |>
  separate(file_name, into = c("cell_type", "suffix"), sep = "_hg38", remove = TRUE) |>
  mutate(disorder = toupper(str_extract(suffix, "scz|bpd|mdd|adhd|ocd"))) |>
  dplyr::select(-Category, -suffix) |>
  filter(!is.na(disorder))                        # drops PTSD (not in the manuscript)

# Keep only config cell types (pipeline labels, so filter before relabelling)
not_in_config <- setdiff(unique(ldsr_tbl$cell_type), cell_types)
if (length(not_in_config) > 0) {
  message("Dropping cell types not in config['cell_types']: ", paste(not_in_config, collapse = ", "))
}
not_in_file <- setdiff(cell_types, unique(ldsr_tbl$cell_type))
if (length(not_in_file) > 0) {
  message("WARNING: config cell types absent from the summary file: ", paste(not_in_file, collapse = ", "))
}

ldsr_tbl <- ldsr_tbl |>
  filter(cell_type %in% cell_types) |>
  mutate(cell_type = relabel_cell_type(cell_type))

# Trajectory rows are flagged before the digit-based L1/L2 split, since bin
# names contain digits and would otherwise be classed as L2 subclusters
ldsr_tbl <- ldsr_tbl |>
  mutate(
    level = case_when(
      str_detect(cell_type, "-to-") ~ "Pseudotime",
      str_detect(cell_type, "\\d")  ~ "L2",
      TRUE                          ~ "L1"),
    main_type = case_when(
      level == "Pseudotime" & str_detect(cell_type, "to-Glu-A") ~ "Trajectory-to-Glu-A",
      level == "Pseudotime" & str_detect(cell_type, "to-Glu-B") ~ "Trajectory-to-Glu-B",
      str_detect(cell_type, "Glu-A")     ~ "Glu-A",
      str_detect(cell_type, "Glu-B")     ~ "Glu-B",
      str_detect(cell_type, "GABA")      ~ "GABA",
      str_detect(cell_type, "NPC")       ~ "NPC",
      str_detect(cell_type, "OPC")       ~ "OPC",
      str_detect(cell_type, "MG")        ~ "MG",
      str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
      TRUE ~ cell_type),
    P_coef_one_sided = pnorm(`Coefficient_z-score`, lower.tail = FALSE),
    neg_log10_P_plot = if_else(`Coefficient_z-score` > 0, -log10(P_coef_one_sided), 0)
  )

## Bonferroni thresholds per level  ----------------------------------------------------
n_tests <- ldsr_tbl |>
  distinct(level, cell_type) |>
  count(level, name = "n_tests")

message("\nCell types tested per level (Bonferroni denominator):")
n_tests |> knitr::kable(format = "simple", align = "l") |> print()

ldsr_tbl <- ldsr_tbl |>
  left_join(n_tests, by = "level") |>
  mutate(
    P_bonferroni_thresh = alpha / n_tests,
    sig_bonferroni      = P_coef_one_sided < P_bonferroni_thresh,
    sig_nominal         = P_coef_one_sided < alpha
  )

## Checks  ------------------------------------------------------------------------------
dups <- ldsr_tbl |> count(cell_type, disorder) |> filter(n > 1)
if (nrow(dups) > 0) {
  stop("Duplicate cell type x disorder rows: ",
       paste(dups$cell_type, dups$disorder, sep = "/", collapse = ", "))
}

n_expected <- n_distinct(ldsr_tbl$cell_type) * n_distinct(ldsr_tbl$disorder)
if (nrow(ldsr_tbl) != n_expected) {
  message("WARNING: ", n_expected - nrow(ldsr_tbl),
          " cell type x disorder combinations missing from the summary file")
}

## Order rows as in the plot  -----------------------------------------------------------
cell_order <- c(
  "Glu-A", sort(grep("^Glu-A-[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  "Glu-B", sort(grep("^Glu-B-[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  "GABA",  sort(grep("^GABA-[0-9]+$",  ldsr_tbl$cell_type, value = TRUE)),
  "NPC",   sort(grep("^NPC-[0-9]+$",   ldsr_tbl$cell_type, value = TRUE)),
  "OPC", "MG", "Endo-Peri",
  sort(grep("^NPC-to-Glu-A-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  sort(grep("^NPC-to-Glu-B-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE))
) |> unique()

unordered <- setdiff(unique(ldsr_tbl$cell_type), cell_order)
if (length(unordered) > 0) {
  stop("Cell types not covered by the plot ordering: ", paste(unordered, collapse = ", "))
}

final_tbl <- ldsr_tbl |>
  mutate(
    disorder = factor(disorder, levels = names(disorder_labels)),
    level    = factor(level, levels = c("L1", "L2", "Pseudotime")),
    order    = match(cell_type, cell_order)
  ) |>
  arrange(disorder, level, order) |>
  dplyr::select(-order) |>
  relocate(disorder, cell_type, level, main_type, annot) |>
  mutate(across(c(disorder, level), as.character))

## Write output  ------------------------------------------------------------------------
sheets <- final_tbl |>
  split(final_tbl$disorder) |>
  map(\(x) dplyr::select(x, -disorder))
sheets <- sheets[intersect(names(disorder_labels), names(sheets))]
names(sheets) <- disorder_labels[names(sheets)]

message("Writing: ", xlsx_out, " (sheets: ", paste(names(sheets), collapse = ", "), ")")
write.xlsx(
  sheets,
  file = xlsx_out,
  overwrite = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
