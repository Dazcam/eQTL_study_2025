#--------------------------------------------------------------------------------------
#
#    pseudotime_split_covariates.R
#
#    Subset the base covariate matrix to a specific combination of genotype
#    and expression PCs for tensorQTL input.
#
#    Adapted from tensorqtl_split_covariates.R for the pseudotime trajectory
#    analysis. Genotype PCs are fixed at 4; expression PCs are tunable via
#    the exp_pc wildcard.
#
#    Covariates retained:
#      - PCW
#      - sex_code (renamed to Sex)
#      - genPC1 to genPC4 (fixed)
#      - expPC1 to expPC{exp_pc} (tunable)
#
#    Usage: called via Snakemake rule pseudotime_split_covariates
#
#--------------------------------------------------------------------------------------

## Set up logging  -------------------------------------------------------------------

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

message("\n\nSplitting covariates for pseudotime tensorQTL ...")

## Load packages  --------------------------------------------------------------------

library(tidyverse)

## Parameters  ------------------------------------------------------------------------

cov_input    <- as.character(snakemake@input)
cov_output   <- as.character(snakemake@output)
trajectory   <- snakemake@wildcards[["trajectory"]]
quantile_n   <- snakemake@wildcards[["quantile_n"]]
bin          <- snakemake@wildcards[["bin"]]
exp_pc       <- as.integer(snakemake@wildcards[["exp_pc"]])
geno_pc      <- 4L

message("============================")
tibble(
  variable = c("cov_input", "cov_output", "trajectory",
               "quantile_n", "bin", "geno_pc", "exp_pc"),
  value    = c(cov_input, cov_output, trajectory,
               quantile_n, bin, geno_pc, exp_pc)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

## Load and subset covariate matrix  -------------------------------------------------

cov_tbl <- read_tsv(cov_input, show_col_types = FALSE) |>
  rename(covariate = `...1`)

message("Input covariate matrix for ", trajectory,
        " Q", quantile_n, " bin ", bin, ":\n")
print(cov_tbl)

cov_final_tbl <- cov_tbl |>
  filter(covariate %in% c(
    "PCW",
    "sex_code",
    paste0("genPC", seq_len(geno_pc)),
    paste0("expPC", seq_len(exp_pc))
  )) |>
  mutate(covariate = if_else(covariate == "sex_code", "Sex", covariate))

message("\nOutput covariate matrix (", nrow(cov_final_tbl),
        " covariates x ", ncol(cov_final_tbl) - 1, " samples):\n")
print(cov_final_tbl)

## Write output  ----------------------------------------------------------------------

cov_matrix           <- as.matrix(cov_final_tbl[, -1])
rownames(cov_matrix) <- cov_final_tbl$covariate
colnames(cov_matrix) <- colnames(cov_final_tbl)[-1]

message("\nCovariate matrix (first 5 rows and columns):\n")
print(cov_matrix[1:min(5, nrow(cov_matrix)), 1:min(5, ncol(cov_matrix))])

dir.create(dirname(cov_output), recursive = TRUE, showWarnings = FALSE)

write.table(cov_matrix,
            file      = cov_output,
            sep       = "\t",
            quote     = FALSE,
            col.names = NA,
            row.names = TRUE,
            append    = FALSE)

message("\nWritten: ", cov_output)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
