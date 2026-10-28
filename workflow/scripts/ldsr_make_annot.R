#--------------------------------------------------------------------------------------
#
#    Make annotation files (maxCPP and CS95) for SLDSR
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09SLDSR | Rule: make_annot (one job per cell type and chromosome)
#            Upstream:   08SUSIE sort_susie (credible sets); get_hg38_refs (.bim files)
#            Downstream: ldsr_ld_scores_hg38
#
# Purpose:   - Build two annotations from SuSiE credible set variants, one row per
#              SNP in the hg38 1000G .bim file, in .bim order:
#              - maxCPP: max PIP across eGenes per variant (0 if not in any CS)
#              - CS95:   1 if the variant is in any 95% credible set, else 0
#
# Inputs:    cred_set  SuSiE credible set variants for one cell type (all chr)
#            bim       hg38 1000G .bim file for one chromosome
#
# Outputs:   maxcpp    <cell_type>.maxCPP.<chr>.annot.gz (CHR, BP, SNP, CM, maxCPP)
#            cs95      <cell_type>.CS95.<chr>.annot.gz   (CHR, BP, SNP, CM, CS95)
#
# Notes:     Variants are matched to the .bim on chr and position, not rsID.
#            ldsc --l2 needs the annotation rows to match the .bim exactly.
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

message("\n\nCreating SLDSR annotation files from SuSiE output ...")

##  Load packages and variables  ------------------------------------------------------
library(tidyverse)

# Input and output paths
input_cs   <- snakemake@input[["cred_set"]]
input_bim  <- snakemake@input[["bim"]]
out_maxcpp <- snakemake@output[["maxcpp"]]
out_cs95   <- snakemake@output[["cs95"]]

# Wildcards
cell_type <- snakemake@wildcards[["cell_type"]]
chr_num   <- snakemake@wildcards[["chr"]]

# Genome build of the reference, for the log only
build <- str_extract(basename(input_bim), "hg\\d+")

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("input_cs", "input_bim", "out_maxcpp", "out_cs95",
               "cell_type", "chr_num", "build"),
  value    = c(input_cs, input_bim, out_maxcpp, out_cs95,
               cell_type, chr_num, build)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

message("Creating maxCPP and CS95 ", build, " annot files for ",
        cell_type, " on chr ", chr_num, " ...")

##  Read data  ------------------------------------------------------------------------
susie_tbl <- read_tsv(input_cs, show_col_types = FALSE)
bim_file  <- read_tsv(input_bim, col_names = FALSE, show_col_types = FALSE)

message(nrow(susie_tbl), " variants in credible sets before annotation calculation ...")

##  Compute annotations  --------------------------------------------------------------
susie_tbl_chr <- susie_tbl |> dplyr::filter(chr == chr_num)

# maxCPP: max PIP across eGenes per position (one row per chr:pos, so the
# join to the .bim can't duplicate rows at multi-allelic sites)
susie_maxcpp <- susie_tbl_chr |>
  summarise(maxCPP = max(pip), .by = c(chr, pos))

# CS95: 1 if the variant is in any credible set
susie_cs95 <- susie_tbl_chr |>
  distinct(chr, pos) |>
  mutate(CS95 = 1L)

message(nrow(susie_maxcpp), " unique variants for maxCPP on chr ", chr_num, " ...")
message(nrow(susie_cs95), " unique variants for CS95 on chr ", chr_num, " ...")

##  Match annotations to .bim  --------------------------------------------------------
# left_join() keeps .bim row order; variants not in a credible set get 0
bim_ref <- bim_file |>
  dplyr::select(CHR = X1, BP = X4, SNP = X2, CM = X3)

bim_maxcpp <- bim_ref |>
  left_join(susie_maxcpp, by = c("CHR" = "chr", "BP" = "pos")) |>
  mutate(maxCPP = replace_na(maxCPP, 0))

bim_cs95 <- bim_ref |>
  left_join(susie_cs95, by = c("CHR" = "chr", "BP" = "pos")) |>
  mutate(CS95 = replace_na(CS95, 0L))

# ldsc --l2 needs the annotation rows to match the .bim exactly
stopifnot(
  "maxCPP annotation does not contain the same SNPs as the .bim file" =
    identical(bim_maxcpp$SNP, bim_ref$SNP),
  "CS95 annotation does not contain the same SNPs as the .bim file" =
    identical(bim_cs95$SNP, bim_ref$SNP)
)

message(sum(bim_maxcpp$maxCPP > 0), " variants with maxCPP > 0 on chr ", chr_num, " ...")
message(sum(bim_cs95$CS95 > 0), " variants in any 95% CS on chr ", chr_num, " ...")

##  Write outputs  --------------------------------------------------------------------
write.table(bim_maxcpp, gzfile(out_maxcpp), col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")
write.table(bim_cs95, gzfile(out_cs95), col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
