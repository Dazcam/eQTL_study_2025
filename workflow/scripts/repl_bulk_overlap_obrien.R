#--------------------------------------------------------------------------------------
#
#    Overlap fetal eGenes with O'Brien 2018 bulk brain eQTL
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08REPLICATION-BULK | Rule: overlap_obrien (single job)
#            Upstream:   extract_unique_egenes
#            Downstream: repl_bulk_report; 17 bulk_egene_overlaps_fig (Supp Fig 9)
#
# Purpose:   - Take the unique fetal eGene set across all cell types
#            - Read O'Brien et al. 2018 supplementary tables:
#                TblS2  gene-level significant eQTL (one row per eGene)
#                TblS3  transcript-level significant eQTL, collapsed to genes
#                       (a gene counts if >= 1 of its transcripts is significant)
#            - Flag each fetal eGene for gene- and transcript-level overlap
#
# Inputs:    fetal_egenes  extract_unique_egenes output
#            obrien_supp   O'Brien 2018 supplementary tables (xlsx)
#
# Outputs:   output        TSV, one row per unique fetal eGene: phenotype_id,
#                          in_gene_bulk, in_transcript_bulk
#
# Notes:     Both O'Brien sheets hold significant results only, so no further
#            threshold is applied. Sheets are read with skip = 2 (title rows).
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

message("\n\nOverlapping fetal eGenes with O'Brien 2018 bulk eQTL ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(readxl)

# Input and output paths
fetal_egenes_file <- snakemake@input[["fetal_egenes"]]
obrien_supp_file  <- snakemake@input[["obrien_supp"]]
output_file       <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("fetal_egenes_file", "obrien_supp_file", "output_file"),
  value    = c(fetal_egenes_file, obrien_supp_file, output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load unique fetal eGenes  ----------------------------------------------------------
fetal_egenes <- read_tsv(fetal_egenes_file, show_col_types = FALSE) |>
  distinct(phenotype_id) |>
  pull(phenotype_id)

message("Unique fetal eGenes: ", length(fetal_egenes))

##  O'Brien TblS2: gene-level significant eQTL  ----------------------------------------
message("\nLoading O'Brien TblS2 (gene-level significant eQTL) ...")
tbls2 <- read_excel(obrien_supp_file, sheet = "TblS2", skip = 2)

gene_bulk_egenes <- tbls2 |>
  distinct(ENSEMBL_geneID) |>
  pull(ENSEMBL_geneID)

message("  TblS2 rows: ", nrow(tbls2))
message("  TblS2 unique ENSEMBL_geneID: ", length(gene_bulk_egenes),
        " (matches row count: ", length(gene_bulk_egenes) == nrow(tbls2), ")")

##  O'Brien TblS3: transcript-level significant eQTL  ----------------------------------
message("\nLoading O'Brien TblS3 (transcript-level significant eQTL) ...")
tbls3 <- read_excel(obrien_supp_file, sheet = "TblS3", skip = 2)

transcript_bulk_egenes <- tbls3 |>
  distinct(ENSEMBL_geneID) |>
  pull(ENSEMBL_geneID)

message("  TblS3 rows: ", nrow(tbls3))
message("  TblS3 unique ENSEMBL_transcriptID: ", n_distinct(tbls3$ENSEMBL_transcriptID))
message("  TblS3 unique ENSEMBL_geneID: ", length(transcript_bulk_egenes))

##  Overlap  ---------------------------------------------------------------------------
overlap_gene_bulk       <- intersect(fetal_egenes, gene_bulk_egenes)
overlap_transcript_bulk <- intersect(fetal_egenes, transcript_bulk_egenes)

message("\nFetal eGenes in O'Brien gene-level bulk (TblS2): ",
        length(overlap_gene_bulk), " / ", length(fetal_egenes))
message("Fetal eGenes in O'Brien transcript-level bulk (TblS3): ",
        length(overlap_transcript_bulk), " / ", length(fetal_egenes))
message("Fetal eGenes in both: ",
        length(intersect(overlap_gene_bulk, overlap_transcript_bulk)))

##  Write output  ----------------------------------------------------------------------
overlap_tbl <- tibble(phenotype_id = fetal_egenes) |>
  mutate(
    in_gene_bulk       = phenotype_id %in% gene_bulk_egenes,
    in_transcript_bulk = phenotype_id %in% transcript_bulk_egenes
  )

write_tsv(overlap_tbl, output_file)
message("\nOverlap table written to: ", output_file)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
