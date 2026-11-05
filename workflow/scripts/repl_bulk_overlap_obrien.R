#--------------------------------------------------------------------------------------
#
#    Replication (bulk): overlap fetal eGenes with O'Brien 2018 bulk eQTL
#
#--------------------------------------------------------------------------------------

#  Reads the unique fetal eGene set (all L1+L2+pseudotime cell types, from
#  extract_unique_egenes) and overlaps it against O'Brien et al. 2018 bulk brain
#  eQTL supplementary tables:
#    - TblS2: gene-level significant eQTL (one row per significant eGene)
#    - TblS3: transcript-level significant eQTL (one row per significant
#      transcript; multiple transcripts can map to the same gene, so this is
#      deduplicated to ENSEMBL_geneID before overlap -- a gene counts as
#      "transcript-bulk significant" if >=1 of its transcripts is significant)
#
#  Both O'Brien sheets are already pre-filtered to significant results only (no
#  further qval threshold applied here).
#
#  Output is one row per unique fetal eGene, flagged for whether it overlaps
#  the gene-level and/or transcript-level bulk significant sets -- scoped
#  deliberately to just this overlap; the per-cell-type breakdown and plotting
#  are a separate downstream rule.
#
#  Step 2 of the replication_bulk pipeline.
#
#--------------------------------------------------------------------------------------

# Set up logging for Snakemake
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

message("\n\nOverlapping unique fetal eGenes with O'Brien 2018 bulk eQTL (gene + transcript level) ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# Input and output paths ---------------------------------------------------------------
fetal_egenes_file <- snakemake@input[["fetal_egenes"]]
obrien_supp_file  <- snakemake@input[["obrien_supp"]]
output_file       <- snakemake@output[[1]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("fetal_egenes_file", "obrien_supp_file", "output_file"),
  value    = c(fetal_egenes_file, obrien_supp_file, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# ========================================================================================
# 1. Load unique fetal eGene set
# ========================================================================================

message("Loading fetal eGenes (per-cell-type long table, deduplicating to unique) ...")
fetal_long <- read_tsv(fetal_egenes_file, show_col_types = FALSE)

fetal_egenes <- fetal_long |>
  distinct(phenotype_id) |>
  pull(phenotype_id)

message("  Unique fetal eGenes: ", length(fetal_egenes))

# ========================================================================================
# 2. Load O'Brien TblS2 (gene-level significant eQTL)
# ========================================================================================

message("\nLoading O'Brien TblS2 (gene-level significant eQTL), skip=2 ...")
tbls2 <- read_excel(obrien_supp_file, sheet = "TblS2", skip = 2)

n_tbls2_rows <- nrow(tbls2)
message("  TblS2 rows: ", n_tbls2_rows)

gene_bulk_egenes <- tbls2 |>
  distinct(ENSEMBL_geneID) |>
  pull(ENSEMBL_geneID)

n_gene_bulk_unique <- length(gene_bulk_egenes)
message("  TblS2 unique ENSEMBL_geneID: ", n_gene_bulk_unique,
        " (matches row count: ", n_gene_bulk_unique == n_tbls2_rows, ")")

# ========================================================================================
# 3. Load O'Brien TblS3 (transcript-level significant eQTL), dedup to gene
# ========================================================================================

message("\nLoading O'Brien TblS3 (transcript-level significant eQTL), skip=2 ...")
tbls3 <- read_excel(obrien_supp_file, sheet = "TblS3", skip = 2)

n_tbls3_rows <- nrow(tbls3)
message("  TblS3 rows: ", n_tbls3_rows)

n_tbls3_transcripts_unique <- n_distinct(tbls3$ENSEMBL_transcriptID)
message("  TblS3 unique ENSEMBL_transcriptID: ", n_tbls3_transcripts_unique)

transcript_bulk_egenes <- tbls3 |>
  distinct(ENSEMBL_geneID) |>
  pull(ENSEMBL_geneID)

n_transcript_bulk_unique <- length(transcript_bulk_egenes)
message("  TblS3 unique ENSEMBL_geneID (after dedup from transcript to gene level): ",
        n_transcript_bulk_unique)

# ========================================================================================
# 4. Overlap
# ========================================================================================

message("\n--- Overlap ---\n")

overlap_gene_bulk       <- intersect(fetal_egenes, gene_bulk_egenes)
overlap_transcript_bulk <- intersect(fetal_egenes, transcript_bulk_egenes)

message("Fetal eGenes overlapping O'Brien gene-level bulk (TblS2): ",
        length(overlap_gene_bulk), " / ", length(fetal_egenes))
message("Fetal eGenes overlapping O'Brien transcript-level bulk (TblS3, gene-deduped): ",
        length(overlap_transcript_bulk), " / ", length(fetal_egenes))

overlap_both <- intersect(overlap_gene_bulk, overlap_transcript_bulk)
message("Fetal eGenes overlapping BOTH gene- and transcript-level bulk: ",
        length(overlap_both))

# ========================================================================================
# 5. Write single tracked output
# ========================================================================================

overlap_tbl <- tibble(phenotype_id = fetal_egenes) |>
  mutate(
    in_gene_bulk       = phenotype_id %in% gene_bulk_egenes,
    in_transcript_bulk = phenotype_id %in% transcript_bulk_egenes
  )

write_tsv(overlap_tbl, output_file)
message("\nFetal-eGene x O'Brien bulk overlap table written to: ", output_file)
