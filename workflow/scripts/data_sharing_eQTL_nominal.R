#--------------------------------------------------------------------------------------
#
#    Generate eQTL nominal files for data sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  15DATA_SHARING | Rule: data_nominal
#            Upstream:   10SMR smr_input (per-cell-type tensorQTL nominal pairs)
#            Downstream: mk_eqtl_tar
#
# Purpose:   For each cell type, take all nominal cis gene-SNP associations, then:
#            - Drop unresolvable SNP IDs ('.') and SNPs absent from the pvar
#            - Add CHROM, POS, REF, ALT from the pvar and gene symbols from BioMart
#            - Flag (but keep) rows with no gene symbol and duplicate gene-SNP pairs
#            - Relabel cell types to manuscript names (Glu-UL -> Glu-A, Glu-DL -> Glu-B,
#              drop "-Q4-")
#
# Inputs:    in_dir           {cell_type}/{cell_type}_nom.cis_qtl_pairs.tsv
#            allele_file      Genotype pvar (chrALL_final.filt.pvar)
#            gene_lookup_file BioMart gene annotation (hg38)
#
# Outputs:   {label}_cis_eQTL_nominal.tsv.gz        One per cell type (manuscript label)
#            eqtl_nominal_row_drop_summary.tsv      Rows in/removed/flagged per cell type
#            duplicated_pvar_snp_counts_nominal.tsv Duplicate SNP IDs in pvar
#            Note: snakemake output is a sentinel; dirname() sets the output dir
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

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

message("\n\neQTL nominal table for data sharing ...")

# -------------------------------------------------------------------------------------
suppressMessages(library(tidyverse))

# --- Set variables
allele_file      <- snakemake@params[["allele_file"]]
in_dir           <- snakemake@params[["in_dir"]]
cell_types       <- snakemake@params[["cell_types"]]
gene_lookup_file <- snakemake@params[["gene_lookup_file"]]
out_file         <- snakemake@output[[1]]
out_dir          <- paste0(dirname(out_file), '/')

# allele_file      <- "../results/05TENSORQTL/prep_input/chrALL_final.filt.pvar"
# in_dir           <- "../results/10SMR/smr_input/"
# cell_types       <- c("Glu-UL", "Glu-DL", "NPC", ...)   # config['cell_types']
# gene_lookup_file <- "../resources/sheets/gene_lookup_hg38.tsv"

# Make a tibble showing what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("allele_file", "in_dir", "cell_types", "gene_lookup_file", "out_file", "out_dir"),
  value    = c(allele_file, in_dir, paste(cell_types, collapse = ", "),
               gene_lookup_file, out_file, out_dir)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")


message("=== eQTL nominal export: ", length(cell_types), " cell types ===")

# --- Relabel helper (display only, output filenames): Glu-UL -> Glu-A, Glu-DL -> Glu-B,
# and strip "-Q4-" from pseudotime bin names, e.g. "NPC-to-Glu-DL-Q4-Bin1" ->
# "NPC-to-Glu-B-Bin1". Same convention as data_sharing_susie_cred_sets.R.
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

cell_type_labels <- relabel_cell_type(cell_types)
if (any(str_detect(cell_type_labels, "UL|DL|Q4"))) {
  stop("Old cell type naming remains in output labels: ",
       paste(cell_type_labels[str_detect(cell_type_labels, "UL|DL|Q4")], collapse = ", "))
}
names(cell_type_labels) <- cell_types  # pipeline label -> manuscript label

gene_lookup_tbl <- suppressMessages(read_tsv(gene_lookup_file, show_col_types = FALSE)) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ----- 1. Load genotype metadata (pvar) for alleles -----
message("Loading pvar file...")
pvar <- read_tsv(allele_file, comment = "#",
                 col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")) |>
  dplyr::select(CHROM, POS, ID, REF, ALT) |>
  filter(ID != '.')

dup_snps <- pvar |>
  count(ID) |>
  filter(n > 1)
write_tsv(dup_snps, paste0(out_dir, "duplicated_pvar_snp_counts_nominal.tsv"))
message("Duplicate SNP IDs in pvar: ", nrow(dup_snps))

message("Dropping SNP duplicates in pvar ...")
pvar <- pvar |>
  distinct(ID, .keep_all = TRUE)

# ----- 2. Iterate through cell types and build a TSV for each -----
drop_summary <- tibble(cell_type = character(), rows_in = integer(), dot_id = integer(),
                       no_pvar_match = integer(), no_gene_symbol = integer(),
                       duplicate_gene_snp = integer(), rows_out = integer())

for (cell_type in cell_types) {

  out_label <- cell_type_labels[[cell_type]]
  file_path <- paste0(in_dir, cell_type, '/', cell_type, '_nom.cis_qtl_pairs.tsv')

  if (!file.exists(file_path)) {
    stop("File not found for ", cell_type, ": ", file_path)
  }

  message('\n--- ', cell_type, ' -> ', out_label, ' ---')
  message('Loading eQTL nominal file: ', file_path)
  eqtl_raw <- read_tsv(file_path, show_col_types = FALSE) |>
    dplyr::rename(ensembl_id = phenotype_id, SNP = variant_id)

  n_in <- nrow(eqtl_raw)
  if (n_in == 0) {
    message("No rows for ", cell_type, " -- skipping")
    drop_summary <- drop_summary |> add_row(cell_type = out_label, rows_in = 0, dot_id = 0,
                                            no_pvar_match = 0, no_gene_symbol = 0,
                                            duplicate_gene_snp = 0, rows_out = 0)
    next
  }

  # --- Drop unresolvable SNP IDs ('.') -- not real IDs, not worth sharing
  n_dot <- sum(eqtl_raw$SNP == '.')
  eqtl_tbl <- eqtl_raw |> filter(SNP != '.')
  if (n_dot > 0) message("SNP ID = '.': ", n_dot, " rows removed")

  # --- Add alleles from pvar
  message('Adding alleles ...')
  eqtl_with_pvar <- eqtl_tbl |> left_join(pvar, by = c("SNP" = "ID"))
  n_no_pvar <- sum(is.na(eqtl_with_pvar$CHROM))
  if (n_no_pvar > 0) message("No match in pvar file: ", n_no_pvar, " rows removed")
  eqtl_with_pvar <- eqtl_with_pvar |>
    filter(!is.na(CHROM)) |>
    mutate(CHROM = ifelse(!str_detect(CHROM, "chr"), paste0("chr", CHROM), CHROM))

  # --- Add gene symbol
  message('Adding gene symbols ...')
  eqtl_enriched <- eqtl_with_pvar |>
    left_join(gene_lookup_tbl, by = join_by(ensembl_id == ensembl_gene_id)) |>
    mutate(CHROM = str_remove(CHROM, "^chr")) |>
    dplyr::select(ensembl_id, symbol = external_gene_name, CHROM,
                  SNP, POS, REF, ALT, AF = af, slope, slope_se, pval_nominal)

  n_no_gene <- sum(is.na(eqtl_enriched$symbol))
  if (n_no_gene > 0) message("No gene symbol match: ", n_no_gene, " rows kept with symbol = NA")

  n_out <- nrow(eqtl_enriched)

  # --- Duplicate gene-SNP pairs are flagged, not removed 
  dup_pairs <- eqtl_enriched |> count(ensembl_id, SNP) |> filter(n > 1)
  n_dup_extra <- if (nrow(dup_pairs) > 0) sum(dup_pairs$n - 1) else 0L
  if (nrow(dup_pairs) > 0) {
    message("Duplicate gene-SNP pairs: ", nrow(dup_pairs), " pairs, ", n_dup_extra,
            " extra rows (kept in output). Example: ", dup_pairs$ensembl_id[1], " / ", dup_pairs$SNP[1])
  }

  message("Processed ", cell_type, ": ", n_in, " in -> ", n_out, " out  (",
          n_dot, " dot-ID, ", n_no_pvar, " no pvar match, ", n_no_gene, " no gene symbol)")

  message('Writing file ...')
  write_tsv(eqtl_enriched, paste0(out_dir, out_label, '_cis_eQTL_nominal.tsv.gz'))

  drop_summary <- drop_summary |> add_row(
    cell_type = out_label, rows_in = n_in, dot_id = n_dot, no_pvar_match = n_no_pvar,
    no_gene_symbol = n_no_gene, duplicate_gene_snp = n_dup_extra, rows_out = n_out
  )
}

# --- Summary: rows removed vs. rows flagged-but-kept, by reason, per cell type and overall
message("\n=== Summary: row counts by reason (dot_id/no_pvar_match = removed; no_gene/dup_pairs = flagged, kept) ===")
totals <- drop_summary |>
  summarise(across(c(rows_in, dot_id, no_pvar_match, no_gene_symbol, duplicate_gene_snp, rows_out), sum)) |>
  mutate(cell_type = "TOTAL", .before = 1)
drop_summary_full <- bind_rows(drop_summary, totals)

pwalk(drop_summary_full, function(cell_type, rows_in, dot_id, no_pvar_match, no_gene_symbol,
                                  duplicate_gene_snp, rows_out) {
  message(sprintf("  %-25s in = %8d  dot_id = %6d  no_pvar = %6d  no_gene(kept) = %6d  dup_pairs(kept) = %6d  out = %8d",
                   cell_type, rows_in, dot_id, no_pvar_match, no_gene_symbol, duplicate_gene_snp, rows_out))
})

drop_summary_file <- paste0(out_dir, "eqtl_nominal_row_drop_summary.tsv")
write_tsv(drop_summary_full, drop_summary_file)
message("Wrote row-drop summary: ", drop_summary_file)

unaccounted <- totals$rows_in - (totals$rows_out + totals$dot_id + totals$no_pvar_match)
if (unaccounted != 0) {
  message("Note: ", unaccounted, " rows unaccounted for between rows_in and rows_out + actual removals ",
          "(no_gene_symbol and duplicate_gene_snp are informational only -- those rows remain in rows_out)")
}

message("\nExport complete: ", nrow(drop_summary), " cell types processed")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
