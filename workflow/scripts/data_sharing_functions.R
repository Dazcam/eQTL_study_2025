#--------------------------------------------------------------------------------------
#
#    Data sharing helper functions
#
#--------------------------------------------------------------------------------------
#
# Functions: read_gene_lookup()     Ensembl ID -> gene symbol, one row per Ensembl ID
#            read_pvar()            Genotype pvar alleles, duplicate SNP IDs dropped
#            annotate_eqtl()        Add alleles and gene symbols to an eQTL table
#            summarise_row_drops()  Per-cell-type drop counts plus a TOTAL row
#
# Packages:  dplyr, readr, stringr, tibble (tidyverse in the calling script)
#
# Used in:   data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R (all);
#            data_sharing_susie_cred_sets.R (read_gene_lookup)
#
#--------------------------------------------------------------------------------------

#' Read the BioMart gene lookup: one row per Ensembl ID, first listed symbol kept.
#' Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R,
#' data_sharing_susie_cred_sets.R
#'
#' @param gene_lookup_file TSV with ensembl_gene_id and external_gene_name columns.
#' @return Tibble with ensembl_gene_id and external_gene_name (blank symbols as NA).
read_gene_lookup <- function(gene_lookup_file) {

  gene_lookup_raw <- readr::read_tsv(gene_lookup_file, show_col_types = FALSE) |>
    dplyr::select(ensembl_gene_id, external_gene_name) |>
    dplyr::distinct()

  n_multi <- gene_lookup_raw |>
    dplyr::count(ensembl_gene_id) |>
    dplyr::filter(n > 1) |>
    nrow()
  message("Gene lookup: ", n_multi, " Ensembl IDs with more than one distinct symbol; ",
          "keeping the first listed for each")

  gene_lookup_raw |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) |>
    dplyr::mutate(external_gene_name = dplyr::na_if(external_gene_name, ""))
}

#' Read pvar alleles, dropping SNP IDs of '.' and keeping the first row per SNP ID.
#' Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R
#'
#' @param pvar_file plink2 pvar (#CHROM POS ID REF ALT INFO).
#' @return List: pvar (CHROM, POS, ID, REF, ALT; one row per ID) and dup_snps
#'   (ID and n for SNP IDs listed more than once).
read_pvar <- function(pvar_file) {

  pvar <- readr::read_tsv(pvar_file, comment = "#", show_col_types = FALSE,
                          col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")) |>
    dplyr::select(CHROM, POS, ID, REF, ALT) |>
    dplyr::filter(ID != ".")

  dup_snps <- pvar |>
    dplyr::count(ID) |>
    dplyr::filter(n > 1)
  message("Duplicate SNP IDs in pvar: ", nrow(dup_snps), " (first row kept)")

  list(pvar = dplyr::distinct(pvar, ID, .keep_all = TRUE), dup_snps = dup_snps)
}

#' Add alleles and gene symbols to a tensorQTL eQTL table.
#' Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R
#'
#' Rows with SNP ID '.' or no pvar match are removed; rows with no gene symbol
#' and duplicate gene-SNP pairs are counted but kept.
#'
#' @param eqtl_tbl Tibble with ensembl_id, SNP, af, slope, slope_se, pval_nominal
#'   and any columns named in extra_cols.
#' @param pvar pvar element of read_pvar().
#' @param gene_lookup Output of read_gene_lookup().
#' @param extra_cols Columns appended after pval_nominal (e.g. pval_beta, qval).
#' @return List: tbl (ensembl_id, symbol, CHROM, SNP, POS, REF, ALT, AF, slope,
#'   slope_se, pval_nominal, extra_cols) and counts (rows_in, dot_id,
#'   no_pvar_match, no_gene_symbol, duplicate_gene_snp, rows_out).
annotate_eqtl <- function(eqtl_tbl, pvar, gene_lookup, extra_cols = character()) {

  n_in  <- nrow(eqtl_tbl)
  n_dot <- sum(eqtl_tbl$SNP == ".")
  if (n_dot > 0) message("SNP ID = '.': ", n_dot, " rows removed")

  with_pvar <- eqtl_tbl |>
    dplyr::filter(SNP != ".") |>
    dplyr::left_join(pvar, by = c("SNP" = "ID"))
  n_no_pvar <- sum(is.na(with_pvar$CHROM))
  if (n_no_pvar > 0) message("No match in pvar file: ", n_no_pvar, " rows removed")

  out_tbl <- with_pvar |>
    dplyr::filter(!is.na(CHROM)) |>
    dplyr::mutate(CHROM = ifelse(!stringr::str_detect(CHROM, "chr"), paste0("chr", CHROM), CHROM)) |>
    dplyr::left_join(gene_lookup, by = dplyr::join_by(ensembl_id == ensembl_gene_id)) |>
    dplyr::mutate(CHROM = stringr::str_remove(CHROM, "^chr")) |>
    dplyr::select(ensembl_id, symbol = external_gene_name, CHROM, SNP, POS, REF, ALT,
                  AF = af, slope, slope_se, pval_nominal, dplyr::all_of(extra_cols))

  n_no_gene <- sum(is.na(out_tbl$symbol))
  if (n_no_gene > 0) message("No gene symbol match: ", n_no_gene, " rows kept with symbol = NA")

  dup_pairs   <- out_tbl |> dplyr::count(ensembl_id, SNP) |> dplyr::filter(n > 1)
  n_dup_extra <- sum(dup_pairs$n - 1L)
  if (nrow(dup_pairs) > 0) {
    message("Duplicate gene-SNP pairs: ", nrow(dup_pairs), " pairs, ", n_dup_extra,
            " extra rows (kept in output). Example: ", dup_pairs$ensembl_id[1], " / ",
            dup_pairs$SNP[1])
  }

  counts <- tibble::tibble(rows_in = n_in, dot_id = n_dot, no_pvar_match = n_no_pvar,
                           no_gene_symbol = n_no_gene, duplicate_gene_snp = n_dup_extra,
                           rows_out = nrow(out_tbl))
  list(tbl = out_tbl, counts = counts)
}

#' Add a TOTAL row to the per-cell-type drop counts and log one line per row.
#' Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R
#'
#' dot_id and no_pvar_match rows are removed; no_gene_symbol and
#' duplicate_gene_snp rows are kept in rows_out.
#'
#' @param drop_summary Tibble: cell_type plus the counts columns of annotate_eqtl().
#' @return drop_summary with a TOTAL row appended.
summarise_row_drops <- function(drop_summary) {

  count_cols <- c("rows_in", "dot_id", "no_pvar_match", "no_gene_symbol",
                  "duplicate_gene_snp", "rows_out")
  totals <- drop_summary |>
    dplyr::summarise(dplyr::across(dplyr::all_of(count_cols), sum)) |>
    dplyr::mutate(cell_type = "TOTAL", .before = 1)
  drop_summary_full <- dplyr::bind_rows(drop_summary, totals)

  message("\nRow counts by reason (dot_id, no_pvar removed; no_gene, dup_pairs kept)")
  for (i in seq_len(nrow(drop_summary_full))) {
    r <- drop_summary_full[i, ]
    message(sprintf(paste0("  %-25s in = %8d  dot_id = %6d  no_pvar = %6d  ",
                           "no_gene = %6d  dup_pairs = %6d  out = %8d"),
                    r$cell_type, r$rows_in, r$dot_id, r$no_pvar_match,
                    r$no_gene_symbol, r$duplicate_gene_snp, r$rows_out))
  }

  unaccounted <- totals$rows_in - (totals$rows_out + totals$dot_id + totals$no_pvar_match)
  if (unaccounted != 0) {
    message("Note: ", unaccounted, " rows unaccounted for between rows_in and ",
            "rows_out + removals")
  }

  drop_summary_full
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
