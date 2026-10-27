#--------------------------------------------------------------------------------------
#
#    General helper functions (shared across pipelines)
#
#--------------------------------------------------------------------------------------
#
# Functions: normalise_counts()  Normalise / filter a CPM matrix by one of four methods
#            get_mhc_genes()     Gene IDs overlapping the extended MHC (or any region)
#
# Packages:  Base R only; limma must be installed (not loaded) for
#            normalise_counts(norm_method = "quantile")
#
# Used in:   tensorqtl_prep_input_files.R (normalise_counts);
#            smr_create_query.R, ctwas_run.R (get_mhc_genes)
#
#--------------------------------------------------------------------------------------

#' Apply one of four normalisation/filtering strategies to a CPM matrix.
#' Used in tensorqtl_prep_input_files.R
#'
#' @param normalised_cnts A genes x samples matrix of TMM-normalised CPM values
#'   (output of edgeR::cpm() upstream in the pipeline).
#' @param norm_method One of "bryois", "fujita", "quantile", "xue".
#' @return The filtered/transformed matrix for the chosen method.
normalise_counts <- function(normalised_cnts, norm_method, cell_type = NA) {

  if (norm_method == "bryois") {
    # Bryois: filter genes with mean CPM < 1; PMID:35915177
    mean_cpm <- rowMeans(normalised_cnts)
    low_expr_genes <- mean_cpm < 1
    if (sum(low_expr_genes) > 0) {
      message('Removed ', sum(low_expr_genes), ' genes with mean CPM < 1 for ', cell_type)
      normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
    }
    print(normalised_cnts[1:5, 1:5])

  } else if (norm_method == "fujita") {
    # Fujita: log2(CPM + 1), filter genes with log2 CPM < 2.0 in all samples; PMID:38514782
    normalised_cnts <- log2(normalised_cnts + 1)
    low_expr_genes <- rowSums(normalised_cnts >= 2.0) == 0
    if (sum(low_expr_genes) > 0) {
      message('Removed ', sum(low_expr_genes), ' genes with log2 CPM < 2.0 in all samples for ', cell_type)
      normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
    }
    print(normalised_cnts[1:5, 1:5])

  } else if (norm_method == "quantile") {
    # Fujita steps + quantile normalisation
    normalised_cnts <- log2(normalised_cnts + 1)
    low_expr_genes <- rowSums(normalised_cnts >= 2.0) == 0
    if (sum(low_expr_genes) > 0) {
      message('Removed ', sum(low_expr_genes), ' genes with log2 CPM < 2.0 for ', cell_type)
      normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
    }
    normalised_cnts <- limma::normalizeQuantiles(normalised_cnts)
    print(normalised_cnts[1:5, 1:5])

  } else if (norm_method == "xue") {
    # Xue Option #11: QC (>=90% zeros), log1p, z-score (STD); PMID:36823676
    zero_prop <- rowMeans(normalised_cnts == 0)
    low_expr_genes <- zero_prop >= 0.9
    if (sum(low_expr_genes) > 0) {
      message('Removed ', sum(low_expr_genes), ' genes with >=90% zeros for ', cell_type)
      normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
    }
    normalised_cnts <- log(normalised_cnts + 1)
    normalised_cnts <- t(scale(t(normalised_cnts)))
    print(normalised_cnts[1:5, 1:5])

  } else {
    stop("Unknown norm_method: ", norm_method)
  }

  normalised_cnts
}


## Extended MHC genes  ------------------------------------------------------------------
# Returns the unique gene IDs in gene_lookup that overlap chr:start-end. A gene counts
# if any part of it overlaps (inclusive)
get_mhc_genes <- function(gene_lookup, chr, start, end,
                          id_col = "ensembl_gene_id", chr_col = "chromosome_name",
                          start_col = "start_position", end_col = "end_position") {

  missing_cols <- setdiff(c(id_col, chr_col, start_col, end_col), names(gene_lookup))
  if (length(missing_cols) > 0) {
    stop("Gene lookup is missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  gene_chr <- sub("^chr", "", as.character(gene_lookup[[chr_col]]))
  gene_st  <- as.numeric(gene_lookup[[start_col]])
  gene_en  <- as.numeric(gene_lookup[[end_col]])
  lo <- pmin(gene_st, gene_en)
  hi <- pmax(gene_st, gene_en)

  in_region <- gene_chr == sub("^chr", "", as.character(chr)) &
    hi >= as.numeric(start) & lo <= as.numeric(end)

  unique(as.character(gene_lookup[[id_col]][in_region %in% TRUE]))
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

