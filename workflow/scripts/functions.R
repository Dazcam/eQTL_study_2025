#--------------------------------------------------------------------------------------
#
#    General helper functions (shared across pipelines)
#
#--------------------------------------------------------------------------------------
#
# Functions: normalise_counts()  Normalise / filter a CPM matrix by one of four methods
#            get_mhc_genes()     Gene IDs overlapping the extended MHC (or any region)
#            relabel_cell_types() Pipeline -> manuscript cell type labels
#            cell_type_levels()   Cell types present, in config order (family or level layout)
#
# Packages:  Base R only; limma must be installed (not loaded) for
#            normalise_counts(norm_method = "quantile")
#
# Used in:   tensorqtl_prep_input_files.R (normalise_counts);
#            smr_create_query.R, ctwas_run.R (get_mhc_genes);
#            16 tables and 17 figures scripts (relabel_cell_types, cell_type_levels);
#            18 data sharing R scripts (relabel_cell_types)
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

#' Relabel pipeline cell type names with manuscript names.
#' Used in 16 tables, 17 figures and 18 data sharing scripts
#'
#' Applies each substitution in config manuscript_labels, in order, to every
#' occurrence in each string, so it works on plain cell types
#' ("Glu-UL-0" -> "Glu-A-0") and on strings that contain them
#' ("NPC-to-Glu-UL-Q4-Bin1" -> "NPC-to-Glu-A-Bin1").
#'
#' @param x Character (or factor) vector of labels.
#' @param labels Named list or vector: pipeline text -> manuscript text.
#' @return Character vector of relabelled strings (NA kept as NA).
relabel_cell_types <- function(x, labels) {
  labels <- unlist(labels)
  if (is.null(names(labels)) || any(names(labels) == "")) {
    stop("relabel_cell_types: labels must be a named list or vector")
  }
  x <- as.character(x)
  for (old in names(labels)) {
    x <- gsub(old, labels[[old]], x, fixed = TRUE)
  }
  x
}

#' Order cell types as in config cell_types.
#' Used in 16 tables and 17 figures scripts
#'
#' L2 subclusters are named "<L1>-<n>"; pseudotime bins contain "-to-".
#' layout = "family": each L1 type followed by its L2 subclusters, in config
#'   order, then pseudotime bins (as on the Fig 2 y-axis).
#' layout = "level": all L1, then all L2, then pseudotime bins, each in config
#'   order.
#'
#' @param x Cell types to order (character or factor; duplicates allowed).
#' @param cell_types Config cell_types (pipeline labels).
#' @param labels Optional named list (config manuscript_labels). When given,
#'   cell_types are relabelled first, so x can use manuscript labels.
#' @param layout "family" or "level".
#' @return Unique cell types of x in order, for use as factor levels. Stops if
#'   x holds a cell type not in cell_types.
cell_type_levels <- function(x, cell_types, labels = NULL, layout = c("family", "level")) {
  layout <- match.arg(layout)
  ref <- as.character(unlist(cell_types))
  if (!is.null(labels)) ref <- relabel_cell_types(ref, labels)

  is_pt  <- grepl("-to-", ref, fixed = TRUE)
  parent <- ifelse(is_pt, NA_character_, sub("-[0-9]+$", "", ref))
  is_l2  <- !is_pt & parent != ref & parent %in% ref
  l1 <- ref[!is_pt & !is_l2]
  l2 <- ref[is_l2]
  pt <- ref[is_pt]

  ordered <- if (layout == "level") {
    c(l1, l2, pt)
  } else {
    c(unlist(lapply(l1, function(p) c(p, l2[parent[is_l2] == p])), use.names = FALSE), pt)
  }

  x <- unique(as.character(x))
  unknown <- setdiff(x, ordered)
  if (length(unknown) > 0) {
    stop("cell_type_levels: not in config cell_types: ", paste(unknown, collapse = ", "))
  }
  ordered[ordered %in% x]
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
