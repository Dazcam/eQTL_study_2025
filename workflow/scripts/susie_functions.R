#--------------------------------------------------------------------------------------
#
#    SuSiE fine-mapping functions
#
#--------------------------------------------------------------------------------------
#
# Functions: importQtlmapCovariates()           read covariates into a sample x covariate matrix
#            splitIntoBatches()                 batch IDs for n items in batches of a given size
#            splitIntoChunks()                  logical selector for one chunk of the gene list
#            finemapPhenotype()                 residualise expression and genotypes, fit SuSiE
#            susie_get_posterior_mean_custom()  posterior mean (fix for a single retained effect)
#            susie_get_posterior_sd_custom()    posterior SD (fix for a single retained effect)
#            extractResults()                   credible set and variant tables from a SuSiE fit
#
# Packages:  dplyr, purrr (and susieR, eQTLUtils for finemapPhenotype())
#
# Used in:   susie_run.R
#
# Notes:     Adapted from the eQTL Catalogue qtlmap run_susie.R. Numerical code is kept
#            as in the original so outputs are unchanged.
#
#--------------------------------------------------------------------------------------

## Read a QTLtools-style covariates file (covariates in rows, samples in columns, first
## column holds covariate names) and return a sample x covariate matrix. Covariates
## with zero variance across samples are dropped.
importQtlmapCovariates <- function(covariates_path) {
  cov_tbl <- read.table(covariates_path, check.names = FALSE, header = TRUE,
                        stringsAsFactors = FALSE)
  cov_tbl <- cov_tbl[apply(cov_tbl[-1], 1, var) != 0, ]
  cov_mat <- t(as.matrix(cov_tbl[, -1]))
  colnames(cov_mat) <- cov_tbl[[1]]
  return(cov_mat)
}

## Batch IDs for n items split into batches of batch_size
splitIntoBatches <- function(n, batch_size) {
  n_batches <- ceiling(n / batch_size)
  batch_ids <- rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

## Logical vector selecting chunk chunk_number of n_chunks from n_total items. Any
## remainder is added to the last chunk.
splitIntoChunks <- function(chunk_number, n_chunks, n_total) {
  chunk_size <- floor(n_total / (n_chunks))
  batches <- splitIntoBatches(n_total, chunk_size)
  batches[batches > n_chunks] <- n_chunks
  selected_batch <- batches == chunk_number
  return(selected_batch)
}

## Fine-map one gene. Expression (already normalised upstream) and genotypes are
## residualised on the covariates plus an intercept; genotypes are centred, not scaled.
##   expr_mat:  gene x sample expression matrix (rownames = phenotype_id)
##   gene_meta: data frame with phenotype_id, chromosome and phenotype_pos (TSS)
finemapPhenotype <- function(phenotype_id, expr_mat, gene_meta, genotype_file,
                             covariates, cis_distance) {
  message("Processing phenotype: ", phenotype_id)

  expression <- expr_mat[phenotype_id, ]
  selected_phenotype <- phenotype_id
  gene_row <- dplyr::filter(gene_meta, phenotype_id == selected_phenotype)

  # Covariates in the same sample order as the expression, plus an intercept
  covariates_matrix <- cbind(covariates[names(expression), ], 1)

  # Genotypes within cis_distance of the TSS
  genotype_matrix <- eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_row$chromosome,
    start = gene_row$phenotype_pos - cis_distance,
    end = gene_row$phenotype_pos + cis_distance,
    dosage_file = genotype_file)

  # Residualise expression on the covariates
  hat <- diag(nrow(covariates_matrix)) - covariates_matrix %*%
    solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)
  expression_vector <- hat %*% expression
  names(expression_vector) <- names(expression)

  gt_matrix <- genotype_matrix[, names(expression_vector)]

  # Exclude variants with no alternative alleles
  gt_matrix <- gt_matrix[rowSums(round(gt_matrix, 0)) != 0, ]

  # Centre and residualise genotypes
  gt_std <- t(gt_matrix - apply(gt_matrix, 1, mean))
  gt_hat <- hat %*% gt_std

  fitted <- susieR::susie(gt_hat, expression_vector,
                          L = 10,
                          estimate_residual_variance = TRUE,
                          estimate_prior_variance = TRUE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE,
                          compute_univariate_zscore = TRUE,
                          min_abs_corr = 0)
  fitted$variant_id <- rownames(gt_matrix)
  return(fitted)
}

## susieR::susie_get_posterior_mean/sd fail when only one effect is retained
## (include_idx of length 1); these versions keep the matrix with drop = FALSE.
susie_get_posterior_mean_custom <- function(res, prior_tol = 1e-9) {
  if (is.numeric(res$V))
    include_idx <- which(res$V > prior_tol)
  else
    include_idx <- 1:nrow(res$alpha)

  if (length(include_idx) > 0)
    return(colSums((res$alpha * res$mu)[include_idx, , drop = FALSE]) /
             res$X_column_scale_factors)
  else
    return(numeric(ncol(res$mu)))
}

susie_get_posterior_sd_custom <- function(res, prior_tol = 1e-9) {
  if (is.numeric(res$V))
    include_idx <- which(res$V > prior_tol)
  else
    include_idx <- 1:nrow(res$alpha)

  if (length(include_idx) > 0)
    return(sqrt(colSums((res$alpha * res$mu2 -
                           (res$alpha * res$mu)^2)[include_idx, , drop = FALSE])) /
             (res$X_column_scale_factors))
  else
    return(numeric(ncol(res$mu)))
}

## Credible set (cs_df) and variant (variant_df) tables from one SuSiE fit.
## Credible sets are taken in SuSiE's order (highest purity first); a set sharing any
## variant with an earlier kept set is dropped. low_purity flags min.abs.corr < 0.5.
extractResults <- function(susie_object) {
  credible_sets <- susie_object$sets$cs
  cs_list <- list()
  susie_object$sets$purity <- dplyr::as_tibble(susie_object$sets$purity) |>
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity),
      cs_size = NA,
      cs_log10bf = NA,
      overlapped = NA)
  added_variants <- c()

  for (index in seq_along(credible_sets)) {
    cs_variants <- credible_sets[[index]]
    cs_id <- susie_object$sets$cs_index[[index]]

    is_overlapped <- any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] <- is_overlapped
    susie_object$sets$purity$cs_size[index] <- length(cs_variants)
    susie_object$sets$purity$cs_log10bf[index] <- log10(exp(susie_object$lbf[cs_id]))
    if (!is_overlapped) {
      cs_list[[index]] <- dplyr::tibble(cs_id = paste0("L", cs_id),
                                        variant_id = susie_object$variant_id[cs_variants])
      added_variants <- append(added_variants, cs_variants)
    }
  }
  df <- purrr::map_df(cs_list, identity)

  # Purity summaries for kept credible sets (none when all PIPs are 0)
  purity_res <- susie_object$sets$purity
  if (nrow(purity_res) > 0) {
    purity_df <- dplyr::as_tibble(purity_res) |>
      dplyr::filter(!overlapped) |>
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2,
        cs_min_r2 = min.abs.corr^2,
        low_purity = min.abs.corr < 0.5) |>
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity)
  } else {
    purity_df <- dplyr::tibble()
  }

  # Posterior means and SDs, overall and per effect
  mean_vec <- susie_get_posterior_mean_custom(susie_object)
  sd_vec <- susie_get_posterior_sd_custom(susie_object)
  alpha_mat <- t(susie_object$alpha)
  colnames(alpha_mat) <- paste0("alpha", seq(ncol(alpha_mat)))
  mean_mat <- t(susie_object$alpha * susie_object$mu) / susie_object$X_column_scale_factors
  colnames(mean_mat) <- paste0("mean", seq(ncol(mean_mat)))
  sd_mat <- sqrt(t(susie_object$alpha * susie_object$mu2 -
                     (susie_object$alpha * susie_object$mu)^2)) /
    (susie_object$X_column_scale_factors)
  colnames(sd_mat) <- paste0("sd", seq(ncol(sd_mat)))
  posterior_df <- dplyr::tibble(variant_id = names(mean_vec),
                                pip = susie_object$pip,
                                z = susie_object$z,
                                posterior_mean = mean_vec,
                                posterior_sd = sd_vec) |>
    dplyr::bind_cols(purrr::map(list(alpha_mat, mean_mat, sd_mat), dplyr::as_tibble))

  if (nrow(df) > 0 & nrow(purity_df) > 0) {
    cs_df <- purity_df
    variant_df <- dplyr::left_join(posterior_df, df, by = "variant_id") |>
      dplyr::left_join(cs_df, by = "cs_id")
  } else {
    cs_df <- NULL
    variant_df <- NULL
  }

  return(list(cs_df = cs_df, variant_df = variant_df))
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
