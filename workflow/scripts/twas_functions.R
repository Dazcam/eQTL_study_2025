#--------------------------------------------------------------------------------------
#
#    Functions for the TWAS weights pipeline (11TWAS)
#
#--------------------------------------------------------------------------------------
#
# Functions: to_plink_pheno()        Transpose a feature x sample table to PLINK FID/IID
#                                    format
#            classify_weights_log()  Per-gene outcome from the compute_weights log
#            read_fusion_wgt()       One-row summary of a FUSION .wgt.RDat file
#
# Packages:  dplyr, tidyr, tibble, stringr (called with ::, so they must be installed
#            but not loaded)
#
# Used in:   twas_prep_expr_for_FUSION.R, twas_prep_covar_for_FUSION.R,
#            twas_weights_summary.Rmd
#
#--------------------------------------------------------------------------------------

## to_plink_pheno  --------------------------------------------------------------------
# Transpose a table with one row per feature (gene or covariate) and one column per
# sample into PLINK phenotype/covariate format: one row per sample, with FID (set
# to 0) and IID first, then one column per feature.
#
# tbl     Tibble: an ID column plus one column per sample, nothing else
# id_col  Name of the ID column (e.g. "TargetID", "covariate")
#
# Returns a tibble with columns FID, IID, <feature IDs>, samples in input column order
to_plink_pheno <- function(tbl, id_col) {
  tbl |>
    tidyr::pivot_longer(cols = -dplyr::all_of(id_col), names_to = "IID", values_to = "value") |>
    tidyr::pivot_wider(names_from = dplyr::all_of(id_col), values_from = "value") |>
    dplyr::mutate(FID = 0) |>
    dplyr::select(FID, IID, dplyr::everything())
}

## classify_weights_log  --------------------------------------------------------------
# Split the compute_weights log into one block per gene (each starts with the rule's
# "Processing <gene>" line) and classify what FUSION did for that gene:
#   hsq_skip          FUSION skipped the gene: heritability < 0 or LRT p > --hsq_p
#   gcta_failed       FUSION skipped the gene: GCTA did not converge
#   no_fusion_output  No FUSION output in the block. The rule only runs FUSION when
#                     plink2 finds cis SNPs, so this means no cis SNPs (inferred)
#   fusion_ran        FUSION wrote output and no skip message was found
# The rule's own "... removed" / "WARNING: Failed to remove" lines are ignored.
#
# log_lines  Character vector, the lines of one compute_weights log
#
# Returns a tibble with columns gene, log_status (one row per gene in the log)
classify_weights_log <- function(log_lines) {
  gene_start <- stringr::str_detect(log_lines, "^Processing \\S+$")
  if (!any(gene_start)) {
    return(tibble::tibble(gene = character(), log_status = character()))
  }

  block_id <- cumsum(gene_start)
  keep <- block_id > 0 &
    !stringr::str_detect(log_lines, "^\\S+ removed$|^WARNING: Failed to remove \\S+$|^All genes completed$")

  tibble::tibble(line = log_lines[keep], block = block_id[keep]) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      gene = stringr::str_remove(dplyr::first(line), "^Processing "),
      fusion_lines = list(line[-1]),
      .groups = "drop"
    ) |>
    dplyr::mutate(log_status = vapply(fusion_lines, function(x) {
      if (length(x) == 0) return("no_fusion_output")
      if (any(stringr::str_detect(x, ": heritability .*skipping gene"))) return("hsq_skip")
      if (any(stringr::str_detect(x, "GCTA could not converge"))) return("gcta_failed")
      "fusion_ran"
    }, character(1))) |>
    dplyr::select(gene, log_status)
}

## read_fusion_wgt  -------------------------------------------------------------------
# Load one FUSION .wgt.RDat file (objects saved by FUSION.compute_weights.R:
# wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot) and summarise it.
#
# file    Path to a .wgt.RDat file
# models  Model names to report CV R2 and p-value for (e.g. c("top1", "lasso", "enet"))
#
# Returns a one-row tibble: gene, n_tot, hsq, hsq_se, hsq_p, n_snps_cis,
#   n_snps_lasso (non-zero lasso weights, NA if no lasso column),
#   best_model (smallest CV p-value), and rsq_<model> / pval_<model> per model
#   (NA if that model was not run)
read_fusion_wgt <- function(file, models) {
  env <- new.env()
  load(file, envir = env)

  wgt <- env$wgt.matrix
  cv <- env$cv.performance

  get_cv <- function(stat, model) {
    if (model %in% colnames(cv)) unname(cv[stat, model]) else NA_real_
  }

  pvals <- cv["pval", , drop = TRUE]
  best_model <- if (all(is.na(pvals))) NA_character_ else colnames(cv)[which.min(pvals)]

  n_snps_lasso <- if ("lasso" %in% colnames(wgt)) {
    lasso <- wgt[, "lasso"]
    if (all(is.na(lasso))) NA_real_ else sum(lasso != 0, na.rm = TRUE)
  } else {
    NA_real_
  }

  out <- tibble::tibble(
    gene = stringr::str_remove(basename(file), "\\.wgt\\.RDat$"),
    n_tot = env$N.tot,
    hsq = env$hsq[1],
    hsq_se = env$hsq[2],
    hsq_p = env$hsq.pv,
    n_snps_cis = nrow(wgt),
    n_snps_lasso = n_snps_lasso,
    best_model = best_model
  )
  for (m in models) {
    out[[paste0("rsq_", m)]] <- get_cv("rsq", m)
    out[[paste0("pval_", m)]] <- get_cv("pval", m)
  }
  out
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
