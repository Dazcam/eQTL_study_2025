#--------------------------------------------------------------------------------------
#
#    Causal TWAS - helper functions
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  12CTWAS | Sourced by: ctwas_run.R
#
# Purpose:   Helper functions for preparing cTWAS inputs and checking weights:
#            - clean_gwas_for_ctwas()  Drop non-rsIDs, duplicate rsIDs, alt contigs
#            - align_gwas_to_ref()     Align GWAS alleles to the reference .bim and
#                                      return the cTWAS z_snp table
#            - make_ld_metatable()     LD matrix / SNP info file paths per LD region
#            - qc_ctwas_weights()      Gene-by-gene weight QC; used to identify the
#                                      failing cell type / GWAS combinations
#
# Notes:     Requires tidyverse to be loaded by the calling script.
#
#--------------------------------------------------------------------------------------

## GWAS cleaning  -----------------------------------------------------------------------
# Removes non-rsID SNPs, duplicate rsIDs (keeps first) and SNPs on 'alt' contigs.
# Returns the cleaned GWAS table.
clean_gwas_for_ctwas <- function(gwas_tbl) {

  # non-rsIDs
  non_rsids <- gwas_tbl |>
    filter(!str_starts(SNP, "rs")) |>
    distinct(SNP) |>
    nrow()
  if (non_rsids > 0) {
    message(str_glue("Found {non_rsids} non-rsID SNPs. Removing them ..."))
    gwas_tbl <- gwas_tbl |> filter(str_starts(SNP, "rs"))
  } else {
    message("No non-rsID SNPs found.")
  }

  # duplicate rsIDs
  duplicate_rsids <- gwas_tbl |>
    group_by(SNP) |>
    filter(n() > 1) |>
    distinct(SNP) |>
    nrow()
  if (duplicate_rsids > 0) {
    message(str_glue("Found {duplicate_rsids} duplicate rsIDs. Keeping first occurrence ..."))
    gwas_tbl <- gwas_tbl |> distinct(SNP, .keep_all = TRUE)
  } else {
    message("No duplicate rsIDs found.")
  }

  # Alt Chrs
  message("Checking GWAS chromosomes ...\n")
  gwas_tbl |>
    group_by(CHR) |>
    count() |>
    print(n = Inf)

  alt_chroms <- gwas_tbl |>
    filter(str_ends(CHR, "alt")) |>
    distinct(CHR) |>
    pull(CHR)

  if (length(alt_chroms) > 0) {
    message("\nChromosomes ending with 'alt' found: ", paste(alt_chroms, collapse = ", "))
    message("\nCount of SNPs with 'alt' chromosomes: ",
            nrow(gwas_tbl |> filter(str_ends(CHR, "alt"))))
    message("Removing SNPs with 'alt' chromosomes ... ")
    gwas_tbl <- gwas_tbl |> filter(!str_ends(CHR, "alt"))
  } else {
    message("No chromosomes ending with 'alt' found.")
  }

  gwas_tbl
}


## GWAS allele alignment  ---------------------------------------------------------------
# Joins the reference .bim alleles onto the GWAS, classifies each SNP as aligned /
# flipped / mismatched, flips Z (and alleles) for flipped SNPs, and returns the
# cTWAS z_snp table (id, A1, A2, z) with incomplete rows removed.
align_gwas_to_ref <- function(gwas_tbl, snp_tbl) {

  message("\nAppending bim ref allele info to GWAS for alignment check...\n")
  mrg_tbl <- gwas_tbl |>
    left_join(snp_tbl, by = c("SNP"), suffix = c("", ".ref"))

  mrg_tbl <- mrg_tbl |>
    mutate(
    # classify allele relationships
    allele_status = case_when(
      A1 == A1.ref & A2 == A2.ref ~ "aligned",
      A1 == A2.ref & A2 == A1.ref ~ "flipped",
      TRUE                       ~ "mismatched"
    ),

    Z = case_when(
      allele_status == "flipped" ~ -Z,
      TRUE                       ~ Z
    ),

    # Swap alleles via new columns: mutate() evaluates in order, so assigning A1
    # first and then A2 = A1 would copy the already-swapped A1 into A2
    A1_new = if_else(allele_status == "flipped", A2, A1),
    A2_new = if_else(allele_status == "flipped", A1, A2)
  ) |>
    dplyr::select(-A1, -A2) |>
    dplyr::rename(A1 = A1_new, A2 = A2_new)

  # Report allele flips
  message("\n=== GWAS allele alignment summary ===")
  mrg_tbl |>
    count(allele_status) |>
    mutate(percent = round(100 * n / sum(n), 1)) |>
    arrange(desc(n)) |>
    print(n = Inf)
  message("========================================\n")

  z_snp <- mrg_tbl |>
    select(id = SNP, A1, A2, z = Z) |>
    filter(if_all(everything(), ~ !is.na(.) & . != ""))

  # Check for any remaining missing or blank values
  if (anyNA(z_snp) || any(z_snp == "", na.rm = TRUE)) {
    stop("Missing values found in GWAS data.")
  } else {
    message("No missing values in GWAS data ...")
  }

  message("Final number of SNPs in GWAS: ", nrow(z_snp))

  z_snp
}


## LD region metadata  ------------------------------------------------------------------
# One row per LD region, with paths to the LD matrix (.RDS) and SNP info (.Rvar)
# files written by ctwas_create_ld_matrices.R (outname = prefix).
make_ld_metatable <- function(region_info, ld_dir, prefix) {

  message('Creating metadata table for LD files ...')
  region_metatable <- region_info
  LD_filestem <- sprintf("%s_chr%s.R_snp.%s_%s",
                         prefix,
                         region_metatable$chrom,
                         region_metatable$start,
                         region_metatable$stop)
  region_metatable$LD_file <- file.path(ld_dir, paste0(LD_filestem, ".RDS"))
  region_metatable$SNP_file <- file.path(ld_dir, paste0(LD_filestem, ".Rvar"))

  region_metatable
}


## Weight QC  ---------------------------------------------------------------------------
# Gene-by-gene check that each processed weight can be used by cTWAS: structure,
# dimensions, SNP presence in z_snp, perfect LD, and the compute_gene_z calculation.
# Returns a tibble: gene, status (PASS / FAIL), reason.
qc_ctwas_weights <- function(weights, z_snp, ld_threshold = 0.9999) {

  results <- vector("list", length(weights))
  names(results) <- names(weights)
  counter <- 1
  total   <- length(weights)

  for (gene_id in names(weights)) {
    g <- weights[[gene_id]]
    cat("Checking gene", counter, "of", total, ":", gene_id, "...\n")
    counter <- counter + 1

    # ----------------------------------------------------
    # 1. Basic structure checks
    # ----------------------------------------------------
    if (is.null(g$wgt) || !is.matrix(g$wgt)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "wgt_not_matrix")
      next
    }
    if (is.null(g$R_wgt)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "R_wgt_NULL")
      next
    }
    if (!is.matrix(g$R_wgt)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "R_wgt_not_matrix")
      next
    }

    wgt <- as.matrix(g$wgt)
    R   <- as.matrix(g$R_wgt)

    # Must be numeric
    if (!is.numeric(wgt) || !is.numeric(R)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "non_numeric")
      next
    }

    # ----------------------------------------------------
    # 2. Dimensionality checks
    # ----------------------------------------------------
    if (nrow(R) != ncol(R)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "R_wgt_not_square")
      next
    }
    if (nrow(wgt) != nrow(R)) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "dim_mismatch")
      next
    }

    # ----------------------------------------------------
    # 3. SNP matching
    # ----------------------------------------------------
    snps <- rownames(wgt)
    idx  <- match(snps, z_snp$id)
    if (any(is.na(idx))) {
      results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "SNP_missing_in_z")
      next
    }
    zvec <- z_snp$z[idx]

    # ----------------------------------------------------
    # 4. Perfect LD detection
    # ----------------------------------------------------
    if (nrow(R) > 1) {
      offdiag <- abs(R[upper.tri(R)])
      if (any(offdiag > ld_threshold)) {
        results[[gene_id]] <- list(id = gene_id, status = "FAIL", reason = "perfect_LD")
        next
      }
    }

    # ----------------------------------------------------
    # 5. Run the exact compute_gene_z math
    # ----------------------------------------------------
    res <- tryCatch({
      storage.mode(wgt) <- "double"
      storage.mode(R) <- "double"
      storage.mode(zvec) <- "double"

      denom <- sqrt(as.numeric(t(wgt) %*% R %*% wgt))
      num   <- as.numeric(crossprod(wgt, zvec))
      z     <- num / denom

      if (!is.finite(z)) stop("computed z is not finite")

      list(id = gene_id, status = "PASS", z = z)
    }, error = function(e) {
      list(id = gene_id, status = "FAIL", reason = paste("compute_error:", conditionMessage(e)))
    })

    results[[gene_id]] <- res
  }

  # Convert list → tibble
  qc_tbl <- tibble::tibble(
    gene   = names(results),
    status = sapply(results, function(x) x$status),
    reason = sapply(results, function(x) if ("reason" %in% names(x)) x$reason else NA_character_)
  )

  return(qc_tbl)
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
