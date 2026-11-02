#--------------------------------------------------------------------------------------
#
#    Functions for the single-cell replication pipeline (also used by eQTL sharing)
#
#--------------------------------------------------------------------------------------
#
# Functions: compute_pi1()              pi1 replication statistic (qvalue)
#            filter_jang_artifacts()    drop Jang placeholder-cohort rows
#            align_jang_beta()          Jang beta relative to Jang's alt allele
#            get_allele_status()        compare our REF/ALT with a reference's
#            harmonise_beta()           reference beta relative to our ALT
#            make_snp_gr()              SNP coordinates -> GRanges
#            perform_permutation_test() block-bootstrap peak overlap enrichment
#
# Packages:  the caller must have qvalue (compute_pi1) and GenomicRanges,
#            GenomeInfoDb and nullranges (make_snp_gr, perform_permutation_test)
#            installed; the other functions use base R only.
#
# Used in:   repl_sc_pi1_enrichment.R    compute_pi1, filter_jang_artifacts,
#                                        align_jang_beta
#            repl_sc_beta_correlation.R  filter_jang_artifacts, align_jang_beta,
#                                        get_allele_status, harmonise_beta
#            repl_sc_atac_enrichment.R   make_snp_gr, perform_permutation_test
#            eqtl_sharing_classify.R  filter_jang_artifacts, align_jang_beta
#            eqtl_sharing_effect_sizes.R     filter_jang_artifacts, align_jang_beta
#            eqtl_sharing_effect_replication.R    filter_jang_artifacts, align_jang_beta,
#                                            get_allele_status, harmonise_beta
#            eqtl_sharing_heterogeneity.R    get_allele_status, harmonise_beta
#            eqtl_sharing_ocr_enrichment.R        make_snp_gr, perform_permutation_test
#            eqtl_sharing_glm_enrichment.R        make_snp_gr
#
#--------------------------------------------------------------------------------------

# pi1 = 1 - pi0 for the reference p-values of the query's significant pairs.
# query_eqtl: variant_id, phenotype_id, slope_my
# full_eqtl:  variant_id, phenotype_id, pval, slope_ref
compute_pi1 <- function(query_eqtl, full_eqtl, ref_name, min_overlap = 50) {
  message("Computing pi1 statistic for ", ref_name)

  merged <- dplyr::left_join(query_eqtl, full_eqtl, by = c("variant_id", "phenotype_id"))
  valid_merged <- dplyr::filter(merged, !is.na(pval), is.finite(pval), pval >= 0, pval <= 1)
  valid_pvals <- valid_merged$pval
  overlap_count <- length(valid_pvals)

  message("Number of overlapping eQTLs: ", overlap_count)

  if (overlap_count < min_overlap) {
    message("Too few overlapping eQTLs (<", min_overlap, ") for pi1 calculation.")
    return(list(pi1 = NA, prop_replicating = NA,
                prop_same_direction = NA,
                overlap_count = overlap_count,
                pvals = valid_pvals))
  }

  qobj <- tryCatch(
    qvalue::qvalue(p = valid_pvals, pi0.method = "smoother"),
    error = function(e) {
      message("Smoother failed, trying bootstrap ...")
      tryCatch(
        qvalue::qvalue(p = valid_pvals, pi0.method = "bootstrap"),
        error = function(e2) {
          message("Bootstrap also failed: ", e2$message)
          NULL
        }
      )
    }
  )

  if (!is.null(qobj)) {
    pi1 <- 1 - qobj$pi0
    prop_replicating <- mean(qobj$qvalues < 0.05)
    rep_idx <- which(qobj$qvalues < 0.05)
    replicating <- valid_merged[rep_idx, ]
    prop_same_direction <- if (nrow(replicating) > 0) {
      mean(sign(replicating$slope_my) == sign(replicating$slope_ref), na.rm = TRUE)
    } else {
      NA
    }
    message("pi1 calculated: ", round(pi1, 3))
  } else {
    pi1 <- NA
    prop_replicating <- NA
    prop_same_direction <- NA
  }

  list(pi1 = pi1,
       prop_replicating = prop_replicating,
       prop_same_direction = prop_same_direction,
       overlap_count = overlap_count,
       pvals = valid_pvals)
}

# Jang et al. 2026 placeholder cohorts encode "not tested" as beta = 0 / sd = 0,
# which forces the meta-analysis to Fixed_P == Random_P == 1 exactly. Those rows
# are removed (rows where either P is NA are dropped too, as dplyr::filter() and
# data.table's i do). Works on data frames, tibbles and data.tables.
filter_jang_artifacts <- function(df, label = "Jang") {
  keep <- !(df$Fixed_P == 1 & df$Random_P == 1)
  keep[is.na(keep)] <- FALSE
  n_before <- nrow(df)
  out <- df[keep, ]
  n_removed <- n_before - nrow(out)
  message("  ", label, ": removed ", n_removed,
          " placeholder-cohort artifact rows (Fixed_P==1 & Random_P==1 exactly) out of ",
          n_before, " (", round(100 * n_removed / n_before, 3), "%)")
  out
}

# Jang betas are relative to the `Allele` column, not always `alt`.
# ref = NULL: beta where Allele == alt, otherwise -beta.
# ref given:  beta where Allele == alt, -beta where Allele == ref, otherwise NA.
align_jang_beta <- function(beta, allele, alt, ref = NULL) {
  if (is.null(ref)) {
    return(ifelse(allele == alt, beta, -beta))
  }
  ifelse(allele == alt, beta, ifelse(allele == ref, -beta, NA_real_))
}

# Compare our REF/ALT (pvar) with a reference dataset's ref/alt. Palindromic
# (A/T, C/G) SNPs are never resolved, as the reference has no allele
# frequencies. If beta_ref is given, a missing beta also counts as missing.
get_allele_status <- function(ref_my, alt_my, ref_ref, alt_ref, beta_ref = NULL) {
  missing <- is.na(ref_my) | is.na(alt_my) | is.na(ref_ref) | is.na(alt_ref)
  if (!is.null(beta_ref)) missing <- missing | is.na(beta_ref)
  palindromic <- (ref_my == "A" & alt_my == "T") | (ref_my == "T" & alt_my == "A") |
                 (ref_my == "C" & alt_my == "G") | (ref_my == "G" & alt_my == "C")
  aligned <- ref_my == ref_ref & alt_my == alt_ref
  flipped <- ref_my == alt_ref & alt_my == ref_ref

  status <- rep("mismatched", length(ref_my))
  status[which(flipped)]     <- "flipped"
  status[which(aligned)]     <- "aligned"
  status[which(palindromic)] <- "palindromic_unresolved"
  status[which(missing)]     <- "missing_allele_info"
  status
}

# Reference beta relative to our ALT: kept if aligned, negated if flipped,
# NA otherwise.
harmonise_beta <- function(beta_ref, allele_status) {
  ifelse(allele_status == "aligned", beta_ref,
         ifelse(allele_status == "flipped", -beta_ref, NA_real_))
}

# SNP coordinates (biomaRt lookup: refsnp_id, chr_name, chrom_start, chrom_end)
# for the given SNP IDs, one range per SNP.
make_snp_gr <- function(snp_ids, coords) {
  coords <- coords[coords$refsnp_id %in% snp_ids, ]
  coords <- coords[!duplicated(coords$refsnp_id), ]
  GenomicRanges::makeGRangesFromDataFrame(
    coords,
    seqnames.field = "chr_name",
    start.field    = "chrom_start",
    end.field      = "chrom_end",
    keep.extra.columns = TRUE
  )
}

# Overlap of SNPs with peaks against a block-bootstrap null (nullranges::bootRanges,
# 500 kb blocks within chromosomes). hg38 chromosome sizes come from UCSC, with
# hard-coded sizes as a fallback if the download fails.
perform_permutation_test <- function(sig_snps_gr, peaks_gr, n_perm = 1000) {
  message("Starting permutation test for ATAC-seq peak enrichment...")

  if (!exists("bootRanges", envir = asNamespace("nullranges"))) {
    stop("bootRanges function not found. Ensure nullranges is installed and loaded.")
  }

  hg38_seqinfo <- tryCatch({
    message("Fetching hg38 seqinfo from UCSC ...")
    hg38_info <- GenomeInfoDb::getChromInfoFromUCSC("hg38")
    GenomeInfoDb::Seqinfo(seqnames = hg38_info$chrom,
                          seqlengths = hg38_info$size,
                          genome = "hg38")
  }, error = function(e) {
    message("UCSC fetch failed (", e$message, "), using hardcoded hg38 chromosome sizes.")
    hg38_seqlengths <- c(
      chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
      chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
      chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
      chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
      chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
      chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415
    )
    GenomeInfoDb::Seqinfo(seqnames = names(hg38_seqlengths),
                          seqlengths = hg38_seqlengths,
                          genome = "hg38")
  })

  GenomeInfoDb::seqinfo(peaks_gr) <- hg38_seqinfo[GenomeInfoDb::seqlevels(peaks_gr)]

  observed_overlap <- length(IRanges::findOverlaps(sig_snps_gr, peaks_gr))
  message("Observed overlaps: ", observed_overlap)

  block_length <- 500000

  perm_overlaps <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    if (i %% 100 == 0) message("Permutation ", i, " of ", n_perm)
    boot_peaks <- nullranges::bootRanges(
      y = peaks_gr,
      blockLength = block_length,
      R = 1,
      withinChrom = TRUE
    )
    perm_overlaps[i] <- length(IRanges::findOverlaps(sig_snps_gr, boot_peaks))
  }

  p_value <- (sum(perm_overlaps >= observed_overlap) + 1) / (n_perm + 1)
  mean_permuted_overlap <- mean(perm_overlaps)
  fold_enrichment <- observed_overlap / mean_permuted_overlap

  message("Permutation p-value: ", signif(p_value, 3))
  message("Fold enrichment: ", signif(fold_enrichment, 3))

  list(
    p_value = p_value,
    observed_overlap = observed_overlap,
    mean_permuted_overlap = mean_permuted_overlap,
    fold_enrichment = fold_enrichment,
    perm_overlaps = perm_overlaps
  )
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
