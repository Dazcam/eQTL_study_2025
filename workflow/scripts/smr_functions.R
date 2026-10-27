#--------------------------------------------------------------------------------------
#
#    SMR helper functions
#
#--------------------------------------------------------------------------------------
#
# Functions: align_freq_to_ref()  Classify alleles against the reference .frq and set
#                                 the allele frequency so it refers to A1
#
# Packages:  dplyr must be loaded by the calling script
#
# Used in:   smr_create_query.R, smr_format_gwas.R
#
#--------------------------------------------------------------------------------------

## Allele frequency alignment  ---------------------------------------------------------
# Joins the reference .frq (A1, A2, MAF) onto tbl by the `by` columns and classifies
# each SNP as aligned / flipped / mismatched against the reference alleles. Alleles
# and effects in tbl are NOT changed; only the frequency column (`freq_col`) is set,
# so that it is the frequency of tbl's A1:
#   aligned     MAF
#   flipped     1 - MAF
#   mismatched  existing freq_col value if keep_unmatched = TRUE, otherwise NA
# Prints the status counts under `label`. Returns tbl with freq_col and
# allele_status added or updated.
align_freq_to_ref <- function(tbl, frq_tbl, by, freq_col, keep_unmatched = FALSE,
                              label = "") {

  ref_tbl <- frq_tbl |>
    dplyr::select(dplyr::all_of(by), A1.ref = A1, A2.ref = A2, MAF)

  out_tbl <- tbl |>
    left_join(ref_tbl, by = by) |>
    mutate(
      allele_status = case_when(
        A1 == A1.ref & A2 == A2.ref ~ "aligned",
        A1 == A2.ref & A2 == A1.ref ~ "flipped",
        TRUE ~ "mismatched"
      ),
      !!freq_col := case_when(
        allele_status == "aligned" ~ MAF,
        allele_status == "flipped" ~ 1 - MAF,
        TRUE ~ if (keep_unmatched) .data[[freq_col]] else NA_real_
      )
    ) |>
    dplyr::select(-A1.ref, -A2.ref, -MAF)

  # Report allele flips
  message("\n=== ", label, " frequency alignment summary ===")
  out_tbl %>%
    count(allele_status) %>%
    mutate(percent = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n)) %>%
    print(n = Inf)
  message("========================================\n")

  out_tbl
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
