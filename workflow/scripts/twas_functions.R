#--------------------------------------------------------------------------------------
#
#    Functions for the TWAS weights pipeline (11TWAS)
#
#--------------------------------------------------------------------------------------
#
# Functions: to_plink_pheno()  Transpose a feature x sample table to PLINK FID/IID format
#
# Packages:  dplyr, tidyr (called with ::, so they must be installed but not loaded)
#
# Used in:   twas_prep_expr_for_FUSION.R, twas_prep_covar_for_FUSION.R
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

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
