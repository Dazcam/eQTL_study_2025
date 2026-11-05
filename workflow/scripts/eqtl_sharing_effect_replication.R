#--------------------------------------------------------------------------------------
#
#    Discovery -> replication eQTL effect-size correlations
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: eqtl_effect_replication (single job)
#            Upstream:   05TENSORQTL tensorqtl_perm; 10SMR nominal export
#                        ({ct}_nom.cis_qtl_pairs); Jang et al. 2026 L1 full-association
#                        files; 05TENSORQTL convert_genotypes (pvar)
#            Downstream: report (effect replication tab);
#                        14 effect_size_heatmaps_plt (Supp Figs 7-8)
#
# Purpose:   For each fetal discovery cell type, take its sig. lead eGene-eSNP pairs
#            (qval < 0.05) and look up the same pairs' effect sizes in each
#            replication cell type, with no significance cutoff there. Pearson r,
#            Spearman r and sign concordance per combination:
#            - Fetal -> fetal: all fetal cell types (L1 + L2); A -> B differs from B -> A
#            - Fetal L1 -> adult: the 7 Jang L1 types
#
# Inputs:    qtl_perm_files  tensorQTL permutation output, one per fetal cell type
#            nom_files       fetal nominal pairs, same order as qtl_perm_files
#            jang_files      Jang L1 full-association files
#            pvar_file       our genotype pvar (REF/ALT per rsID)
#            cell_types      fetal cell types, same order as the files (params)
#            fetal_l1        fetal L1 cell types (params)
#            jang_order      Jang L1 order for the fetal -> adult matrix (params)
#            fdr_thresh      qval threshold for discovery sets (params)
#
# Outputs:   output_file  rds list: fetal_vs_fetal, fetal_vs_adult (incl.
#                         allele_status), discovery_counts, settings
#
# Notes:     - Discovery slope is from the perm file; replication slope from the
#              nominal file (fetal) or Jang fixed_beta. The diagonal (perm vs own
#              nominal) is a sanity check (r ~ 1).
#            - Jang betas are aligned to Jang's alt, then harmonised to our alt;
#              palindromic, mismatched and missing-allele pairs are dropped.
#            - Pairs not tested in a replication cell type are absent; n_tested
#              records how many pairs contributed to each r.
#
#--------------------------------------------------------------------------------------

## Set up logging for Snakemake
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

message("\n\nComputing discovery -> replication eQTL effect-size correlations ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
qtl_perm_files <- snakemake@input[["qtl_perm_files"]]   # perm files, all fetal cell types
nom_files      <- snakemake@input[["nom_files"]]        # nominal files, same order
jang_files     <- snakemake@input[["jang_files"]]       # Jang full-assoc, 7 L1-pooled types
pvar_file      <- snakemake@input[["pvar_file"]]
output_file    <- snakemake@output[[1]]
cell_types     <- snakemake@params[["cell_types"]]      # all fetal (L1+L2), no pseudotime bins
fetal_l1       <- snakemake@params[["fetal_l1"]]        # 7 fetal L1 types
FDR_THRESH     <- snakemake@params[["fdr_thresh"]]

# Preserve cell_type <-> file correspondence and verify it from the file names,
# so an ordering mismatch between helpers fails loudly instead of silently
stopifnot(length(qtl_perm_files) == length(cell_types),
          length(nom_files)      == length(cell_types))
names(qtl_perm_files) <- cell_types
names(nom_files)      <- cell_types
stopifnot(all(startsWith(basename(qtl_perm_files), paste0(cell_types, "_"))),
          all(startsWith(basename(nom_files),      paste0(cell_types, "_"))))

jang_types <- sub("_eqtl_full_assoc\\.tsv\\.gz$", "", basename(jang_files))
names(jang_files) <- jang_types
jang_order <- intersect(unlist(snakemake@params[["jang_order"]]), jang_types)
stopifnot(length(jang_order) == length(jang_types))

stopifnot(all(fetal_l1 %in% cell_types))

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("n_fetal_cell_types", "n_fetal_l1", "jang_types", "pvar_file",
               "fdr_thresh", "output_file"),
  value    = c(length(cell_types), length(fetal_l1), paste(jang_types, collapse = ","),
               pvar_file, FDR_THRESH, output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  ------------------------------------------------------------------------

# One row of results for a discovery/replication combination
cor_row <- function(disc_ct, rep_ct, n_discovery, x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  can_cor <- n >= 3 && sd(x) > 0 && sd(y) > 0
  data.table(
    discovery      = disc_ct,
    replication    = rep_ct,
    n_discovery    = n_discovery,
    n_tested       = n,
    pearson_r      = if (can_cor) cor(x, y) else NA_real_,
    spearman_r     = if (can_cor) cor(x, y, method = "spearman") else NA_real_,
    prop_same_sign = if (n > 0) mean(sign(x) == sign(y)) else NA_real_
  )
}

# Long table -> discovery (rows) x replication (columns) matrix
to_matrix <- function(dt, value, rows, cols) {
  m <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  m[cbind(match(dt$discovery, rows), match(dt$replication, cols))] <- dt[[value]]
  m
}

##  Discovery sets (sig. lead pairs per fetal cell type)  -----------------------------

message("--- Reading discovery sets (qval < ", FDR_THRESH, ") ---\n")

disc <- rbindlist(lapply(cell_types, function(ct) {
  message("Reading perm file for: ", ct, " (", qtl_perm_files[[ct]], ")")
  dt <- fread(qtl_perm_files[[ct]],
              select = c("phenotype_id", "variant_id", "slope", "qval"),
              showProgress = FALSE)
  dt <- dt[!is.na(qval) & qval < FDR_THRESH]
  dt <- unique(dt, by = "phenotype_id")   # one lead SNP per gene per cell type
  message("  Discovery pairs: ", nrow(dt))
  dt[, .(discovery = ct, phenotype_id, variant_id, slope_disc = slope)]
}))

discovery_counts <- disc[, .(n_discovery = .N), by = discovery]
discovery_counts <- discovery_counts[match(cell_types, discovery)]
message("\nTotal discovery rows: ", nrow(disc),
        " | unique (gene, SNP) pairs: ", uniqueN(disc, by = c("phenotype_id", "variant_id")))

##  Fetal -> fetal  -------------------------------------------------------------------

message("\n--- Fetal -> fetal replication slopes (nominal files) ---\n")

pairs_all <- unique(disc[, .(phenotype_id, variant_id)])
setkey(pairs_all, phenotype_id, variant_id)

fetal_rep <- rbindlist(lapply(cell_types, function(ct) {
  message("Reading nominal file for: ", ct, " (", nom_files[[ct]], ")")
  nom <- fread(nom_files[[ct]],
               select = c("phenotype_id", "variant_id", "slope"),
               showProgress = FALSE)
  setkey(nom, phenotype_id, variant_id)
  m <- nom[pairs_all, nomatch = 0L]
  m <- unique(m, by = c("phenotype_id", "variant_id"))
  message("  Universe pairs found in nominal file: ", nrow(m), " / ", nrow(pairs_all))
  rm(nom); gc()
  m[, .(replication = ct, phenotype_id, variant_id, slope_rep = slope)]
}))

fetal_rep_split <- split(fetal_rep, by = "replication", keep.by = FALSE)

message("\nComputing fetal -> fetal correlations ...")
fetal_pairs <- list()
fetal_rows  <- list()

for (a in cell_types) {
  da <- disc[discovery == a]
  for (b in cell_types) {
    m <- merge(da, fetal_rep_split[[b]], by = c("phenotype_id", "variant_id"))
    fetal_rows[[paste(a, b)]]  <- cor_row(a, b, nrow(da), m$slope_disc, m$slope_rep)
    fetal_pairs[[paste(a, b)]] <- m[, .(discovery = a, replication = b,
                                        phenotype_id, variant_id, slope_disc, slope_rep)]
  }
}
fetal_long  <- rbindlist(fetal_rows)
fetal_pairs <- rbindlist(fetal_pairs)

message("Fetal -> fetal diagonal (expect r ~= 1): ",
        paste(round(fetal_long[discovery == replication, pearson_r], 3), collapse = " "))

##  Fetal L1 -> adult (Jang L1)  ------------------------------------------------------

message("\n--- Fetal L1 -> Jang adult replication slopes ---\n")

# Discovery restricted to fetal L1; rsID/ENSG only so IDs match Jang's
disc_l1 <- disc[discovery %in% fetal_l1]
n_l1_total <- nrow(disc_l1)
disc_l1 <- disc_l1[grepl("^rs", variant_id) & grepl("^ENSG", phenotype_id)]
message("Fetal L1 discovery rows: ", n_l1_total, " | with rsID + ENSG ids: ", nrow(disc_l1))

pairs_l1 <- unique(disc_l1[, .(phenotype_id, variant_id)])
setkey(pairs_l1, phenotype_id, variant_id)

# Our REF/ALT (pvar), restricted to SNPs we need
message("\nLoading REF/ALT from pvar ...")
pvar_dt <- fread(pvar_file, select = c("ID", "REF", "ALT"),
                 col.names = c("variant_id", "ref_my", "alt_my"), showProgress = FALSE)
pvar_dt <- pvar_dt[variant_id %in% pairs_l1$variant_id]
dup_ids <- pvar_dt[duplicated(variant_id), variant_id]
if (length(dup_ids) > 0) {
  message("  Dropping ", length(unique(dup_ids)), " SNP IDs with >1 pvar record (ambiguous alleles)")
  pvar_dt <- pvar_dt[!variant_id %in% dup_ids]
}
setkey(pvar_dt, variant_id)
message("  Needed SNPs found in pvar: ", nrow(pvar_dt), " / ", uniqueN(pairs_l1$variant_id))

jang_rep <- list()
allele_summary <- list()

for (jct in jang_order) {
  message("\nReading Jang full-association file for: ", jct, " (", jang_files[[jct]], ")")
  jd <- fread(jang_files[[jct]],
              select = c("feature", "variant_id", "ref", "alt", "Allele",
                         "fixed_beta", "Fixed_P", "Random_P"),
              showProgress = FALSE)

  jd <- filter_jang_artifacts(jd, "Jang full-association")
  jd[, c("Fixed_P", "Random_P") := NULL]

  jd[, phenotype_id := sub("\\..*", "", feature)]
  jd <- jd[grepl("^rs", variant_id) & grepl("^ENSG", phenotype_id)]
  jd <- unique(jd, by = c("phenotype_id", "variant_id"))
  setkey(jd, phenotype_id, variant_id)

  m <- jd[pairs_l1, nomatch = 0L]
  rm(jd); gc()
  message("  Universe pairs found in Jang: ", nrow(m), " / ", nrow(pairs_l1))

  # Sign-normalise to Jang alt; NA if Allele is neither ref nor alt
  m[, beta_jang_alt := align_jang_beta(fixed_beta, Allele, alt, ref)]

  m <- merge(m, pvar_dt, by = "variant_id", all.x = TRUE)
  m[, allele_status := get_allele_status(ref_my, alt_my, ref, alt, beta_jang_alt)]
  m[, slope_rep := harmonise_beta(beta_jang_alt, allele_status)]

  status_tbl <- m[, .N, by = allele_status][, jang_cell_type := jct]
  allele_summary[[jct]] <- status_tbl
  message("  Allele status: ",
          paste(status_tbl$allele_status, status_tbl$N, sep = "=", collapse = ", "))

  jang_rep[[jct]] <- m[!is.na(slope_rep), .(phenotype_id, variant_id, slope_rep)]
  message("  Usable (harmonised) pairs: ", nrow(jang_rep[[jct]]))
}
allele_summary <- rbindlist(allele_summary, use.names = TRUE)

message("\nComputing fetal L1 -> Jang correlations ...")
adult_pairs <- list()
adult_rows  <- list()

for (a in fetal_l1) {
  da <- disc_l1[discovery == a]
  n_disc_raw <- discovery_counts[discovery == a, n_discovery]
  for (j in jang_order) {
    m <- merge(da, jang_rep[[j]], by = c("phenotype_id", "variant_id"))
    adult_rows[[paste(a, j)]]  <- cor_row(a, j, n_disc_raw, m$slope_disc, m$slope_rep)
    adult_pairs[[paste(a, j)]] <- m[, .(discovery = a, replication = j,
                                        phenotype_id, variant_id, slope_disc, slope_rep)]
  }
}
adult_long  <- rbindlist(adult_rows)
adult_pairs <- rbindlist(adult_pairs)

##  Write output  ---------------------------------------------------------------------

make_block <- function(long, pairs, rows, cols) {
  list(
    long           = as_tibble(long),
    pairs          = as_tibble(pairs),
    pearson        = to_matrix(long, "pearson_r",      rows, cols),
    spearman       = to_matrix(long, "spearman_r",     rows, cols),
    prop_same_sign = to_matrix(long, "prop_same_sign", rows, cols),
    n_tested       = to_matrix(long, "n_tested",       rows, cols)
  )
}

result <- list(
  fetal_vs_fetal   = make_block(fetal_long, fetal_pairs, cell_types, cell_types),
  fetal_vs_adult   = c(make_block(adult_long, adult_pairs, fetal_l1, jang_order),
                       list(allele_status = as_tibble(allele_summary))),
  discovery_counts = as_tibble(discovery_counts),
  settings         = list(fdr_thresh = FDR_THRESH, fetal_cell_types = cell_types,
                          fetal_l1 = fetal_l1, jang_cell_types = jang_order)
)

message("\n--- Summary ---\n")
message("Fetal -> fetal Pearson r (rows = discovery, cols = replication):")
print(round(result$fetal_vs_fetal$pearson, 2))
message("\nFetal L1 -> Jang Pearson r:")
print(round(result$fetal_vs_adult$pearson, 2))
message("\nFetal L1 -> Jang n_tested:")
print(result$fetal_vs_adult$n_tested)

write_rds(result, output_file)
message("\nEffect-size replication correlations written to: ", output_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
