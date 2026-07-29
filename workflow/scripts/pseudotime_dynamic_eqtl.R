#--------------------------------------------------------------------------------------
#
#    pseudotime_dynamic_eqtl.R
#
#    Dynamic eQTL mapping across pseudotime bins using linear mixed models.
#
#    Adapted from Yazar et al. (Science 2022) dyn-eqtl_mapping.R.
#    Tests whether eQTL effect sizes change linearly or quadratically across
#    six pseudotime quantile bins using lmer with individual random effects.
#
#    Analysis design:
#      - Input eGene set: cell-type level eQTLs (Option A, matching Yazar)
#      - Filter: eGenes expressed in >= min_bins Q6 bins
#      - Expression: log(raw_counts + 1) from pseudobulk CSVs per bin
#      - Covariates: PCW, sex, genPC1-4 only (no expression PCs)
#        See report methods note for rationale
#      - Models: linear (snp*quantile) and quadratic (snp*quantile + snp*quantile^2)
#        vs respective null models via likelihood ratio test
#      - FDR: qvalue on LRT p-values
#      - Parallelisation: mclapply across SNP-gene pairs
#
#    Outputs:
#      - dynamic_eqtl_results.rds: full LMM results per SNP-gene pair
#      - dynamic_eqtl_significant.tsv: significant dynamic pairs only
#      - dynamic_eqtl_candidates.tsv: stacked per-individual data for all
#        dynamic pairs — used for Panel E style plots in the report
#
#    Usage: called via Snakemake rule pseudotime_dynamic_eqtl
#
#--------------------------------------------------------------------------------------

## Logging  --------------------------------------------------------------------------

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, append = TRUE)
sink(con, append = TRUE, type = "message")

message("\n\nRunning pseudotime dynamic eQTL analysis ...")
message("Timestamp: ", Sys.time())

## Libraries  ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(qvalue)
  library(parallel)
})

## Parameters  -----------------------------------------------------------------------

trajectory         <- as.character(snakemake@wildcards[["trajectory"]])
dosage_sub         <- as.character(snakemake@input[["dosage_sub"]])
pvar_file          <- as.character(snakemake@input[["pvar"]])
candidates_tsv     <- as.character(snakemake@output[["candidates_tsv"]])
n_threads          <- as.integer(snakemake@threads)

# Paths passed as params — wildcards already resolved by Snakemake
perm_base          <- as.character(snakemake@params[["perm_dir"]])
pseudobulk_dir     <- as.character(snakemake@params[["pseudobulk_dir"]])
cov_dir            <- as.character(snakemake@params[["cov_dir"]])
results_rds        <- as.character(snakemake@params[["results_rds"]])
results_tsv        <- as.character(snakemake@params[["results_tsv"]])
cell_type_perm_dir <- as.character(snakemake@params[["cell_type_perm_dir"]])

# Derive cell type from trajectory name (NPC-to-Glu-DL -> Glu-DL)
cell_type <- str_split(trajectory, "-to-")[[1]][2]

# Config parameters
min_bins      <- as.integer(snakemake@config[["pseudotime"]][["dynamic_eqtl"]][["min_bins"]])
fdr_threshold <- as.numeric(snakemake@config[["pseudotime"]][["dynamic_eqtl"]][["fdr_threshold"]])
quantile_n    <- as.integer(snakemake@config[["pseudotime"]][["dynamic_eqtl"]][["quantile_n"]])
exp_pc        <- as.integer(snakemake@config[["pseudotime"]][["dynamic_eqtl"]][["exp_pc"]])
norm_method   <- as.character(snakemake@config[["tensorQTL"]][["norm_methods"]][[1]])
geno_pcs      <- as.integer(snakemake@config[["tensorQTL"]][["geno_pcs"]])

# Cell type expPC for cell type perm file
cell_type_exp_pc_map <- snakemake@config[["pseudotime"]][["cell_type_exp_pc_map"]]
cell_type_exp_pc     <- as.integer(cell_type_exp_pc_map[[trajectory]])

# Output directory
out_dir <- dirname(candidates_tsv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("\n============================")
message("Parameters:")
message("  trajectory:           ", trajectory)
message("  cell_type:            ", cell_type)
message("  cell_type_exp_pc:     ", cell_type_exp_pc)
message("  quantile_n:           ", quantile_n)
message("  exp_pc:               ", exp_pc)
message("  min_bins:             ", min_bins)
message("  fdr_threshold:        ", fdr_threshold)
message("  n_threads:            ", n_threads)
message("  perm_base:            ", perm_base)
message("  pseudobulk_dir:       ", pseudobulk_dir)
message("  cov_dir:              ", cov_dir)
message("  cell_type_perm_dir:   ", cell_type_perm_dir)
message("  dosage_sub:           ", dosage_sub)
message("  pvar_file:            ", pvar_file)
message("  results_rds:          ", results_rds)
message("  results_tsv:          ", results_tsv)
message("  candidates_tsv:       ", candidates_tsv)
message("============================\n")

## 1. Load cell type eGenes  ---------------------------------------------------------
#
#  Load eGenes from the cell type level analysis (Glu-DL or Glu-UL).
#  These form the starting set for dynamic eQTL testing, matching Yazar et al.
#  who tested their 2,343 B cell eQTLs for dynamic effects.

cell_type_perm_file <- file.path(
  cell_type_perm_dir,
  sprintf("%s_%s_genPC_%d_expPC_%d", cell_type, norm_method, geno_pcs, cell_type_exp_pc),
  sprintf("%s_%s_perm.cis_qtl.txt.gz", cell_type, norm_method)
)

message("Loading cell type eGenes from: ", cell_type_perm_file)
if (!file.exists(cell_type_perm_file)) {
  stop("Cell type perm file not found: ", cell_type_perm_file)
}

cell_type_eqtl <- fread(cell_type_perm_file) |>
  filter(!is.na(qval) & qval < fdr_threshold) |>
  select(phenotype_id, variant_id, slope, slope_se, qval) |>
  distinct(phenotype_id, .keep_all = TRUE)

message("Cell type eGenes at FDR < ", fdr_threshold, ": ", nrow(cell_type_eqtl))

## 2. Load Q6 perm files and check bin presence  ------------------------------------
#
#  For each Q6 bin, load perm output and record which eGenes are detected.
#  Filter to eGenes present in >= min_bins bins, matching Yazar filter.

message("\nLoading Q6 permutation files ...")
bin_egenes <- list()

for (b in seq_len(quantile_n)) {
  perm_file <- file.path(
    perm_base,
    sprintf("Q%d_bin%d_genPC%d_expPC%d", quantile_n, b, geno_pcs, exp_pc),
    sprintf("Q%d_bin%d_perm.cis_qtl.txt.gz", quantile_n, b)
  )
  if (!file.exists(perm_file)) {
    message("  Bin ", b, ": file not found — ", perm_file)
    bin_egenes[[b]] <- character(0)
    next
  }
  df  <- fread(perm_file)
  sig <- df[!is.na(qval) & qval < fdr_threshold, phenotype_id]
  bin_egenes[[b]] <- sig
  message("  Bin ", b, ": ", length(sig), " eGenes at FDR < ", fdr_threshold)
}

# Count bins per eGene across Q6
all_genes  <- unique(cell_type_eqtl$phenotype_id)
bin_counts <- map_int(all_genes, function(gene) {
  sum(map_lgl(bin_egenes, ~ gene %in% .x))
})
names(bin_counts) <- all_genes

# Apply min_bins filter
genes_pass <- names(bin_counts)[bin_counts >= min_bins]
genes_fail <- names(bin_counts)[bin_counts <  min_bins]

message("\n============================")
message("eGene filtering:")
message("  Cell type eGenes total:             ", length(all_genes))
message("  eGenes in >= ", min_bins, " bins (retained): ", length(genes_pass))
message("  eGenes in <  ", min_bins, " bins (excluded): ", length(genes_fail))
message("============================\n")

if (length(genes_pass) == 0) {
  stop("No eGenes passed the minimum bin filter. Check perm file paths and FDR threshold.")
}

eqtl_test <- cell_type_eqtl |> filter(phenotype_id %in% genes_pass)
message("SNP-gene pairs to test: ", nrow(eqtl_test))

## 3. Load pseudobulk expression (raw counts, log transformed)  ---------------------
#
#  Load raw pseudobulk counts per Q6 bin per sample.
#  Apply log(x + 1) transformation matching Yazar et al.
#  Input CSVs: samples as rows, genes as columns.

message("\nLoading pseudobulk expression ...")
exp_list <- list()

for (b in seq_len(quantile_n)) {
  csv_file <- file.path(
    pseudobulk_dir,
    sprintf("Q%d", quantile_n),
    sprintf("Q%d_bin%d_%s_pseudobulk.csv", quantile_n, b, trajectory)
  )
  if (!file.exists(csv_file)) stop("Pseudobulk CSV not found: ", csv_file)

  mat <- fread(csv_file, data.table = FALSE)
  rownames(mat) <- mat[[1]]
  mat <- mat[, -1, drop = FALSE]

  mat_t   <- t(mat)
  mat_log <- log(mat_t + 1)

  genes_present <- intersect(genes_pass, rownames(mat_log))
  exp_list[[b]] <- mat_log[genes_present, , drop = FALSE]

  message("  Bin ", b, ": ", nrow(mat_log), " genes x ", ncol(mat_log),
          " samples loaded | ", length(genes_present), " test genes present")
}

## 4. Load covariates  ---------------------------------------------------------------
#
#  Load base covariate file for each Q6 bin.
#  Retain only: PCW, sex_code, genPC1-genPC{geno_pcs}.
#  See report methods note — expression PCs excluded for LMM comparability.

message("\nLoading covariates ...")
cov_list <- list()

for (b in seq_len(quantile_n)) {
  cov_file <- file.path(
    cov_dir,
    sprintf("Q%d_bin%d_quantile_base_covariates.txt", quantile_n, b)
  )
  if (!file.exists(cov_file)) stop("Covariate file not found: ", cov_file)

  cov_mat <- fread(cov_file, data.table = FALSE)
  rownames(cov_mat) <- cov_mat[[1]]
  cov_mat <- cov_mat[, -1, drop = FALSE]

  keep_rows <- c("PCW", "sex_code", paste0("genPC", seq_len(geno_pcs)))
  keep_rows <- intersect(keep_rows, rownames(cov_mat))
  cov_sub   <- t(cov_mat[keep_rows, , drop = FALSE]) |> as.data.frame()
  cov_sub$sampleid <- rownames(cov_sub)

  cov_list[[b]] <- cov_sub
  message("  Bin ", b, ": ", nrow(cov_sub), " samples x ",
          ncol(cov_sub) - 1, " covariates")
}

## 5. Load genotype dosages  ---------------------------------------------------------
#
#  Dosages pre-extracted from full dosage file via pseudotime_extract_dosages rule.
#  Plain TSV with header: CHROM POS REF ALT sample1 sample2 ...
#  rsIDs joined from pvar file on CHROM + POS.

message("\nLoading genotype dosages ...")
snp_ids <- unique(eqtl_test$variant_id)
message("  Unique lead eSNPs needed: ", length(snp_ids))

pvar <- fread(cmd = sprintf("grep -v '^##' %s", pvar_file), data.table = FALSE)
colnames(pvar)[1] <- "CHROM"
pvar$CHROM <- as.character(pvar$CHROM)
pvar_sub <- pvar |> filter(ID %in% snp_ids)
message("  SNPs found in pvar: ", nrow(pvar_sub), " / ", length(snp_ids))
if (nrow(pvar_sub) == 0) stop("No SNPs matched in pvar file")

dose_mat <- fread(dosage_sub, data.table = FALSE) |>
  mutate(CHROM = as.character(CHROM)) |>
  left_join(
    pvar_sub |> select(CHROM, POS, ID, pvar_REF = REF, pvar_ALT = ALT),
    by = c("CHROM", "POS")
  ) |>
  filter(!is.na(ID), ID %in% snp_ids) |>
  distinct(ID, .keep_all = TRUE)

rownames(dose_mat) <- dose_mat$ID

snp_alleles <- dose_mat |>
  rownames_to_column("variant_id") |>
  select(variant_id, REF = pvar_REF, ALT = pvar_ALT)

sample_cols <- setdiff(colnames(dose_mat),
                       c("CHROM", "POS", "REF", "ALT", "ID",
                         "pvar_REF", "pvar_ALT"))

message("  Dosages loaded: ", nrow(dose_mat), " SNPs x ",
        length(sample_cols), " samples")
message("  SNPs matched: ", sum(snp_ids %in% rownames(dose_mat)),
        " / ", length(snp_ids))

## 6. Build stacked LMM input and fit models  ----------------------------------------
#
#  For each SNP-gene pair:
#    - Stack expression + genotype + covariates across all 6 bins
#    - Fit linear and quadratic LMMs
#    - LRT via anova() matching Yazar dyn-eqtl_mapping.R

build_pair_data <- function(gene_id, snp_id) {

  if (!snp_id %in% rownames(dose_mat)) return(NULL)

  dose_vec        <- as.numeric(dose_mat[snp_id, sample_cols])
  names(dose_vec) <- sample_cols

  rows <- list()
  for (b in seq_len(quantile_n)) {
    exp_mat <- exp_list[[b]]
    cov_df  <- cov_list[[b]]

    if (!gene_id %in% rownames(exp_mat)) next

    samples_b <- intersect(colnames(exp_mat), names(dose_vec))
    samples_b <- intersect(samples_b, cov_df$sampleid)
    if (length(samples_b) < 10) next

    gene_expr <- exp_mat[gene_id, samples_b]
    geno_dose <- dose_vec[samples_b]
    cov_sub   <- cov_df[cov_df$sampleid %in% samples_b, ]
    rownames(cov_sub) <- cov_sub$sampleid
    cov_sub   <- cov_sub[samples_b, ]

    df <- data.frame(
      gene     = as.numeric(gene_expr),
      snp      = as.numeric(geno_dose),
      sampleid = samples_b,
      quantile = b,
      stringsAsFactors = FALSE
    )
    df <- cbind(df,
                cov_sub[, setdiff(colnames(cov_sub), "sampleid"), drop = FALSE])
    rows[[b]] <- df
  }

  if (length(rows) == 0) return(NULL)
  bind_rows(rows)
}

fit_lmm_pair <- function(pair_idx) {

  gene_id <- eqtl_test$phenotype_id[pair_idx]
  snp_id  <- eqtl_test$variant_id[pair_idx]

  dat <- tryCatch(build_pair_data(gene_id, snp_id), error = function(e) NULL)

  if (is.null(dat) || nrow(dat) < 20) {
    return(list(
      phenotype_id       = gene_id,
      variant_id         = snp_id,
      p_linear           = NA_real_,
      p_quadratic        = NA_real_,
      n_obs              = if (is.null(dat)) 0L else nrow(dat),
      n_bins             = bin_counts[gene_id],
      singular_linear    = NA,
      singular_quadratic = NA,
      error              = "insufficient data"
    ))
  }

  cov_cols    <- setdiff(colnames(dat),
                         c("gene", "snp", "sampleid", "quantile", "quantile_2"))
  cov_formula <- if (length(cov_cols) > 0) paste(cov_cols, collapse = " + ") else "1"

  # ── Linear model ─────────────────────────────────────────────────────────
  f0_lin <- as.formula(sprintf(
    "gene ~ snp + %s + (1 | sampleid) + quantile", cov_formula))
  f1_lin <- as.formula(sprintf(
    "gene ~ snp + %s + (1 | sampleid) + snp * quantile", cov_formula))

  fit0_lin <- tryCatch(lmer(f0_lin, data = dat, REML = FALSE), error = function(e) NULL)
  fit1_lin <- tryCatch(lmer(f1_lin, data = dat, REML = FALSE), error = function(e) NULL)

  p_linear <- tryCatch(
    anova(fit1_lin, fit0_lin)[2, "Pr(>Chisq)"],
    error = function(e) NA_real_
  )

  # ── Quadratic model ───────────────────────────────────────────────────────
  dat$quantile_2 <- dat$quantile^2

  f0_quad <- as.formula(sprintf(
    "gene ~ snp + %s + (1 | sampleid) + quantile + quantile_2", cov_formula))
  f1_quad <- as.formula(sprintf(
    "gene ~ snp + %s + (1 | sampleid) + snp * quantile + snp * quantile_2",
    cov_formula))

  fit0_quad <- tryCatch(lmer(f0_quad, data = dat, REML = FALSE), error = function(e) NULL)
  fit1_quad <- tryCatch(lmer(f1_quad, data = dat, REML = FALSE), error = function(e) NULL)

  p_quadratic <- tryCatch(
    anova(fit1_quad, fit0_quad)[2, "Pr(>Chisq)"],
    error = function(e) NA_real_
  )

  list(
    phenotype_id       = gene_id,
    variant_id         = snp_id,
    p_linear           = p_linear,
    p_quadratic        = p_quadratic,
    n_obs              = nrow(dat),
    n_bins             = bin_counts[gene_id],
    singular_linear    = if (!is.null(fit1_lin))  isSingular(fit1_lin)  else NA,
    singular_quadratic = if (!is.null(fit1_quad)) isSingular(fit1_quad) else NA,
    error              = NA_character_
  )
}

message("Fitting LMMs for ", nrow(eqtl_test),
        " SNP-gene pairs using ", n_threads, " threads ...")

results_list <- mclapply(
  seq_len(nrow(eqtl_test)),
  fit_lmm_pair,
  mc.cores = n_threads
)

## 7. Aggregate results and apply FDR  -----------------------------------------------

message("\nAggregating results ...")
results_tbl <- bind_rows(map(results_list, as.data.frame)) |> as_tibble()

valid_linear <- !is.na(results_tbl$p_linear)
valid_quad   <- !is.na(results_tbl$p_quadratic)

if (sum(valid_linear) > 10) {
  results_tbl$fdr_linear               <- NA_real_
  results_tbl$fdr_linear[valid_linear] <-
    qvalue(results_tbl$p_linear[valid_linear])$qvalues
} else {
  results_tbl$fdr_linear <- NA_real_
  message("WARNING: too few valid linear p-values for qvalue correction")
}

if (sum(valid_quad) > 10) {
  results_tbl$fdr_quadratic             <- NA_real_
  results_tbl$fdr_quadratic[valid_quad] <-
    qvalue(results_tbl$p_quadratic[valid_quad])$qvalues
} else {
  results_tbl$fdr_quadratic <- NA_real_
  message("WARNING: too few valid quadratic p-values for qvalue correction")
}

results_tbl <- results_tbl |>
  mutate(
    dynamic_linear    = !is.na(fdr_linear)    & fdr_linear    < fdr_threshold,
    dynamic_quadratic = !is.na(fdr_quadratic) & fdr_quadratic < fdr_threshold,
    dynamic           = dynamic_linear | dynamic_quadratic,
    dynamic_type      = case_when(
      dynamic_linear & dynamic_quadratic ~ "linear_and_quadratic",
      dynamic_linear                     ~ "linear",
      dynamic_quadratic                  ~ "quadratic",
      TRUE                               ~ "none"
    )
  )

## 8. Summary logging  ---------------------------------------------------------------

message("\n============================")
message("Results summary — ", trajectory)
message("  SNP-gene pairs tested:           ", nrow(results_tbl))
message("  Pairs with valid linear p:       ", sum(valid_linear))
message("  Pairs with valid quadratic p:    ", sum(valid_quad))
message("  Linear dynamic (FDR<",  fdr_threshold, "):        ",
        sum(results_tbl$dynamic_linear,    na.rm = TRUE))
message("  Quadratic dynamic (FDR<", fdr_threshold, "):    ",
        sum(results_tbl$dynamic_quadratic, na.rm = TRUE))
message("  Linear only:                     ",
        sum(results_tbl$dynamic_type == "linear",               na.rm = TRUE))
message("  Quadratic only:                  ",
        sum(results_tbl$dynamic_type == "quadratic",            na.rm = TRUE))
message("  Linear and quadratic:            ",
        sum(results_tbl$dynamic_type == "linear_and_quadratic", na.rm = TRUE))
message("  Total dynamic eGenes:            ",
        sum(results_tbl$dynamic, na.rm = TRUE))
message("  Singular linear fits:            ",
        sum(results_tbl$singular_linear,    na.rm = TRUE))
message("  Singular quadratic fits:         ",
        sum(results_tbl$singular_quadratic, na.rm = TRUE))
message("============================\n")

## 9. Write results RDS and TSV  -----------------------------------------------------

message("Writing full results RDS: ", results_rds)
saveRDS(results_tbl, results_rds)

message("Writing significant results TSV: ", results_tsv)
results_tbl |>
  filter(dynamic) |>
  arrange(dynamic_type, fdr_linear) |>
  write_tsv(results_tsv)

## 10. Extract candidate data for Panel E plots  -------------------------------------
#
#  For all dynamic SNP-gene pairs, extract the per-individual stacked expression
#  and genotype data used in the LMM. Includes genotype calls derived from
#  rounded dosage values and REF/ALT alleles from the dosage file.
#  The report filters to top candidates from this complete dataset.

message("\nExtracting candidate data for Panel E plots ...")

dynamic_pairs <- results_tbl |>
  filter(dynamic) |>
  mutate(min_fdr = pmin(fdr_linear, fdr_quadratic, na.rm = TRUE)) |>
  arrange(min_fdr)

message("  Dynamic pairs to extract: ", nrow(dynamic_pairs))

candidate_data <- map_dfr(seq_len(nrow(dynamic_pairs)), function(i) {
  gene_id <- dynamic_pairs$phenotype_id[i]
  snp_id  <- dynamic_pairs$variant_id[i]

  dat <- tryCatch(build_pair_data(gene_id, snp_id), error = function(e) NULL)
  if (is.null(dat)) return(NULL)

  alleles <- snp_alleles |> filter(variant_id == snp_id)
  ref     <- if (nrow(alleles) > 0) alleles$REF[1] else "R"
  alt     <- if (nrow(alleles) > 0) alleles$ALT[1] else "A"

  dat |>
    mutate(
      phenotype_id   = gene_id,
      variant_id     = snp_id,
      dosage_rounded = round(snp),
      genotype_call  = case_when(
        dosage_rounded == 0 ~ paste0(ref, ref),
        dosage_rounded == 1 ~ paste0(ref, alt),
        dosage_rounded == 2 ~ paste0(alt, alt),
        TRUE                ~ NA_character_
      ),
      dynamic_type = dynamic_pairs$dynamic_type[i],
      min_fdr      = dynamic_pairs$min_fdr[i],
      fdr_linear   = dynamic_pairs$fdr_linear[i],
      fdr_quadratic= dynamic_pairs$fdr_quadratic[i]
    ) |>
    select(phenotype_id, variant_id, sampleid, quantile,
           gene_expr = gene, dosage = snp, dosage_rounded,
           genotype_call, dynamic_type, min_fdr,
           fdr_linear, fdr_quadratic)
})

message("  Candidate data rows: ", nrow(candidate_data))
message("  Writing candidates TSV: ", candidates_tsv)
write_tsv(candidate_data, candidates_tsv)

message("Done.")
sink()
sink(type = "message")
