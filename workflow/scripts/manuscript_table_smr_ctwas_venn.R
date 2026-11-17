#--------------------------------------------------------------------------------------
#
#    SMR vs cTWAS gene overlap tables, one per GWAS trait (Fig 6)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: smr_ctwas_venn_tbl
#            Upstream:   smr_sig_tbl (16); 14 cTWAS results and FUSION weights
#            Downstream: 17 smr_ctwas_venn_fig (Fig 6)
#
# Purpose:   For each GWAS trait:
#            - SMR genes: union across cell types of the smr_sig_tbl genes
#            - cTWAS genes: union across cell types of fine-mapped genes
#              (PIP > ctwas_pip_thresh, in a credible set), optionally with
#              Bonferroni P < ctwas_bf_thresh (factor = testable genes per cell type)
#            - Write a membership table (both / SMR_only / cTWAS_only) and log
#              the gene lists per category
#
# Inputs:    smr_tbl_dir       Directory with {gwas}_smr.tsv (smr_sig_tbl)
#            ctwas_in_dir      14 cTWAS results directory (output/ctwas_{ct}_{gwas}_ctwas.rds)
#            weights_dir       14 FUSION weights directory ({ct}/*.wgt.RDat)
#            gene_lookup_file  Ensembl ID to gene symbol lookup
#
# Outputs:   venn_{gwas}_genes.tsv, one per trait, same order as traits
#            (gene, symbol, in_SMR, in_cTWAS, category)
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

message("\n\nBuilding SMR vs cTWAS gene overlap tables ...")

##  Load packages and variables  ------------------------------------------------------
suppressMessages({
  library(tidyverse)
})

# Inputs, outputs and settings
traits           <- snakemake@params[["traits"]]
smr_tbl_dir      <- snakemake@params[["smr_tbl_dir"]]
ctwas_in_dir     <- snakemake@params[["ctwas_in_dir"]]
weights_dir      <- snakemake@params[["weights_dir"]]
cell_types       <- snakemake@params[["cell_types"]]
gene_lookup_file <- snakemake@params[["gene_lookup_file"]]
ctwas_pip_thresh <- snakemake@params[["ctwas_pip_thresh"]]
ctwas_use_bf     <- snakemake@params[["ctwas_use_bf"]]
ctwas_bf_thresh  <- snakemake@params[["ctwas_bf_thresh"]]
use_symbols      <- snakemake@params[["use_symbols"]]
genes_out_files  <- snakemake@output

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("traits", "smr_tbl_dir", "ctwas_in_dir", "weights_dir", "cell_types",
               "gene_lookup_file", "ctwas_pip_thresh", "ctwas_use_bf", "ctwas_bf_thresh",
               "use_symbols", "genes_out_files"),
  value    = c(paste(traits, collapse = ", "), smr_tbl_dir, ctwas_in_dir, weights_dir,
               paste(cell_types, collapse = ", "), gene_lookup_file, ctwas_pip_thresh,
               ctwas_use_bf, ctwas_bf_thresh, use_symbols,
               paste(unlist(genes_out_files), collapse = ", "))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  ------------------------------------------------------------------------
# SMR significant genes: {smr_tbl_dir}/{gwas}_smr.tsv from smr_sig_tbl
get_smr_genes <- function(gwas, smr_tbl_dir) {
  f <- file.path(smr_tbl_dir, paste0(gwas, "_smr.tsv"))
  if (!file.exists(f)) {
    message("No SMR sig table for ", gwas, " (", f, ") - treating as empty set")
    return(character(0))
  }
  message("Reading SMR table: ", f)
  tbl <- tryCatch(
    read_tsv(f, show_col_types = FALSE),
    error = function(e) stop("Failed to read SMR table (", f, "): ", conditionMessage(e))
  )
  if (!"Ensembl ID" %in% names(tbl)) {
    stop("SMR table ", f, " has no 'Ensembl ID' column - columns present: ",
         paste(names(tbl), collapse = ", "))
  }
  unique(tbl[["Ensembl ID"]])
}

# cTWAS significant genes, as in ctwas_report.Rmd: group != 'SNP',
# susie_pip > pip_thresh and in a credible set (!is.na(cs))
`%||%` <- function(a, b) if (is.null(a)) b else a

# Bonferroni factor, as in the ctwas_report.Rmd Correction tab: the number of
# non-empty .wgt.RDat FUSION weight files for a cell type, i.e. the genes that
# were testable there (hsq_p < 0.01), not the number that came out significant
get_bf_factor <- function(ct, weights_dir) {
  ct_dir <- file.path(weights_dir, ct)
  if (!dir.exists(ct_dir)) {
    message("  [", ct, "] weights dir not found for BF factor: ", ct_dir)
    return(NA_integer_)
  }
  files <- list.files(ct_dir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
  files <- files[file.info(files)$size > 0]
  n <- length(files)
  if (n == 0) {
    message("  [", ct, "] zero non-empty .wgt.RDat files found in ", ct_dir)
    return(NA_integer_)
  }
  n
}

get_ctwas_genes <- function(gwas, ctwas_in_dir, weights_dir, cell_types,
                             pip_thresh, use_bf, bf_thresh) {
  bf_skipped_cts <- character(0)  # cell types skipped for missing/empty weights

  genes <- map(cell_types, function(ct) {
    f <- file.path(ctwas_in_dir, "output", paste0("ctwas_", ct, "_", gwas, "_ctwas.rds"))
    if (!file.exists(f)) {
      message("  [", ct, "] file not found: ", f)
      return(character(0))
    }
    message("  [", ct, "] reading ", f)
    res <- tryCatch(
      read_rds(f),
      error = function(e) stop("Failed to read RDS for cell type '", ct, "' (", f, "): ",
                                conditionMessage(e))
    )
    if (!is.null(res$status) && res$status == "skipped") {
      message("  [", ct, "] skipped: ", res$reason %||% "no reason given")
      return(character(0))
    }
    if (is.null(res$finemap_res) || nrow(res$finemap_res) == 0) {
      message("  [", ct, "] no finemap_res rows")
      return(character(0))
    }

    fm <- tryCatch({
      as_tibble(res$finemap_res) |>
        filter(group != "SNP", susie_pip > pip_thresh, !is.na(cs)) |>
        separate(id, into = c("gene", "id"), sep = "\\|") |>
        separate(id, into = c("cell_type", "id"), sep = "_")
    }, error = function(e) {
      stop("Failed to parse finemap_res for cell type '", ct, "' (", f, "): ",
           conditionMessage(e))
    })

    if (use_bf) {
      if (nrow(fm) == 0) {
        # Nothing passed the PIP / credible-set filter in this cell type
        return(character(0))
      }
      bf_factor <- get_bf_factor(ct, weights_dir)
      if (is.na(bf_factor)) {
        # A cell type can have no FUSION weights (e.g. small subclusters or
        # pseudotime bins): it contributes zero genes
        message("  [", ct, "] skipping BF correction (see reason above) - contributes 0 genes")
        bf_skipped_cts <<- c(bf_skipped_cts, ct)
        return(character(0))
      }
      fm <- fm |>
        mutate(pval = 2 * pnorm(-abs(z)),
               pval_bf = pmin(pval * bf_factor, 1)) |>
        filter(pval_bf < bf_thresh)
    }
    message("  [", ct, "] ", length(unique(fm$gene)), " significant genes")
    unique(fm$gene)
  })

  if (use_bf && length(bf_skipped_cts) == length(cell_types) && length(cell_types) > 0) {
    message("WARNING: ALL ", length(cell_types), " cell types were skipped for BF ",
            "correction (missing/empty weights) for ", gwas, " - this looks like a ",
            "'weights_dir' path problem rather than normal per-cell-type variation. ",
            "Check the path: ", weights_dir)
  } else if (use_bf && length(bf_skipped_cts) > 0) {
    message("  BF correction skipped for ", length(bf_skipped_cts), "/", length(cell_types),
            " cell types (no weights found): ", paste(bf_skipped_cts, collapse = ", "))
  }

  unique(unlist(genes))
}

# Ensembl -> symbol lookup, falling back to the Ensembl ID
to_display <- function(ids, use_symbols, gene_lookup_file) {
  if (length(ids) == 0) return(character(0))
  if (!use_symbols || is.null(gene_lookup_file) || !file.exists(gene_lookup_file)) {
    if (use_symbols) message("Gene lookup file not found (", gene_lookup_file,
                              ") - falling back to Ensembl IDs")
    return(ids)
  }
  lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
    dplyr::select(ensembl_gene_id, external_gene_name) |>
    distinct()
  tibble(ensembl_gene_id = ids) |>
    left_join(lookup, by = "ensembl_gene_id") |>
    mutate(external_gene_name = case_when(
      is.na(external_gene_name) | external_gene_name == "" ~ ensembl_gene_id,
      TRUE ~ external_gene_name
    )) |>
    pull(external_gene_name)
}

# Per-trait gene sets -> membership table, logging the full gene lists
log_gene_list <- function(label, symbols) {
  if (length(symbols) == 0) {
    message("  ", label, ": (none)")
  } else {
    message("  ", label, " (n=", length(symbols), "): ", paste(sort(symbols), collapse = ", "))
  }
}

build_membership <- function(gwas, smr_tbl_dir, ctwas_in_dir, weights_dir, cell_types,
                              ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh,
                              use_symbols, gene_lookup_file) {

  message("--- ", gwas, " ---")
  message("Getting SMR genes...")
  smr_genes   <- get_smr_genes(gwas, smr_tbl_dir)
  message("Getting cTWAS genes...")
  ctwas_genes <- get_ctwas_genes(gwas, ctwas_in_dir, weights_dir, cell_types,
                                  ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh)
  message("SMR genes: ", length(smr_genes), " | cTWAS genes: ", length(ctwas_genes))

  membership <- tibble(gene = union(smr_genes, ctwas_genes)) |>
    mutate(
      symbol   = to_display(gene, use_symbols, gene_lookup_file),
      in_SMR   = gene %in% smr_genes,
      in_cTWAS = gene %in% ctwas_genes,
      category = case_when(
        in_SMR  & in_cTWAS  ~ "both",
        in_SMR  & !in_cTWAS ~ "SMR_only",
        !in_SMR & in_cTWAS  ~ "cTWAS_only"
      )
    ) |>
    arrange(category, gene)

  message(sprintf(
    "%s: SMR=%d, cTWAS=%d, both=%d, SMR_only=%d, cTWAS_only=%d",
    gwas, length(smr_genes), length(ctwas_genes),
    sum(membership$category == "both"),
    sum(membership$category == "SMR_only"),
    sum(membership$category == "cTWAS_only")
  ))
  # Full gene lists per category, not just counts
  log_gene_list("Both (SMR & cTWAS)", membership$symbol[membership$category == "both"])
  log_gene_list("SMR only",           membership$symbol[membership$category == "SMR_only"])
  log_gene_list("cTWAS only",         membership$symbol[membership$category == "cTWAS_only"])

  membership
}

##  Membership table per trait  -------------------------------------------------------
if (length(genes_out_files) != length(traits)) {
  stop("Number of Snakemake outputs (", length(genes_out_files),
       ") doesn't match number of traits (", length(traits),
       ") -- check the expand(..., gwas = ...) list matches `traits` exactly.")
}

for (i in seq_along(traits)) {
  gwas <- traits[i]
  membership <- build_membership(gwas, smr_tbl_dir, ctwas_in_dir, weights_dir, cell_types,
                                  ctwas_pip_thresh, ctwas_use_bf, ctwas_bf_thresh,
                                  use_symbols, gene_lookup_file)
  write_tsv(membership, genes_out_files[[i]])
  message("Wrote membership table: ", genes_out_files[[i]])
}

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
