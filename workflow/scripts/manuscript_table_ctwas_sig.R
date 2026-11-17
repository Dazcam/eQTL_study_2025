#--------------------------------------------------------------------------------------
#
#    Significant cTWAS genes per GWAS (Supp Table 10)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: ctwas_sig_tbl
#            Upstream:   14 cTWAS results and FUSION weights
#            Downstream: Supp Table 10
#
# Purpose:   - Bonferroni factor per cell type: the number of testable genes
#              (non-empty FUSION .wgt.RDat files), as in ctwas_report.Rmd;
#              cell types without weights are skipped
#            - For each GWAS trait and cell type, fine-mapped genes (PIP >
#              pip_thresh, in a credible set) with Bonferroni P < bf_p_thresh,
#              using that cell type's factor (P = 2 * pnorm(-|z|))
#            - Stack the cell types (one row per gene and cell type); one xlsx
#              sheet per disorder, manuscript labels
#
# Inputs:    in_dir            14 cTWAS results directory (output/ctwas_{ct}_{gwas}_ctwas.rds)
#            weights_dir       14 FUSION weights directory ({ct}/*.wgt.RDat)
#            gene_lookup_file  Ensembl ID to gene symbol, position and biotype
#
# Outputs:   correction        Bonferroni factor per cell type, plus the number of
#                              unique testable genes across all cell types (TSV;
#                              the all-cell-types count is not used for correction)
#            sig_out           Significant genes, one sheet per disorder (xlsx)
#
# Notes:     ctwas_report.Rmd (14) computes the same factors independently,
#            because it runs before this pipeline.
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

message("\n\nBuilding significant cTWAS gene tables ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
traits            <- snakemake@params[["traits"]]
cell_types        <- snakemake@params[["cell_types"]]
in_dir            <- snakemake@params[["in_dir"]]
weights_dir       <- snakemake@params[["weights_dir"]]
gene_lookup_file  <- snakemake@params[["gene_lookup_file"]]
manuscript_labels <- snakemake@params[["labels"]]
pip_thresh        <- as.numeric(snakemake@params[["pip_thresh"]])
bf_p_thresh       <- as.numeric(snakemake@params[["bf_p_thresh"]])
correction_out    <- snakemake@output[["correction"]]
sig_out           <- snakemake@output[["sig_out"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("traits", "cell_types", "in_dir", "weights_dir", "gene_lookup_file",
               "labels", "pip_thresh", "bf_p_thresh", "correction_out", "sig_out"),
  value    = c(paste(traits, collapse = ", "), paste(cell_types, collapse = ", "),
               in_dir, weights_dir, gene_lookup_file,
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               pip_thresh, bf_p_thresh, correction_out, sig_out)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

# Full disorder names for sheet titles - order must match `traits`
disorder_labels <- c(
  scz  = "Schizophrenia",
  bpd  = "Bipolar disorder",
  mdd  = "Major depressive disorder",
  adhd = "ADHD",
  ocd  = "Obsessive-compulsive disorder"
)
missing_labels <- setdiff(traits, names(disorder_labels))
if (length(missing_labels) > 0) {
  stop("No sheet-name label defined for trait(s): ", paste(missing_labels, collapse = ", "),
       " -- add them to disorder_labels above.")
}

# gene_lookup_hg38.tsv has >1 row for some gene IDs; without distinct() the
# left_join below duplicates significant rows
gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  dplyr::select(chr = chromosome_name, start = start_position, end = end_position,
         gene_biotype, ensembl_id = ensembl_gene_id, symbol = external_gene_name) |>
  distinct(ensembl_id, .keep_all = TRUE)

##  Bonferroni factors per cell type  -------------------------------------------------
# Non-empty FUSION weight files (.wgt.RDat) per cell type: genes that were
# testable (hsq_p < 0.01), i.e. the number of tests
message("Computing BF correction factors (testable gene counts) per cell type...")

fusion_tbl <- map_dfr(cell_types, function(ct) {
  ct_dir <- file.path(weights_dir, ct)
  if (!dir.exists(ct_dir)) {
    message("  [", ct, "] weights dir not found: ", ct_dir)
    return(tibble())
  }
  files <- list.files(ct_dir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
  files <- files[file.info(files)$size > 0]
  tibble(cell_type = ct, gene_id = basename(files) |> str_remove("\\.wgt\\.RDat$"))
})

per_cell_counts <- fusion_tbl |>
  group_by(cell_type) |>
  summarise(n_sig_genes = n_distinct(gene_id), .groups = "drop")

unique_genes_all <- fusion_tbl |> distinct(gene_id) |> nrow()

correction_tbl <- per_cell_counts |>
  bind_rows(tibble(cell_type = "Unique across all cell types", n_sig_genes = unique_genes_all)) |>
  arrange(desc(n_sig_genes))

message("Testable genes per cell type (BF factors):")
for (i in seq_len(nrow(per_cell_counts))) {
  message("  [", per_cell_counts$cell_type[i], "] ", per_cell_counts$n_sig_genes[i])
}
message("Unique across all cell types: ", unique_genes_all)

write_tsv(correction_tbl, correction_out)
message("Wrote correction factors: ", correction_out)

##  Significant fine-mapped genes per trait  ------------------------------------------
get_ctwas_sig <- function(ct, gwas, in_dir, gene_lookup, per_cell_counts,
                           pip_thresh, bf_p_thresh) {
  file <- file.path(in_dir, "output", paste0("ctwas_", ct, "_", gwas, "_ctwas.rds"))
  if (!file.exists(file)) {
    message("  [", ct, "] file not found: ", file)
    return(NULL)
  }

  ctw_res <- read_rds(file)

  if (!is.null(ctw_res$status) && ctw_res$status == "skipped") {
    message("  [", ct, "] skipped: ", ctw_res$reason %||% "no reason given")
    return(NULL)
  }
  if (is.null(ctw_res$finemap_res) || nrow(ctw_res$finemap_res) == 0) {
    message("  [", ct, "] no finemapped genes (no credible sets identified)")
    return(NULL)
  }

  bf_factor <- per_cell_counts |> filter(cell_type == ct) |> pull(n_sig_genes)
  if (length(bf_factor) == 0) {
    message("  [", ct, "] no BF correction factor available (no FUSION weights) - skipping")
    return(NULL)
  }

  finemap_res <- as_tibble(ctw_res$finemap_res) |>
    # A cell type with no credible sets has an all-NA logical cs column,
    # which cannot be bound to character cs IDs from other cell types
    mutate(cs = as.character(cs)) |>
    filter(group != 'SNP' & susie_pip > pip_thresh & !is.na(cs)) |>
    separate(id, into = c('gene', 'id'), sep = "\\|") |>
    separate(id, into = c('cell_type', 'id'), sep = "_") |>
    left_join(gene_lookup, by = join_by(gene == ensembl_id)) |>
    mutate(symbol = na_if(symbol, "")) |>  # blank symbol -> NA, not ""
    mutate(pval = 2 * pnorm(-abs(z)),
           pval_bf = pmin(pval * bf_factor, 1)) |>
    filter(pval_bf < bf_p_thresh) |>
    arrange(pval_bf) |>
    mutate(across(c(pval, pval_bf, susie_pip), ~ signif(.x, 2)), z = round(z, 2)) |>
    select(cell_type, symbol, gene, z, pval, pval_bf, susie_pip, cs, gene_biotype, chr)

  message("  [", ct, "] ", nrow(finemap_res), " significant genes")
  finemap_res
}

`%||%` <- function(a, b) if (is.null(a)) b else a

empty_sig_tbl <- tibble(
  cell_type = character(), symbol = character(), gene = character(),
  z = double(), pval = double(), pval_bf = double(), susie_pip = double(),
  cs = character(), gene_biotype = character(), chr = character()
)

# Final column names for the xlsx, applied to every sheet (including empty
# ones, so every sheet has the same header)
format_sig_tbl <- function(tbl) {
  tbl |>
    mutate(cell_type = relabel_cell_types(cell_type, manuscript_labels)) |>
    rename(
      `Cell type`        = cell_type,
      `Gene symbol`      = symbol,
      `Ensembl ID`       = gene,
      `Z-score`          = z,
      `P (uncorrected)`  = pval,
      `Bonferroni P`     = pval_bf,
      PIP                = susie_pip,
      CS                 = cs,
      Type               = gene_biotype,
      Chr                = chr
    )
}

sig_sheets <- list()

for (i in seq_along(traits)) {
  gwas <- traits[i]
  message("--- ", gwas, " ---")

  gwas_sig_results <- map_dfr(cell_types, get_ctwas_sig, gwas = gwas, in_dir = in_dir,
                               gene_lookup = gene_lookup, per_cell_counts = per_cell_counts,
                               pip_thresh = pip_thresh, bf_p_thresh = bf_p_thresh)

  if (is.null(gwas_sig_results) || ncol(gwas_sig_results) == 0) {
    message(gwas, ": no readable cTWAS results for any cell type")
    gwas_sig_results <- empty_sig_tbl
  } else {
    message(gwas, ": ", nrow(gwas_sig_results), " significant (gene, cell type) rows across ",
            n_distinct(gwas_sig_results$cell_type), " cell types, ",
            n_distinct(gwas_sig_results$gene), " unique genes")
  }

  sig_sheets[[disorder_labels[[gwas]]]] <- format_sig_tbl(gwas_sig_results)
}

##  Write combined xlsx  --------------------------------------------------------------
message("Writing combined cTWAS significant gene table: ", sig_out)
write.xlsx(sig_sheets, file = sig_out, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Wrote: ", sig_out)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
