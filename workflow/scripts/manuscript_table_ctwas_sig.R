#--------------------------------------------------------------------------------------
#
#    Generate cTWAS significant gene tables for manuscript
#
#--------------------------------------------------------------------------------------

# Two outputs:
#   A. Per-cell-type BF correction factors (count of cTWAS-testable genes,
#      i.e. non-empty FUSION weight files) - ported from ctwas_report.Rmd's
#      "Correction" tab (Table 1). Needed as an input to B, not just useful
#      on its own. Still written as its own TSV.
#   B. Significant fine-mapped genes (PIP > 0.8, in a credible set,
#      BF-corrected P < 0.05), unioned across all cell types - ported from
#      the "Single (Fine-mapped only)" tab's per-(gwas,cell_type) loop. That
#      tab built the same results_single_lst this script builds, but its
#      write.xlsx() (one sheet per gwas x cell_type) was disabled ("Not
#      working yet"). Replaced here with a single xlsx, one sheet per GWAS
#      trait (named by full disorder name), formatted column names, and
#      cell types relabelled Glu-UL/Glu-DL -> Glu-A/Glu-B.
#
# Independent of ctwas_report.Rmd (same reasoning as smr_sig_tbl vs
# smr_report.Rmd): that report renders as part of 12ctwas.smk, which runs
# BEFORE the tables pipeline where this script lives, so it can't depend on
# this script's output without inverting that ordering. Both compute the
# same correction-factor / significance logic independently.

## Info  ------------------------------------------------------------------------------

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

message("\n\nGenerating cTWAS significant gene tables for the manuscript ...")

# -------------------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# --- Set variables
traits           <- snakemake@params[["traits"]]
cell_types       <- snakemake@params[["cell_types"]]
in_dir           <- snakemake@params[["in_dir"]]
weights_dir      <- snakemake@params[["weights_dir"]]
gene_lookup_file <- snakemake@params[["gene_lookup_file"]]
pip_thresh       <- as.numeric(snakemake@params[["pip_thresh"]])
bf_p_thresh      <- as.numeric(snakemake@params[["bf_p_thresh"]])

correction_out <- snakemake@output[["correction"]]
sig_out        <- snakemake@output[["sig_out"]]  # single xlsx, one sheet per trait

# traits           <- c("scz", "bpd", "mdd", "adhd", "ocd")
# cell_types       <- c("Glu-UL", "Glu-DL", ..., "NPC-to-Glu-UL-Q4-Bin4")
# in_dir           <- "../results/12CTWAS/"
# weights_dir      <- "../results/12CTWAS/weights"
# gene_lookup_file <- "../resources/sheets/gene_lookup_hg38.tsv"
# pip_thresh       <- 0.8
# bf_p_thresh      <- 0.05

message("Cell types (", length(cell_types), "): ", paste(cell_types, collapse = ", "))
message("Traits (", length(traits), "): ", paste(traits, collapse = ", "))

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

# Relabel cell type names for display: Glu-UL -> Glu-A, Glu-DL -> Glu-B
# (covers L2 subtypes too, e.g. Glu-UL-0 -> Glu-A-0)
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

gene_lookup <- suppressMessages(read_tsv(gene_lookup_file, show_col_types = FALSE)) |>
  select(chr = chromosome_name, start = start_position, end = end_position,
         gene_biotype, ensembl_id = ensembl_gene_id, symbol = external_gene_name) |>
  # gene_lookup_hg38.tsv has >1 row for some gene IDs; without this the
  # left_join below fans out and duplicates significant rows (same fix
  # applied in manuscript_table_smr_sig.R for the same underlying file)
  distinct(ensembl_id, .keep_all = TRUE)

# ---- A. Per-cell-type BF correction factors ------------------------------------
# Count of non-empty FUSION weight files (.wgt.RDat) per cell type - genes
# that were actually testable (hsq_p < 0.01 passed), i.e. the number of
# tests, not the number that came out significant.
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

# ---- B. Per-GWAS significant fine-mapped genes, unioned across cell types -----
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
    mutate(cs = as.character(cs)) |>  # a cell type with 0 real cs values gets an
                                       # all-NA logical column from R by default,
                                       # which then can't bind_rows() against a
                                       # cell type with genuine character cs IDs
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

# Final column names/order for the xlsx, applied uniformly to every sheet
# (including the empty-result case, so every sheet has the same header)
format_sig_tbl <- function(tbl) {
  tbl |>
    mutate(cell_type = relabel_cell_type(cell_type)) |>
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

message("Writing combined cTWAS significant gene table: ", sig_out)
write.xlsx(sig_sheets, file = sig_out, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Wrote: ", sig_out)

message("Done: ", paste(traits, collapse = ", "))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
