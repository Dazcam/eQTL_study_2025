#!/usr/bin/env Rscript
# --------------------------------------------------------------------------------------
#
#    singler_classify.R
#
#    Runs SingleR label transfer from Qian et al. 2025 MERFISH H2 reference
#    onto snRNA-seq query data.
#
#    - Loads reference and query h5ad files via anndataR
#    - Normalises both with scuttle::logNormCounts()
#    - Trains SingleR reference from Qian H2 annotations
#    - Classifies all query cells
#    - Outputs:
#        - CSV: per-cell labels, deltas, and full score matrix
#        - RDS: full pred object for diagnostic plots in Rmd
#        - RDS: subsampled query SCE (200 cells/L1) for score heatmap
#        - RDS: cell type metadata for concordance analysis
#
#    Usage: called via Snakemake rule singler_classify
#
# --------------------------------------------------------------------------------------

# ── Logging setup ─────────────────────────────────────────────────────────────────────

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, append = TRUE)
sink(con, append = TRUE, type = "message")

log_msg <- function(msg) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"))
}

log_msg("Starting SingleR classification")
log_msg(paste("R version:", R.version$version.string))

# ── Libraries ─────────────────────────────────────────────────────────────────────────

log_msg("Loading libraries...")
suppressPackageStartupMessages({
    library(SingleR)
    library(SingleCellExperiment)
    library(scuttle)
    library(anndataR)
    library(BiocParallel)
    library(Matrix)
})

log_msg(paste("SingleR version:",  as.character(packageVersion("SingleR"))))
log_msg(paste("anndataR version:", as.character(packageVersion("anndataR"))))

# ── Inputs / outputs ──────────────────────────────────────────────────────────────────

ref_h5ad       <- as.character(snakemake@input[["ref_h5ad"]])
query_h5ad     <- as.character(snakemake@input[["query_h5ad"]])
output_csv     <- as.character(snakemake@output[["csv"]])
output_pred    <- as.character(snakemake@output[["pred_rds"]])
output_sce     <- as.character(snakemake@output[["sce_rds"]])
n_threads      <- snakemake@resources[["threads"]]
n_subsample    <- 200L  # cells per L1 population for score heatmap

log_msg(paste("Reference h5ad:  ", ref_h5ad))
log_msg(paste("Query h5ad:      ", query_h5ad))
log_msg(paste("Output CSV:      ", output_csv))
log_msg(paste("Output pred RDS: ", output_pred))
log_msg(paste("Output SCE RDS:  ", output_sce))
log_msg(paste("Threads:         ", n_threads))
log_msg(paste("Subsample N/pop: ", n_subsample))

BPPARAM <- MulticoreParam(n_threads)

# ── Load reference ─────────────────────────────────────────────────────────────────────

log_msg("Loading Qian MERFISH reference...")
ref_adata <- anndataR::read_h5ad(ref_h5ad)
ref_sce   <- ref_adata$as_SingleCellExperiment()
log_msg(paste("Reference shape:", nrow(ref_sce), "genes x", ncol(ref_sce), "cells"))
log_msg(paste("Reference assay names:", paste(assayNames(ref_sce), collapse = ", ")))

# Sanity check counts
ref_counts <- assay(ref_sce, "X")
log_msg(paste("Reference max value:", max(ref_counts)))
log_msg(paste("Reference counts are whole numbers:",
              all(ref_counts@x == floor(ref_counts@x))))

if (!"H2_annotation" %in% colnames(colData(ref_sce))) {
    stop("H2_annotation column not found in reference colData")
}
ref_labels <- ref_sce$H2_annotation
log_msg(paste("Reference H2 clusters:", length(unique(ref_labels))))
log_msg("H2 cluster cell counts:")
print(sort(table(ref_labels)))

log_msg("Normalising reference with logNormCounts...")
ref_sce <- scuttle::logNormCounts(ref_sce, assay.type = "X")

# ── Load query ────────────────────────────────────────────────────────────────────────

log_msg("Loading snRNA-seq query...")
query_adata <- anndataR::read_h5ad(query_h5ad)
query_sce   <- query_adata$as_SingleCellExperiment()
log_msg(paste("Query shape:", nrow(query_sce), "genes x", ncol(query_sce), "cells"))
log_msg(paste("Query assay names:", paste(assayNames(query_sce), collapse = ", ")))

# Sanity check counts
query_counts <- assay(query_sce, "X")
log_msg(paste("Query max value:", max(query_counts)))
log_msg(paste("Query counts are whole numbers:",
              all(query_counts@x == floor(query_counts@x))))

if (!"cell_type" %in% colnames(colData(query_sce))) {
    stop("cell_type column not found in query colData")
}
log_msg("Query L1 cell type distribution:")
print(sort(table(query_sce$cell_type)))

log_msg("Normalising query with logNormCounts...")
query_sce <- scuttle::logNormCounts(query_sce, assay.type = "X")

# ── Shared genes ──────────────────────────────────────────────────────────────────────

shared_genes <- intersect(rownames(ref_sce), rownames(query_sce))
log_msg(paste("Genes in reference:         ", nrow(ref_sce)))
log_msg(paste("Genes in query:             ", nrow(query_sce)))
log_msg(paste("Shared genes for transfer:  ", length(shared_genes),
              sprintf("(%.1f%% of MERFISH panel)",
                      100 * length(shared_genes) / nrow(ref_sce))))

if (length(shared_genes) < 200) {
    warning(paste("Only", length(shared_genes),
                  "shared genes — label transfer quality may be reduced"))
}

log_msg("Subsetting reference and query to shared genes...")
ref_sce   <- ref_sce[shared_genes, ]
query_sce <- query_sce[shared_genes, ]
log_msg(paste("Reference shape after subsetting:", nrow(ref_sce), "x", ncol(ref_sce)))
log_msg(paste("Query shape after subsetting:    ", nrow(query_sce), "x", ncol(query_sce)))

# ── Train ─────────────────────────────────────────────────────────────────────────────

log_msg("Training SingleR reference...")
trained_ref <- SingleR::trainSingleR(
    ref        = ref_sce,
    labels     = ref_labels,
    assay.type = "logcounts",
    BPPARAM    = BPPARAM
)
log_msg("Training complete")

# ── Classify ──────────────────────────────────────────────────────────────────────────

log_msg(paste("Classifying", ncol(query_sce), "cells..."))
pred <- SingleR::classifySingleR(
    test       = query_sce,
    trained    = trained_ref,
    assay.type = "logcounts",
    BPPARAM    = BPPARAM
)
log_msg("Classification complete")

# ── Summary ───────────────────────────────────────────────────────────────────────────

log_msg("Predicted H2 label distribution:")
print(sort(table(pred$labels)))

n_pruned <- sum(is.na(pred$pruned.labels))
log_msg(paste("Cells pruned (low confidence):", n_pruned,
              sprintf("(%.1f%% of total)", 100 * n_pruned / ncol(query_sce))))

log_msg("Pruning rate per L1 population:")
pruning_df <- data.frame(
    cell_type = query_sce$cell_type,
    pruned    = is.na(pred$pruned.labels)
)
pruning_summary <- aggregate(pruned ~ cell_type, data = pruning_df,
                              FUN = function(x) sprintf("%.1f%%", 100 * mean(x)))
print(pruning_summary)

# ── Compute per-cell deltas ───────────────────────────────────────────────────────────
# Delta = score for assigned label minus median score across all labels
# Low delta = ambiguous assignment

log_msg("Computing per-cell deltas...")
score_mat    <- pred$scores
assigned_idx <- match(pred$labels, colnames(score_mat))
assigned_scores <- score_mat[cbind(seq_len(nrow(score_mat)), assigned_idx)]
median_scores   <- matrixStats::rowMedians(score_mat)
deltas          <- assigned_scores - median_scores

log_msg(paste("Delta summary:"))
log_msg(paste("  Mean:  ", round(mean(deltas, na.rm = TRUE), 4)))
log_msg(paste("  Median:", round(median(deltas, na.rm = TRUE), 4)))
log_msg(paste("  Min:   ", round(min(deltas,  na.rm = TRUE), 4)))
log_msg(paste("  Max:   ", round(max(deltas,  na.rm = TRUE), 4)))

# ── Build CSV ─────────────────────────────────────────────────────────────────────────

log_msg("Building results dataframe...")
results_df <- data.frame(
    cell_id    = colnames(query_sce),
    cell_type  = query_sce$cell_type,
    H2_label   = pred$labels,
    H2_pruned  = pred$pruned.labels,
    delta      = deltas,
    stringsAsFactors = FALSE
)
scores_df  <- as.data.frame(pred$scores)
results_df <- cbind(results_df, scores_df)
log_msg(paste("Results dataframe dimensions:", nrow(results_df), "x", ncol(results_df)))

log_msg(paste("Writing CSV to:", output_csv))
write.csv(results_df, output_csv, row.names = FALSE)
log_msg(paste("Written", nrow(results_df), "rows"))

# ── Save full pred object ─────────────────────────────────────────────────────────────
# Required for plotScoreHeatmap() and plotDeltaDistribution() in the Rmd report

log_msg(paste("Saving pred RDS to:", output_pred))
saveRDS(pred, output_pred)
log_msg("pred RDS saved")

# ── Save subsampled query SCE ─────────────────────────────────────────────────────────
# 200 cells per L1 population for plotScoreHeatmap() - full 679k is too large to render

log_msg(paste("Subsampling query SCE:", n_subsample, "cells per L1 population..."))
set.seed(42L)
cell_types  <- unique(query_sce$cell_type)
keep_idx    <- unlist(lapply(cell_types, function(ct) {
    idx <- which(query_sce$cell_type == ct)
    sample(idx, min(n_subsample, length(idx)))
}))
keep_idx    <- sort(keep_idx)

query_sub   <- query_sce[, keep_idx]
pred_sub    <- pred[keep_idx, ]

# Store cell_type and H2 label in SCE colData for use in diagnostic plots
query_sub$H2_label    <- pred_sub$labels
query_sub$H2_pruned   <- pred_sub$pruned.labels
query_sub$delta       <- deltas[keep_idx]

log_msg(paste("Subsampled SCE shape:", nrow(query_sub), "genes x", ncol(query_sub), "cells"))
log_msg(paste("Subsampled cell type distribution:"))
print(sort(table(query_sub$cell_type)))

log_msg(paste("Saving subsampled SCE RDS to:", output_sce))
saveRDS(list(query_sub = query_sub, pred_sub = pred_sub), output_sce)
log_msg("SCE RDS saved")

log_msg("SingleR classification complete.")

sink(type = "message")
sink()
