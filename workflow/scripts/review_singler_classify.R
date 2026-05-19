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
#    - Outputs per-cell labels and full score matrix as CSV
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

log_msg(paste("SingleR version:", as.character(packageVersion("SingleR"))))
log_msg(paste("anndataR version:", as.character(packageVersion("anndataR"))))

# ── Inputs / outputs ──────────────────────────────────────────────────────────────────

ref_h5ad    <- snakemake@input[["ref_h5ad"]]
query_h5ad  <- snakemake@input[["query_h5ad"]]
output_csv  <- snakemake@output[[1]]
n_threads   <- snakemake@resources[["threads"]]

log_msg(paste("Reference h5ad:", ref_h5ad))
log_msg(paste("Query h5ad:    ", query_h5ad))
log_msg(paste("Output CSV:    ", output_csv))
log_msg(paste("Threads:       ", n_threads))

# Set up parallel backend
BPPARAM <- MulticoreParam(n_threads)

# ── Load reference ─────────────────────────────────────────────────────────────────────

log_msg("Loading Qian MERFISH reference...")
ref_adata <- anndataR::read_h5ad(ref_h5ad)
ref_sce   <- ref_adata$as_SingleCellExperiment()
log_msg(paste("Reference shape:", nrow(ref_sce), "genes x", ncol(ref_sce), "cells"))
log_msg(paste("Reference assay names:", paste(assayNames(ref_sce), collapse = ", ")))
ref_counts <- assay(ref_sce, "X")  # replace "X" with actual assay name once confirmed
log_msg(paste("Reference count matrix class:", class(ref_counts)))
log_msg(paste("Reference count dtype:", typeof(ref_counts)))
log_msg(paste("Reference max value:", max(ref_counts)))
log_msg(paste("Reference counts are whole numbers:", 
              all(ref_counts@x == floor(ref_counts@x))))

# Check H2 annotation column
if (!"H2_annotation" %in% colnames(colData(ref_sce))) {
    stop("H2_annotation column not found in reference colData")
}
ref_labels <- ref_sce$H2_annotation
log_msg(paste("Reference H2 clusters:", length(unique(ref_labels))))
log_msg("H2 cluster cell counts:")
print(sort(table(ref_labels)))

# Normalise reference
log_msg("Normalising reference with logNormCounts...")
ref_sce <- scuttle::logNormCounts(ref_sce, assay.type = "X")

# ── Load query ────────────────────────────────────────────────────────────────────────

log_msg("Loading snRNA-seq query...")
query_adata <- anndataR::read_h5ad(query_h5ad)
query_sce   <- query_adata$as_SingleCellExperiment()
log_msg(paste("Query shape:", nrow(query_sce), "genes x", ncol(query_sce), "cells"))
log_msg(paste("Query assay names:", paste(assayNames(query_sce), collapse = ", ")))
query_counts <- assay(query_sce, "X")
log_msg(paste("Query count matrix class:", class(query_counts)))
log_msg(paste("Query count dtype:", typeof(query_counts)))
log_msg(paste("Query max value:", max(query_counts)))
log_msg(paste("Query counts are whole numbers:", 
              all(query_counts@x == floor(query_counts@x))))

# Check cell_type column
if (!"cell_type" %in% colnames(colData(query_sce))) {
    stop("cell_type column not found in query colData")
}
log_msg("Query L1 cell type distribution:")
print(sort(table(query_sce$cell_type)))

# Normalise query
log_msg("Normalising query with logNormCounts...")
query_sce <- scuttle::logNormCounts(query_sce, assay.type = "X")

# ── Confirm shared genes ───────────────────────────────────────────────────────────────

shared_genes <- intersect(rownames(ref_sce), rownames(query_sce))
log_msg(paste("Genes in reference:          ", nrow(ref_sce)))
log_msg(paste("Genes in query:              ", nrow(query_sce)))
log_msg(paste("Shared genes for transfer:   ", length(shared_genes),
              sprintf("(%.1f%% of MERFISH panel)", 
                      100 * length(shared_genes) / nrow(ref_sce))))

if (length(shared_genes) < 200) {
    warning(paste("Only", length(shared_genes), 
                  "shared genes — label transfer quality may be reduced"))
}

# Subset both reference and query to shared genes only
log_msg("Subsetting reference and query to shared genes...")
ref_sce   <- ref_sce[shared_genes, ]
query_sce <- query_sce[shared_genes, ]
log_msg(paste("Reference shape after subsetting:", nrow(ref_sce), "genes x", ncol(ref_sce), "cells"))
log_msg(paste("Query shape after subsetting:", nrow(query_sce), "genes x", ncol(query_sce), "cells"))

# ── Train SingleR reference ────────────────────────────────────────────────────────────

log_msg("Training SingleR reference...")
trained_ref <- SingleR::trainSingleR(
    ref        = ref_sce,
    labels     = ref_labels,
    assay.type = "logcounts",
    BPPARAM    = BPPARAM
)
log_msg("Training complete")

# ── Classify query cells ───────────────────────────────────────────────────────────────

log_msg("Classifying query cells with SingleR...")
log_msg(paste("Classifying", ncol(query_sce), "cells..."))

pred <- SingleR::classifySingleR(
    test       = query_sce,
    trained    = trained_ref,
    assay.type = "logcounts",
    BPPARAM    = BPPARAM
)
log_msg("Classification complete")

# ── Classification summary ─────────────────────────────────────────────────────────────

log_msg("Predicted H2 label distribution:")
print(sort(table(pred$labels)))

n_pruned <- sum(is.na(pred$pruned.labels))
log_msg(paste("Cells pruned (low confidence):", n_pruned,
              sprintf("(%.1f%% of total)", 100 * n_pruned / ncol(query_sce))))

# Pruning rate per L1 population
log_msg("Pruning rate per L1 population:")
pruning_df <- data.frame(
    cell_type = query_sce$cell_type,
    pruned    = is.na(pred$pruned.labels)
)
pruning_summary <- aggregate(pruned ~ cell_type, data = pruning_df, 
                              FUN = function(x) sprintf("%.1f%%", 100 * mean(x)))
print(pruning_summary)

# ── Build results dataframe ────────────────────────────────────────────────────────────

log_msg("Building results dataframe...")
results_df <- data.frame(
    cell_id      = colnames(query_sce),
    cell_type    = query_sce$cell_type,
    H2_label     = pred$labels,
    H2_pruned    = pred$pruned.labels,
    stringsAsFactors = FALSE
)

# Append full score matrix - cells x H2 clusters
scores_df <- as.data.frame(pred$scores)
results_df <- cbind(results_df, scores_df)

log_msg(paste("Results dataframe dimensions:", nrow(results_df), "x", ncol(results_df)))

# ── Write output ───────────────────────────────────────────────────────────────────────

log_msg(paste("Writing results to:", output_csv))
write.csv(results_df, output_csv, row.names = FALSE)
log_msg(paste("Written", nrow(results_df), "rows"))

log_msg("SingleR classification complete.")

# Close log connection
sink(type = "message")
sink()

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
