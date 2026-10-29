# --------------------------------------------------------------------------------------
#
#    Pseudotime helper functions
#
# --------------------------------------------------------------------------------------
#
# Functions: get_gene_expr()       One gene's expression from adata.X as a flat array
#            select_anchor_cell()  Max-marker-expression cell in a cell type, with a
#                                  multi-scale-space centroid fallback (root/terminal)
#            is_integer_valued()   Check all non-zero values in a matrix are integers
#
# Packages:  numpy, scipy (no scanpy/palantir needed, so tests stay light)
#
# Used in:   pseudotime_trajectory.py (get_gene_expr, select_anchor_cell);
#            pseudotime_pseudobulk.py (is_integer_valued)
#
# --------------------------------------------------------------------------------------

import logging

import numpy as np
import scipy.sparse as sp

log = logging.getLogger()


def get_gene_expr(adata, gene, cells=None):
    """Return expression of `gene` from adata.X as a 1D dense array.

    cells: optional integer index array to subset rows before densifying.
    Returns None if the gene is not in adata.var_names.
    """
    if gene not in adata.var_names:
        return None
    gene_col = adata.var_names.get_loc(gene)
    expr = adata.X[:, gene_col] if cells is None else adata.X[cells, gene_col]
    if sp.issparse(expr):
        return np.array(expr.todense()).flatten()
    return np.array(expr).flatten()


def select_anchor_cell(adata, cell_type, gene, ms_key="X_palantir_ms"):
    """Pick the cell of `cell_type` with the highest `gene` expression.

    Falls back to the cell closest to the cell-type centroid in adata.obsm[ms_key]
    if the gene is absent or has zero expression in every cell of that type.

    Returns (cell_id, method), where method is a string for the log and uns.
    """
    mask = adata.obs["cell_type"] == cell_type
    if mask.sum() == 0:
        raise ValueError(f"No cells found with cell_type == '{cell_type}'")
    idx = np.where(mask)[0]

    expr = get_gene_expr(adata, gene, cells=idx)
    if expr is None:
        log.warning(f"{gene} not found in adata.var_names — falling back to centroid")
    else:
        max_expr = expr.max()
        if max_expr > 0:
            global_idx = idx[np.argmax(expr)]
            return adata.obs_names[global_idx], f"{gene} expression (max={max_expr:.3f})"
        log.warning(
            f"{gene} expression is zero for all {cell_type} cells "
            "— falling back to centroid"
        )

    ms       = adata.obsm[ms_key][idx]
    centroid = ms.mean(axis=0)
    dists    = np.linalg.norm(ms - centroid, axis=1)
    local    = np.argmin(dists)
    return (adata.obs_names[idx[local]],
            f"{cell_type} multi-scale centroid (dist={dists[local]:.4f})")


def is_integer_valued(x):
    """True if every non-zero value in x (sparse or dense) is a whole number."""
    vals = x.data if sp.issparse(x) else np.asarray(x).ravel()
    return bool(np.all(np.mod(vals[vals != 0], 1) == 0))
