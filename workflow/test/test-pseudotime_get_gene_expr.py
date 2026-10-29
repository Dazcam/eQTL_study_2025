# --------------------------------------------------------------------------------------
#
#    test-pseudotime_get_gene_expr.py
#    - Used in pseudotime_trajectory.py
#
# --------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from pseudotime_functions import get_gene_expr

X = np.array([[0, 1, 2],
              [3, 0, 5],
              [6, 7, 0],
              [9, 0, 1]], dtype=float)
VAR = pd.DataFrame(index=["PAX6", "TLE4", "SATB2"])


def make(sparse):
    return ad.AnnData(X=sp.csr_matrix(X) if sparse else X.copy(), var=VAR)


def test_dense_and_sparse_give_same_flat_array():
    for sparse in (False, True):
        out = get_gene_expr(make(sparse), "TLE4")
        assert out.ndim == 1
        np.testing.assert_array_equal(out, [1, 0, 7, 0])


def test_cell_subset():
    for sparse in (False, True):
        out = get_gene_expr(make(sparse), "SATB2", cells=np.array([1, 3]))
        np.testing.assert_array_equal(out, [5, 1])


def test_missing_gene_returns_none():
    assert get_gene_expr(make(False), "FOXP2") is None
