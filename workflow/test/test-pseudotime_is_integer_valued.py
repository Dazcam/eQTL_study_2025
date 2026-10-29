# --------------------------------------------------------------------------------------
#
#    test-pseudotime_is_integer_valued.py
#    - Used in pseudotime_pseudobulk.py
#
# --------------------------------------------------------------------------------------

import numpy as np
import scipy.sparse as sp

from pseudotime_functions import is_integer_valued

INT = np.array([[0, 2, 0], [5, 0, 1]], dtype=float)
FLT = np.array([[0, 2.5, 0], [5, 0, 1]], dtype=float)


def test_integer_dense_and_sparse():
    assert is_integer_valued(INT) is True
    assert is_integer_valued(sp.csr_matrix(INT)) is True


def test_non_integer_dense_and_sparse():
    assert is_integer_valued(FLT) is False
    assert is_integer_valued(sp.csr_matrix(FLT)) is False


def test_all_zero_counts_as_integer():
    assert is_integer_valued(np.zeros((2, 2))) is True
