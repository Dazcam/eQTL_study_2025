# --------------------------------------------------------------------------------------
#
#    test-pseudotime_select_anchor_cell.py
#    - Used in pseudotime_trajectory.py
#
# --------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import anndata as ad

from pseudotime_functions import select_anchor_cell

# 5 cells: c0-c2 NPC, c3-c4 Glu-DL. Gene columns: PAX6, TLE4, ZERO
X = np.array([[2.0, 0.0, 0.0],
              [5.0, 0.0, 0.0],   # NPC PAX6 max
              [1.0, 9.0, 0.0],   # NPC with highest TLE4 overall: must be ignored for Glu-DL
              [0.0, 3.0, 0.0],
              [0.0, 4.0, 0.0]])  # Glu-DL TLE4 max
OBS = pd.DataFrame({"cell_type": pd.Categorical(
    ["NPC", "NPC", "NPC", "Glu-DL", "Glu-DL"])},
    index=[f"c{i}" for i in range(5)])
VAR = pd.DataFrame(index=["PAX6", "TLE4", "ZERO"])
# Multi-scale space: NPC centroid (1, 1) is closest to c1; Glu-DL centroid tie -> first (c3)
MS = np.array([[0.0, 0.0],
               [1.0, 1.0],
               [2.0, 2.0],
               [5.0, 5.0],
               [6.0, 6.0]])


def make(sparse=False):
    a = ad.AnnData(X=sp.csr_matrix(X) if sparse else X.copy(), obs=OBS.copy(), var=VAR)
    a.obsm["X_palantir_ms"] = MS
    return a


@pytest.mark.parametrize("sparse", [False, True])
def test_max_expression_within_cell_type(sparse):
    assert select_anchor_cell(make(sparse), "NPC", "PAX6") == (
        "c1", "PAX6 expression (max=5.000)")
    assert select_anchor_cell(make(sparse), "Glu-DL", "TLE4") == (
        "c4", "TLE4 expression (max=4.000)")


def test_missing_gene_falls_back_to_centroid():
    cell, method = select_anchor_cell(make(), "NPC", "FOXP2")
    assert cell == "c1"
    assert method == "NPC multi-scale centroid (dist=0.0000)"


def test_zero_expression_falls_back_to_centroid():
    cell, method = select_anchor_cell(make(), "Glu-DL", "ZERO")
    assert cell == "c3"
    assert method.startswith("Glu-DL multi-scale centroid (dist=0.7071)")


def test_absent_cell_type_raises():
    with pytest.raises(ValueError, match="No cells found"):
        select_anchor_cell(make(), "OPC", "PAX6")
