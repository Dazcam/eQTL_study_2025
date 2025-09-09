#!/usr/bin/env python3
import pandas as pd
import glob
import os
import sys

out_dir = sys.argv[1]

cell_types = [
    "Astrocytes",
    "Endothelial.cells",
    "Excitatory.neurons",
    "Inhibitory.neurons",
    "Microglia",
    "Oligodendrocytes",
    "OPCs...COPs",
    "Pericytes",
]

for cell in cell_types:
    files = sorted(glob.glob(f"{out_dir}/{cell}.[0-9]*.gz"))
    if files:
        out_cell = "OPCs" if cell == "OPCs...COPs" else cell.replace(".", "_")
        df = pd.concat(
            [
                pd.read_table(
                    f,
                    compression="gzip",
                    header=None,
                    names=["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"],
                )
                for f in files
            ]
        )
        df.to_csv(
            f"{out_dir}/{out_cell}_nom.cis_qtl.txt.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        for f in files:
            os.remove(f)
