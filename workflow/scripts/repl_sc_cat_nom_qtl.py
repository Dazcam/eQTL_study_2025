#!/usr/bin/env python3
#--------------------------------------------------------------------------------------
#
#    Concatenate per-chromosome tensorQTL nominal outputs
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  06QTL-REPLICATION | Rule: cat_nom_qtl (one job per cell type)
#            Upstream:   05TENSORQTL tensorqtl_nominal (chr 1-22 parquet files)
#            Downstream: pi1_enrich_internal; pi1_enrich_jang (qtl_all)
#
# Purpose:   - Read all per-chromosome nominal parquet files for a cell type
#            - Rename start_distance -> tss_distance and af -> maf
#            - Write the selected columns to one gzipped TSV
#
# Usage:     repl_sc_cat_nom_qtl.py <output> <parquet1> <parquet2> ...
#
# Inputs:    parquets   tensorQTL nominal parquet files, one per chromosome
#
# Outputs:   output     gzipped TSV: phenotype_id, variant_id, tss_distance,
#                       ma_samples, ma_count, maf, pval_nominal, slope, slope_se
#
# Notes:     Called from a shell rule; stderr is redirected to the rule log.
#
#--------------------------------------------------------------------------------------

import logging
import os
import sys

import pandas as pd

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

COLS = [
    "phenotype_id", "variant_id", "tss_distance", "ma_samples",
    "ma_count", "maf", "pval_nominal", "slope", "slope_se",
]


def main(parquets, output):
    log.info("Reading %d files ...", len(parquets))
    df = pd.concat([pd.read_parquet(f) for f in parquets], ignore_index=True)

    # Rename columns to standardise names
    df = df.rename(columns={"start_distance": "tss_distance", "af": "maf"})

    df[COLS].to_csv(output, sep="\t", index=False, compression="gzip")
    log.info("Wrote %d rows to %s", len(df), output)


if __name__ == "__main__":
    output = sys.argv[1]
    parquets = sys.argv[2:]

    log.info("Variables")
    log.info("============================")
    log.info("output:     %s", output)
    log.info("n_parquets: %d", len(parquets))
    log.info("============================")

    os.makedirs(os.path.dirname(output), exist_ok=True)
    main(parquets, output)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
