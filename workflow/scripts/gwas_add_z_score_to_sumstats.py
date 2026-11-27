#--------------------------------------------------------------------------------------
#
#    Add signed Z scores to GWAS summary statistics
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09PREP-GWAS | Rule: add_z_score (one job per GWAS)
#            Upstream:   standardise_sumstats
#            Downstream: add_N
#
# Purpose:   Keep an existing Z column if it has any non-missing values; otherwise
#            compute Z from BETA / SE, log(OR) / SE or sign(log(OR)) x |Z| from P
#
# Inputs:    --sumstats  standardised sumstats (tsv)
#
# Outputs:   --out       sumstats with a Z column (tsv)
#
#--------------------------------------------------------------------------------------

import argparse
import logging
import sys

import numpy as np
import pandas as pd

logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Adding Z scores to GWAS sumstats ...")

##  Arguments  -------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--sumstats", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

variables = {"sumstats": args.sumstats, "out": args.out}
width = max(len(k) for k in variables)
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  Add Z  -----------------------------------------------------------------------------
df = pd.read_csv(args.sumstats, sep='\t', dtype=str, low_memory=False)
log.info(f"Read {len(df):,} rows from {args.sumstats}")

for c in ['BETA', 'SE', 'OR', 'PVAL', 'Z']:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')

if 'Z' in df.columns and df['Z'].notna().any():
    log.info("Z column already present, keeping it")
else:
    if 'BETA' in df.columns and 'SE' in df.columns:
        df['Z'] = df['BETA'] / df['SE']
        log.info("Computed Z = BETA / SE")
    elif 'OR' in df.columns and 'SE' in df.columns:
        df['Z'] = np.log(df['OR']) / df['SE']
        log.info("Computed Z = log(OR) / SE")
    elif 'OR' in df.columns and 'PVAL' in df.columns:
        from scipy.stats import norm
        df['Z'] = np.sign(np.log(df['OR'])) * norm.ppf(df['PVAL']/2).abs()
        log.info("Computed Z = sign(log(OR)) x |Z| from P")
    else:
        sys.exit("ERROR: Cannot compute Z: need BETA + SE, OR + SE or OR + P")

log.info(f"Z missing in {df['Z'].isna().sum():,} rows")

df.to_csv(args.out, sep='\t', index=False, float_format='%.6g',
          compression='gzip' if args.out.endswith(('.gz', '.bgz')) else None)
log.info(f"Wrote {len(df):,} rows to {args.out}")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
