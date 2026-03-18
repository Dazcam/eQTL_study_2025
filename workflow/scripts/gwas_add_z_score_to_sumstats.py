#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--sumstats", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

# Read (auto-detects .gz)
df = pd.read_csv(args.sumstats, sep='\t', dtype=str, low_memory=False)
print(f"Read {len(df):,} rows from {args.sumstats}")

# Convert relevant columns to numeric
for c in ['BETA','SE','OR','PVAL','Z']:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')

# If Z is already good → do nothing
if 'Z' in df.columns and df['Z'].notna().any():
    print("Z column already present and valid → keeping it")
else:
    print("Computing signed Z...")
    if 'BETA' in df.columns and 'SE' in df.columns:
        df['Z'] = df['BETA'] / df['SE']
        print("   → Z = BETA / SE")
    elif 'OR' in df.columns and 'SE' in df.columns:
        df['Z'] = np.log(df['OR']) / df['SE']
        print("   → Z = log(OR) / SE")
    elif 'OR' in df.columns and 'PVAL' in df.columns:
        from scipy.stats import norm
        df['Z'] = np.sign(np.log(df['OR'])) * norm.ppf(df['PVAL']/2).abs()
        print("   → Z = sign(log(OR)) × |Z| from P")
    else:
        sys.exit("ERROR: Cannot compute Z — need BETA+SE or OR+SE or OR+P")

# Write output (preserves .gz if filename ends with .gz)
df.to_csv(args.out, sep='\t', index=False, float_format='%.6g',
          compression='gzip' if args.out.endswith(('.gz','.bgz')) else None)

print(f"Wrote {len(df):,} rows with Z column to {args.out}")
