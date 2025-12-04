#!/usr/bin/env python3
import pandas as pd
import sys
import os

def main(parquets, output):
    print(f"Reading {len(parquets)} files...", file=sys.stderr)
    df = pd.concat([pd.read_parquet(f) for f in parquets], ignore_index=True)

    # Rename columns to standardize names
    df = df.rename(columns={"start_distance": "tss_distance", "af": "maf"})

    # Select and save
    cols = [
        "phenotype_id", "variant_id", "tss_distance", "ma_samples",
        "ma_count", "maf", "pval_nominal", "slope", "slope_se"
    ]
    df[cols].to_csv(output, sep="\t", index=False, compression="gzip")
    print(f"Wrote {output}", file=sys.stderr)

if __name__ == "__main__":
    # Arguments: script.py <output> <parquet1> <parquet2> ...
    output = sys.argv[1]
    parquets = sys.argv[2:]
    os.makedirs(os.path.dirname(output), exist_ok=True)
    main(parquets, output)

