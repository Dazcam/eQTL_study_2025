#!/usr/bin/env python3
"""
generate_sra_manifest.py

Reads the per-plate sample JSON files and produces a rename manifest
(sra_rename_manifest.json) mapping each fastq's current relative path
(run_folder/original_filename) to its new SRA-compatible unique filename.

Rename scheme:
    {plate}_sublib{sublib_id}_run{run_index}_L{lane}_R{read}_001.fastq.gz

run_index is assigned by sorting run-folder names alphabetically (which is
also chronological given the YYMMDD prefix in Illumina run folder names).

Usage:
    python scripts/generate_sra_manifest.py \
        --json-dir ../config \
        --plates plate1 plate2 plate3 \
        --output ../config/sra_rename_manifest.json

Output JSON structure:
    {
      "plate1": {
        "241009_A00748_0603_BHVF5CDSXC_fastq/10_S9_L001_R1_001.fastq.gz":
            "plate1_sublib10_run2_L001_R1_001.fastq.gz",
        ...
      },
      "plate2": { ... },
      "plate3": { ... }
    }
"""

import argparse
import json
import os
import sys


def build_plate_manifest(plate: str, json_path: str) -> dict:
    """
    Build the rename mapping for one plate.

    Returns a dict of {old_rel_path: new_filename}.
    old_rel_path is relative to the plate's DATA_SHARING directory,
    e.g. '241009_A00748_0603_BHVF5CDSXC_fastq/10_S9_L001_R1_001.fastq.gz'
    """
    with open(json_path) as f:
        data = json.load(f)

    # Collect all run folders for this plate and assign chronological indices
    run_folders = sorted(
        set(path.split("/")[-2] for sample in data.values() for path in sample["R1"])
    )
    run_idx_map = {rf: i + 1 for i, rf in enumerate(run_folders)}

    mapping = {}
    for key, sample in data.items():
        sublib_id = key.split("_")[0]
        biosample = f"{plate}_sublib{sublib_id}"

        r1_paths = sample.get("R1", [])
        r2_paths = sample.get("R2", [])

        if len(r1_paths) != len(r2_paths):
            print(
                f"WARNING: R1/R2 count mismatch for {key} "
                f"(R1={len(r1_paths)}, R2={len(r2_paths)})",
                file=sys.stderr,
            )

        for r1, r2 in zip(r1_paths, r2_paths):
            run_folder = r1.split("/")[-2]
            lane = r1.split("_L")[1][:3]          # e.g. '001'
            run_idx = run_idx_map[run_folder]

            old_r1 = f"{run_folder}/{r1.split('/')[-1]}"
            old_r2 = f"{run_folder}/{r2.split('/')[-1]}"
            new_r1 = f"{biosample}_run{run_idx}_L{lane}_R1_001.fastq.gz"
            new_r2 = f"{biosample}_run{run_idx}_L{lane}_R2_001.fastq.gz"

            mapping[old_r1] = new_r1
            mapping[old_r2] = new_r2

    return mapping


def validate_manifest(manifest: dict) -> bool:
    """Check for duplicate new filenames within each plate (indicates a bug)."""
    ok = True
    for plate, mapping in manifest.items():
        new_names = list(mapping.values())
        seen = {}
        for new in new_names:
            seen[new] = seen.get(new, 0) + 1
        dupes = {n: c for n, c in seen.items() if c > 1}
        if dupes:
            print(f"ERROR: collisions in {plate}:", file=sys.stderr)
            for name, count in dupes.items():
                print(f"  {name}  (x{count})", file=sys.stderr)
            ok = False
        else:
            print(f"{plate}: {len(mapping)} renames — no collisions")
    return ok


def main():
    parser = argparse.ArgumentParser(
        description="Generate SRA rename manifest from plate JSON files."
    )
    parser.add_argument(
        "--json-dir",
        default="../config",
        help="Directory containing samples_{plate}.json files (default: ../config)",
    )
    parser.add_argument(
        "--plates",
        nargs="+",
        default=["plate1", "plate2", "plate3"],
        help="Plate names to process (default: plate1 plate2 plate3)",
    )
    parser.add_argument(
        "--output",
        default="../config/sra_rename_manifest.json",
        help="Output manifest path (default: ../config/sra_rename_manifest.json)",
    )
    args = parser.parse_args()

    manifest = {}
    for plate in args.plates:
        json_path = os.path.join(args.json_dir, f"samples_{plate}.json")
        if not os.path.exists(json_path):
            print(f"ERROR: {json_path} not found", file=sys.stderr)
            sys.exit(1)
        manifest[plate] = build_plate_manifest(plate, json_path)

    if not validate_manifest(manifest):
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nManifest written to {args.output}")

    # Print run-folder -> run_index legend for each plate
    print("\nRun index legend:")
    for plate in args.plates:
        json_path = os.path.join(args.json_dir, f"samples_{plate}.json")
        with open(json_path) as f:
            data = json.load(f)
        run_folders = sorted(
            set(p.split("/")[-2] for s in data.values() for p in s["R1"])
        )
        print(f"  {plate}:")
        for i, rf in enumerate(run_folders, 1):
            print(f"    run{i}: {rf}")


if __name__ == "__main__":
    main()
