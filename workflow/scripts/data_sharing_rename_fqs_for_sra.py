import json
import os
import sys

plate = snakemake.wildcards.plate
plate_dir = snakemake.params.plate_dir
manifest_path = snakemake.input.json
log_path = snakemake.log[0]

with open(manifest_path) as f:
    manifest = json.load(f)

plate_map = manifest.get(plate, {})
if not plate_map:
    print(f"ERROR: no entries for {plate} in manifest", file=sys.stderr)
    sys.exit(1)

renamed = 0
skipped = 0
errors = []

with open(log_path, 'a') as log:
    log.write(f"Starting rename for {plate} at {os.popen('date').read()}")
    
    for old_rel, new_name in plate_map.items():
        src = os.path.join(plate_dir, old_rel)
        dest = os.path.join(plate_dir, new_name)

        if not os.path.exists(src):
            msg = f"MISSING src: {src}"
            log.write(msg + "\n")
            errors.append(msg)
            continue

        if os.path.exists(dest):
            log.write(f"SKIP (already exists): {dest}\n")
            skipped += 1
            continue

        # Hard-link to preserve originals without doubling disk space
        os.link(src, dest)
        log.write(f"LINKED {old_rel} -> {new_name}\n")
        renamed += 1

    log.write(f"\nDone. renamed={renamed} skipped={skipped} errors={len(errors)}\n")

if errors:
    print(f"{len(errors)} files missing — see {log_path}", file=sys.stderr)
    sys.exit(1)
