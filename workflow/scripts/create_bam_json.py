import os
import json
import sys
from datetime import datetime

# Input directories
input_dirs = [
    "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1",
    "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2",
    "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate3"
]

# Output JSON file and log file
output_json = "bam_paths.json"
log_file = "bam_paths.log"

# Ensure paths end with /
input_dirs = [os.path.normpath(d) + "/" for d in input_dirs]

# Dictionary to store BAM paths
bam_paths = {}

# Open log file for appending
with open(log_file, "a") as log:
    # Log start time
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log.write(f"==== Starting BAM path generation at {start_time} ====\n")
    print(f"==== Starting BAM path generation at {start_time} ====")

    # Iterate through each input directory
    for input_dir in input_dirs:
        plate = os.path.basename(os.path.dirname(input_dir))  # e.g., "plate1"
        log.write(f"==== Processing directory: {input_dir} ====\n")
        print(f"==== Processing directory: {input_dir} ====")

        # Walk through the directory
        try:
            for root, dirs, files in os.walk(input_dir):
                # Check for barcode_headAligned_anno.bam in process/ subdirectory
                if "barcode_headAligned_anno.bam" in files and os.path.basename(root) == "process":
                    # Get the sublibrary directory name (e.g., 13_plate1, 2_plate2)
                    sublib_dir = os.path.basename(os.path.dirname(root))
                    # Skip directories like DGE_filtered
                    if sublib_dir == "DGE_filtered":
                        continue
                    # Construct full BAM path
                    bam_path = os.path.join(root, "barcode_headAligned_anno.bam")
                    # Normalize path to ensure consistency
                    bam_path = os.path.normpath(bam_path)
                    # Use sublibrary name as key
                    sublib_name = sublib_dir
                    # Check for duplicates
                    if sublib_name in bam_paths:
                        log.write(f"Warning: Duplicate sublibrary {sublib_name} found at {bam_path}\n")
                        print(f"Warning: Duplicate sublibrary {sublib_name} found at {bam_path}")
                    else:
                        bam_paths[sublib_name] = bam_path
                        log.write(f"Found BAM: {sublib_name} -> {bam_path}\n")
                        print(f"Found BAM: {sublib_name} -> {bam_path}")
        except Exception as e:
            log.write(f"Error processing {input_dir}: {str(e)}\n")
            print(f"Error processing {input_dir}: {str(e)}")

    # Write BAM paths to JSON
    try:
        with open(output_json, "w") as f:
            json.dump(bam_paths, f, indent=4)
        log.write(f"==== Generated JSON file: {output_json} ====\n")
        print(f"==== Generated JSON file: {output_json} ====")
    except Exception as e:
        log.write(f"Error writing JSON file: {str(e)}\n")
        print(f"Error writing JSON file: {str(e)}")
        sys.exit(1)

    # Log completion
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log.write(f"==== Completed BAM path generation at {end_time} ====\n\n")
    print(f"==== Completed BAM path generation at {end_time} ====\n")
