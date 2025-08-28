import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

# Modify parser to accept multiple directories
parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dirs", nargs='+', help="Required. FULL paths to the fastq folders, space-separated for multiple directories.")
parser.add_argument("--plate", required=True, help="Plate identifier to append to sample names (e.g., plate1, plate2).")
args = parser.parse_args()

assert args.fastq_dirs is not None, "please provide at least one path to the fastq folder"

# Initialize dictionary
FILES = defaultdict(lambda: defaultdict(list))

# Iterate over multiple directories
for fastq_dir in args.fastq_dirs:
    for root, dirs, files in os.walk(fastq_dir):
        for file in files:
            if file.endswith("fastq.gz"):
                full_path = join(root, file)
                print(full_path)
                print(file)
                #R1 will be forward reads, R2 will be reverse reads
                #Filename structure: 10_S6_L001_R1_001.fastq.gz
                m = re.search(r'(\d+)_S\d+_(L\d{3})_(R\d)', file) 
                print(m)
                if m:
                    sample = m.group(1)
#                    if sample == '1':  # Exclude sublibrary 1 # We excluded sublibrary 1 from plate1
#                        continue
                    lane = m.group(2)
                    reads = m.group(3)
                    FILES[sample][reads].append(full_path)

# Sort R1 and R2 reads across lanes
FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
    for read in FILES[sample]:
        FILES_sorted[sample][read] = sorted(FILES[sample][read])

# Add plate identifier to the sample names
FILES_with_plate = defaultdict(lambda: defaultdict(list))

plate = args.plate

for sample in FILES_sorted.keys():
    sample_with_plate = f"{sample}_{plate}"
    for read in FILES_sorted[sample]:
        FILES_with_plate[sample_with_plate][read] = FILES_sorted[sample][read]

# Output summary information
print()
print("total {} unique samples will be processed".format(len(FILES_with_plate.keys())))
print("------------------------------------------")
for sample in FILES_with_plate.keys():
    for read in sorted(FILES_with_plate[sample]):
        print("{sample} {read} has {n} fastq".format(sample=sample, read=read, n=len(FILES_with_plate[sample][read])))
print("------------------------------------------")
print(f"check the samples_{plate}.json file for fastqs belonging to each sample")
print()

# Write the JSON output file with sample names including plate identifier
js = json.dumps(FILES_with_plate, indent=4, sort_keys=True)
with open(f'samples_{plate}.json', 'w') as json_file:
    json_file.write(js)
