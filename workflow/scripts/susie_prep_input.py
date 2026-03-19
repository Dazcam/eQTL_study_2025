import pandas as pd

# Access input and output files from Snakemake
input_exp_bed = snakemake.input.exp_bed
input_covar = snakemake.input.covar
output_exp_bed = snakemake.output.exp_bed
output_covar = snakemake.output.covar
output_samp_lst = snakemake.output.samp_lst

# Process exp_bed: Extract columns from 4th onward and rename 'TargetID' to 'phenotype_id'
with open(input_exp_bed, 'r') as f:
    lines = f.readlines()
header = lines[0].strip().split('\t')
data = [line.strip().split('\t') for line in lines[1:]]
columns = header[3:]  # Columns from 4th onward (0-based index)
if 'TargetID' in columns:
    columns[columns.index('TargetID')] = 'phenotype_id'
with open(output_exp_bed, 'w') as f:
    f.write('\t'.join(columns) + '\n')
    for row in data:
        f.write('\t'.join(row[3:]) + '\n')

# Process covar: Add 'SampleID' to the first line if it starts with a tab
with open(input_covar, 'r') as f:
    lines = f.readlines()
if lines[0].startswith('\t'):
    lines[0] = 'SampleID\t' + lines[0].lstrip('\t')
with open(output_covar, 'w') as f:
    f.writelines(lines)

# Process samp_lst: Extract sample IDs from exp_bed header (columns 5 onward)
sample_ids = header[4:]  # Columns from 5th onward (0-based index)
with open(output_samp_lst, 'w') as f:
    f.write("sample_id\n")
    f.write('\n'.join(sample_ids) + '\n')
