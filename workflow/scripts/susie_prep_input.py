#--------------------------------------------------------------------------------------
#
#    Prepare expression, covariate and sample files for SuSiE fine-mapping
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08SUSIE | Rule: prep_susie_input (one job per cell type)
#            Upstream:   05TENSORQTL (prep_input expression bed and split covariates)
#            Downstream: run_susie
#
# Purpose:   - Expression: drop the first three bed columns (chr, start, end) and
#              rename TargetID to phenotype_id
#            - Covariates: add a 'SampleID' label to the blank first header field
#            - Samples: write the sample IDs from the bed header (column 5 onward)
#
# Inputs:    exp_bed   tensorQTL expression bed
#            covar     tensorQTL split covariates file
#
# Outputs:   exp_bed   expression matrix (phenotype_id + one column per sample)
#            covar     covariates file with labelled first header field
#            samp_lst  one-column sample list (header: sample_id)
#
#--------------------------------------------------------------------------------------

import logging
import sys

## Set up logging for Snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Preparing SuSiE input files ...")

##  Variables  --------------------------------------------------------------------------
exp_bed_in   = snakemake.input.exp_bed
covar_in     = snakemake.input.covar
exp_bed_out  = snakemake.output.exp_bed
covar_out    = snakemake.output.covar
samp_lst_out = snakemake.output.samp_lst

variables = {
    "exp_bed_in": exp_bed_in,
    "covar_in": covar_in,
    "exp_bed_out": exp_bed_out,
    "covar_out": covar_out,
    "samp_lst_out": samp_lst_out,
}
width = max(len(k) for k in variables)
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  Expression matrix  ------------------------------------------------------------------
log.info("Writing expression matrix ...")
with open(exp_bed_in, 'r') as f:
    lines = f.readlines()
header = lines[0].strip().split('\t')
data = [line.strip().split('\t') for line in lines[1:]]
columns = header[3:]  # Columns from 4th onward (0-based index)
if 'TargetID' in columns:
    columns[columns.index('TargetID')] = 'phenotype_id'
with open(exp_bed_out, 'w') as f:
    f.write('\t'.join(columns) + '\n')
    for row in data:
        f.write('\t'.join(row[3:]) + '\n')
log.info(f"{len(data)} genes written.")

##  Covariates  -------------------------------------------------------------------------
log.info("Writing covariates file ...")
with open(covar_in, 'r') as f:
    lines = f.readlines()
if lines[0].startswith('\t'):
    lines[0] = 'SampleID\t' + lines[0].lstrip('\t')
with open(covar_out, 'w') as f:
    f.writelines(lines)
log.info(f"{len(lines) - 1} covariates written.")

##  Sample list  ------------------------------------------------------------------------
log.info("Writing sample list ...")
sample_ids = header[4:]  # Columns from 5th onward (0-based index)
with open(samp_lst_out, 'w') as f:
    f.write("sample_id\n")
    f.write('\n'.join(sample_ids) + '\n')
log.info(f"{len(sample_ids)} samples written.")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
