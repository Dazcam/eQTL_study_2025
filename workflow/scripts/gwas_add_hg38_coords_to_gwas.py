#--------------------------------------------------------------------------------------
#
#    Replace hg19 coordinates in GWAS summary statistics with lifted hg38 coordinates
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09PREP-GWAS | Rule: add_hg38_coords_to_gwas (one job per GWAS)
#            Upstream:   liftover_gwas_to_hg38, add_N
#            Downstream: munge_sumstats; 12SMR, 14CTWAS (hg38 sumstats)
#
# Purpose:   - Convert the lifted bed file to 1-based hg38 CHR and BP per SNP
#            - Inner join with the sumstats on SNP (SNPs that failed liftOver are dropped)
#            - Replace CHR and BP, keep the required columns and any optional ones present
#
# Inputs:    lifted    liftOver output bed (hg38)
#            sumstats  sumstats with Z and N (hg19)
#
# Outputs:   sumstats with hg38 coordinates (tsv)
#
#--------------------------------------------------------------------------------------

import logging
import sys

import pandas as pd

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()
sys.excepthook = lambda *exc: log.error("Uncaught exception", exc_info=exc)

log.info("Adding hg38 coordinates to GWAS sumstats ...")

##  Variables  -------------------------------------------------------------------------
lifted_file = str(snakemake.input["lifted"])
sumstats_file = str(snakemake.input["sumstats"])
output_file = str(snakemake.output[0])

variables = {"lifted_file": lifted_file, "sumstats_file": sumstats_file,
             "output_file": output_file}
width = max(len(k) for k in variables)
log.info("Variables")
log.info("============================")
for name, value in variables.items():
    log.info(f"{name:<{width}}  {value}")
log.info("============================")

##  Read lifted coordinates  -----------------------------------------------------------
lifted = pd.read_csv(lifted_file, sep='\t', header=None, comment='#')
lifted.columns = ['CHR', 'START', 'END', 'SNP']
lifted['CHR'] = lifted['CHR'].str.replace('chr', '', regex=False)
lifted['BP'] = lifted['START'] + 1
lifted = lifted[['SNP', 'CHR', 'BP']]
log.info(f"Lifted SNPs: {lifted.shape[0]}")
log.info(f"Duplicated SNP IDs in lifted file: {lifted['SNP'].duplicated().sum()}")

##  Merge with sumstats  ---------------------------------------------------------------
sumstats = pd.read_csv(sumstats_file, sep='\t')
log.info(f"Sumstats SNPs: {sumstats.shape[0]}")
log.info(f"Duplicated SNP IDs in sumstats: {sumstats['SNP'].duplicated().sum()}")

merged = pd.merge(lifted, sumstats, on='SNP', how='inner')
log.info(f"Merged SNPs: {merged.shape[0]}")

merged['CHR'] = merged['CHR_x']
merged['BP'] = merged['BP_x']
merged = merged.drop(columns=['CHR_x', 'CHR_y', 'BP_x', 'BP_y'])

##  Select columns and write  ----------------------------------------------------------
required_columns = ['SNP', 'CHR', 'BP', 'PVAL', 'A1', 'A2', 'Z', 'N']
optional_columns = ['BETA', 'OR', 'SE', 'INFO', 'DIRECTION', 'NCAS', 'NCON', 'NCA', 'NCO']
output_columns = required_columns + [col for col in optional_columns if col in merged.columns]

missing_required = [col for col in required_columns if col not in merged.columns]
if missing_required:
    raise ValueError(f"Missing required columns in merged DataFrame: {missing_required}")

merged = merged[output_columns]
merged.to_csv(output_file, sep='\t', index=False)
log.info(f"Wrote {merged.shape[0]} rows to {output_file}")

log.info("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
