# --------------------------------------------------------------------------------------
#
#    repl_bulk_isoform_eqtl.py
#
#    Compares snRNA-seq gene-level eQTL (all L1+L2+pseudotime cell types) against
#    bulk transcript-level eQTL (O'Brien et al. 2018) to assess whether identified
#    eGenes show evidence of isoform-specific regulation.
#
#    Strategy:
#      1. Parse GTF to build ENST -> ENSG mapping (Ensembl 113)
#      2. Load unique fetal eGenes from the shared extract_unique_egenes output
#         (repl_bulk_extract_unique_egenes.R) -- NOT re-derived from raw perm
#         files here (see MIGRATION NOTE below)
#      3. Pre-filter the ~12G bulk file using awk (subprocess) to extract only
#         rows matching our target SNP + transcript pairs -- avoids loading
#         the full file into memory
#      4. For each eGene, classify the bulk transcript pattern:
#           - n_transcripts_tested   : transcripts of this gene in bulk data
#           - n_transcripts_nominal  : transcripts with p < nominal_threshold
#           - isoform_specific       : only 1 (of >1) transcript significant (flag)
#           - direction_concordant   : snRNA-seq and bulk slope same sign
#           - pattern                : one of [all_sig, isoform_specific,
#                                       none_sig, not_in_bulk]
#      5. Add gene symbols and write the per-eGene summary table (single
#         tracked Snakemake output); per-transcript detail, the bulk-hits
#         extract, and the ENST->ENSG map are written as untracked side
#         files for caching/debugging
#
#    Usage: called via Snakemake rule isoform_eqtl_comparison
#
# --------------------------------------------------------------------------------------

import sys
import logging
import subprocess
import gzip
import io
import os
from pathlib import Path

import numpy as np
import pandas as pd

# ── Logging ───────────────────────────────────────────────────────────────────────────

log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# ── Parameters ────────────────────────────────────────────────────────────────────────

gtf_path          = str(snakemake.input.gtf)
bulk_eqtl_path    = str(snakemake.input.bulk_eqtl)
fetal_egenes_path = str(snakemake.input.fetal_egenes)
gene_lookup_path  = str(snakemake.input.gene_lookup)

out_summary  = str(snakemake.output[0])
out_dir      = Path(out_summary).parent
out_enst_map  = str(out_dir / "enst_ensg_map.tsv")
out_bulk_hits = str(out_dir / "bulk_transcript_hits.tsv")
out_detail    = str(out_dir / "isoform_comparison_detail.tsv")

nominal_threshold = float(snakemake.params.nominal_threshold)

logger.info(f"Fetal eGenes file: {fetal_egenes_path}")
logger.info(f"Bulk nominal threshold: {nominal_threshold}")

# ── 1. Build ENST -> ENSG mapping from GTF ───────────────────────────────────────────
#
#  Parse transcript-level records only. Each transcript line contains both
#  gene_id (ENSG) and transcript_id (ENST) attributes.

logger.info(f"Parsing GTF for ENST -> ENSG mapping: {gtf_path}")

enst_ensg = {}

opener = gzip.open if gtf_path.endswith('.gz') else open

with opener(gtf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9 or fields[2] != 'transcript':
            continue
        attrs = fields[8]

        # Extract gene_id and transcript_id from attribute string
        gene_id, transcript_id = None, None
        for part in attrs.split(';'):
            part = part.strip()
            if part.startswith('gene_id'):
                gene_id = part.split('"')[1] if '"' in part else part.split()[1]
            elif part.startswith('transcript_id'):
                transcript_id = part.split('"')[1] if '"' in part else part.split()[1]

        if gene_id and transcript_id:
            # Strip version suffixes (e.g. ENSG00000001.5 -> ENSG00000001)
            enst_clean = transcript_id.split('.')[0]
            ensg_clean = gene_id.split('.')[0]
            enst_ensg[enst_clean] = ensg_clean

logger.info(f"ENST -> ENSG mapping built: {len(enst_ensg):,} transcripts")

enst_ensg_df = pd.DataFrame(
    list(enst_ensg.items()),
    columns=['transcript_id', 'gene_id']
)
enst_ensg_df.to_csv(out_enst_map, sep='\t', index=False)
logger.info(f"ENST -> ENSG map written: {out_enst_map}")

# ── 2. Load unique fetal eGenes from the shared extract_unique_egenes output ─────────

logger.info("Loading fetal eGenes from shared extract_unique_egenes output...")

all_egenes = pd.read_csv(fetal_egenes_path, sep='\t')
all_egenes['gene_id'] = all_egenes['phenotype_id'].str.split('.').str[0]

logger.info(f"Rows loaded (eGene x cell_type pairs): {len(all_egenes):,}")
logger.info(f"Unique eGenes: {all_egenes['gene_id'].nunique():,}")
logger.info(f"Unique SNPs:   {all_egenes['variant_id'].nunique():,}")
if 'level' in all_egenes.columns:
    logger.info(f"Cell types by level: "
                f"{all_egenes[['cell_type','level']].drop_duplicates()['level'].value_counts().to_dict()}")

# ── 3. Build lookup sets for bulk file filtering ──────────────────────────────────────
#
#  Strategy: for each eGene we want to find all transcripts of that gene in bulk,
#  tested against the same top SNP from snRNA-seq.
#
#  We need:
#    - The set of all target SNPs (variant_ids from our eGenes)
#    - The ENSG -> ENST mapping (reverse of the GTF mapping above)
#    - Then: for each bulk row, keep if variant_id in target_snps
#            AND gene_id (mapped ENST->ENSG) in target_genes
#
#  This is done via awk for speed on the ~12G file.

target_snps  = set(all_egenes['variant_id'].unique())
target_genes = set(all_egenes['gene_id'].unique())

# Build reverse map: ENSG -> set of ENSTs
ensg_to_enst = {}
for enst, ensg in enst_ensg.items():
    if ensg not in ensg_to_enst:
        ensg_to_enst[ensg] = set()
    ensg_to_enst[ensg].add(enst)

# Get all ENSTs belonging to our target eGenes
target_ensts = set()
for ensg in target_genes:
    target_ensts.update(ensg_to_enst.get(ensg, set()))

logger.info(f"Target SNPs:        {len(target_snps):,}")
logger.info(f"Target eGenes:      {len(target_genes):,}")
logger.info(f"Target transcripts: {len(target_ensts):,}")

# ── 4. Pre-filter bulk file using awk ─────────────────────────────────────────────────
#
#  Write lookup files for awk, then use FNR==NR pattern to build in-memory
#  hash tables of target SNPs and transcripts, filtering the bulk file in one pass.
#  This avoids loading ~12G into Python memory.

tmp_dir = out_dir / "tmp"
tmp_dir.mkdir(parents=True, exist_ok=True)

snp_list_path  = tmp_dir / "target_snps.txt"
enst_list_path = tmp_dir / "target_ensts.txt"

logger.info(f"Writing target SNP list ({len(target_snps):,} SNPs)...")
with open(snp_list_path, 'w') as f:
    f.write('\n'.join(target_snps) + '\n')

logger.info(f"Writing target ENST list ({len(target_ensts):,} transcripts)...")
with open(enst_list_path, 'w') as f:
    f.write('\n'.join(target_ensts) + '\n')

# awk command:
#   File 1 (snps.txt):   build snp hash
#   File 2 (ensts.txt):  build enst hash
#   File 3 (bulk file):  print header + matching rows
#   Column indices (1-based): gene_id=1, variant_id=2

awk_cmd = (
    f"awk 'BEGIN{{OFS=\"\\t\"}} "
    f"FILENAME==\"{snp_list_path}\" {{snps[$1]=1; next}} "
    f"FILENAME==\"{enst_list_path}\" {{ensts[$1]=1; next}} "
    f"FNR==1 {{print; next}} "
    f"($2 in snps) && ($1 in ensts) {{print}}' "
    f"{snp_list_path} {enst_list_path} "
    f"<(zcat {bulk_eqtl_path})"
)

logger.info("Pre-filtering bulk file with awk (this may take a few minutes)...")
logger.info(f"awk command: {awk_cmd}")

result = subprocess.run(
    awk_cmd,
    shell=True,
    executable='/bin/bash',       # needed for process substitution <(zcat ...)
    capture_output=True,
    text=True
)

if result.returncode != 0:
    logger.error(f"awk filtering failed:\n{result.stderr}")
    raise RuntimeError("awk bulk file filtering failed")

logger.info(f"awk filtering complete. Output size: {len(result.stdout):,} characters")

bulk_filtered = pd.read_csv(
    io.StringIO(result.stdout),
    sep='\t',
    dtype={'gene_id': str, 'variant_id': str}
)

# Strip version suffixes from ENST IDs if present
bulk_filtered['transcript_id'] = bulk_filtered['gene_id'].str.split('.').str[0]

# Map transcript -> gene
bulk_filtered['gene_id_ensg'] = bulk_filtered['transcript_id'].map(enst_ensg)

n_missing_map = bulk_filtered['gene_id_ensg'].isna().sum()
if n_missing_map > 0:
    logger.warning(f"{n_missing_map:,} bulk rows could not be mapped to ENSG "
                   f"(transcript not in GTF -- may be from older Ensembl release)")

bulk_filtered = bulk_filtered.dropna(subset=['gene_id_ensg'])

logger.info(f"Bulk rows after filtering and mapping: {len(bulk_filtered):,}")
logger.info(f"Unique transcripts in bulk hits:  {bulk_filtered['transcript_id'].nunique():,}")
logger.info(f"Unique genes in bulk hits:        {bulk_filtered['gene_id_ensg'].nunique():,}")
logger.info(f"Unique SNPs in bulk hits:         {bulk_filtered['variant_id'].nunique():,}")

bulk_filtered.to_csv(out_bulk_hits, sep='\t', index=False)
logger.info(f"Bulk hits written: {out_bulk_hits}")

# ── 5. Classify isoform patterns per eGene ────────────────────────────────────────────
#
#  For each unique ENSG + variant_id pair across all cell types, characterise
#  the bulk transcript pattern.
#
#  Classifications:
#    not_in_bulk       : no transcripts of this gene tested for this SNP in bulk
#    none_sig          : transcripts present in bulk but none significant
#    all_sig           : all tested transcripts significant (gene-wide effect)
#    isoform_specific  : only a subset (>=1 but not all) transcripts significant
#
#  isoform_specific is the key flag for the reviewer response

logger.info("Classifying isoform patterns for each eGene x SNP pair...")

# Unique eGene x SNP pairs (collapse across cell types for the bulk lookup)
egene_snp_pairs = (all_egenes[['gene_id', 'variant_id', 'slope']]
                   .drop_duplicates(subset=['gene_id', 'variant_id']))

detail_rows  = []
summary_rows = []

for _, row in egene_snp_pairs.iterrows():
    ensg     = row['gene_id']
    snp      = row['variant_id']
    slope_sn = row['slope']    # snRNA-seq slope for direction comparison

    bulk_sub = bulk_filtered[
        (bulk_filtered['gene_id_ensg'] == ensg) &
        (bulk_filtered['variant_id']   == snp)
    ].copy()

    n_transcripts = len(bulk_sub)

    if n_transcripts == 0:
        pattern = 'not_in_bulk'
        n_sig   = 0
        isoform_specific   = False
        direction_majority = np.nan
        most_sig_transcript = np.nan
        most_sig_pval       = np.nan
        most_sig_slope      = np.nan
    else:
        bulk_sub['significant'] = bulk_sub['pval_nominal'] < nominal_threshold

        n_sig = bulk_sub['significant'].sum()

        if n_sig == 0:
            pattern = 'none_sig'
            isoform_specific = False
        elif n_sig == n_transcripts:
            pattern = 'all_sig'
            isoform_specific = False
        else:
            pattern = 'isoform_specific'
            isoform_specific = True

        sig_slopes = bulk_sub.loc[bulk_sub['significant'], 'slope']
        if len(sig_slopes) > 0:
            direction_majority = (
                (sig_slopes > 0).sum() > (sig_slopes < 0).sum()
            )
            direction_concordant = (
                (slope_sn > 0) == direction_majority
            )
        else:
            direction_concordant = np.nan

        idx_min = bulk_sub['pval_nominal'].idxmin()
        most_sig_transcript = bulk_sub.loc[idx_min, 'transcript_id']
        most_sig_pval       = bulk_sub.loc[idx_min, 'pval_nominal']
        most_sig_slope      = bulk_sub.loc[idx_min, 'slope']

        for _, br in bulk_sub.iterrows():
            detail_rows.append({
                'gene_id'            : ensg,
                'variant_id'         : snp,
                'transcript_id'      : br['transcript_id'],
                'bulk_pval_nominal'  : br['pval_nominal'],
                'bulk_slope'         : br['slope'],
                'bulk_slope_se'      : br['slope_se'],
                'bulk_maf'           : br['maf'],
                'bulk_significant'   : br['significant'],
                'pattern'            : pattern,
                'isoform_specific'   : isoform_specific,
                'snrna_slope'        : slope_sn,
            })

    ct_list = '|'.join(sorted(
        all_egenes.loc[
            (all_egenes['gene_id'] == ensg) &
            (all_egenes['variant_id'] == snp),
            'cell_type'
        ].unique()
    ))

    sn_row = all_egenes[
        (all_egenes['gene_id'] == ensg) &
        (all_egenes['variant_id'] == snp)
    ].iloc[0]

    summary_rows.append({
        'gene_id'                    : ensg,
        'variant_id'                 : snp,
        'cell_types'                 : ct_list,
        'n_cell_types'               : len(ct_list.split('|')),
        'snrna_pval_nominal'         : sn_row['pval_nominal'],
        'snrna_slope'                : sn_row['slope'],
        'snrna_slope_se'             : sn_row['slope_se'],
        'n_bulk_transcripts_tested'  : n_transcripts,
        'n_bulk_transcripts_nominal' : n_sig,
        'pattern'                    : pattern,
        'isoform_specific'           : isoform_specific,
        'direction_concordant'       : direction_concordant if n_transcripts > 0 else np.nan,
        'most_sig_transcript'        : most_sig_transcript if n_transcripts > 0 else np.nan,
        'most_sig_bulk_pval'         : most_sig_pval if n_transcripts > 0 else np.nan,
        'most_sig_bulk_slope'        : most_sig_slope if n_transcripts > 0 else np.nan,
    })

summary_df = pd.DataFrame(summary_rows)
detail_df  = pd.DataFrame(detail_rows)

# ── 6. Summary statistics ─────────────────────────────────────────────────────────────

logger.info("\n── Pattern classification summary ──────────────────────────")
pattern_counts = summary_df['pattern'].value_counts()
logger.info(f"\n{pattern_counts.to_string()}")
logger.info(f"\nTotal eGene x SNP pairs:      {len(summary_df):,}")
logger.info(f"Isoform-specific candidates:  "
            f"{summary_df['isoform_specific'].sum():,} "
            f"({100*summary_df['isoform_specific'].mean():.1f}%)")
logger.info(f"Not found in bulk:            "
            f"{(summary_df['pattern']=='not_in_bulk').sum():,} "
            f"({100*(summary_df['pattern']=='not_in_bulk').mean():.1f}%)")
logger.info(f"Direction concordant (of those in bulk): "
            f"{summary_df['direction_concordant'].sum():.0f} / "
            f"{summary_df['direction_concordant'].notna().sum():.0f}")

logger.info("\n── Per cell type eGene counts and isoform-specific rates ────")
ct_stats = (all_egenes[['cell_type','gene_id','variant_id']]
            .merge(summary_df[['gene_id','variant_id','pattern','isoform_specific']],
                   on=['gene_id','variant_id'], how='left')
            .groupby('cell_type')
            .agg(
                n_egenes        = ('gene_id', 'count'),
                n_isoform_spec  = ('isoform_specific', 'sum'),
                n_not_in_bulk   = ('pattern', lambda x: (x=='not_in_bulk').sum()),
                n_none_sig      = ('pattern', lambda x: (x=='none_sig').sum()),
                n_all_sig       = ('pattern', lambda x: (x=='all_sig').sum()),
            )
            .reset_index())
ct_stats['pct_isoform_specific'] = (
    ct_stats['n_isoform_spec'] / ct_stats['n_egenes'] * 100
).round(1)
logger.info(f"\n{ct_stats.to_string(index=False)}")

# ── 7. Add gene symbols and write outputs ─────────────────────────────────────────────

logger.info(f"\nLoading gene symbol lookup: {gene_lookup_path}")
gene_lookup = pd.read_csv(gene_lookup_path, sep='\t', usecols=['ensembl_gene_id', 'external_gene_name'])
gene_lookup = gene_lookup.drop_duplicates(subset='ensembl_gene_id')

summary_df = summary_df.merge(
    gene_lookup.rename(columns={'ensembl_gene_id': 'gene_id', 'external_gene_name': 'gene_symbol'}),
    on='gene_id', how='left'
)
n_missing_symbol = summary_df['gene_symbol'].isna().sum()
if n_missing_symbol > 0:
    logger.warning(f"{n_missing_symbol:,} eGenes had no matching gene symbol in the lookup")

# Move gene_symbol next to gene_id for readability
cols = summary_df.columns.tolist()
cols.insert(cols.index('gene_id') + 1, cols.pop(cols.index('gene_symbol')))
summary_df = summary_df[cols]

summary_df.to_csv(out_summary, sep='\t', index=False)
logger.info(f"\nSummary (tracked output) written: {out_summary} ({summary_df.shape})")

if detail_df.empty:
    logger.warning("No transcript-level detail rows generated -- "
                   "check that bulk SNPs match snRNA-seq variant IDs")
    detail_df = pd.DataFrame(columns=[
        'gene_id','variant_id','transcript_id',
        'bulk_pval_nominal','bulk_slope','bulk_slope_se','bulk_maf',
        'bulk_significant','pattern','isoform_specific','snrna_slope'
    ])
detail_df.to_csv(out_detail, sep='\t', index=False)
logger.info(f"Detail (untracked) written:  {out_detail} ({detail_df.shape})")

for f in [snp_list_path, enst_list_path]:
    try:
        os.remove(f)
    except Exception:
        pass

logger.info("Done.")
