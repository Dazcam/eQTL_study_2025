# Retired from 16tables.smk (Nov 2026): superseded by smr_sig_tbl (Supp Table 9);
# script moved to scripts/draft_scripts/, config block tables.smr_tbl removed.

rule smr_tbl:
    input:  smr_sig_tbls = smr_sig_tbls
    output: config['tables']['smr_tbl']['out_file']
    params: smr_dir      = config['tables']['smr_tbl']['smr_dir'],
            eqtl_nom_dir = config['tables']['smr_tbl']['eqtl_nom_dir'],
            cell_types   = config['cell_types'],
            disorders    = GWAS_TRAITS
    singularity: config["containers"]["r_eqtl"]
    resources: time="2:00:00",threads = 10, mem_mb = 80000
    log:  config['tables']['smr_tbl']['log']
    script: "../scripts/manuscript_table_smr.R"
