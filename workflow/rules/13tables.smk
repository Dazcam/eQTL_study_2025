configfile: '../config/config.yaml'

GWAS_TRAITS = list(config['gwas'].keys())
venn_gene_tbls = expand(config['tables']['smr_ctwas_venn_tbl']['genes_out'], gwas=GWAS_TRAITS)
smr_sig_tbls = expand(config['tables']['smr_sig_tbl']['out_file'], gwas = GWAS_TRAITS)

smr_raw_files = expand(
    "../results/10SMR/smr/{cell_type}/{cell_type}_{gwas}.smr",
    cell_type = config['cell_types'], gwas = GWAS_TRAITS
)

rule all:
    input:
#        config['tables']['eqtl_tbl']['out_file'],
#        config['tables']['smr_tbl']['out_file'],
#        config['tables']['ctwas_sig_tbl']['sig_out'],
#        config["tables"]["diff_exp_tbl"]["html_out"],
#        config["tables"]["diff_exp_tbl"]["xlsx_out"],
#        venn_gene_tbls
#        config['tables']['subcluster_specific_tbl']['out_file'],
#       config['tables']['smr_sig_tbl']['xlsx_out'],
       config['tables']['susie_diag_tbl']['out_file'], 
        

#rule eqtl_tbl:
#    output: config['tables']['eqtl_tbl']['out_file']
#    params: in_dir = config['tables']['eqtl_tbl']['in_dir'],
#            allele_file = config['tables']['eqtl_tbl']['allele_file'],
#            peak_dir = config['tables']['eqtl_tbl']['peak_dir'],
#            cell_types = config['cell_types'],
#            exp_pc_map = config['exp_pc_map']
#    singularity: config["containers"]["r_eqtl"]
#    resources: time="2:00:00"
#    log:  config['tables']['eqtl_tbl']['log']
#    script: "../scripts/manuscript_table_smr.R"

rule diff_exp_tbl:
    input:
        nb      = config["tables"]["diff_exp_tbl"]["nb"],
        pt_h5ad = expand(
            "../results/17PSEUDOTIME/{trajectory}/adata_pseudotime.h5ad",
            trajectory=config["trajectories"]
        )
    output:  html = config["tables"]["diff_exp_tbl"]["html_out"],
             xlsx = config["tables"]["diff_exp_tbl"]["xlsx_out"]
    conda:     config["scanpy"]["env"]
    resources: threads = 16, mem_mb = 360000, time = "2:00:00"
    params:    nb_out = config["tables"]["diff_exp_tbl"]["nb_out"]
    message: "Computing L1/L2/pseudotime bin differential expression tables"
    log:     config["tables"]["diff_exp_tbl"]["log"]
    shell:
        "papermill {input.nb} {params.nb_out} -p plate extra >> {log} 2>&1 && "
        "jupyter nbconvert --to html {params.nb_out} "
        "--output-dir=$(dirname {output.html}) "
        "--output=$(basename {output.html}) >> {log} 2>&1"

rule smr_sig_tbl:
    # ST 7
    input:  smr_raw_files 
    output: tsv  = smr_sig_tbls,
            xlsx = config['tables']['smr_sig_tbl']['xlsx_out']
    params: traits           = GWAS_TRAITS,
            cell_types       = config['cell_types'],
            in_dir           = config['tables']['smr_sig_tbl']['in_dir'],
            gene_lookup_file = config['dev_specificity']['report']['gene_lookup'],
            p_smr            = config['p_smr'],
            p_heidi          = config['p_heidi']
    singularity: config["containers"]["r_eqtl"]
    resources: time="0:45:00", threads=4, mem_mb=32000
    log:    config['tables']['smr_sig_tbl']['log']
    script: "../scripts/manuscript_table_smr_sig.R"

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

rule smr_ctwas_venn_tbl:
    input:  smr_sig_tbls = smr_sig_tbls    # <-- same
    output: venn_gene_tbls
    params: traits           = GWAS_TRAITS,
            smr_tbl_dir      = config['tables']['smr_ctwas_venn_tbl']['smr_tbl_dir'],
            ctwas_in_dir     = config['tables']['smr_ctwas_venn_tbl']['ctwas_in_dir'],
            weights_dir      = config['tables']['smr_ctwas_venn_tbl']['weights_dir'],
            cell_types       = config['cell_types'],
            gene_lookup_file = config['dev_specificity']['report']['gene_lookup'],
            ctwas_pip_thresh = 0.8,
            ctwas_use_bf     = True,
            ctwas_bf_thresh  = 0.05,
            use_symbols      = True
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00", threads=4, mem_mb=32000
    log:    config['tables']['smr_ctwas_venn_tbl']['log']
    script: "../scripts/manuscript_table_smr_ctwas_venn.R"

rule ctwas_sig_tbl:
    output: correction = config['tables']['ctwas_sig_tbl']['correction_out'],
            sig_out     = config['tables']['ctwas_sig_tbl']['sig_out']
    params: traits           = GWAS_TRAITS,
            cell_types       = config['cell_types'],
            in_dir           = config['tables']['ctwas_sig_tbl']['in_dir'],
            weights_dir      = config['tables']['ctwas_sig_tbl']['weights_dir'],
            gene_lookup_file = config['dev_specificity']['report']['gene_lookup'],
            pip_thresh       = 0.8,
            bf_p_thresh      = 0.05
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00", threads=4, mem_mb=32000
    log:    config['tables']['ctwas_sig_tbl']['log']
    script: "../scripts/manuscript_table_ctwas_sig.R"

rule subcluster_specific_tbl:
    # Supp Table 3
    output: config['tables']['subcluster_specific_tbl']['out_file']
    params: eqtl_effects = config['dev_specificity']['eqtl_effects']['out_file']
    singularity: config["containers"]["r_eqtl"]
    resources: time="0:20:00"
    log:    config['tables']['subcluster_specific_tbl']['log']
    script: "../scripts/manuscript_table_subcluster_specific.R"

rule susie_diag_tbl:
    # ST 9
    input:  
        susie_files = expand(config["susie"]["sort_susie"]["output"],
                             cell_type = config["cell_types"]),
        gene_meta_files = expand(config["susie"]["prep_susie_gene_meta"]["output"],
                                 cell_type = config["cell_types"])
    output: xlsx = config['tables']['susie_diag_tbl']['out_file']
    params: cell_types = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time="0:30:00", threads=2, mem_mb=16000
    message: "SuSiE fine-mapping diagnostics table (Supp Table 9)"
    log:    config['tables']['susie_diag_tbl']['log']
    script: "../scripts/manuscript_table_susie_diag.R"
