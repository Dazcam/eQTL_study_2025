configfile: '../config/config.yaml'

import re

GWAS_TRAITS = list(config['gwas'].keys())
TBL = config['tables']

smr_sig_tbls = expand(TBL['smr_sig_tbl']['out_file'], gwas = GWAS_TRAITS)
venn_gene_tbls = expand(TBL['smr_ctwas_venn_tbl']['genes_out'], gwas = GWAS_TRAITS)

# L1 cell types with subclusters, in config cell_types order
L2_POPS = [ct for ct in config['cell_types']
           if '-to-' not in ct and any(re.fullmatch(re.escape(ct) + r'-\d+', sub) for sub in config['cell_types'])]


rule diff_exp_tbl:
    # Supp Table 1
    input:
        nb      = TBL["diff_exp_tbl"]["nb"],
        pt_h5ad = expand(config["pseudotime"]["palantir"]["h5ad"],
                         trajectory = config["trajectories"])
    output:  html = TBL["diff_exp_tbl"]["html_out"],
             xlsx = TBL["diff_exp_tbl"]["xlsx_out"]
    params:  nb_out         = TBL["diff_exp_tbl"]["nb_out"],
             pseudotime_dir = TBL["diff_exp_tbl"]["pseudotime_dir"],
             trajectories   = ",".join(config["trajectories"]),
             cell_types     = ",".join(config["cell_types"]),
             labels         = ",".join(f"{k}:{v}" for k, v in config["manuscript_labels"].items())
    conda:     config["scanpy"]["env"]
    resources: threads = 16, mem_mb = 160000, time = "2:00:00"
    message:   "Computing L1/L2/pseudotime bin differential expression tables"
    benchmark: TBL["diff_exp_tbl"]["benchmark"]
    log:       TBL["diff_exp_tbl"]["log"]
    shell:
        "papermill {input.nb} {params.nb_out} -p plate extra "
        "-p xlsx_out {output.xlsx} -p pseudotime_dir {params.pseudotime_dir} "
        "-p trajectories {params.trajectories} -p cell_types {params.cell_types} "
        "-p labels {params.labels} > {log} 2>&1 && "
        "jupyter nbconvert --to html {params.nb_out} "
        "--output-dir=$(dirname {output.html}) "
        "--output=$(basename {output.html}) >> {log} 2>&1"

rule eqtl_tbl:
    # Supp Table 2
    input:  allele_file = config['tensorQTL']['convert_genotypes']['pvar'],
            peak_file   = config['replication_sc']['dwnld_ziffra']['output']
    output: TBL['eqtl_tbl']['out_file']
    params: in_dir           = TBL['eqtl_tbl']['in_dir'],
            labels           = config['manuscript_labels'],
            cell_types       = config['cell_types'],
            exp_pc_map       = config['exp_pc_map'],
            geno_pc          = config['tensorQTL']['geno_pcs'],
            norm_method      = config['tensorQTL']['norm_methods'][0],
            gene_lookup_file = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00", mem_mb = 8000
    benchmark: TBL['eqtl_tbl']['benchmark']
    log:       TBL['eqtl_tbl']['log']
    script:    TBL['eqtl_tbl']['script']

rule subcluster_specific_tbl:
    # Supp Table 3
    output: TBL['subcluster_specific_tbl']['out_file']
    params: eqtl_effects = config['eqtl_sharing']['eqtl_effects']['out_file'],
            labels       = config['manuscript_labels']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['subcluster_specific_tbl']['benchmark']
    log:       TBL['subcluster_specific_tbl']['log']
    script:    TBL['subcluster_specific_tbl']['script']

rule eqtl_specificity_tbl:
    # Supp Table 4
    input:  classification = config['eqtl_sharing']['classify_specificity']['out_file'],
            gene_lookup    = config['gene_lookup']
    output: TBL['eqtl_specificity_tbl']['out_file']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['eqtl_specificity_tbl']['benchmark']
    log:       TBL['eqtl_specificity_tbl']['log']
    script:    TBL['eqtl_specificity_tbl']['script']

rule eqtl_heterogeneity_tbl:
    # Supp Table 5
    input:  config['eqtl_sharing']['heterogeneity_test']['out_file']
    output: TBL['eqtl_heterogeneity_tbl']['out_file']
    params: labels = config['manuscript_labels']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['eqtl_heterogeneity_tbl']['benchmark']
    log:       TBL['eqtl_heterogeneity_tbl']['log']
    script:    TBL['eqtl_heterogeneity_tbl']['script']

# Supp Table 6: made outside the pipeline

rule isoform_tbl:
    # Supp Table 7
    input:  config['replication_bulk']['isoform_eqtl']['summary']
    output: TBL['isoform_tbl']['out_file']
    params: cell_types        = config['cell_types'],
            labels            = config['manuscript_labels'],
            nominal_threshold = config['replication_bulk']['isoform_eqtl']['nominal_threshold'],
            bulk_p_thresh     = TBL['isoform_tbl']['bulk_p_thresh']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['isoform_tbl']['benchmark']
    log:       TBL['isoform_tbl']['log']
    script:    TBL['isoform_tbl']['script']

rule ldsr_tbl:
    # Supp Table 8
    input:  config['sldsr']['ldsr_strat_summary']['output'].format(annot_type = TBL['ldsr_tbl']['annot_type'])
    output: TBL['ldsr_tbl']['out_file']
    params: cell_types = config['cell_types'],
            labels     = config['manuscript_labels']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['ldsr_tbl']['benchmark']
    log:       TBL['ldsr_tbl']['log']
    script:    TBL['ldsr_tbl']['script']

rule smr_sig_tbl:
    # Supp Table 9
    input:  expand(config['smr']['smr']['smr'], cell_type = config['cell_types'], gwas = GWAS_TRAITS)
    output: tsv  = smr_sig_tbls,
            xlsx = TBL['smr_sig_tbl']['xlsx_out']
    params: traits           = GWAS_TRAITS,
            labels           = config['manuscript_labels'],
            cell_types       = config['cell_types'],
            in_dir           = TBL['smr_sig_tbl']['in_dir'],
            gene_lookup_file = config['gene_lookup'],
            p_smr            = config['p_smr'],
            p_heidi          = config['p_heidi']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['smr_sig_tbl']['benchmark']
    log:       TBL['smr_sig_tbl']['log']
    script:    TBL['smr_sig_tbl']['script']

rule ctwas_sig_tbl:
    # Supp Table 10
    output: correction = TBL['ctwas_sig_tbl']['correction_out'],
            sig_out    = TBL['ctwas_sig_tbl']['sig_out']
    params: traits           = GWAS_TRAITS,
            labels           = config['manuscript_labels'],
            cell_types       = config['cell_types'],
            in_dir           = TBL['ctwas_sig_tbl']['in_dir'],
            weights_dir      = TBL['ctwas_sig_tbl']['weights_dir'],
            gene_lookup_file = config['gene_lookup'],
            pip_thresh       = config['ctwas']['sig_pip_thresh'],
            bf_p_thresh      = config['ctwas']['sig_bf_thresh']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "2:00:00", mem_mb = 32000
    benchmark: TBL['ctwas_sig_tbl']['benchmark']
    log:       TBL['ctwas_sig_tbl']['log']
    script:    TBL['ctwas_sig_tbl']['script']

rule susie_diag_tbl:
    # Supp Table 11
    input:
        susie_files     = expand(config["susie"]["sort_susie"]["output"],
                                 cell_type = config["cell_types"]),
        gene_meta_files = expand(config["susie"]["prep_susie_gene_meta"]["output"],
                                 cell_type = config["cell_types"])
    output: xlsx = TBL['susie_diag_tbl']['out_file']
    params: cell_types = config['cell_types'],
            labels     = config['manuscript_labels']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: TBL['susie_diag_tbl']['benchmark']
    log:       TBL['susie_diag_tbl']['log']
    script:    TBL['susie_diag_tbl']['script']

# Supp Table 12: made outside the pipeline

rule smr_ctwas_venn_tbl:
    # Gene lists for Fig 6
    input:  smr_sig_tbls = smr_sig_tbls
    output: venn_gene_tbls
    params: traits           = GWAS_TRAITS,
            smr_tbl_dir      = TBL['smr_ctwas_venn_tbl']['smr_tbl_dir'],
            ctwas_in_dir     = TBL['smr_ctwas_venn_tbl']['ctwas_in_dir'],
            weights_dir      = TBL['smr_ctwas_venn_tbl']['weights_dir'],
            cell_types       = config['cell_types'],
            gene_lookup_file = config['gene_lookup'],
            ctwas_pip_thresh = config['ctwas']['sig_pip_thresh'],
            ctwas_use_bf     = True,
            ctwas_bf_thresh  = config['ctwas']['sig_bf_thresh'],
            use_symbols      = True
    singularity: config["containers"]["r_eqtl"]
    resources: time = "4:00:00", mem_mb = 16000
    benchmark: TBL['smr_ctwas_venn_tbl']['benchmark']
    log:       TBL['smr_ctwas_venn_tbl']['log']
    script:    TBL['smr_ctwas_venn_tbl']['script']

rule cell_counts_tbl:
    # Nuclei counts for Fig 2B
    input:  l1 = config['scanpy']['clusters_h5ad'],
            l2 = expand(config['scanpy']['subclusters_h5ad'], l1_pop = L2_POPS),
            pt = expand(config['pseudotime']['palantir']['h5ad'], trajectory = config['trajectories'])
    output: TBL['cell_counts_tbl']['out_file']
    params: cell_types = config['cell_types'],
            labels     = config['manuscript_labels']
    conda:  config["scanpy"]["env"]
    resources: time = "0:30:00", mem_mb = 8000
    benchmark: TBL['cell_counts_tbl']['benchmark']
    log:       TBL['cell_counts_tbl']['log']
    script:    TBL['cell_counts_tbl']['script']


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["16tables"] = [rules.diff_exp_tbl.output,
                       rules.eqtl_tbl.output,
                       rules.subcluster_specific_tbl.output,
                       rules.eqtl_specificity_tbl.output,
                       rules.eqtl_heterogeneity_tbl.output,
                       rules.isoform_tbl.output,
                       rules.ldsr_tbl.output,
                       rules.smr_sig_tbl.output,
                       rules.ctwas_sig_tbl.output,
                       rules.susie_diag_tbl.output,
                       rules.smr_ctwas_venn_tbl.output,
                       rules.cell_counts_tbl.output]
