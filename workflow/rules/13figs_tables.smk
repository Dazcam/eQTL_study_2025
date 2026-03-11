configfile: '../config/config.yaml'

rule all:
    input:
#        config['plotting']['eqtl_tbl']['out_file'],
#        config['plotting']['smr_tbl']['out_file'],
#        config['plotting']['eqtl_boxplots_py']['output'],
#        config['plotting']['eqtl_qc_plt']['out_file'],
#        config['plotting']['ldsr_plt']['out_file'],
#        config['plotting']['rep_plt']['out_file'],
#        config['plotting']['supp_plt']['out_file'],
        config['plotting']['data_perm']['out_file'],
        config['plotting']['data_nominal']['out_file'],
        config['plotting']['data_mk_eqtl_tar']['out_file']
#        config['plotting']['data_weights']['out_file']

rule eqtl_tbl:
    output: config['plotting']['eqtl_tbl']['out_file']
    params: in_dir = config['plotting']['eqtl_tbl']['in_dir'],
            allele_file = config['plotting']['eqtl_tbl']['allele_file'],
            peak_dir = config['plotting']['eqtl_tbl']['peak_dir']    
    singularity: config["containers"]["r_eqtl"]
    resources: time="2:00:00"
    log:  config['plotting']['eqtl_tbl']['log']
    script: "../scripts/manuscript_eQTL_table.R"

rule smr_tbl:
    output: config['plotting']['smr_tbl']['out_file']
    params: smr_dir = config['plotting']['smr_tbl']['smr_dir'],
            eqtl_nom_dir = config['plotting']['smr_tbl']['eqtl_nom_dir']
    singularity: config["containers"]["r_eqtl"]
    resources: time="2:00:00",threads = 10, mem_mb = 80000
    log:  config['plotting']['smr_tbl']['log']
    script: "../scripts/manuscript_smr_table.R"


rule eqtl_qc_plt:
    output: config['plotting']['eqtl_qc_plt']['out_file']
    params: in_dir = config['plotting']['eqtl_qc_plt']['in_dir'],
            ziffra_dir = config['plotting']['eqtl_qc_plt']['ziffra_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['eqtl_qc_plt']['log']
    script: "../scripts/manuscript_plot_eqtl_QC.R"

rule ldsr_plt:
    output: config['plotting']['ldsr_plt']['out_file']
    params: in_dir = config['plotting']['ldsr_plt']['in_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['ldsr_plt']['log']
    script: "../scripts/manuscript_plot_ldsr.R"

rule replication_plt:
    output: config['plotting']['rep_plt']['out_file']
    params: in_dir = config['plotting']['rep_plt']['in_dir'],
            internal_dir = config['plotting']['rep_plt']['internal_dir'],
            fugita_dir = config['plotting']['rep_plt']['fugita_dir'],
            beta_dir = config['plotting']['rep_plt']['beta_dir']
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['rep_plt']['log']
    script: "../scripts/manuscript_plot_replication.R"

rule supplementary_plt:
    output: config['plotting']['supp_plt']['out_file']
    params: geno_dir = config['plotting']['supp_plt']['geno_dir'],
            expr_dir = config['plotting']['supp_plt']['expr_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    resources: threads = 4, mem_mb = 20000
    log:  config['plotting']['supp_plt']['log']
    script: "../scripts/manuscript_plot_supplementary.R"

rule data_perm:
    output: config['plotting']['data_perm']['out_file']
    params: in_dir = config['plotting']['data_perm']['in_dir'],
            allele_file = config['plotting']['data_perm']['allele_file'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['data_perm']['log']
    script: "../scripts/data_sharing_eQTL_perm.R"

rule data_nominal:
    output: config['plotting']['data_nominal']['out_file']
    params: in_dir = config['plotting']['data_nominal']['in_dir'],
            allele_file = config['plotting']['data_nominal']['allele_file'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="4:00:00", threads = 10, mem_mb = 80000
    log:  config['plotting']['data_nominal']['log']
    script: "../scripts/data_sharing_eQTL_nominal.R"

rule mk_eqtl_tar:
    input:  perm = config['plotting']['data_perm']['out_file'],
            nom = config['plotting']['data_nominal']['out_file']
    output: config['plotting']['data_mk_eqtl_tar']['out_file']
    params: out_dir = config['plotting']['data_mk_eqtl_tar']['out_dir'],
            cell_types = config['cell_types']
    envmodules: "pigz"
    resources: time="5:00:00", threads = 8, mem_mb = 16000
    log:  config['plotting']['data_mk_eqtl_tar']['log']
    shell: """
           scripts/data_sharing_mk_eqtl_tar.sh \
           {params.out_dir} \
           {threads} \
           "{params.cell_types}" \
           >> {log} 2>&1
           """

rule data_weights:
    output: config['plotting']['data_weights']['out_file']
    params: out_dir = config['plotting']['data_weights']['out_dir'],
            file_name = config['plotting']['data_weights']['file_name'],
            cell_types = config['cell_types'] 
    envmodules: "pigz"
    resources: time="3:00:00", threads = 10, mem_mb = 80000
    log:  config['plotting']['data_weights']['log']
    shell: 
           """
           scripts/data_sharing_twas_weights.sh \
             {params.out_dir} \
             {params.file_name} \
             "{params.cell_types}" \
             {threads} \
            >> {log} 2>&1
           """

#rule manuscript_tables_report:
#    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
#    input:  ctwas_multi = expand(../results/12CTWAS/multi/ctwas_multi_{gwas}_ctwas.rds, gwas = config['gwas']),
#            
#            rmd_script = "scripts/ctwas_report.Rmd"
#    output: "reports/12CTWAS/12ctwas_report.html"
#    params: in_dir = "../../results/12CTWAS/multi/",
#            bmark_dir = "../reports/benchmarks/",
#            lookup_dir = "../../resources/sheets/",
#            output_file = "../reports/12CTWAS/12ctwas_report.html"
#    singularity: config["containers"]["r_eqtl"] # Need to add ctwas to r_eqtl conatiner to print locus plot
#    message: "Generate cTWAS report"
#    benchmark: "reports/benchmarks/12ctwas.ctwas_report.benchmark.txt"
#    log:     "../results/00LOG/12CTWAS/ctwas_report.log"
#    shell:
#        """
#        Rscript -e "rmarkdown::render('{input.rmd_script}', \
#            output_file = '{params.output_file}', \
#            params = list(in_dir = '{params.in_dir}', \
#            bmark_dir = '{params.bmark_dir}', \
#            lookup_dir = '{params.lookup_dir}'))" > {log} 2>&1
#        """

#rule eqtl_boxplots:
#    output: config['plotting']['eqtl_boxplots']['output'] 
#    params: exp_dir = config['plotting']['eqtl_boxplots']['exp_dir'],
#            pval_dir = config['plotting']['eqtl_boxplots']['pval_dir'],
#            geno_prefix = config['plotting']['eqtl_boxplots']['geno_prefix'],
#            gene_id = config['plotting']['eqtl_boxplots']['gene_id'],
#            snp_id = config['plotting']['eqtl_boxplots']['snp_id']
#    singularity: config["containers"]["r_eqtl"]
#    resources: threads = 4, mem_mb = 20000
#    envmodules: "plink/2.0"
#    log:  config['plotting']['eqtl_boxplots']['log']
#    script: "../scripts/plot_eQTL_boxplots.R"

#rule eqtl_boxplots_py:
#    input:  pairs_csv = config['plotting']['eqtl_boxplots_py']['pairs_csv'],
#            geno = config['geno_post_impute']['exclude_SNPs']['output'],
#    output: config['plotting']['eqtl_boxplots_py']['output']
#    params: exp_dir = config['plotting']['eqtl_boxplots_py']['exp_dir'],
#            out_dir = config['plotting']['eqtl_boxplots_py']['out_dir'],
##    envmodules: "bcftools"
#    conda:  config["scanpy"]["env"]
#    resources: threads = 4, mem_mb = 20000
#    log:    config['plotting']['eqtl_boxplots_py']['log']
#    shell:  """
#            python3 scripts/eqtl_plot.py --pairs_file {input.pairs_csv} \
#              --genotype_file {input.geno} \
#              --expression_dir {params.exp_dir} \
#              --output_dir {params.out_dir} >> {log} 2>&1
#            """
