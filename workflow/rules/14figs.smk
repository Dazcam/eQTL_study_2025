configfile: '../config/config.yaml'

rule all:
    input:
#        config['plotting']['eqtl_boxplots_py']['output'],
#        config['plotting']['eqtl_qc_plt']['out_file'],
#        config['plotting']['ldsr_plt']['out_file'],
#        config['plotting']['rep_plt']['out_file'],
#        config['plotting']['supp_plt']['out_file'],
#        config['plotting']['data_perm']['out_file'],
#        config['plotting']['data_nominal']['out_file'],
#        config['plotting']['data_mk_eqtl_tar']['out_file']
#        config['plotting']['data_weights']['out_file']
#        "../results/13MANUSCRIPT_PLOTS_TABLES/data_sharing/eqtl_atlas_genotypes.vcf.gz"
#       config['plotting']['smr_ctwas_venn']['out_file_grid']
#       config["plotting"]["diff_exp_tbl"]["html_out"],
#       config["plotting"]["diff_exp_tbl"]["xlsx_out"] 
#       config['plotting']['eqtl_boxplots_py']['output'],
#       config['plotting']['bulk_egene_overlaps_plt']['out_file'],
#       config['plotting']['family_heatmaps_plt']['glu_a'],
#       config['plotting']['family_heatmaps_plt']['glu_b'],
#       config['plotting']['family_heatmaps_plt']['gaba'],
#       config['plotting']['family_heatmaps_plt']['npc'],
#       config['plotting']['effect_size_heatmaps_plt']['fetal_vs_fetal'],
#       config['plotting']['effect_size_heatmaps_plt']['fetal_vs_jang'],
#       config['plotting']['fig4_plt']['out_file'],
#       config['plotting']['tss_density_plt']['out_file']
       config['plotting']['snrna_qc_plt']['out_file_sample'],
       config['plotting']['snrna_qc_plt']['out_file_celltype'],


rule eqtl_qc_plt:
    # Fig 2
    output: config['plotting']['eqtl_qc_plt']['out_file']
    params: in_dir = config['plotting']['eqtl_qc_plt']['in_dir'],
            ziffra_dir = config['plotting']['eqtl_qc_plt']['ziffra_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['eqtl_qc_plt']['log']
    script: "../scripts/manuscript_plot_eqtl_QC.R"

rule eqtl_boxplots_py:
    # Prep script for eQTL plots - Feeds Fig 3
    input:  pairs_csv = config['plotting']['eqtl_boxplots_py']['pairs_csv'],
            geno = config['geno_post_impute']['exclude_SNPs']['output'],
    output: config['plotting']['eqtl_boxplots_py']['output']
    params: exp_dir = config['plotting']['eqtl_boxplots_py']['exp_dir'],
            out_dir = config['plotting']['eqtl_boxplots_py']['out_dir'],
    envmodules: "BCFtools"
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 20000
    log:    config['plotting']['eqtl_boxplots_py']['log']
    shell:  """
            python3 scripts/manuscript_eqtl_plot.py --pairs_file {input.pairs_csv} \
              --genotype_file {input.geno} \
              --expression_dir {params.exp_dir} \
              --output_dir {params.out_dir} >> {log} 2>&1
            """

rule replication_plt:
    # Fig 3
    input:  beta_file = config['plotting']['rep_plt']['jang_beta_file'],
            eqtl_boxplots_output  = config['plotting']['eqtl_boxplots_py']['output']
    output: config['plotting']['rep_plt']['out_file']
    params: in_dir = config['plotting']['rep_plt']['in_dir'],
            internal_dir = config['plotting']['rep_plt']['internal_dir'],
            jang_pi1_dir = config['plotting']['rep_plt']['jang_pi1_dir']
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['rep_plt']['log']
    script: "../scripts/manuscript_plot_replication.R"

rule fig4_plt:
    # Fig 4 -- full figure: A/B (pseudotime UMAP + features) + C/D (overlap barplots)
    input:  glu_a_h5ad = config['pseudotime']['palantair']['h5ad'].format(trajectory="NPC-to-Glu-UL"),
            glu_b_h5ad = config['pseudotime']['palantair']['h5ad'].format(trajectory="NPC-to-Glu-DL")
    output: config['plotting']['fig4_plt']['out_file']
    params: perm_dir       = config['dev_specificity']['extract_universe']['qtl_dir'],
            pseudotime_dir = config['plotting']['fig4_plt']['pseudotime_dir'],
            cell_types     = config['cell_types'],
            exp_pc_map     = config['exp_pc_map']
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 120000, time = "1:00:00"
    log:    config['plotting']['fig4_plt']['log']
    script: "../scripts/manuscript_plot_fig4.py"

rule ldsr_plt:
    # Fig 5
    output: config['plotting']['ldsr_plt']['out_file']
    params: in_dir = config['plotting']['ldsr_plt']['in_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['plotting']['ldsr_plt']['log']
    script: "../scripts/manuscript_plot_ldsr.R"

rule compare_smr_ctwas:
    # Fig 6 and Supp Fig X
    output: plot_grid = config['plotting']['smr_ctwas_venn']['out_file_grid']
    singularity: config["containers"]["r_eqtl"]
    resources: time="3:30:00"
    log: config['plotting']['smr_ctwas_venn']['log']
    script: "../scripts/manuscript_plot_compare_smr_ctwas.R"

rule supplementary_plt:
    # Supp Fig 1
    output: config['plotting']['supp_plt']['out_file']
    params: geno_dir = config['plotting']['supp_plt']['geno_dir'],
            expr_dir = config['plotting']['supp_plt']['expr_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    resources: threads = 4, mem_mb = 20000
    log:  config['plotting']['supp_plt']['log']
    script: "../scripts/manuscript_plot_supplementary.R"

rule family_heatmaps_plt:
    # SF 3-6 -- subcluster-specific family heatmaps
    input:  eqtl_effects = config['dev_specificity']['eqtl_effects']['out_file'],
            gene_lookup  = config['dev_specificity']['report']['gene_lookup']
    output: glu_a = config['plotting']['family_heatmaps_plt']['glu_a'],
            glu_b = config['plotting']['family_heatmaps_plt']['glu_b'],
            gaba  = config['plotting']['family_heatmaps_plt']['gaba'],
            npc   = config['plotting']['family_heatmaps_plt']['npc']
    params: eqtl_effects = config['dev_specificity']['eqtl_effects']['out_file'],
            gene_lookup  = config['dev_specificity']['report']['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['plotting']['family_heatmaps_plt']['log']
    script:      "../scripts/manuscript_plot_family_heatmaps.R"

rule effect_size_heatmaps_plt:
    # SF 7-8 -- effect-size correlation heatmaps (discovery -> replication,
    input:  eqtl_effect_replication = config['dev_specificity']['eqtl_effect_replication']['out_file']
    output: fetal_vs_fetal = config['plotting']['effect_size_heatmaps_plt']['fetal_vs_fetal'],
            fetal_vs_jang  = config['plotting']['effect_size_heatmaps_plt']['fetal_vs_jang']
    params: eqtl_effect_replication = config['dev_specificity']['eqtl_effect_replication']['out_file']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['plotting']['effect_size_heatmaps_plt']['log']
    script:      "../scripts/manuscript_plot_effect_size_heatmaps.R"

rule pseud_tss_density_plt:
    # SF 9 -- combined TSS-distance density plots per trajectory
    output: config['plotting']['tss_density_plt']['out_file']
    params: pseudotime_dir = config['plotting']['tss_density_plt']['pseudotime_dir'],
            exp_pc_map     = config['exp_pc_map'],
            trajectories   = config['trajectories']
    singularity: config["containers"]["r_eqtl"]
    resources: time="0:30:00"
    log:    config['plotting']['tss_density_plt']['log']
    script: "../scripts/manuscript_plot_tss_density.R"

rule bulk_egene_overlaps_plt:
    # SF 10
    input:  egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file']
    output: config['plotting']['bulk_egene_overlaps_plt']['out_file']
    params: egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['plotting']['bulk_egene_overlaps_plt']['log']
    script:	 "../scripts/manuscript_plot_bulk_egene_overlaps.R"

rule snrna_qc_plt:
    # Supp Fig X -- snRNA-seq QC per sample and per cell type (reviewer, line 473)
    input:  h5ad = config['plotting']['snrna_qc_plt']['in_h5ad']
    output: per_sample   = config['plotting']['snrna_qc_plt']['out_file_sample'],
            per_celltype = config['plotting']['snrna_qc_plt']['out_file_celltype']
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 64000, time = "1:00:00"
    log:    config['plotting']['snrna_qc_plt']['log']
    script: "../scripts/manuscript_plot_snrna_qc.py"

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
#    envmodules: "PLINK"
#    log:  config['plotting']['eqtl_boxplots']['log']
#    script: "../scripts/plot_eQTL_boxplots.R"
