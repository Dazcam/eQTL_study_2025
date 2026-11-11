configfile: '../config/config.yaml'


# F1 is the infographic

rule eqtl_qc_fig:
    # Fig 2
    input:  ziffra_tbl = config['replication_sc']['ziffra_summary']['out_file']
    output: config['figures']['eqtl_qc_fig']['out_file']
    params: in_dir = config['figures']['eqtl_qc_fig']['in_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['figures']['eqtl_qc_fig']['log']
    script: "../scripts/manuscript_figure_eqtl_QC.R"

rule eqtl_boxplots_py:
    # Prep script for eQTL plots - Feeds Fig 3
    input:  pairs_csv = config['figures']['eqtl_boxplots_py']['pairs_csv'],
            geno = config['geno_post_impute']['exclude_SNPs']['output'],
    output: config['figures']['eqtl_boxplots_py']['output']
    params: exp_dir = config['figures']['eqtl_boxplots_py']['exp_dir'],
            out_dir = config['figures']['eqtl_boxplots_py']['out_dir'],
    envmodules: "BCFtools"
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 20000
    log:    config['figures']['eqtl_boxplots_py']['log']
    shell:  """
            python3 scripts/manuscript_figure_eqtl_boxplots.py --pairs_file {input.pairs_csv} \
              --genotype_file {input.geno} \
              --expression_dir {params.exp_dir} \
              --output_dir {params.out_dir} >> {log} 2>&1
            """

rule replication_fig:
    # Fig 3
    input:  beta_file = config['figures']['rep_fig']['jang_beta_file'],
            eqtl_boxplots_output  = config['figures']['eqtl_boxplots_py']['output']
    output: config['figures']['rep_fig']['out_file']
    params: in_dir = config['figures']['rep_fig']['in_dir'],
            internal_dir = config['figures']['rep_fig']['internal_dir'],
            jang_pi1_dir = config['figures']['rep_fig']['jang_pi1_dir'],
            gene_lookup  = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['figures']['rep_fig']['log']
    script: "../scripts/manuscript_figure_replication.R"

rule fig4:
    # Fig 4 -- full figure: A/B (pseudotime UMAP + features) + C/D (overlap barplots)
    input:  h5ad = config['pseudotime']['palantir']['h5ad'].format(trajectory=config['figures']['fig4']['trajectory'])
    output: config['figures']['fig4']['out_file']
    params: trajectory       = config['figures']['fig4']['trajectory'],
            root_gene        = config['pseudotime']['palantir']['root_gene'],
            terminal_markers = config['pseudotime']['palantir']['terminal_markers'],
            perm_dir         = config['tensorQTL']['tensorqtl_perm']['output'].split('{cell_type}')[0],
            pseudotime_dir   = config['figures']['fig4']['pseudotime_dir'],
            cell_types       = config['cell_types'],
            exp_pc_map       = config['exp_pc_map']
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 120000, time = "1:00:00"
    log:    config['figures']['fig4']['log']
    script: "../scripts/manuscript_figure_fig4.py"

rule ldsr_fig:
    # Fig 5
    output: config['figures']['ldsr_fig']['out_file']
    params: in_dir = config['figures']['ldsr_fig']['in_dir'],
            cell_types = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['figures']['ldsr_fig']['log']
    script: "../scripts/manuscript_figure_ldsr.R"

rule smr_ctwas_venn_fig:
    # Fig 6 and Supp Fig X
    input:  smr_sig_tbls = expand(config['tables']['smr_sig_tbl']['out_file'], gwas = list(config['gwas'].keys()))
    output: plot_scz  = config['figures']['smr_ctwas_venn']['out_file_scz'],
            plot_grid = config['figures']['smr_ctwas_venn']['out_file_grid']
    params: traits       = list(config['gwas'].keys()),
            smr_tbl_dir  = config['figures']['smr_ctwas_venn']['smr_tbl_dir'],
            ctwas_in_dir = config['figures']['smr_ctwas_venn']['ctwas_in_dir'],
            weights_dir  = config['figures']['smr_ctwas_venn']['weights_dir'],
            genes_out    = lambda wc: config['figures']['smr_ctwas_venn']['genes_out'],
            cell_types   = config['cell_types'],
            gene_lookup  = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time="3:30:00"
    log: config['figures']['smr_ctwas_venn']['log']
    script: "../scripts/manuscript_figure_smr_ctwas_venn.R"

rule supplementary_fig:
    # Supp Fig 2 & 12
    output: config['figures']['supp_fig']['out_file']
    params: geno_dir = config['geno_post_impute']['calc_genotype_pcs']['params'].rsplit('/', 1)[0] + '/',
            expr_dir = config['figures']['supp_fig']['expr_dir'],
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    resources: threads = 4, mem_mb = 20000
    log:  config['figures']['supp_fig']['log']
    script: "../scripts/manuscript_figure_supplementary.R"

rule family_heatmaps_fig:
    # SF 3-6 -- subcluster-specific family heatmaps
    input:  eqtl_effects = config['eqtl_sharing']['eqtl_effects']['out_file'],
            gene_lookup  = config['gene_lookup']
    output: glu_a = config['figures']['family_heatmaps_fig']['glu_a'],
            glu_b = config['figures']['family_heatmaps_fig']['glu_b'],
            gaba  = config['figures']['family_heatmaps_fig']['gaba'],
            npc   = config['figures']['family_heatmaps_fig']['npc']
    params: eqtl_effects = config['eqtl_sharing']['eqtl_effects']['out_file'],
            gene_lookup  = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['figures']['family_heatmaps_fig']['log']
    script:      "../scripts/manuscript_figure_family_heatmaps.R"

rule effect_size_heatmaps_fig:
    # SF 7-8 -- effect-size correlation heatmaps (discovery -> replication,
    input:  eqtl_effect_replication = config['eqtl_sharing']['eqtl_effect_replication']['out_file']
    output: fetal_vs_fetal = config['figures']['effect_size_heatmaps_fig']['fetal_vs_fetal'],
            fetal_vs_jang  = config['figures']['effect_size_heatmaps_fig']['fetal_vs_jang']
    params: eqtl_effect_replication = config['eqtl_sharing']['eqtl_effect_replication']['out_file']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['figures']['effect_size_heatmaps_fig']['log']
    script:      "../scripts/manuscript_figure_effect_size_heatmaps.R"

rule bulk_egene_overlaps_fig:
    # SF 10
    input:  egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file']
    output: config['figures']['bulk_egene_overlaps_fig']['out_file']
    params: egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file'],
            cell_types          = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources:   time="0:30:00"
    log:         config['figures']['bulk_egene_overlaps_fig']['log']
    script:	 "../scripts/manuscript_figure_bulk_egene_overlaps.R"

rule snrna_qc_fig:
    # Supp Fig X -- snRNA-seq QC per sample and per cell type (reviewer, line 473)
    input:  h5ad = config['scanpy']['clusters_h5ad']
    output: per_sample   = config['figures']['snrna_qc_fig']['out_file_sample'],
            per_celltype = config['figures']['snrna_qc_fig']['out_file_celltype']
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 64000, time = "1:00:00"
    log:    config['figures']['snrna_qc_fig']['log']
    script: "../scripts/manuscript_figure_snrna_qc.py"


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["17figures"] = [rules.eqtl_qc_fig.output,
                        rules.replication_fig.output,
                        rules.fig4.output,
                        rules.ldsr_fig.output,
                        rules.smr_ctwas_venn_fig.output,
                        rules.supplementary_fig.output,
                        rules.snrna_qc_fig.output,
                        rules.eqtl_boxplots_py.output,
                        rules.bulk_egene_overlaps_fig.output,
                        rules.family_heatmaps_fig.output,
                        rules.effect_size_heatmaps_fig.output]
