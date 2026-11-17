configfile: '../config/config.yaml'

FIG = config['figures']

BOXPLOT_PAIRS = FIG['eqtl_boxplots_py']['pairs']
BOXPLOT_EXP = {f"exp_{panel}": config['tensorQTL']['prep_tensorQTL_input']['exp_out'].format(
                   cell_type = pair['cell_type'], norm_method = config['tensorQTL']['norm_methods'][0])
               for panel, pair in BOXPLOT_PAIRS.items()}
BOXPLOT_CSV = {panel: FIG['eqtl_boxplots_py']['out_csv'].format(**pair)
               for panel, pair in BOXPLOT_PAIRS.items()}


rule eqtl_qc_fig:
    # Fig 2
    input:  ziffra_tbl    = config['replication_sc']['ziffra_summary']['out_file'],
            nuclei_counts = config['tables']['cell_counts_tbl']['out_file']
    output: FIG['eqtl_qc_fig']['out_file']
    params: in_dir      = FIG['eqtl_qc_fig']['in_dir'],
            labels      = config['manuscript_labels'],
            cell_types  = config['cell_types'],
            l1_types    = config['fetal_l1_types'],
            exp_pc_map  = config['exp_pc_map'],
            geno_pc     = config['tensorQTL']['geno_pcs'],
            norm_method = config['tensorQTL']['norm_methods'][0],
            eqtl_fdr    = config['eqtl_fdr']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['eqtl_qc_fig']['benchmark']
    log:       FIG['eqtl_qc_fig']['log']
    script:    FIG['eqtl_qc_fig']['script']

rule eqtl_boxplots_py:
    # Fig 3E and 3G boxplot data
    input:  vcf = config['geno_post_impute']['exclude_SNPs']['output'],
            **BOXPLOT_EXP
    output: **BOXPLOT_CSV
    params: pairs       = BOXPLOT_PAIRS,
            plot_checks = FIG['eqtl_boxplots_py']['plot_checks']
    envmodules: "BCFtools"
    conda:  config["scanpy"]["env"]
    resources: time = "0:30:00", mem_mb = 20000
    benchmark: FIG['eqtl_boxplots_py']['benchmark']
    log:       FIG['eqtl_boxplots_py']['log']
    script:    FIG['eqtl_boxplots_py']['script']

rule replication_fig:
    # Fig 3
    input:  beta_file = FIG['replication_fig']['jang_beta_file'],
            boxplot_e = rules.eqtl_boxplots_py.output.fig3e,
            boxplot_g = rules.eqtl_boxplots_py.output.fig3g
    output: FIG['replication_fig']['out_file']
    params: in_dir        = FIG['replication_fig']['in_dir'],
            labels        = config['manuscript_labels'],
            internal_dir  = FIG['replication_fig']['internal_dir'],
            jang_pi1_dir  = FIG['replication_fig']['jang_pi1_dir'],
            gene_lookup   = config['gene_lookup'],
            boxplot_pairs = BOXPLOT_PAIRS,
            cell_types    = config['cell_types'],
            l1_types      = config['fetal_l1_types'],
            jang_types    = config['jang_l1_types'],
            exp_pc_map    = config['exp_pc_map'],
            geno_pc       = config['tensorQTL']['geno_pcs'],
            norm_method   = config['tensorQTL']['norm_methods'][0],
            eqtl_fdr      = config['eqtl_fdr']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['replication_fig']['benchmark']
    log:       FIG['replication_fig']['log']
    script:    FIG['replication_fig']['script']

rule fig4:
    # Fig 4
    input:  h5ad = config['pseudotime']['palantir']['h5ad'].format(trajectory = FIG['fig4']['trajectory'])
    output: FIG['fig4']['out_file']
    params: trajectory       = FIG['fig4']['trajectory'],
            labels           = config['manuscript_labels'],
            root_gene        = config['pseudotime']['palantir']['root_gene'],
            terminal_markers = config['pseudotime']['palantir']['terminal_markers'],
            perm_dir         = config['tensorQTL']['tensorqtl_perm']['output'].split('{cell_type}')[0],
            pseudotime_dir   = FIG['fig4']['pseudotime_dir'],
            cell_types       = config['cell_types'],
            exp_pc_map       = config['exp_pc_map'],
            geno_pc          = config['tensorQTL']['geno_pcs'],
            norm_method      = config['tensorQTL']['norm_methods'][0]
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 100000, time = "1:00:00"
    benchmark: FIG['fig4']['benchmark']
    log:       FIG['fig4']['log']
    script:    FIG['fig4']['script']

rule ldsr_fig:
    # Fig 5
    output: FIG['ldsr_fig']['out_file']
    params: in_dir     = FIG['ldsr_fig']['in_dir'],
            labels     = config['manuscript_labels'],
            cell_types = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['ldsr_fig']['benchmark']
    log:       FIG['ldsr_fig']['log']
    script:    FIG['ldsr_fig']['script']

rule smr_ctwas_venn_fig:
    # Fig 6
    input:  genes_tbls = expand(config['tables']['smr_ctwas_venn_tbl']['genes_out'], gwas = list(config['gwas'].keys()))
    output: plot_scz  = FIG['smr_ctwas_venn_fig']['out_file_scz'],
            plot_grid = FIG['smr_ctwas_venn_fig']['out_file_grid']
    params: traits = list(config['gwas'].keys())
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['smr_ctwas_venn_fig']['benchmark']
    log:       FIG['smr_ctwas_venn_fig']['log']
    script:    FIG['smr_ctwas_venn_fig']['script']

rule supplementary_fig:
    # Supp Figs 2 and 11
    output: FIG['supplementary_fig']['out_file']
    params: geno_dir    = config['geno_post_impute']['calc_genotype_pcs']['params'].rsplit('/', 1)[0] + '/',
            labels      = config['manuscript_labels'],
            expr_dir    = FIG['supplementary_fig']['expr_dir'],
            cell_types  = config['cell_types'],
            l1_types    = config['fetal_l1_types'],
            exp_pcs     = config['tensorQTL']['exp_pcs'],
            geno_pc     = config['tensorQTL']['geno_pcs'],
            norm_method = config['tensorQTL']['norm_methods'][0]
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['supplementary_fig']['benchmark']
    log:       FIG['supplementary_fig']['log']
    script:    FIG['supplementary_fig']['script']

rule family_heatmaps_fig:
    # Supp Figs 3-6
    input:  eqtl_effects = config['eqtl_sharing']['eqtl_effects']['out_file'],
            gene_lookup  = config['gene_lookup']
    output: glu_a = FIG['family_heatmaps_fig']['glu_a'],
            glu_b = FIG['family_heatmaps_fig']['glu_b'],
            gaba  = FIG['family_heatmaps_fig']['gaba'],
            npc   = FIG['family_heatmaps_fig']['npc']
    params: eqtl_effects = config['eqtl_sharing']['eqtl_effects']['out_file'],
            labels       = config['manuscript_labels'],
            gene_lookup  = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['family_heatmaps_fig']['benchmark']
    log:       FIG['family_heatmaps_fig']['log']
    script:    FIG['family_heatmaps_fig']['script']

rule effect_size_heatmaps_fig:
    # Supp Figs 7 and 8
    input:  eqtl_effect_replication = config['eqtl_sharing']['eqtl_effect_replication']['out_file']
    output: fetal_vs_fetal = FIG['effect_size_heatmaps_fig']['fetal_vs_fetal'],
            fetal_vs_jang  = FIG['effect_size_heatmaps_fig']['fetal_vs_jang']
    params: eqtl_effect_replication = config['eqtl_sharing']['eqtl_effect_replication']['out_file'],
            labels                  = config['manuscript_labels']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['effect_size_heatmaps_fig']['benchmark']
    log:       FIG['effect_size_heatmaps_fig']['log']
    script:    FIG['effect_size_heatmaps_fig']['script']

rule bulk_egene_overlaps_fig:
    # Supp Fig 9A
    input:  egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file']
    output: FIG['bulk_egene_overlaps_fig']['out_file']
    params: egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            labels              = config['manuscript_labels'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file'],
            cell_types          = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['bulk_egene_overlaps_fig']['benchmark']
    log:       FIG['bulk_egene_overlaps_fig']['log']
    script:    FIG['bulk_egene_overlaps_fig']['script']

rule snrna_qc_fig:
    # Supp Fig 12
    input:  h5ad = config['scanpy']['clusters_h5ad']
    output: per_sample   = FIG['snrna_qc_fig']['out_file_sample'],
            per_celltype = FIG['snrna_qc_fig']['out_file_celltype']
    params: labels = config['manuscript_labels']
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 20000, time = "1:00:00"
    benchmark: FIG['snrna_qc_fig']['benchmark']
    log:       FIG['snrna_qc_fig']['log']
    script:    FIG['snrna_qc_fig']['script']


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["17figures"] = [rules.eqtl_qc_fig.output,
                        rules.eqtl_boxplots_py.output,
                        rules.replication_fig.output,
                        rules.fig4.output,
                        rules.ldsr_fig.output,
                        rules.smr_ctwas_venn_fig.output,
                        rules.supplementary_fig.output,
                        rules.family_heatmaps_fig.output,
                        rules.effect_size_heatmaps_fig.output,
                        rules.bulk_egene_overlaps_fig.output,
                        rules.snrna_qc_fig.output]
