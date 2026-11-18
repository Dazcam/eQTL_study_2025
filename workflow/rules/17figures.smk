configfile: '../config/config.yaml'

FIG = config['figures']

BOXPLOT_PAIRS = FIG['eqtl_boxplots_py']['pairs']
BOXPLOT_EXP = {f"exp_{panel}": config['tensorQTL']['prep_tensorQTL_input']['exp_out'].format(
                   cell_type = pair['cell_type'], norm_method = config['tensorQTL']['norm_methods'][0])
               for panel, pair in BOXPLOT_PAIRS.items()}
BOXPLOT_CSV = {panel: FIG['eqtl_boxplots_py']['out_csv'].format(**pair)
               for panel, pair in BOXPLOT_PAIRS.items()}


# Fig 1A: made outside the pipeline

rule umap_fig:
    # Fig 1B-D
    input:  h5ad = config['scanpy']['clusters_h5ad']
    output: bc_pdf = FIG['umap_fig']['bc_pdf'],
            bc_svg = FIG['umap_fig']['bc_svg'],
            d_pdf  = FIG['umap_fig']['d_pdf'],
            d_svg  = FIG['umap_fig']['d_svg']
    params: labels   = config['manuscript_labels'],
            l1_types = config['fetal_l1_types']
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 150000, time = "2:00:00"
    benchmark: FIG['umap_fig']['benchmark']
    log:       FIG['umap_fig']['log']
    script:    FIG['umap_fig']['script']

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

rule pseudotime_fig:
    # Fig 4
    input:  h5ad = config['pseudotime']['palantir']['h5ad'].format(trajectory = FIG['pseudotime_fig']['trajectory'])
    output: FIG['pseudotime_fig']['out_file']
    params: trajectory       = FIG['pseudotime_fig']['trajectory'],
            labels           = config['manuscript_labels'],
            root_gene        = config['pseudotime']['palantir']['root_gene'],
            terminal_markers = config['pseudotime']['palantir']['terminal_markers'],
            perm_dir         = config['tensorQTL']['tensorqtl_perm']['output'].split('{cell_type}')[0],
            pseudotime_dir   = FIG['pseudotime_fig']['pseudotime_dir'],
            cell_types       = config['cell_types'],
            exp_pc_map       = config['exp_pc_map'],
            geno_pc          = config['tensorQTL']['geno_pcs'],
            norm_method      = config['tensorQTL']['norm_methods'][0]
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 100000, time = "1:00:00"
    benchmark: FIG['pseudotime_fig']['benchmark']
    log:       FIG['pseudotime_fig']['log']
    script:    FIG['pseudotime_fig']['script']

rule ldsr_fig:
    # Fig 5
    input:  config['sldsr']['ldsr_strat_summary']['output'].format(annot_type = config['tables']['ldsr_tbl']['annot_type'])
    output: FIG['ldsr_fig']['out_file']
    params: labels     = config['manuscript_labels'],
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

rule subcluster_umap_fig:
    # Supp Fig 1
    input:  h5ad = expand(config['scanpy']['subclusters_h5ad'], l1_pop = FIG['subcluster_umap_fig']['panel_order'])
    output: pdf = FIG['subcluster_umap_fig']['out_file']
    params: panel_order = FIG['subcluster_umap_fig']['panel_order'],
            labels      = config['manuscript_labels']
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 80000, time = "1:00:00"
    benchmark: FIG['subcluster_umap_fig']['benchmark']
    log:       FIG['subcluster_umap_fig']['log']
    script:    FIG['subcluster_umap_fig']['script']

rule supplementary_fig:
    # Supp Figs 2 and 11
    input:  eigenval = config['geno_post_impute']['calc_genotype_pcs']['params'] + '.eigenval'
    output: expr_pcs = FIG['supplementary_fig']['out_file_expr'],
            geno_pcs = FIG['supplementary_fig']['out_file_geno']
    params: labels      = config['manuscript_labels'],
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
    params: labels     = config['manuscript_labels'],
            cell_types = config['cell_types']
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
    params: labels     = config['manuscript_labels'],
            cell_types = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['effect_size_heatmaps_fig']['benchmark']
    log:       FIG['effect_size_heatmaps_fig']['log']
    script:    FIG['effect_size_heatmaps_fig']['script']

rule bulk_egene_overlaps_fig:
    # Supp Fig 9
    input:  egenes_per_celltype = config['replication_bulk']['extract_unique_egenes']['out_file'],
            overlap             = config['replication_bulk']['overlap_obrien']['out_file']
    output: FIG['bulk_egene_overlaps_fig']['out_file']
    params: labels     = config['manuscript_labels'],
            cell_types = config['cell_types']
    singularity: config["containers"]["r_eqtl"]
    resources: time = "0:30:00"
    benchmark: FIG['bulk_egene_overlaps_fig']['benchmark']
    log:       FIG['bulk_egene_overlaps_fig']['log']
    script:    FIG['bulk_egene_overlaps_fig']['script']

# Supp Fig 10: made outside the pipeline

rule snrna_qc_fig:
    # Supp Fig 12
    input:  h5ad = config['scanpy']['clusters_h5ad']
    output: per_sample = FIG['snrna_qc_fig']['out_file'],
            summary    = FIG['snrna_qc_fig']['summary']
    params: labels   = config['manuscript_labels'],
            l1_types = config['fetal_l1_types']
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 20000, time = "1:00:00"
    benchmark: FIG['snrna_qc_fig']['benchmark']
    log:       FIG['snrna_qc_fig']['log']
    script:    FIG['snrna_qc_fig']['script']


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["17figures"] = [rules.umap_fig.output,
                        rules.eqtl_qc_fig.output,
                        rules.eqtl_boxplots_py.output,
                        rules.replication_fig.output,
                        rules.pseudotime_fig.output,
                        rules.ldsr_fig.output,
                        rules.smr_ctwas_venn_fig.output,
                        rules.subcluster_umap_fig.output,
                        rules.supplementary_fig.output,
                        rules.family_heatmaps_fig.output,
                        rules.effect_size_heatmaps_fig.output,
                        rules.bulk_egene_overlaps_fig.output,
                        rules.snrna_qc_fig.output]
