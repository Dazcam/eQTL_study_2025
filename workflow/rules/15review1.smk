configfile: "../config/config.yaml" 

rule all:
    input: 
#        config["review"]["singler_report"]["html"],
#        config["review"]["batch_corr_report"]["html"],
#        config["review"]["pseudotime_paga"]["h5ad"],
#        config["review"]["pseudotime_paga"]["paga_plot"],
#        config["review"]["pseudotime_paga"]["dpt_plot"],
#        config["review"]["pseudotime_paga"]["fate_plot"],
#        config["review"]["pseudotime_paga"]["summary_csv"],
#        config["review"]["isoform_eqtl"]["enst_ensg_map"],
#        config["review"]["isoform_eqtl"]["egenes_combined"],
#        config["review"]["isoform_eqtl"]["bulk_hits"],
#        config["review"]["isoform_eqtl"]["summary"],
#        config["review"]["isoform_eqtl"]["detail"]      
#        config["review"]["singler_snrnaseq"]["sce_cameron"],
#        config["review"]["singler_snrnaseq"]["sce_polioudakis"],
#        config["review"]["singler_snrnaseq"]["sce_nowakowski"],
#       config["review"]["mashr"]["sentinel"],
#        config["review"]["mashr_report"]["html"]
#        config["review"]["palantir"]["fig_corr"],
        config["review"]["palantir_report"]["html"]

# Move to set up / or ref download rule
rule dwnld_qian: 
    output: config["review"]["dwnld_qian"]["output"]
    params: web_link = config["review"]["dwnld_qian"]["web_link"],
    message: "Download snRNAseq spatial data from Qian 2025, PMID:40369074"
    benchmark: "reports/benchmarks/15review.dwnld_qian.txt"
    log:    config["review"]["dwnld_qian"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule prep_qian_ref:
    input:  config["review"]["dwnld_qian"]["output"]
    output: config["review"]["prep_qian_ref"]["output"]
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 64000, time = "2:00:00"
    params:    n_cells_per_cluster = 500,
               random_seed         = 1234
    log:    config["review"]["prep_qian_ref"]["log"]
    message: "Preparing Qian 2025 MERFISH reference: filtering V1/V2, subsampling H2 clusters"
    script: config["review"]["prep_qian_ref"]["script"]

rule prep_singler_query:
    input:   query_h5ad = config["review"]["prep_singler_query"]["query_h5ad"],
             ref_h5ad   = config["review"]["prep_qian_ref"]["output"]
    output:  h5ad = config["review"]["prep_singler_query"]["h5ad"],
             umap_csv = config["review"]["prep_singler_query"]["umap_csv"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    log:     config["review"]["prep_singler_query"]["log"]
    message: "Preparing snRNA-seq query: assigning L1 labels, subsetting to MERFISH panel"
    script:  config["review"]["prep_singler_query"]["script"]

rule singler_classify:
    input:   ref_h5ad   = config["review"]["prep_qian_ref"]["output"],
             query_h5ad = config["review"]["prep_singler_query"]["h5ad"]
    output:  csv = config["review"]["singler_classify"]["csv"],
             pred_rds = config["review"]["singler_classify"]["pred_rds"],
             sce_rds  = config["review"]["singler_classify"]["sce_rds"]             
    singularity: config["containers"]["singler"]
    resources: threads  = 8, mem_mb = 128000, time = "2:00:00"
    threads: 8
    log:     config["review"]["singler_classify"]["log"]
    message: "Running SingleR classification: Qian H2 reference vs snRNA-seq query"
    script:  config["review"]["singler_classify"]["script"]

rule batch_correction_lisi:
    input:   h5ad = config["review"]["prep_singler_query"]["query_h5ad"]
    output:  csv  = config["review"]["batch_correction_lisi"]["csv"],
             plot = config["review"]["batch_correction_lisi"]["plot"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    threads: 10
    params:  batch_col     = "sample",
             cell_type_col = "cell_type",
             random_seed   = config["review"]["batch_correction_lisi"].get("random_seed", 1234),
             perplexity    = config["review"]["batch_correction_lisi"].get("perplexity", 30)
    log:     config["review"]["batch_correction_lisi"]["log"]
    message: "Running LISI (iLISI / cLISI) batch correction assessment pre- and post-Harmony"
    script:  config["review"]["batch_correction_lisi"]["script"]

rule batch_correction_kbet:
    input:   h5ad = config["review"]["prep_singler_query"]["query_h5ad"]
    output:  csv  = config["review"]["batch_correction_kbet"]["csv"],
             plot = config["review"]["batch_correction_kbet"]["plot"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    threads: 10
    params:  batch_col      = "sample",
             cell_type_col  = "cell_type",
             n_subsample    = config["review"]["batch_correction_kbet"].get("n_subsample", 10000),
             random_seed    = config["review"]["batch_correction_kbet"].get("random_seed", 1234),
             k0             = config["review"]["batch_correction_kbet"].get("k0", 30)
    log:     config["review"]["batch_correction_kbet"]["log"]
    message: "Running kBET batch correction assessment (pre- and post-Harmony) on subsampled cells"
    script:  config["review"]["batch_correction_kbet"]["script"]

rule batch_corr_report:
    input:   kbet_csv  = config["review"]["batch_correction_kbet"]["csv"],
             lisi_csv  = config["review"]["batch_correction_lisi"]["csv"].replace('.csv', '_summary.csv'),
             kbet_plot = config["review"]["batch_correction_kbet"]["plot"],
             lisi_plot = config["review"]["batch_correction_lisi"]["plot"]
    output:  html = config["review"]["batch_corr_report"]["html"]
    conda:   config["scanpy"]["env"]
    resources:  mem_mb  = 8000,
                time    = "0:30:00"
    params:     n_cells      = 679738,
                n_batches    = 136,
                n_cell_types = 7,
                harmony_vars = "sample"
    log:     config["review"]["batch_corr_report"]["log"]
    message: "Generating batch correction QC report"
    script:  config["review"]["batch_corr_report"]["script"]

rule singler_report:
    input:   singler_csv = config["review"]["singler_classify"]["csv"],
             umap_csv    = config["review"]["prep_singler_query"]["umap_csv"],
             pred_rds    = config["review"]["singler_classify"]["pred_rds"],
             sce_rds     = config["review"]["singler_classify"]["sce_rds"],             
             rmd_script  = "scripts/review_singler_report.Rmd"
    output:  config["review"]["singler_report"]["html"]
    params:  singler_csv = "../../results/15REVIEW/singler_results.csv", # Paths relative to the Rmd script location (scripts/)
             umap_csv    = "../../results/15REVIEW/umap_coords.csv",
             pred_rds    = "../../results/15REVIEW/singler_pred.rds",
             sce_rds     = "../../results/15REVIEW/singler_sce_sub.rds",
             output_file = "../reports/15REVIEW/15singler_report.html"
    resources:   threads = 4, mem_mb = 32000, time = "1:00:00"
    singularity: config["containers"]["singler"]
    log:     config["review"]["singler_report"]["log"]
    message: "Generating SingleR spatial validation report"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(
                singler_csv = '{params.singler_csv}',
                umap_csv    = '{params.umap_csv}',
                pred_rds    = '{params.pred_rds}',
                sce_rds     = '{params.sce_rds}'))" > {log} 2>&1
        """

rule pseudotime_paga:
    input:   h5ad = "../results/02SCANPY/adata_clusters_v3.h5ad"
    output:  h5ad        = config["review"]["pseudotime_paga"]["h5ad"],
             paga_plot   = config["review"]["pseudotime_paga"]["paga_plot"],
             dpt_plot    = config["review"]["pseudotime_paga"]["dpt_plot"],
             fate_plot   = config["review"]["pseudotime_paga"]["fate_plot"],
             summary_csv = config["review"]["pseudotime_paga"]["summary_csv"]
    conda:   config["review"]["pseudotime_paga"]["env"]
    resources: threads = 16, mem_mb = 120000, time = "2-00:00"
    threads: 16
    params:  root_cell_type = config["review"]["pseudotime_paga"].get("root_cell_type", "NPC"),
             n_dcs          = config["review"]["pseudotime_paga"].get("n_dcs", 15),
             paga_threshold = config["review"]["pseudotime_paga"].get("paga_threshold", 0.05),
             random_seed    = config["review"]["pseudotime_paga"].get("random_seed", 1234)
    log:     config["review"]["pseudotime_paga"]["log"]
    message: "Running PAGA + DPT + CellRank pseudotime analysis on L1 cell populations"
    script:  config["review"]["pseudotime_paga"]["script"]

rule palantir_pseudotime:
    input:   h5ad = config["review"]["pseudotime_paga"]["h5ad"]
    output:  h5ad         = config["review"]["palantir"]["h5ad"],
             summary      = config["review"]["palantir"]["summary"],
             fig_auto     = config["review"]["palantir"]["fig_auto"],
             fig_manual   = config["review"]["palantir"]["fig_manual"],
             fig_corr     = config["review"]["palantir"]["fig_corr"]
    conda:   config["review"]["palantir"]["env"]
    resources: threads = 8, mem_mb = 120000, time  = "6:00:00"
    threads: 8
    params:  random_seed = config["review"]["palantir"].get("random_seed", 1234),
             n_dcs       = config["review"]["palantir"].get("n_dcs", 15)
    log:     config["review"]["palantir"]["log"]
    message: "Running Palantir pseudotime and fate probability analysis"
    script:  config["review"]["palantir"]["script"]

rule palantir_report:
    input:
        adata = config["review"]["palantir"]["h5ad"],
        nb    = config["review"]["palantir_report"]["notebook"]
    output:
        config["review"]["palantir_report"]["html"]
    conda:   config["scanpy"]["env"]
    resources:
        threads = 4,
        mem_mb  = 80000,
        time    = "1:00:00"
    params:
        nb_out         = lambda wildcards: os.path.abspath(config["review"]["palantir_report"]["nb_out"]),
        html_dir       = lambda wildcards: os.path.abspath(config["review"]["palantir_report"]["html_dir"]),
        adata_abs      = lambda wildcards, input: os.path.abspath(input.adata),
        paga_threshold = config["review"]["pseudotime_paga"].get("paga_threshold", 0.05),
        n_cells        = 679738,
        age_range      = "12-21 PCW",
        n_samples      = 134
    log:     config["review"]["palantir_report"]["log"]
    message: "Rendering Palantir trajectory report"
    shell:
        "papermill {input.nb} {params.nb_out} "
        "--execution-timeout -1 --request-save-on-cell-execute "
        "-p adata_path {params.adata_abs} "
        "-p paga_threshold {params.paga_threshold} "
        "-p n_cells {params.n_cells} "
        "-p age_range '{params.age_range}' "
        "-p n_samples {params.n_samples} "
        ">> {log} 2>&1 && "
        "jupyter nbconvert --to html {params.nb_out} "
        "--no-input --no-prompt --template lab "
        "--output-dir {params.html_dir} "
        ">> {log} 2>&1"

rule isoform_eqtl_comparison:
    input:
        tensorqtl_done = expand(
            "../results/05TENSORQTL/tensorqtl_perm/{ct}_quantile_genPC_4_expPC_{npc}/{ct}_quantile_perm.cis_qtl.txt.gz",
            zip,
            ct  = config["cell_types"],
            npc = [config["exp_pc_map"][ct] for ct in config["cell_types"]]
        ),
        gtf = config["parse"]["get_refs"]["gtf_out"],
        bulk_eqtl = config["review"]["isoform_eqtl"]["bulk_eqtl"]
    output:  enst_ensg_map   = config["review"]["isoform_eqtl"]["enst_ensg_map"],
             egenes_combined = config["review"]["isoform_eqtl"]["egenes_combined"],
             bulk_hits       = config["review"]["isoform_eqtl"]["bulk_hits"],
             summary         = config["review"]["isoform_eqtl"]["summary"],
             detail          = config["review"]["isoform_eqtl"]["detail"]
    conda:   config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 48000, time = "3:00:00"
    threads: 8
    params:  cell_types      = config["cell_types"],
             exp_pc_map      = config["exp_pc_map"],
             tensorqtl_dir   = config["review"]["isoform_eqtl"]["tensorqtl_dir"],
             fdr_threshold   = config["review"]["isoform_eqtl"].get("fdr_threshold", 0.05),
             nominal_threshold = config["review"]["isoform_eqtl"].get("nominal_threshold", 0.05)
    log:     config["review"]["isoform_eqtl"]["log"]
    message: "Comparing snRNA-seq eGenes against bulk transcript-level eQTL (O'Brien et al.)"
    script:  config["review"]["isoform_eqtl"]["script"]

#rule prep_snrnaseq_refs:
#    input:   cameron_rds        = config["review"]["prep_snrnaseq_refs"]["cameron_rds"],
#             polioudakis_counts = config["review"]["prep_snrnaseq_refs"]["polioudakis_counts"],
#             polioudakis_meta   = config["review"]["prep_snrnaseq_refs"]["polioudakis_meta"],
#             nowakowski_expr    = config["review"]["prep_snrnaseq_refs"]["nowakowski_expr"],
#             nowakowski_meta    = config["review"]["prep_snrnaseq_refs"]["nowakowski_meta"]
#    output:  cameron      = config["review"]["prep_snrnaseq_refs"]["cameron_out"],
#             polioudakis  = config["review"]["prep_snrnaseq_refs"]["polioudakis_out"],
#             nowakowski   = config["review"]["prep_snrnaseq_refs"]["nowakowski_out"]
#    singularity: config["containers"]["R"]  # your existing R container
#    resources: threads = 4,
#               mem_mb  = 32000,
#               time    = "1:00:00"
#    log:     config["review"]["prep_snrnaseq_refs"]["log"]
#    message: "Extracting reference datasets for SingleR classification"
#    script:  config["review"]["prep_snrnaseq_refs"]["script"]

#rule prep_singler_query_full:
#    input:
#        query_h5ad      = config["review"]["prep_singler_query"]["query_h5ad"],
#        ref_cameron     = config["review"]["prep_snrnaseq_refs"]["cameron"],
#        ref_polioudakis = config["review"]["prep_snrnaseq_refs"]["polioudakis"],
#        ref_nowakowski  = config["review"]["prep_snrnaseq_refs"]["nowakowski"]
#    output:
#        h5ad_cameron     = config["review"]["prep_singler_query_full"]["h5ad_cameron"],
#        h5ad_polioudakis = config["review"]["prep_singler_query_full"]["h5ad_polioudakis"],
#        h5ad_nowakowski  = config["review"]["prep_singler_query_full"]["h5ad_nowakowski"]
#    conda:
#        config["singler"]["py_env"]
#    resources:
#        threads = 4,
#        mem_mb  = 128000,
#        time    = "2:00:00"
#    log:
#        config["review"]["prep_singler_query_full"]["log"]
#    message:
#        "Preparing gene-subsetted query h5ads for snRNA-seq SingleR classification"
#    script:
#        config["review"]["prep_singler_query_full"]["script"]

#rule singler_snrnaseq:
#    input:
#        query_cameron     = config["review"]["prep_singler_query_full"]["h5ad_cameron"],
#        query_polioudakis = config["review"]["prep_singler_query_full"]["h5ad_polioudakis"],
#        query_nowakowski  = config["review"]["prep_singler_query_full"]["h5ad_nowakowski"],
#        ref_cameron       = config["review"]["prep_snrnaseq_refs"]["cameron"],
#        ref_polioudakis   = config["review"]["prep_snrnaseq_refs"]["polioudakis"],
#        ref_nowakowski    = config["review"]["prep_snrnaseq_refs"]["nowakowski"]
#    output:
#        sce_cameron     = config["review"]["singler_snrnaseq"]["sce_cameron"],
#        sce_polioudakis = config["review"]["singler_snrnaseq"]["sce_polioudakis"],
#        sce_nowakowski  = config["review"]["singler_snrnaseq"]["sce_nowakowski"]
#    singularity:
#        config["containers"]["singler"]
#    resources:
#        threads = 8,
#        mem_mb  = 200000,
#        time    = "6:00:00"
#    threads: 8
#    log:
#        config["review"]["singler_snrnaseq"]["log"]
#    message:
#        "Running SingleR classification against three snRNA-seq references"
#    script:
#        config["review"]["singler_snrnaseq"]["script"]


rule review_mashr:
# May need to add an input flag or this would run automatically 
#    input:   done = config["review"]["mashr"]["input_sentinel"]
#    input:   done = config["review"]["mashr"]["input_sentinel"]
    output:  config["review"]["mashr"]["sentinel"]
    params:  in_dir      = config["review"]["mashr"]["in_dir"],
             nom_dir     = config["review"]["mashr"]["nom_dir"],
             fugita_dir  = config["review"]["mashr"]["fugita_dir"],
             out_dir     = config["review"]["mashr"]["out_dir"],
             lfsr_thresh = config["review"]["mashr"]["lfsr_thresh"],
             n_random    = config["review"]["mashr"]["n_random"]
    singularity: config["containers"]["mashr"]
    resources: threads = 1,
               mem_mb  = 64000,
               time    = "4:00:00"
    threads: 1
    log:     config["review"]["mashr"]["log"]
    message: "Running matched mash models on prenatal vs adult neuronal eQTL"
    script:  config["review"]["mashr"]["script"]

#rule review_mashr_report:
#    input:   classified_inh = config["review"]["mashr"]["classified_inh"]
#    output:  config["review"]["mashr_report"]["html"]
#    params:  pi1_dir   = config["review"]["mashr_report"]["pi1_dir"],
#             mashr_dir = config["review"]["mashr_report"]["mashr_dir"],
#             rmd       = config["review"]["mashr_report"]["rmd"]
#    singularity: config["containers"]["mashr"]
#    resources: threads = 1, mem_mb  = 16000, time = "0:30:00"
#    threads: 1
#    log:     config["review"]["mashr_report"]["log"]
#    message: "Rendering mash and pi1 eQTL sharing report"
#    shell:
#        """
#        Rscript -e "
#          rmarkdown::render(
#            input       = '{params.rmd}',
#            output_file = '../{output}',
#            params      = list(
#              pi1_dir   = '{params.pi1_dir}',
#              mashr_dir = '{params.mashr_dir}'
#            ),
#            quiet = FALSE
#          )
#        " > {log} 2>&1
#        """
