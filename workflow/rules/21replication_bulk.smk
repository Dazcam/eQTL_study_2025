configfile: "../config/config.yaml"

def get_all_egene_perm_files():
    qtl_dir = config["replication_bulk"]["extract_unique_egenes"]["qtl_dir"]
    geno_pc = config["tensorQTL"]["geno_pcs"]
    norm_method = config["tensorQTL"]["norm_methods"][0]
    return [
        f"{qtl_dir}{ct}_{norm_method}_genPC_{geno_pc}_expPC_{config['exp_pc_map'][ct]}/"
        f"{ct}_{norm_method}_perm.cis_qtl.txt.gz"
        for ct in config["cell_types"]
    ]

rule all:
    input:
         config["replication_bulk"]["report"]["html"]

rule extract_unique_egenes:
    input:  perm_files = get_all_egene_perm_files()
    output: config["replication_bulk"]["extract_unique_egenes"]["out_file"]
    params: cell_types = config["cell_types"],
            geno_pc = config["tensorQTL"]["geno_pcs"],
            norm_method = config["tensorQTL"]["norm_methods"][0]
    singularity: config["containers"]["r_eqtl"]
    message: "Extract unique FDR-sig eGenes across all L1+L2+pseudotime cell types"
    benchmark: "reports/benchmarks/replication_bulk.extract_unique_egenes.txt"
    log:    config["replication_bulk"]["extract_unique_egenes"]["log"]
    script: "../scripts/repl_bulk_extract_unique_egenes.R"

rule overlap_obrien:
    input:  fetal_egenes = rules.extract_unique_egenes.output[0],
            obrien_supp  = config["replication_bulk"]["overlap_obrien"]["obrien_supp"]
    output: config["replication_bulk"]["overlap_obrien"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Overlap unique fetal eGenes against O'Brien 2018 bulk gene- and transcript-level significant eQTL"
    benchmark: "reports/benchmarks/replication_bulk.overlap_obrien.txt"
    log:    config["replication_bulk"]["overlap_obrien"]["log"]
    script: "../scripts/repl_bulk_overlap_obrien.R"

rule isoform_eqtl_comparison:
    input:   fetal_egenes = rules.extract_unique_egenes.output[0],
             gtf          = config["parse"]["get_refs"]["gtf_out"],
             bulk_eqtl    = config["replication_bulk"]["isoform_eqtl"]["bulk_eqtl"],
             gene_lookup  = config["dev_specificity"]["report"]["gene_lookup"]
    output:  config["replication_bulk"]["isoform_eqtl"]["summary"]
    params:  nominal_threshold = config["replication_bulk"]["isoform_eqtl"].get("nominal_threshold", 0.05)
    conda:   config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 48000, time = "3:00:00"
    threads: 8
    log:     config["replication_bulk"]["isoform_eqtl"]["log"]
    benchmark: "reports/benchmarks/replication_bulk.isoform_eqtl_comparison.txt"
    message: "Comparing snRNA-seq eGenes against bulk transcript-level eQTL (O'Brien et al.)"
    script:  config["replication_bulk"]["isoform_eqtl"]["script"]


rule repl_bulk_report:
    input:  egenes_per_celltype = rules.extract_unique_egenes.output[0],
            overlap             = rules.overlap_obrien.output[0],
            isoform_summary     = rules.isoform_eqtl_comparison.output[0],
            rmd_script          = config["replication_bulk"]["report"]["script"]
    output: config["replication_bulk"]["report"]["html"]
    params: egenes_per_celltype = lambda wc, input: "../" + input.egenes_per_celltype,
            overlap             = lambda wc, input: "../" + input.overlap,
            isoform_summary     = lambda wc, input: "../" + input.isoform_summary,
            bmark_dir           = "../" + config["replication_bulk"]["report"]["bmark_dir"],
            output_file         = config["replication_bulk"]["report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate replication_bulk report"
    benchmark: "reports/benchmarks/replication_bulk.report.txt"
    log: config["replication_bulk"]["report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(egenes_per_celltype = '{params.egenes_per_celltype}', \
                          overlap = '{params.overlap}', \
                          isoform_summary = '{params.isoform_summary}', \
                          bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """
