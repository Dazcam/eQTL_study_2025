configfile: "../config/config.yaml" 

rule all:
    input: 
        config["singler"]["singler_report"]["html"]

rule extract_reference_data:
    input:
        cameron_rds        = config["singler"]["extract_reference_data"]["cameron_rds"],
        polioudakis_counts = config["singler"]["extract_reference_data"]["polioudakis_counts"],
        polioudakis_meta   = config["singler"]["extract_reference_data"]["polioudakis_meta"],
        nowakowski_expr    = config["singler"]["extract_reference_data"]["nowakowski_expr"],
        nowakowski_meta    = config["singler"]["extract_reference_data"]["nowakowski_meta"]
    output:
        cameron     = config["singler"]["extract_reference_data"]["cameron"],
        polioudakis = config["singler"]["extract_reference_data"]["polioudakis"],
        nowakowski  = config["singler"]["extract_reference_data"]["nowakowski"]
    singularity:
        config["containers"]["R"]
    resources:
        threads = 4, mem_mb = 32000, time = "1:00:00"
    log:
        config["singler"]["extract_reference_data"]["log"]
    message:
        "Extracting raw counts and metadata from source reference datasets"
    script:
        config["singler"]["extract_reference_data"]["script"]


rule prep_snrnaseq_refs:
    input:
        cameron     = config["singler"]["extract_reference_data"]["cameron"],
        polioudakis = config["singler"]["extract_reference_data"]["polioudakis"],
        nowakowski  = config["singler"]["extract_reference_data"]["nowakowski"]
    output:
        cameron     = config["singler"]["prep_snrnaseq_refs"]["cameron"],
        polioudakis = config["singler"]["prep_snrnaseq_refs"]["polioudakis"],
        nowakowski  = config["singler"]["prep_snrnaseq_refs"]["nowakowski"],
        cameron_genes    = config["singler"]["prep_snrnaseq_refs"]["cameron_genes"],
        polioudakis_genes = config["singler"]["prep_snrnaseq_refs"]["polioudakis_genes"],
        nowakowski_genes  = config["singler"]["prep_snrnaseq_refs"]["nowakowski_genes"]
    singularity:
        config["containers"]["singler"]
    resources:
        threads = 4, mem_mb = 64000, time = "2:00:00"
    threads: 4
    params:
        n_cells_per_cluster = 500,
        random_seed         = 1234
    log:
        config["singler"]["prep_snrnaseq_refs"]["log"]
    message:
        "Building normalised SCE reference objects for SingleR"
    script:
        config["singler"]["prep_snrnaseq_refs"]["script"]


rule prep_singler_query_full:
    input:
        query_h5ad            = config["singler"]["prep_singler_query"]["query_h5ad"],
        ref_cameron_genes     = config["singler"]["prep_snrnaseq_refs"]["cameron_genes"],
        ref_polioudakis_genes = config["singler"]["prep_snrnaseq_refs"]["polioudakis_genes"],
        ref_nowakowski_genes  = config["singler"]["prep_snrnaseq_refs"]["nowakowski_genes"]
    output:
        sentinel = config["singler"]["prep_singler_query_full"]["sentinel"]
    conda:
        config["scanpy"]["env"]
    resources:
        threads = 10, mem_mb = 500000, time = "2:00:00"
    threads: 10
    log:
        config["singler"]["prep_singler_query_full"]["log"]
    message:
        "Preparing gene-subsetted query MTX files for snRNA-seq SingleR classification"
    script:
        config["singler"]["prep_singler_query_full"]["script"]


rule singler_snrnaseq:
    input:
        sentinel    = config["singler"]["prep_singler_query_full"]["sentinel"],
        ref_cameron     = config["singler"]["prep_snrnaseq_refs"]["cameron"],
        ref_polioudakis = config["singler"]["prep_snrnaseq_refs"]["polioudakis"],
        ref_nowakowski  = config["singler"]["prep_snrnaseq_refs"]["nowakowski"]
    output:
        sentinel = config["singler"]["singler_snrnaseq"]["sentinel"]
    singularity:
        config["containers"]["singler"]
    resources:
        threads = 8, mem_mb = 500000, time = "4:00:00"
    threads: 8
    log:
        config["singler"]["singler_snrnaseq"]["log"]
    message:
        "Running SingleR classification against three snRNA-seq references"
    script:
        config["singler"]["singler_snrnaseq"]["script"]

rule singler_report:
    input:  sentinel        = config["singler"]["singler_snrnaseq"]["sentinel"],
            umap_csv        = config["singler"]["prep_singler_query"]["umap_csv"],
            rmd_script      = "scripts/review_singler_report.Rmd"
    output: config["singler"]["singler_report"]["html"]
    params:
        cameron_csv     = "../../results/15REVIEW/singler_snrnaseq_cameron.csv",
        polioudakis_csv = "../../results/15REVIEW/singler_snrnaseq_polioudakis.csv",
        nowakowski_csv  = "../../results/15REVIEW/singler_snrnaseq_nowakowski.csv",
        umap_csv        = "../../results/15REVIEW/umap_coords.csv",
        output_file     = "../reports/15REVIEW/singler_report.html"
    resources:
        threads = 4,
        mem_mb  = 32000,
        time    = "1:00:00"
    singularity:
        config["containers"]["singler"]
    log:
        config["singler"]["singler_report"]["log"]
    message:
        "Generating SingleR snRNA-seq reference correspondence report"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(
                cameron_csv     = '{params.cameron_csv}',
                polioudakis_csv = '{params.polioudakis_csv}',
                nowakowski_csv  = '{params.nowakowski_csv}',
                umap_csv        = '{params.umap_csv}'))" > {log} 2>&1
        """
