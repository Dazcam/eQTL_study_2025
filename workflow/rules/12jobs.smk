rule run_jobs:
    input:  bulk_file = config['jobs']['run_jobs']['bulk_file']
    output: config['jobs']['run_jobs']['output']
    params: sc_dir = config['jobs']['run_jobs']['sc_dir'],
            out_dir = config['jobs']['run_jobs']['out_dir']
    resources: threads = 20, mem_mb = 200000, time="24:00:00"
    singularity: config["containers"]["r_eqtl"]
    message: "Running JOBS to boost sc-eQTL signal in TensorQTL nominal output"
    benchmark: "reports/benchmarks/12jobs.run_jobs.txt"    
    log:    config['jobs']['run_jobs']['log']
    script: "../scripts/jobs_eqtl_boost.R"  

rule merge_jobs:
    input:  tensor = config["smr"]["cat_tensorqtl_nom_snps"]["cat_snps"],
            jobs = config['jobs']['run_jobs']['output']
    output: nom_p = config['jobs']['merge_jobs']['nom_p'],
            fdr = config['jobs']['merge_jobs']['fdr'], # For Rmd reporting
    params: jobs_dir = config['jobs']['merge_jobs']['jobs_dir']
    resources: threads = 20, mem_mb = 200000, time="24:00:00"
    singularity: config["containers"]["r_eqtl"]
    message: "Merge JOBS boosted sc-eQTL signal with original TensorQTL nominal output"
    benchmark: "reports/benchmarks/12jobs.merge_jobs_{cell_type}.txt"
    log:    config['jobs']['merge_jobs']['log']
    script: "../scripts/jobs_merge_boosted_eQTL.R"


# Note that I've kept the smr rule names the same as that in 10smr.smk for now
# As still in testing, need to comment that out of Snakefile to avoid rule clashes
rule create_query:
    input:  eqtl = config['jobs']['merge_jobs']['nom_p'],
#            eqtl = config["output_files"]["tensorqtl_perm_output"], # All probes failed HEIDI when using perm
            snps = config['smr']['get_snp_positions']['snp_pos'],
            genes = config["susie"]["prep_susie_gene_meta"]["output"],
            frq = config['smr']['cat_refs']['cat_frq'],
            r_script = config["smr"]["create_query"]["r_script"]
    output: query = config["jobs"]["create_query"]["query"],
            gene_lst = config["jobs"]["create_query"]["gene_lst"]
#    params: qval_thresh = config["eqtl_fdr"]
    singularity: config["containers"]["R"]
    resources: threads = 4, mem_mb = 20000
    message: "Create query file for SMR"
    benchmark: "reports/benchmarks/12jobs.create_query_{cell_type}.txt"
    log: config["jobs"]["create_query"]["log"]
    script: "../scripts/prep_eQTL_for_smr.R"

rule create_besd:
    input: query = config["jobs"]["create_query"]["query"],
           smr = config['smr']['get_smr_binary']['smr'],
    output: config["jobs"]["create_besd"]["besd"]
    params: config["jobs"]["create_besd"]["prefix"]
    resources: threads = 4, mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    message: "Create besd file for SMR"
    benchmark: "reports/benchmarks/12jobs.create_besd_{cell_type}.txt"
    log: config["jobs"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --make-besd --out {params} > {log} 2>&1
        """

rule smr:
    input:  bin = config['smr']['get_smr_binary']['smr'],
            gwas = config['smr']['format_gwas']['ma'],
            besd = config["jobs"]["create_besd"]["besd"]
    output: config["jobs"]["smr"]["smr"]
    params: geno_prefix = config["smr"]["smr"]["geno_prefix"],
            besd_prefix = config["jobs"]["create_besd"]["prefix"],
            out_prefix = config["jobs"]["smr"]["smr_prefix"]
    resources: threads = 4, mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    message: "Run SMR"
    benchmark: "reports/benchmarks/12jobs.smr_{cell_type}_{gwas}.txt"
    log:    config["jobs"]["smr"]["log"]
    shell:  """
            {input.bin} --bfile {params.geno_prefix} \
                --gwas-summary {input.gwas} \
                --beqtl-summary {params.besd_prefix} \
                --out {params.out_prefix} \
                --peqtl-smr 0.01 >> {log} 2>&1 
            """

rule jobs_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  jobs = expand(config['jobs']['smr']['smr'], cell_type = config['cell_types'], gwas = config['gwas']),
            rmd_script = "scripts/jobs_report.Rmd"
    output: "reports/12JOBS/jobs_report.html"
    params: in_dir = "../../results/12JOBS/mrgd/",
            bmark_dir = "../reports/benchmarks/",
            tensor_dir = "../../results/05TENSORQTL/tensorqtl_perm/",
            log_dir = "../../results/00LOG/12JOBS/",
            smr_dir = "../../results/12JOBS/smr/",
            p_smr = config["p_smr"],
            p_heidi =  config["p_heidi"],    
            output_file = "../reports/12JOBS/jobs_report.html"
    resources: threads = 8, mem_mb = 80000, time="24:00:00"
    singularity: config["containers"]["r_eqtl"]
    message: "Generate JOBS report"
    benchmark: "reports/benchmarks/12jobs.jobs_report.benchmark.txt"
    log:     "../results/00LOG/12JOBS/jobs_report.log"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', \
            bmark_dir = '{params.bmark_dir}', \
            log_dir = '{params.log_dir}', \
            smr_dir = '{params.smr_dir}', \
            p_smr = '{params.p_smr}', \
            p_heidi = '{params.p_heidi}', \
            tensor_dir = '{params.tensor_dir}'))" > {log} 2>&1
        """
