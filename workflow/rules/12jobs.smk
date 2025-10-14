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
