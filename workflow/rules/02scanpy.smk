rule scanpy_qc:
    input:   plate = "../results/02PARSE/combine_{plate}/all-sample/DGE_filtered/anndata.h5ad",
             nb = "scripts/scanpy_qc.ipynb"
    output:  "../results/03SCANPY/scanpy_qc_{plate}.html"
    conda:   "../envs/eqtl_study.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = "../results/03SCANPY/scanpy_qc_{plate}_pm.ipynb",
             plate = "{plate}",
             pm_path = "/scratch/c.c1477909/.conda/envs/snakemake/bin/papermill"
    message: "Running Scanpy in Jupyter notebook and producing HTML output"
    log:     "../results/00LOG/03SCANPY/scanpy_qc_{plate}.log"
    shell:
             "echo 'Running on: {params.pm_path}' >> {log} && "
             "echo 'Using PATH: $PATH' >> {log} && "
             "{params.pm_path} {input.nb} {params.nb_out} -p param1 {params.plate} >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {output} >> {log} 2>&1"
