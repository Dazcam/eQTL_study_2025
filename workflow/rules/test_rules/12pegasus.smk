rule pegasus_qc:
    input:   plate = "../results/02PARSE/combine_{plate}/all-sample/DGE_filtered/anndata.h5ad",
             nb = "scripts/pegasus_testing.ipynb"
    output:  "../results/04PEGASUS/pegasus_qc_{plate}.html"
    conda:   "../envs/pegasus.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = "../results/04SCANPY/pegasus_{plate}_pm.ipynb",
             plate = "{plate}",
             html_out = "pegasus_{plate}.html"
    message: "Running Pegasus QC in Jupyter notebook and producing HTML output"
    log:     "../results/00LOG/04SCANPY/pegasus_qc_{plate}.log"
    shell:
             "papermill {input.nb} {params.nb_out} -p plate {params.plate} >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1" 
