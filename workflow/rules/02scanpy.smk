rule scanpy_qc:
    input:   plate = "../results/02PARSE/combine_{plate}/all-sample/DGE_filtered/anndata.h5ad",
             nb = "scripts/scanpy_qc.ipynb"
    output:  "../results/03SCANPY/scanpy_qc_{plate}.html"
    conda:   "../envs/eqtl_study.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = "../results/03SCANPY/scanpy_qc_{plate}_pm.ipynb",
             plate = "{plate}",
             html_out = "scanpy_qc_{plate}.html"
    message: "Running Scanpy in Jupyter notebook and producing HTML output"
    log:     "../results/00LOG/03SCANPY/scanpy_qc_{plate}.log"
    shell:
             "papermill {input.nb} {params.nb_out} -p plate {params.plate} >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_clustering:
    # Note that plate is on requried to run initialize_env: need to add functionality to omit this here
    input:   html = expand("../results/03SCANPY/scanpy_qc_{plate}.html", plate = config['plate']),
             nb = "scripts/scanpy_clustering.ipynb"
    output:  "../results/03SCANPY/scanpy_clustering.html"
    conda:   "../envs/eqtl_study.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = "../results/03SCANPY/scanpy_clustering_pm.ipynb",
             html_out = "scanpy_clustering.html"
    message: "Running Scanpy in Jupyter notebook and producing HTML output"
    log:     "../results/00LOG/03SCANPY/scanpy_clustering.log"
    shell:
             "papermill {input.nb} {params.nb_out} -p plate clustering >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"
