rule scanpy_qc:
    input:   plate = "../results/02PARSE/combine_{plate}/all-sample/DGE_filtered/anndata.h5ad",
             nb = "scripts/scanpy_qc.ipynb"
    output:  "../results/03SCANPY/scanpy_qc_{plate}.html"
    conda:   "../envs/eqtl_study.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    message: "Running Scanpy in Jupyter notebook and producing HTML output"
    log:     notebook="../results/00LOG/03SCANPY/scanpy_qc_{plate}.ipynb"
    notebook:"{input.nb}"
