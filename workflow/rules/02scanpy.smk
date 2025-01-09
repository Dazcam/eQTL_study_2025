rule scanpy_qc:
    input:   plate = "../results/02PARSE/combine_{plate}/all-sample/DGE_filtered/anndata.h5ad",
             nb = "scripts/scanpy_qc.ipynb"
    output:  "../results/03SCANPY/scanpy_qc_{plate}.html"
    conda:   "../envs/eqtl_study.yml"
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    message: "Running Scanpy in Jupyter notebook and producing HTML output"
    log:     "../results/00LOG/03SCANPY/scanpy_qc_{plate}.log"
    shell:
        r"""
        jupyter nbconvert --to html --execute {input.nb} --output {output} \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.extra_arguments='--plate {wildcards.plate}' > {log} 2>&1
        """

