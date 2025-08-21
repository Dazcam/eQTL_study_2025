rule scanpy_qc:
    input:   plate = config["scanpy"]["qc"]["input"]["plate"],
             nb = config["scanpy"]["qc"]["input"]["nb"]
    output:  config["scanpy"]["qc"]["output"]
    conda:   config["scanpy"]["env"]
    resources: threads = 16, mem_mb = 160000, time = "0-3:00:00"
    params:  nb_out = config["scanpy"]["qc"]["nb_out"],
             html_out = config["scanpy"]["qc"]["html_out"]
    benchmark: config["scanpy"]["qc"]["benchmark"]
    message: "Running Scanpy QC in Jupyter notebook and producing HTML output"
    log:     config["scanpy"]["qc"]["log"]
    shell:
        """
        papermill {input.nb} {params.nb_out} -p plate {wildcards.plate} >> {log} 2>&1 && 
        jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1
        """

rule scanpy_clustering:
    # Note that plate is on requried to run initialize_env: need to add functionality to omit this here
    input:   html = expand(config["scanpy"]["qc"]["output"], plate = config['plate']),
             nb = config["scanpy"]["clustering"]["input"]
    output:  config["scanpy"]["clustering"]["output"]
    conda:   config["scanpy"]["env"]
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
#    resources: threads = 10, mem_mb = 100000, time="0-3:00:00"
    params:  nb_out = config["scanpy"]["clustering"]["nb_out"],
             html_out = config["scanpy"]["clustering"]["html_out"]
    benchmark: config["scanpy"]["clustering"]["benchmark"]
    message: "Running Scanpy clustering in Jupyter notebook and producing HTML output"
    log:     config["scanpy"]["clustering"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate clustering >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_annotation:
    # Note that plate is on requried to run initialize_env: need to add functionality to omit this here
    input:   html = config["scanpy"]["clustering"]["output"],
             nb = config["scanpy"]["annotation"]["input"]
    output:  config["scanpy"]["annotation"]["output"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time="3-0:00:00"
    params:  nb_out = config["scanpy"]["annotation"]["nb_out"],
             html_out = config["scanpy"]["annotation"]["html_out"]
    message: "Running Scanpy annotation in Jupyter notebook and producing HTML output"
    log:     config["scanpy"]["annotation"]["log"]
    shell:   
             "papermill {input.nb} {params.nb_out} -p plate annotation >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_cell_label_and_pseudobulk:
    input:   html = config["scanpy"]["clustering"]["output"],
             nb = config["scanpy"]["pseudobulk"]["input"]
    output:  config["scanpy"]["pseudobulk"]["output"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 100000, time="3-0:00:00"
    params:  nb_out = config["scanpy"]["pseudobulk"]["nb_out"],
             html_out = config["scanpy"]["pseudobulk"]["html_out"]
    message: "Running Scanpy cell label and pseudobulk in Jupyter notebook and producing HTML output"
    log:     config["scanpy"]["pseudobulk"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate cell_label_and_pseudobulk >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_extra:
    # Rule to run quick minor tests
    input:   html = config["scanpy"]["clustering"]["output"],
             nb = config["scanpy"]["extra"]["input"]
    output:  config["scanpy"]["extra"]["output"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 100000, time="3-0:00:00"
    params:  nb_out = config["scanpy"]["extra"]["nb_out"],
             html_out = config["scanpy"]["extra"]["html_out"]
    message: "Running Scanpy extra in Jupyter notebook and producing HTML output"
    log:     config["scanpy"]["extra"]["log"]
    shell:   
             "papermill {input.nb} {params.nb_out} -p plate extra >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"
