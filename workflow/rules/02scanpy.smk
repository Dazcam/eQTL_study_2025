##  Settings  ------------------------------------------------------------------------------
SCP = config["scanpy"]


##  Rules  ---------------------------------------------------------------------------------

rule scanpy_qc:
    input:   plate = SCP["qc"]["input"],
             nb = SCP["qc"]["script"]
    output:  SCP["qc"]["output"]
    conda:   SCP["env"]
    resources: threads = 16, mem_mb = 160000, time = "3:00:00"
    params:  nb_out = SCP["qc"]["nb_out"],
             html_out = SCP["qc"]["html_out"]
    benchmark: SCP["qc"]["benchmark"]
    message: "Running Scanpy QC in Jupyter notebook and producing HTML output"
    log:     SCP["qc"]["log"]
    shell:
        """
        papermill {input.nb} {params.nb_out} -p plate {wildcards.plate} >> {log} 2>&1 && 
        jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1
        """

rule scanpy_clustering:
    input:   html = expand(SCP["qc"]["output"], plate = config['plate']),
             nb = SCP["clustering"]["script"]
    output:  SCP["clustering"]["output"]
    conda:   SCP["env"]
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = SCP["clustering"]["nb_out"],
             html_out = SCP["clustering"]["html_out"]
    benchmark: SCP["clustering"]["benchmark"]
    message: "Running Scanpy clustering in Jupyter notebook and producing HTML output"
    log:     SCP["clustering"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate clustering >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_clustering2:
    input:   html = SCP["clustering"]["output"],
             nb = SCP["clustering2"]["script"]
    output:  SCP["clustering2"]["output"]
    conda:   SCP["env"]
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = SCP["clustering2"]["nb_out"],
             html_out = SCP["clustering2"]["html_out"]
    message: "Running Scanpy clustering 2 in Jupyter notebook and producing HTML output"
    benchmark: SCP["clustering2"]["benchmark"]
    log:     SCP["clustering2"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate clustering2 >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_clustering3:
    input:   html = SCP["clustering2"]["output"],
             nb = SCP["clustering3"]["script"]
    output:  SCP["clustering3"]["output"]
    conda:   SCP["env"]
    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    params:  nb_out = SCP["clustering3"]["nb_out"],
             html_out = SCP["clustering3"]["html_out"]
    message: "Running Scanpy clustering 3 in Jupyter notebook and producing HTML output"
    benchmark: SCP["clustering3"]["benchmark"]
    log:     SCP["clustering3"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate clustering2 >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_annotation:
    input:   html = SCP["clustering3"]["output"],
             nb = SCP["annotation"]["script"]
    output:  SCP["annotation"]["output"]
    envmodules: "Miniforge3"
    conda:   SCP["env"]
    resources: threads = 10, mem_mb = 200000, time="1:00:00"
    params:  nb_out = SCP["annotation"]["nb_out"],
             html_out = SCP["annotation"]["html_out"]
    message: "Running Scanpy annotation in Jupyter notebook and producing HTML output"
    benchmark: SCP["annotation"]["benchmark"]
    log:     SCP["annotation"]["log"]
    shell:   
             """
             papermill {input.nb} {params.nb_out} -p plate annotation --kernel eqtl_study >> {log} 2>&1 && 
             jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1
             """

rule scanpy_pseudobulk:
    input:   html = SCP["clustering3"]["output"],
             nb = SCP["pseudobulk"]["script"]
    output:  SCP["pseudobulk"]["output"]
    conda:   SCP["env"]
    resources: threads = 10, mem_mb = 100000, time="1-0:00:00"
    params:  nb_out = SCP["pseudobulk"]["nb_out"],
             html_out = SCP["pseudobulk"]["html_out"]
    message: "Running Scanpy cell label and pseudobulk in Jupyter notebook and producing HTML output"
    benchmark: SCP["pseudobulk"]["benchmark"]
    log:     SCP["pseudobulk"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate pseudobulk >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1 "

rule scanpy_subclustering:
    input:   html = SCP["clustering3"]["output"],
             nb = SCP["subclustering"]["script"]
    output:  SCP["subclustering"]["output"]
    conda:   SCP["env"]
    resources: threads = 8, mem_mb = 80000, time="2:00:00"
    params:  nb_out = SCP["subclustering"]["nb_out"],
             html_out = SCP["subclustering"]["html_out"]
    benchmark: SCP["subclustering"]["benchmark"]
    message: "Running Scanpy subclustering in Jupyter notebook and producing HTML output"
    log:     SCP["subclustering"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} -p plate subclustering >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_extra:
    # Rule to run quick minor tests
    input:   html = SCP["clustering"]["output"],
             nb = SCP["extra"]["script"]
    output:  SCP["extra"]["output"]
    conda:   SCP["env"]
    resources: threads = 8, mem_mb = 80000, time="1:00:00"
    params:  nb_out = SCP["extra"]["nb_out"],
             html_out = SCP["extra"]["html_out"]
    message: "Running Scanpy extra in Jupyter notebook and producing HTML output"
    benchmark: SCP["extra"]["benchmark"]
    log:     SCP["extra"]["log"]
    shell:
             "papermill {input.nb} {params.nb_out} --execution-timeout -1 --request-save-on-cell-execute -p plate extra >> {log} 2>&1 && "
             "jupyter nbconvert --to html {params.nb_out} --output {params.html_out} >> {log} 2>&1"

rule scanpy_lisi:
    input:  h5ad = SCP["clusters_h5ad"]
    output: csv  = SCP["lisi"]["csv"],
            plot = SCP["lisi"]["plot"]
    params: batch_col     = SCP["lisi"]["batch_col"],
            cell_type_col = SCP["lisi"]["cell_type_col"],
            random_seed   = SCP["lisi"]["random_seed"],
            perplexity    = SCP["lisi"]["perplexity"]
    conda:  SCP["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    threads: 10
    message:   "Running LISI (iLISI / cLISI) batch correction assessment pre- and post-Harmony"
    benchmark: SCP["lisi"]["benchmark"]
    log:       SCP["lisi"]["log"]
    script:    SCP["lisi"]["script"]


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["02scanpy"] = [expand(rules.scanpy_qc.output, plate = config["plate"]),
                       rules.scanpy_annotation.output,
                       rules.scanpy_extra.output,
                       rules.scanpy_pseudobulk.output,
                       rules.scanpy_subclustering.output,
                       rules.scanpy_lisi.output]
