rule review_mashr:
# May need to add an input flag or this would run automatically 
#    input:   done = config["review"]["mashr"]["input_sentinel"]
#    input:   done = config["review"]["mashr"]["input_sentinel"]
    output:  config["review"]["mashr"]["sentinel"]
    params:  in_dir      = config["review"]["mashr"]["in_dir"],
             nom_dir     = config["review"]["mashr"]["nom_dir"],
             fugita_dir  = config["review"]["mashr"]["fugita_dir"],
             out_dir     = config["review"]["mashr"]["out_dir"],
             lfsr_thresh = config["review"]["mashr"]["lfsr_thresh"],
             n_random    = config["review"]["mashr"]["n_random"]
    singularity: config["containers"]["mashr"]
    resources: threads = 1,
               mem_mb  = 64000,
               time    = "4:00:00"
    threads: 1
    log:     config["review"]["mashr"]["log"]
    message: "Running matched mash models on prenatal vs adult neuronal eQTL"
    script:  config["review"]["mashr"]["script"]

#rule review_mashr_report:
#    input:   classified_inh = config["review"]["mashr"]["classified_inh"]
#    output:  config["review"]["mashr_report"]["html"]
#    params:  pi1_dir   = config["review"]["mashr_report"]["pi1_dir"],
#             mashr_dir = config["review"]["mashr_report"]["mashr_dir"],
#             rmd       = config["review"]["mashr_report"]["rmd"]
#    singularity: config["containers"]["mashr"]
#    resources: threads = 1, mem_mb  = 16000, time = "0:30:00"
#    threads: 1
#    log:     config["review"]["mashr_report"]["log"]
#    message: "Rendering mash and pi1 eQTL sharing report"
#    shell:
#        """
#        Rscript -e "
#          rmarkdown::render(
#            input       = '{params.rmd}',
#            output_file = '../{output}',
#            params      = list(
#              pi1_dir   = '{params.pi1_dir}',
#              mashr_dir = '{params.mashr_dir}'
#            ),
#            quiet = FALSE
#          )
#        " > {log} 2>&1
#        """
