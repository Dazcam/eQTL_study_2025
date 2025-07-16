localrules: ldsr_strat_summary

rule get_hg38_refs:
    output: baseline = expand(config["slsdr"]["get_hg38_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(config["slsdr"]["get_hg38_refs"]["weights"], chr=range(1, 23)),
            frq = expand(config["slsdr"]["get_hg38_refs"]["frq"], chr=range(1, 23)),
            bim = expand(config["slsdr"]["get_hg38_refs"]["bim"], chr=range(1, 23))
#            snps = config["slsdr"]["get_hg38_refs"]["snps"]
    params: config["slsdr"]["get_hg38_refs"]["params"] 
    log: expand(config["slsdr"]["get_hg38_refs"]["log"], chr=range(1, 23))
    shell: "scripts/get_ldsr_hg38_refs.sh {params} > {log} 2>&1"

rule get_hg19_refs:
    output: baseline = expand(config["slsdr"]["get_hg19_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(config["slsdr"]["get_hg19_refs"]["weights"], chr=range(1, 23)),
            frq = expand(config["slsdr"]["get_hg19_refs"]["frq"], chr=range(1, 23)),
            bim = expand(config["slsdr"]["get_hg19_refs"]["bim"], chr=range(1, 23))
#            snps = config["slsdr"]["get_hg19_refs"]["snps"]
    params: config["slsdr"]["get_hg19_refs"]["params"]
    log: expand(config["slsdr"]["get_hg19_refs"]["log"], chr=range(1, 23))
    shell: "scripts/get_ldsr_hg19_refs.sh {params} > {log} 2>&1"

rule make_annot_maxCPP:
    input:  cred_set = rules.sort_susie.output,
            bim = config["slsdr"]["get_hg38_refs"]["bim"] 
    output: cred_set = config["slsdr"]["make_annot_maxCPP"]["output"]
    log:    config["slsdr"]["make_annot_maxCPP"]["log"] 
    singularity: config["containers"]["R"]
    script: "../scripts/make_annot_maxCPP.R"

rule lift_hapmap3_snps:
   # Lift hapmap3 SNPs to hg38
   input:  bim_hg38 = lambda wildcards: expand(rules.get_hg38_refs.output.bim, chr=range(1,23)),
           bim_hg19 = lambda wildcards: expand(rules.get_hg19_refs.output.bim, chr=range(1,23))
   output: config["slsdr"]["lift_hapmap3_snps"]["output"],
   singularity: config["containers"]["ubuntu"]   
   log: config["slsdr"]["lift_hapmap3_snps"]["log"]
   shell: "/scratch/c.c1477909/eQTL_study_2025/workflow/scripts/liftover_hapmap3_snps.sh {output} > {log} 2>&1"

rule ldsr_ld_scores_hg38:
    input:   annot = rules.make_annot_maxCPP.output,
             bfile = rules.get_hg38_refs.output.bim,
             snps = rules.lift_hapmap3_snps.output
    output:  config["slsdr"]["ldsr_ld_scores_hg38"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    params:  bfile = config["slsdr"]["ldsr_ld_scores_hg38"]["bfile"],
             ldscores = config["slsdr"]["ldsr_ld_scores_hg38"]["ldscores"]
    message: "Generating LD scores on hg38 for {wildcards.cell_type}, chr {wildcards.chr}" 
    log:     config["slsdr"]["ldsr_ld_scores_hg38"]["log"]
    shell:
             """
             eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
             conda activate ldsr
             python ../resources/ldsr/ldsc/ldsc.py --l2 \
               --bfile {params.bfile} --ld-wind-cm 1 \
               --annot {input.annot} \
               --out {params.ldscores} \
               --print-snps {input.snps} 2> {log} 2>&1
             """

rule ldsr_strat_hg38_bl_v12:
    input:   gwas = rules.munge_sumstats.output,
             ldsr = expand(rules.ldsr_ld_scores_hg38.output, cell_type = config["cell_types"], chr = range(1,23))
    output:  config["slsdr"]["ldsr_strat_hg38_bl_v12"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    params:  weights = config["slsdr"]["ldsr_strat_hg38_bl_v12"]["weights"],
             baseline = config["slsdr"]["ldsr_strat_hg38_bl_v12"]["baseline"],
             frq = config["slsdr"]["ldsr_strat_hg38_bl_v12"]["frq"],
             ldscores = config["slsdr"]["ldsr_strat_hg38_bl_v12"]["ldscores"],
             out_prefix = config["slsdr"]["ldsr_strat_hg38_bl_v12"]["out_prefix"]
    message: "Running stratified LDSR with {wildcards.cell_type} using hg38 refs, baseline 1.2 and {wildcards.gwas} GWAS"
    log:     config["slsdr"]["ldsr_strat_hg38_bl_v12"]["log"]
    shell:
             """
             eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
             conda activate ldsr
             echo "Starting Snakemake script..." > {log}
             which python >> {log}

             python ../resources/ldsr/ldsc/ldsc.py --h2 {input.gwas} \
               --w-ld-chr {params.weights} \
               --ref-ld-chr {params.baseline},{params.ldscores} \
               --overlap-annot \
               --frqfile-chr {params.frq} \
               --out {params.out_prefix} \
               --print-coefficients >> {log} 2>&1
              """

rule ldsr_strat_summary:
    input:  expand(rules.ldsr_strat_hg38_bl_v12.output, cell_type = config["cell_types"], gwas = config["gwas"])
    output: config["slsdr"]["ldsr_strat_summary"]["output"]
    log:    config["slsdr"]["ldsr_strat_summary"]["log"]
    shell:  """
            head -1 {input[0]} > {output}
            grep L2_1 {input}  | cut -f 7- -d/ >> {output}
            """
