localrules: prep_susie_gene_meta, get_gwas_sumstats

# This rule is not required for susie but will be useful elsewhere (maybe qtl smk)
rule get_sig_eGenes:
    input:  config["output_files"]["tensorqtl_perm_output"]
    output: config["slsdr"]["get_sig_eGenes"]["output"]
    singularity: config["containers"]["R"]
    log:    config["slsdr"]["get_sig_eGenes"]["log"]
    script: "../scripts/get_sig_egenes.R"

rule prep_susie_gene_meta:
    input:  config["input_files"]["counts"]
    output: config["slsdr"]["prep_susie_gene_meta"]["output"]
    singularity: config["containers"]["R"]
    log:    config["slsdr"]["prep_susie_gene_meta"]["log"]
    script: "../scripts/prep_gene_meta_for_susie.R"

rule vcf_to_dosage:
    input:  config["input_files"]["genotypes"]
    output: dosage = config["slsdr"]["vcf_to_dosage"]["out_dosage"],
            idx = config["slsdr"]["vcf_to_dosage"]["out_idx"]
    envmodules: "bcftools","htslib"
    log:    config["slsdr"]["vcf_to_dosage"]["log"]
    shell: """
           # Extract sample IDs into a single tab-separated line
           bcftools query -l {input} | paste -sd '\t' - > header_samples.tsv

           # Create header with CHROM POS REF ALT followed by sample IDs
           echo -e "CHROM\tPOS\tREF\tALT\t$(cat header_samples.tsv)" > header.tsv

           # Extract dosage matrix, remove 'chr' prefix
           bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input} | sed 's/^chr//' > dose_matrix.tsv

           # Combine header and dosage matrix, compress, and index
           cat header.tsv dose_matrix.tsv | bgzip > {output.dosage}
           tabix -s1 -b2 -e2 -S1 {output.dosage}

           # Clean up
           rm header_samples.tsv header.tsv dose_matrix.tsv
           """

rule prep_susie_input: 
    input:  exp_bed = config["input_files"]["counts"],
            covar = config["input_files"]["covariates"],
    output: exp_bed = config["slsdr"]["prep_susie_input"]["out_exp"],
            covar = config["slsdr"]["prep_susie_input"]["out_covar"],
            samp_lst = config["slsdr"]["prep_susie_input"]["out_samp_lst"]
    log:    config["slsdr"]["prep_susie_input"]["log"]
    shell: """
           cat {input.exp_bed} | cut -f4- | sed 's/TargetID/phenotype_id/g' > {output.exp_bed}
           cat  {input.covar} | sed '1s/^\t/SampleID\t/' > {output.covar}
           (echo "sample_id"; cat {input.exp_bed} | head -1 | cut -f5- | tr '\\t' '\\n') > {output.samp_lst}
           """

#    input:  exp_mat = "testdata/GEUVADIS_cqn.tsv",
#            gene_meta = "testdata/GEUVADIS_phenotype_metadata.tsv",
#            smpl_lst = "testdata/GEUVADIS_sample_metadata.tsv",
#            eqt_to_test = "testdata/susie_debug/GEUVADIS_test_ge.permuted.tsv.gz",
#            covar = "testdata/susie_debug/GEUVADIS_test_ge.covariates.txt",
#            geno_mat  = "testdata/susie_debug/LCL.dose.tsv.gz",


rule run_susie:
    input:  exp_mat = rules.prep_susie_input.output.exp_bed,
            gene_meta = rules.prep_susie_gene_meta.output,
            smpl_lst = rules.prep_susie_input.output.samp_lst,
            eqt_to_test = config["output_files"]["tensorqtl_perm_output"],
            covar = rules.prep_susie_input.output.covar,
            geno_mat  = rules.vcf_to_dosage.output.dosage,
            geno_idx   = rules.vcf_to_dosage.output.idx
    output: cs_hp_out = config["slsdr"]["run_susie"]["cs_hp_out"],
            cs_out = config["slsdr"]["run_susie"]["cs_out"],
            snp_out = config["slsdr"]["run_susie"]["snp_out"]
    params: chunk = lambda wildcards: f"{wildcards.batch_index} {config['susie_batches']}",
            out_prefix = lambda wildcards: f"../results/05SLDSR/susie/{wildcards.cell_type}/{wildcards.cell_type}.{wildcards.batch_index}_{config['susie_batches']}.susie",
            cis_window = config["susie_window"],
            write_full = config["write_full_susie"]
    singularity: config["containers"]["susie"]
    log:    config["slsdr"]["run_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Running run_susie for {wildcards.cell_type}, batch {wildcards.batch_index}" > {log}
            Rscript scripts/run_susie.R \
              --expression_matrix {input.exp_mat} \
              --phenotype_meta {input.gene_meta} \
              --sample_meta {input.smpl_lst} \
              --phenotype_list {input.eqt_to_test} \
              --covariates {input.covar} \
              --genotype_matrix {input.geno_mat} \
              --chunk '{params.chunk}' \
              --cisdistance {params.cis_window} \
              --out_prefix '{params.out_prefix}' \
              --write_full_susie {params.write_full} \
              >> {log} 2>&1
            echo "Completed run_susie for {wildcards.cell_type}, batch {wildcards.batch_index}" >> {log} 
            """

rule merge_susie:
    input:
        lambda wildcards: expand(
            "../results/05SLDSR/susie/{cell_type}/{cell_type}.{batch_index}_{susie_batches}.susie.{susie_suffix}",
            cell_type=wildcards.cell_type,
            batch_index=range(1, config["susie_batches"] + 1),
            susie_batches=config["susie_batches"],
            susie_suffix=wildcards.susie_suffix
        )    
    output: config["slsdr"]["merge_susie"]["output"]
    log:    config["slsdr"]["merge_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Merging SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" > {log}
            awk 'NR==1 || FNR>1{{print}}' {input} | bgzip -c > {output}
            echo "Completed merge of SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" >> {log}
            """
 
rule sort_susie:
    input: 
        expand(config["slsdr"]["merge_susie"]["output"], 
            cell_type=config["cell_types"], 
            susie_suffix=["cred.hp.txt", "cred.txt", "snp.txt"])
    output: config["slsdr"]["sort_susie"]["output"],
    envmodules: "htslib/1.9"
    log:    config["slsdr"]["sort_susie"]["log"]
    shell:
            """
            set -euo pipefail
            echo "Sorting SuSiE files for {wildcards.cell_type}" > {log}
            gunzip -c {input[0]} > {wildcards.cell_type}_temp.txt
            (head -n 1 {wildcards.cell_type}_temp.txt && tail -n +2 {wildcards.cell_type}_temp.txt | sort -k3 -k4n) | bgzip > {output}
            rm {wildcards.cell_type}_temp.txt
            echo "Completed sorting SuSiE files for {wildcards.cell_type}" >> {log}
            """

##### TODO: 
##### There is an issue with the get refs rules as they are tracking
##### {chr} so running 22 instances instead of 1
##### Two options use a touch file and a pause before running 
##### make_annot_maxCPP which tracks 22 chrs or extract these rules
##### to a set up snakefile to download all prerequisite packages
##### containers and files
##### Solution: use expand in output and logs forces one instance of rule to be run with all chrs

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

localrules: get_gwas_sumstats

rule get_gwas_sumstats:
    output:  "../results/GWAS/{GWAS}_hg19_raw.tsv"
    params:  lambda wildcards: config['GWAS'][wildcards.GWAS]
    message: "Download sumstats file"
    log:     "../results/00LOG/get_and_munge_GWAS/{GWAS}_download_sumstats.log"
    run:

             if wildcards.GWAS in ("MDD", "NEUROTICISM"):

                 shell("""

                 cp {params} temp; gunzip -c temp > {output}; rm temp 2> {log}

                 """)
             
             elif wildcards.GWAS in ("SCZ_EUR_ONLY", "BPD"):

                 shell("""

                 wget -O - {params} | gunzip -c | sed '/##/d' > {output} 

                 """)

             else:
                 
                 shell("""

                 wget -O - {params} | gunzip -c > {output} 

                 """)

#rule ldsr_strat_hg38_bl_12:
#    input:   gwas = "../results/03SUMSTATS/{GWAS}_hg19_LDSR_ready.sumstats.gz",
#             ldsr = expand(ldsr_ld_scores_hg38.output, chr = range(1,23))
#    output:  config["slsdr"]["strat_hg38_bl_v12"]["output"]
#    conda:   "../envs/ldsr.yml"
#    params:  weights = config["slsdr"]["strat_hg38_bl_v12"]["weights"]
#             baseline = config["slsdr"]["strat_hg38_bl_v12"]["baseline"]
#             frq = config["slsdr"]["strat_hg38_bl_v12"]["frq"]
#             ldscores = config["slsdr"]["strat_hg38_bl_v12"]["ldscores"]
#             out_prefix = config["slsdr"]["strat_hg38_bl_v12"]["out_prefix"]
#    message: "Running stratified LDSR with {wildcards.cell_type} using hg38 refs, baseline 1.2 and {wildcards.GWAS} GWAS"
#    log:     config["slsdr"]["strat_hg38_bl_v12"]["log"]
#    shell:
#             """
#             python ../resources/ldsr/ldsc.py --h2 {input.gwas} \
#               --w-ld-chr {params.weights} \
#               --ref-ld-chr {params.baseline},{params.ldscores} \ 
#               --overlap-annot \
#               --frqfile-chr {params.frq} \
#               --out {params.out_prefix} \
#               --print-coefficients 2> {log} 2>&1"

