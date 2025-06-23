localrules: prep_susie_gene_meta, get_gwas_sumstats, ldsr_strat_summary

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
    script: "../scripts/prep_susie_input.py"

#    shell: """
#           cat {input.exp_bed} | cut -f4- | sed 's/TargetID/phenotype_id/g' > {output.exp_bed}
#           cat  {input.covar} | sed '1s/^\t/SampleID\t/' > {output.covar}
#           (echo "sample_id"; cat {input.exp_bed} | head -1 | cut -f5- | tr '\\t' '\\n') > {output.samp_lst}
#           """

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
    envmodules: "htslib/1.9"    
    output: config["slsdr"]["merge_susie"]["output"]
    log:    config["slsdr"]["merge_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Merging SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" > {log}
            awk 'NR==1 || FNR>1{{print}}' {input} | bgzip -c > {output}
            echo "Completed merge of SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" >> {log}
            """
 
rule sort_susie:
    input:  config["slsdr"]["merge_susie"]["output"].replace("{susie_suffix}", "cred.hp.txt")
    output: config["slsdr"]["sort_susie"]["output"]
    log:    config["slsdr"]["sort_susie"]["log"]
    envmodules: "htslib/1.9"
    shell:
        """
        set -euo pipefail
        echo "Sorting SuSiE cred.hp.txt for {wildcards.cell_type}" > {log}
        gunzip -c {input} > {wildcards.cell_type}_temp.txt
        (head -n 1 {wildcards.cell_type}_temp.txt && tail -n +2 {wildcards.cell_type}_temp.txt | sort -k3 -k4n) | bgzip > {output}
        rm {wildcards.cell_type}_temp.txt
        echo "Completed sorting SuSiE cred.hp.txt for {wildcards.cell_type}" >> {log}
        """

#rule sort_susie:
#    input: 
#        expand(config["slsdr"]["merge_susie"]["output"], 
#            cell_type=config["cell_types"], 
#            susie_suffix=["cred.hp.txt", "cred.txt", "snp.txt"])
#    output: config["slsdr"]["sort_susie"]["output"],
#    envmodules: "htslib/1.9"
#    log:    config["slsdr"]["sort_susie"]["log"]
#    shell:
#            """
#            set -euo pipefail
#            echo "Sorting SuSiE files for {wildcards.cell_type}" > {log}
#            gunzip -c {input[0]} > {wildcards.cell_type}_temp.txt
#            (head -n 1 {wildcards.cell_type}_temp.txt && tail -n +2 {wildcards.cell_type}_temp.txt | sort -k3 -k4n) | bgzip > {output}
#            rm {wildcards.cell_type}_temp.txt
#            echo "Completed sorting SuSiE files for {wildcards.cell_type}" >> {log}
#            """

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

rule get_gwas_sumstats:
    # Download GWAS sumsatts files
    output:  config["slsdr"]["get_gwas_sumstats"]["output"]
    params:  lambda wildcards: config['gwas'][wildcards.gwas]
    message: "Download {wildcards.gwas} sumstats file"
    log:     config["slsdr"]["get_gwas_sumstats"]["log"]
    run:

             if wildcards.gwas in ("mdd", "neuroticism"):

                 shell("""

                 cp {params} temp; gunzip -c temp > {output}; rm temp 2> {log}

                 """)
             
             elif wildcards.gwas in ("scz", "bd"):

                 shell("""

                 wget -O - {params} | gunzip -c | sed '/##/d' > {output} 

                 """)

             else:
                 
                 shell("""

                 wget -O - {params} | gunzip -c > {output} 

                 """)


rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    # python convert available here: https://github.com/precimed/python_convert/tree/master
    input:   rules.get_gwas_sumstats.output
    output:  config["slsdr"]["standardise_sumstats"]["output"]
    message: "Standardising {input}"
    params: config["slsdr"]["standardise_sumstats"]["temp"]
    log:    config["slsdr"]["standardise_sumstats"]["log"] 
    run:

             if wildcards.gwas in ("bd", "scz"):

                 shell("""

                 cat {input} | sed 's/ID/SNP/g' | sed 's/#CHROM/CHR/g' > {params};
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}

                  """)

             elif "height" in wildcards.gwas:

                 shell("""
                 cat {input} | sed 's/Tested_Allele/A1/g' | sed 's/Other_Allele/A2/g' > {params};
    
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
       	           --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}
             
                  """)

             elif "neuroticism" in wildcards.gwas:

                 shell("""
 
                 cat {input} | sed 's/POS/BP/g' | sed 's/RSID_UKB/SNP/g' | sed 's/REF/A1/g' | sed 's/ALT/A2/g' > {params};

                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}

                 """)

             else:

                 shell("""

                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {input} \
                   --out {output} --force --auto --head 5 \
                   --log {log}

                  """)

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking
    input:   rules.standardise_sumstats.output
    output:  config["slsdr"]["add_z_score"]["output"]
    message: "Adding Z score to {input} if required"
    log:     config["slsdr"]["add_z_score"]["log"]
    shell:
             """
             python ../resources/python_convert/sumstats.py zscore \
             --sumstats {input} \
             --out {output} --force \
             --log {log} \
             --a1-inc
             """

rule add_N:
    # N to GWAS sumstats lacking 
    input:   rules.add_z_score.output
    output:  config["slsdr"]["add_N"]["output"]
    message: "Adding N to {input} if required"
    log:     config["slsdr"]["add_N"]["log"]
    run:

             if "scz" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{{{s=(NR==1)?"N":"130644"; $0=$0 OFS s}}}}1' {input} > {output} 2> {log}

                 """)

             elif "adhd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"225534";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)

             elif "asd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"46350";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)
             
             elif "bd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"413466";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             elif "neuroticism" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"313467";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             else:

                 shell("cp {input} {output}")

rule make_gwas_bed_hg19:
    # Generate gwas hg19 bed file for lift over
    input:   rules.add_N.output
    output:  config["slsdr"]["make_gwas_bed_hg19"]["output"]
    message: "Create bed input file for {input} hg19 to hg38 LiftOver"
    log: config["slsdr"]["make_gwas_bed_hg19"]["log"]
    run:

        shell("""
        
        awk 'BEGIN {{OFS="\t"}} NR > 1 {{print "chr"$2, $3-1, $3, $1}}' {input} > {output} 2> {log}

        """)

rule liftover_gwas_to_hg38:
    # Liftover sumstats from hg19 to hg38
    input: bed = rules.make_gwas_bed_hg19.output,
           chain = config["slsdr"]["liftover_gwas_to_hg38"]["chain"]
    output:  config["slsdr"]["liftover_gwas_to_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params: config["slsdr"]["liftover_gwas_to_hg38"]["unlifted"]
    message: "Lifting over {input} from hg19 to hg38"
    log: config["slsdr"]["liftover_gwas_to_hg38"]["log"]
    shell:
        """

        ../resources/liftover/liftOver {input.bed} {input.chain} {output} {params} 2> {log}
    
        """

rule add_hg38_coords_to_gwas:
    #Restore summary statistics file with hg38 coords
    input: lifted = rules.liftover_gwas_to_hg38.output,
           sumstats = rules.add_N.output
    output: config["slsdr"]["add_hg38_coords_to_gwas"]["output"]
    message: "Adding hg38 coords to {input.sumstats} sumstats"
    log: config["slsdr"]["add_hg38_coords_to_gwas"]["log"]
    script: "../scripts/add_hg38_coords_to_gwas.py"


rule munge_sumstats:
    # Format sumstats for LDSR input
    input:   snps = "../results/05SLDSR/sldsr/1000G.EUR.hg38.w_hm3_test.snplist",
             gwas = rules.add_hg38_coords_to_gwas.output
    output:  config["slsdr"]["munge_sumstats"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    params:  config["slsdr"]["munge_sumstats"]["prefix"]
    log:     config["slsdr"]["munge_sumstats"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        echo "Starting Snakemake script..." > {log}
        which python >> {log}
        python ../resources/ldsr/ldsc/munge_sumstats.py  --sumstats {input.gwas} \
          --merge-alleles {input.snps} \
          --out {params} \
          --a1-inc \
          --p PVAL >> {log} 2>&1
        echo "Finished at $(date)" >> {log}

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
