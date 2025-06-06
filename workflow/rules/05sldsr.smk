localrules: prep_susie_gene_meta

# This rul is not required for susie but will be useful elsewhere (maybe qtl smk)
#rule get_sig_eGenes:
#    input:  config["output_files"]["tensorqtl_perm_output"]
#    output: config["slsdr"]["get_sig_eGenes"]["output"]
#    singularity: config["containers"]["R"]
#    log:    config["slsdr"]["get_sig_eGenes"]["log"]
#    script: "../scripts/get_sig_egenes_with_cis_window.R"

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
           # Extract the static header columns
           printf "CHROM\\tPOS\\tREF\\tALT\\t" > header_static.tsv

           # Extract sample IDs and paste onto the header line
           bcftools query -l {input} | paste -sd '\\t' - >> header_static.tsv

           # Transpose into a single header line
           tr '\\n' '\\t' < header_static.tsv | sed 's/\\t$//' | gzip > header_row.tsv.gz

           # Extract dosage matrix
           bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n" {input} | gzip > dose_matrix.tsv.gz

           # Combine header and dosage, compress, and index
           zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > {output.dosage}
           tabix -s1 -b2 -e2 -S1 {output.dosage}
           rm header_static.tsv header_row.tsv.gz dose_matrix.tsv.gz
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

rule run_susie:
    input:  exp_mat = rules.prep_susie_input.output.exp_bed,
            gene_meta = rules.prep_susie_gene_meta.output,
            samp_lst = rules.prep_susie_input.output.samp_lst,
            eqt_to_test = rules.prep_susie_input.output.exp_bed,
            covar = rules.prep_susie_input.output.covar,
            geno_mat  = rules.vcf_to_dosage.output.dosage,
            geno_idx   = rules.vcf_to_dosage.output.idx
    output: cs_variant = config["slsdr"]["run_susie"]["cs_output"],
            lbf_variable = config["slsdr"]["run_susie"]["lbf_output"],
            full_susie = config["slsdr"]["run_susie"]["full_output"]
    params: chunk           = lambda wildcards: f"{wildcards.batch_index} {config['n_batches']}",
            out_prefix      = lambda wildcards: f"{wildcards.qtl_subset}.{wildcards.batch_index}_{config['n_batches']}",
            cis_window      = config["cis_window"],
            write_full      = config["write_full_susie"]
    singularity: config["containers"]["susie"]
    log:    config["slsdr"]["run_susie"]["log"]
    shell:  """
            Rscript bin/run_susie.R \
              --expression_matrix {input.exp_mat} \
              --phenotype_meta {input.gene_meta} \
              --sample_meta {input.smpl_lst} \
              --phenotype_list {input.eqt_to_test} \
              --covariates {input.covar} \
              --genotype_matrix {input.geno_mat} \
              --chunk '{params.chunk}' \
              --cisdistance {params.cis_window} \
              --out_prefix '{params.out_prefix}' \
              --eqtlutils null \
              --write_full_susie {params.write_full}
            """

#rule prep_susie_input:
#    input:  sig_eGenes = config["slsdr"]["get_sig_eGenes"]["output"],
#            pseudobulk = config["input_files"]["counts"],
#            genotypes = config["output_files"]["genotypes_plink"],
#            covariates = config["input_files"]["covariates"] 
#    output: config["slsdr"]["prep_susie_input"]["output"]
#    envmodules: "plink/2.0", "R/3.5.1"
#    resources: threads = 8, mem_mb = 64000, time="1:00:00"
#    message: "Generating input files for SuSIE for all sig eGenes"
#    log:    config["slsdr"]["prep_susie_input"]["log"]
#    script: "../scripts/prep_susie_input_files.py"
