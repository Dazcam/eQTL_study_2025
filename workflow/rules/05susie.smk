localrules: prep_susie_gene_meta

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
