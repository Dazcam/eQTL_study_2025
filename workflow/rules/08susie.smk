configfile: "../config/config.yaml"
localrules: prep_susie_gene_meta

rule all:
    input:
       "reports/08SUSIE/08susie_report.html"

rule get_sig_eGenes:
    input:  lambda w, norm_method=config['tensorQTL']['norm_methods'][0],
                geno_pc=config['tensorQTL']['geno_pcs']:
                f"../results/05TENSORQTL/tensorqtl_perm/{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{config['exp_pc_map'][w.cell_type]}/{w.cell_type}_{norm_method}_perm.cis_qtl.txt.gz"
    output: config["susie"]["get_sig_eGenes"]["output"]
    singularity: config["containers"]["R"]
    message: "Pull out signifcant eGenes for {wildcards.cell_type} from TensorQTL perm output"
    benchmark: "reports/benchmarks/08susie.get_sig_eGenes_{cell_type}.txt",
    log: config["susie"]["get_sig_eGenes"]["log"]
    script: "../scripts/susie_get_sig_egenes.R"

rule prep_susie_gene_meta:
    input:  lambda w, norm_method=config['tensorQTL']['norm_methods'][0]:
                f"../results/05TENSORQTL/prep_input/{w.cell_type}_{norm_method}.bed"
#config["susie"]["prep_susie_gene_meta"]["pseudblk"]
    output: config["susie"]["prep_susie_gene_meta"]["output"]
    singularity: config["containers"]["R"]
    message: "Create gene metadata file for {wildcards.cell_type} for SuSiE"
    benchmark: "reports/benchmarks/08susie.prep_susie_gene_meta_{cell_type}.txt",
    log:    config["susie"]["prep_susie_gene_meta"]["log"]
    script: "../scripts/susie_prep_gene_meta.R"

rule vcf_to_dosage:
    input:  config["susie"]["vcf_to_dosage"]["vcf"]
    output: dosage = config["susie"]["vcf_to_dosage"]["out_dosage"],
            idx = config["susie"]["vcf_to_dosage"]["out_idx"]
    envmodules: "bcftools","htslib"
    message: "Generate a allele dosage file from genotypes VCF for SuSiE"
    benchmark: "reports/benchmarks/08susie.vcf_to_dosage.txt",
    log:    config["susie"]["vcf_to_dosage"]["log"]
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
    input:
        #sig_eGenes = lambda wc: f"../results/08SUSIE/susie_input/{wc.cell_type}_eGenes_fdr_0.05.tsv",
        exp_bed = lambda w, norm_method=config['tensorQTL']['norm_methods'][0]:
            f"../results/05TENSORQTL/prep_input/{w.cell_type}_{norm_method}.bed",
        #genotypes = config["geno_post_impute"]["vcf_to_plink"]["output"],
        covar = lambda w, norm_method=config['tensorQTL']['norm_methods'][0],
            geno_pc=config['tensorQTL']['geno_pcs']:
            f"../results/05TENSORQTL/prep_input/{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{config['exp_pc_map'][w.cell_type]}_split_covariates.txt"
    output: exp_bed = config["susie"]["prep_susie_input"]["out_exp"],
            covar = config["susie"]["prep_susie_input"]["out_covar"],
            samp_lst = config["susie"]["prep_susie_input"]["out_samp_lst"]
    message: "Extract genotype, expression, and covariate data for {wildcards.cell_type} for SuSiE fine-mapping"
    benchmark: "reports/benchmarks/08susie.prep_susie_input_{cell_type}.txt",
    log:    config["susie"]["prep_susie_input"]["log"]
    script: "../scripts/prep_susie_input.py"

rule run_susie:
    input:  exp_mat = rules.prep_susie_input.output.exp_bed,
            gene_meta = rules.prep_susie_gene_meta.output,
            smpl_lst = rules.prep_susie_input.output.samp_lst,
            eqt_to_test = lambda w, norm_method=config['tensorQTL']['norm_methods'][0],
                geno_pc=config['tensorQTL']['geno_pcs']:
                f"../results/05TENSORQTL/tensorqtl_perm/{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{config['exp_pc_map'][w.cell_type]}/{w.cell_type}_{norm_method}_perm.cis_qtl.txt.gz",
            covar = rules.prep_susie_input.output.covar,
            geno_mat  = rules.vcf_to_dosage.output.dosage,
            geno_idx   = rules.vcf_to_dosage.output.idx
    output: cs_hp_out = config["susie"]["run_susie"]["cs_hp_out"],
            cs_out = config["susie"]["run_susie"]["cs_out"],
            snp_out = config["susie"]["run_susie"]["snp_out"]
    params: chunk = lambda wildcards: f"{wildcards.batch_index} {config['susie_batches']}",
            out_prefix = lambda wildcards: f"../results/08SUSIE/susie/{wildcards.cell_type}/{wildcards.cell_type}.{wildcards.batch_index}_{config['susie_batches']}.susie",
            cis_window = config["susie_window"],
            write_full = config["write_full_susie"]
    singularity: config["containers"]["susie"]
    message: "Run SuSiE fine-mapping for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/08susie.run_susie_{cell_type}_{batch_index}_{susie_batches}.txt",    
    log:    config["susie"]["run_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Running run_susie for {wildcards.cell_type}, batch {wildcards.batch_index}" > {log}
            Rscript scripts/susie_run.R \
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
            "../results/08SUSIE/susie/{cell_type}/{cell_type}.{batch_index}_{susie_batches}.susie.{susie_suffix}",
            cell_type=wildcards.cell_type,
            batch_index=range(1, config["susie_batches"] + 1),
            susie_batches=config["susie_batches"],
            susie_suffix=wildcards.susie_suffix
        )
    envmodules: "htslib/1.9"
    output: config["susie"]["merge_susie"]["output"]
    message: "Merge SuSiE fine-mapping output into a single file for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/08susie.merge_susie_{cell_type}_{susie_suffix}.txt",
    log:    config["susie"]["merge_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Merging SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" > {log}
            awk 'NR==1 || FNR>1{{print}}' {input} | bgzip -c > {output}
            echo "Completed merge of SuSiE {wildcards.susie_suffix} for {wildcards.cell_type}" >> {log}
            """

rule sort_susie:
    input:  config["susie"]["merge_susie"]["output"].replace("{susie_suffix}", "cred.hp.txt")
    output: config["susie"]["sort_susie"]["output"]
    log:    config["susie"]["sort_susie"]["log"]
    message: "Sort SuSiE fine-mapping output for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/08susie.sort_susie_{cell_type}.txt",
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

rule susie_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  susie_files = expand(rules.sort_susie.output, cell_type = config["cell_types"]),
            rmd_script = "scripts/susie_report.Rmd"
    output: config["susie"]["susie_report"]["html"]
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            in_dir = config["susie"]["susie_report"]["in_dir"],
            output_file = "../reports/08SUSIE/08susie_report.html"
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SuSiE report"
    benchmark: "reports/benchmarks/08susie.susie_report.txt"
    log: config["susie"]["susie_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), in_dir = '{params.in_dir}'))" > {log} 2>&1
        """
