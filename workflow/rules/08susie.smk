localrules: prep_susie_gene_meta

rule all:
    input:
        config["susie"]["susie_report"]["html"]


def susie_exp_bed(w):
    norm_method = config['tensorQTL']['norm_methods'][0]
    return f"../results/05TENSORQTL/prep_input/{w.cell_type}_{norm_method}.bed"


def susie_covar(w):
    norm_method = config['tensorQTL']['norm_methods'][0]
    geno_pc     = config['tensorQTL']['geno_pcs']
    exp_pc      = config['exp_pc_map'][w.cell_type]
    return (f"../results/05TENSORQTL/prep_input/"
            f"{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}_split_covariates.txt")


def susie_perm(w):
    norm_method = config['tensorQTL']['norm_methods'][0]
    geno_pc     = config['tensorQTL']['geno_pcs']
    exp_pc      = config['exp_pc_map'][w.cell_type]
    return (f"../results/05TENSORQTL/tensorqtl_perm/"
            f"{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}/"
            f"{w.cell_type}_{norm_method}_perm.cis_qtl.txt.gz")


rule prep_susie_gene_meta:
    input:  bed = susie_exp_bed
    output: config["susie"]["prep_susie_gene_meta"]["output"]
    singularity: config["containers"]["R"]
    message:     "Create gene metadata file for {wildcards.cell_type} for SuSiE"
    benchmark:   "reports/benchmarks/08susie.prep_susie_gene_meta_{cell_type}.benchmark.txt"
    log:         config["susie"]["prep_susie_gene_meta"]["log"]
    script: "../scripts/susie_prep_gene_meta.R"

rule vcf_to_dosage:
    input:  config["susie"]["vcf_to_dosage"]["vcf"]
    output: dosage = config["susie"]["vcf_to_dosage"]["out_dosage"],
            idx    = config["susie"]["vcf_to_dosage"]["out_idx"]
    envmodules: "BCFtools/1.23-GCC-14.3.0", "HTSlib/1.23-GCC-14.3.0"
    message:    "Generate an allele dosage file from the genotype VCF for SuSiE"
    benchmark:  "reports/benchmarks/08susie.vcf_to_dosage.benchmark.txt"
    log:        config["susie"]["vcf_to_dosage"]["log"]
    shell:  """
            set -euo pipefail
            exec 2> {log}

            # Header (CHROM POS REF ALT + sample IDs), then dosages with 'chr' removed
            (
              printf 'CHROM\tPOS\tREF\tALT'
              bcftools query -l {input} | awk '{{printf "\t%s", $0}}'
              printf '\n'
              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input} | sed 's/^chr//'
            ) | bgzip > {output.dosage}
            tabix -s1 -b2 -e2 -S1 {output.dosage}
            """

rule prep_susie_input:
    input:  exp_bed  = susie_exp_bed,
            covar    = susie_covar
    output: exp_bed  = config["susie"]["prep_susie_input"]["out_exp"],
            covar    = config["susie"]["prep_susie_input"]["out_covar"],
            samp_lst = config["susie"]["prep_susie_input"]["out_samp_lst"]
    message:   "Prepare expression, covariate and sample files for {wildcards.cell_type} for SuSiE"
    benchmark: "reports/benchmarks/08susie.prep_susie_input_{cell_type}.benchmark.txt"
    log:       config["susie"]["prep_susie_input"]["log"]
    script: "../scripts/susie_prep_input.py"

rule run_susie:
    input:  exp_mat     = rules.prep_susie_input.output.exp_bed,
            gene_meta   = rules.prep_susie_gene_meta.output,
            smpl_lst    = rules.prep_susie_input.output.samp_lst,
            eqt_to_test = susie_perm,
            covar       = rules.prep_susie_input.output.covar,
            geno_mat    = rules.vcf_to_dosage.output.dosage,
            geno_idx    = rules.vcf_to_dosage.output.idx
    output: cs_hp_out   = config["susie"]["run_susie"]["cs_hp_out"],
            cs_out      = config["susie"]["run_susie"]["cs_out"],
            snp_out     = config["susie"]["run_susie"]["snp_out"]
    params: chunk       = lambda w: f"{w.batch_index} {config['susie_batches']}",
            cis_window  = config["window"],
            fdr         = config["eqtl_fdr"]
    singularity: config["containers"]["susie"]
    message:     "Run SuSiE fine-mapping for {wildcards.cell_type}, batch {wildcards.batch_index}"
    benchmark:   "reports/benchmarks/08susie.run_susie_{cell_type}_{batch_index}_{susie_batches}.benchmark.txt"
    log:         config["susie"]["run_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Running run_susie for {wildcards.cell_type}, batch {wildcards.batch_index}" > {log}
            Rscript scripts/susie_run.R \
              --expression_matrix {input.exp_mat} \
              --phenotype_meta {input.gene_meta} \
              --sample_meta {input.smpl_lst} \
              --phenotype_list {input.eqt_to_test} \
              --fdr {params.fdr} \
              --covariates {input.covar} \
              --genotype_matrix {input.geno_mat} \
              --chunk '{params.chunk}' \
              --cisdistance {params.cis_window} \
              --out_cred_hp {output.cs_hp_out} \
              --out_cred {output.cs_out} \
              --out_snp {output.snp_out} \
              >> {log} 2>&1
            echo "Completed run_susie for {wildcards.cell_type}, batch {wildcards.batch_index}" >> {log}
            """

def susie_batch_files(w):
    return expand(rules.run_susie.output.cs_hp_out,
                  cell_type     = w.cell_type,
                  batch_index   = range(1, config["susie_batches"] + 1),
                  susie_batches = config["susie_batches"])


rule merge_susie:
    input:      susie_batch_files
    output:     config["susie"]["merge_susie"]["output"]
    envmodules: "HTSlib/1.21-GCC-13.3.0"
    message:    "Merge SuSiE credible sets into a single file for {wildcards.cell_type}"
    benchmark:  "reports/benchmarks/08susie.merge_susie_{cell_type}.benchmark.txt"
    log:        config["susie"]["merge_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Merging SuSiE cred.hp.txt for {wildcards.cell_type}" > {log}
            awk 'NR==1 || FNR>1{{print}}' {input} | bgzip -c > {output}
            echo "Completed merge of SuSiE cred.hp.txt for {wildcards.cell_type}" >> {log}
            """

rule sort_susie:
    input:      rules.merge_susie.output
    output:     config["susie"]["sort_susie"]["output"]
    envmodules: "HTSlib/1.21-GCC-13.3.0"
    message:    "Sort SuSiE credible sets by chromosome and position for {wildcards.cell_type}"
    benchmark:  "reports/benchmarks/08susie.sort_susie_{cell_type}.benchmark.txt"
    log:        config["susie"]["sort_susie"]["log"]
    shell:  """
            set -euo pipefail
            echo "Sorting SuSiE cred.hp.txt for {wildcards.cell_type}" > {log}
            gunzip -c {input} | (
              IFS= read -r header
              printf '%s\n' "$header"
              LC_ALL=C sort -t$'\t' -k3,3V -k4,4n
            ) | bgzip > {output}
            echo "Completed sorting SuSiE cred.hp.txt for {wildcards.cell_type}" >> {log}
            """

rule susie_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  susie_files     = expand(rules.sort_susie.output, cell_type = config["cell_types"]),
            gene_meta_files = expand(rules.prep_susie_gene_meta.output, cell_type = config["cell_types"]),
            gene_lookup     = "../resources/sheets/gene_lookup_hg38.tsv",
            rmd_script      = config["susie"]["susie_report"]["script"]
    output: config["susie"]["susie_report"]["html"]
    params: cell_types    = ",".join(config["cell_types"]),
            in_dir        = config["susie"]["susie_report"]["in_dir"],
            gene_meta_dir = config["susie"]["susie_report"]["gene_meta_dir"],
            gene_lookup   = config["susie"]["susie_report"]["gene_lookup"],
            bmark_dir     = config["susie"]["susie_report"]["bmark_dir"],
            cis_window    = config["window"],
            n_batches     = config["susie_batches"],
            eqtl_fdr      = config["eqtl_fdr"],
            output_file   = config["susie"]["susie_report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    message:     "Generate SuSiE report"
    benchmark:   "reports/benchmarks/08susie.susie_report.benchmark.txt"
    log:         config["susie"]["susie_report"]["log"]
    shell:  """
            Rscript -e "rmarkdown::render('{input.rmd_script}', \
                output_file = '{params.output_file}', \
                params = list(cell_types = '{params.cell_types}', \
                              in_dir = '{params.in_dir}', \
                              gene_meta_dir = '{params.gene_meta_dir}', \
                              gene_lookup = '{params.gene_lookup}', \
                              bmark_dir = '{params.bmark_dir}', \
                              cis_window = {params.cis_window}, \
                              n_batches = {params.n_batches}, \
                              eqtl_fdr = {params.eqtl_fdr}))" > {log} 2>&1
            """
