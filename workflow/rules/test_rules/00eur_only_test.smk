#### 

# This is a test script to restrict our genotype analysis and GeX data
# to EUR samples only. I'm running this in a seprate script rather than 
# Rejigging the existing pipeline, in case we decide not to go down this 
# route. To do this requires:

# 1. Resticting geneotypes to EUR only and recaclculating the geno PCs
# 2. Re


rule vcf_to_plink:
    input:  vcf = rules.exclude_SNPs.output,
            idx = rules.idx_vcf.output
    output: config['geno_post_impute']['vcf_to_plink']['output']
    params: config['geno_post_impute']['vcf_to_plink']['params']
    envmodules: "plink/1.9"
    message: "Convert genotypes VCF to plink format"
    benchmark: "reports/benchmarks/geno_post_impute.vcf_to_plink.txt"
    log:    config['geno_post_impute']['vcf_to_plink']['log']
    shell:  "plink --vcf {input.vcf} --double-id --make-bed --out {params} > {log} 2>&1"

rule get_ld_pruned_snps:
    input:  rules.vcf_to_plink.output
    output: config['geno_post_impute']['get_ld_pruned_snps']['output']
    params: input_prefix = config['geno_post_impute']['vcf_to_plink']['params'],
            output_prefix = config['geno_post_impute']['get_ld_pruned_snps']['params']
    envmodules: "plink/1.9"
    message: "LD prune SNPs before running PCA on genotypes"
    benchmark: "reports/benchmarks/geno_post_impute.get_ld_pruned_snps.txt"
    log:    config['geno_post_impute']['get_ld_pruned_snps']['log']
    shell:  """
            plink --bfile {params.input_prefix} \
                  --indep-pairwise 250 5 0.2 \
                  --out {params.output_prefix} > {log} 2>&1
            """

rule prune_genotypes:
    input:  bfile = rules.vcf_to_plink.output,
            included = rules.get_ld_pruned_snps.output
    output: config['geno_post_impute']['prune_genotypes']['output']
    params: input_prefix = config['geno_post_impute']['vcf_to_plink']['params'],
            output_prefix = config['geno_post_impute']['prune_genotypes']['params']
    message: "Create genotypes plink-format file containing only pruned SNPs"
    benchmark: "reports/benchmarks/geno_post_impute.prune_genotypes.txt"
    envmodules: "plink/1.9"
    log:    config['geno_post_impute']['prune_genotypes']['log']
    shell:  """
            plink --bfile {params.input_prefix} \
                  --extract {input.included} \
                  --make-bed \
                  --out {params.output_prefix} > {log} 2>&1
            """ 
