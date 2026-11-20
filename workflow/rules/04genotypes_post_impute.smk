##  Settings  ------------------------------------------------------------------------------
GPO         = config["geno_post_impute"]
CHROMOSOMES = list(range(1, 23)) + ["X"]


##  Rules  ---------------------------------------------------------------------------------

rule impute_check:
    input:  GPO['impute_check']['input']
    output: tsv = GPO['impute_check']['tsv'],
            vcf = GPO['impute_check']['vcf']
    params: dir = GPO['impute_check']['in_dir'],
            pwd = GPO['impute_check']['pwd']
    envmodules: config['envmodules']['bcftools']
    message: "Evaluate imputation quality for each chromosome"
    benchmark: GPO['impute_check']['benchmark']
    log:    GPO['impute_check']['log']
    shell:
        """
        exec > {log} 2>&1

        # Unzip the TOPMed download for this chromosome
        if [ ! -f "{output.vcf}.tbi" ]; then
            pwd=$(cat {params.pwd})
            unzip -P "$pwd" -n {input} -d {params.dir}
        fi

        # Count total SNPs from dose.vcf.gz
        snp_count=$(bcftools view -H {output.vcf} | wc -l)

        # Create index for vcf file
        tabix -p vcf {output.vcf}

        # Count MAF < 0.05 and >= 0.05 from info file
        maf_counts=$(gzip -d -c {params.dir}chr{wildcards.chr}.info.gz | grep -v "^#" | awk -F'\t' 'BEGIN {{common=0; rare=0}} {{split($8,a,";"); for(i in a) if(a[i]~/^MAF=/) {{maf=substr(a[i],5); if(maf>=0.05) common++; else rare++}}}} END {{print rare " " common}}')
        maf_lt_05=$(echo $maf_counts | cut -d' ' -f1)
        maf_gteq_05=$(echo $maf_counts | cut -d' ' -f2)

        # Count R² < 0.8 and >= 0.8 from info file
        rsq_counts=$(gzip -d -c {params.dir}chr{wildcards.chr}.info.gz | grep -v "^#" | awk -F'\t' 'BEGIN {{low=0; high=0}} {{split($8,a,";"); for(i in a) if(a[i]~/^R2=/) {{rsq=substr(a[i],4); if(rsq<0.8) low++; else high++}}}} END {{print low " " high}}')
        rsq_lt_08=$(echo $rsq_counts | cut -d' ' -f1)
        rsq_gteq_08=$(echo $rsq_counts | cut -d' ' -f2)

        # Count MAF >= 0.05 AND R² >= 0.8
        maf_rsq_count=$(gzip -d -c {params.dir}chr{wildcards.chr}.info.gz | grep -v "^#" | awk -F'\t' 'BEGIN {{count=0}} {{split($8,a,";"); maf=0; rsq=0; for(i in a) {{if(a[i]~/^MAF=/) maf=substr(a[i],5); if(a[i]~/^R2=/) rsq=substr(a[i],4)}} if(maf>=0.05 && rsq>=0.8) count++}} END {{print count}}')

        # Write to summary table
        echo -e "chr\tsnp_count\tmaf_lt_0.05\tmaf_gteq_0.05\trsq_lt_0.08\trsq_gteq_0.08\tmaf_gteq_0.05_rsq_gteq_0.08" > {output.tsv}
        echo -e "{wildcards.chr}\t$snp_count\t$maf_lt_05\t$maf_gteq_05\t$rsq_lt_08\t$rsq_gteq_08\t$maf_rsq_count" >> {output.tsv}
        """

rule impute_check_cat:
    input:  expand(rules.impute_check.output.tsv, chr=CHROMOSOMES)
    output: GPO['impute_check_cat']['output']
    message: "Cat impute check counts for each chr into single file"
    benchmark: GPO['impute_check_cat']['benchmark']
    log:    GPO['impute_check_cat']['log']
    shell:
        """
        exec > {log} 2>&1
        # Write header
        echo -e "chr\tsnp_count\tmaf_lt_0.05\tmaf_gteq_0.05\trsq_lt_0.08\trsq_gteq_0.08\tmaf_gteq_0.05_rsq_gteq_0.08" > {output}

        # Combine per-chromosome summaries
        cat {input} | grep -v "^chr" >> {output}

        # Add overall summary
        awk 'NR>1 {{snp+=$2; mlt+=$3; mgteq+=$4; rlt+=$5; rgteq+=$6; mr+=$7}} END {{print "Overall\t" snp "\t" mlt "\t" mgteq "\t" rlt "\t" rgteq "\t" mr}}' {output} >> {output}
        """

rule dwnld_dbsnp_ref:
    output: GPO['dwnld_dbsnp_ref']['output']
    params: web_link = GPO['dwnld_dbsnp_ref']['web_link']
    message: "Download dbSNP reference"
    benchmark: GPO['dwnld_dbsnp_ref']['benchmark']
    log:    GPO['dwnld_dbsnp_ref']['log']
    shell:
        """
        exec > {log} 2>&1
        wget {params.web_link} -O {output}
        wget {params.web_link}.tbi -O {output}.tbi
        """

rule add_rsID:
    input:  vcf = rules.impute_check.output.vcf,
            dbsnp = rules.dwnld_dbsnp_ref.output
    output: GPO['add_rsID']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Add dbSNP rsIDs to each imputed chr-specific VCF file"
    benchmark: GPO['add_rsID']['benchmark']
    log:    GPO['add_rsID']['log']
    shell:
        """
        bcftools annotate -a {input.dbsnp} -c ID {input.vcf} -O z -o {output} > {log} 2>&1
        """

rule vcf_cat:
    input:  expand(rules.add_rsID.output, chr=CHROMOSOMES)
    output: GPO['vcf_cat']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Cat imputed chr-specific VCF files into single file"
    benchmark: GPO['vcf_cat']['benchmark']
    log:    GPO['vcf_cat']['log']
    shell:
        """
        bcftools concat {input} -O z -o {output} > {log} 2>&1
        """

# Compute HWE and filter the concatenated VCF; Keep SNPs passing MAF, R2 and HWE thresholds
rule filter_tags:
    input:  rules.vcf_cat.output
    output: GPO['filter_tags']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Compute HWE and filt concated VCF; Keep SNPs passing MAF, R^2 and HWE threshs"
    benchmark: GPO['filter_tags']['benchmark']
    params: hwe = GPO['filter_tags']['hwe'],
            maf = GPO['filter_tags']['maf'],
            rsq = GPO['filter_tags']['rsq']
    log:    GPO['filter_tags']['log']
    shell:
        """
        exec > {log} 2>&1
        bcftools +fill-tags {input} -O z -- -t HWE |\
        bcftools view -e 'HWE<{params.hwe} |\
                          INFO/MAF<{params.maf} |\
                          INFO/R2<{params.rsq}' -O z -o {output}
        """

rule check_VCF:
    input:  vcf = rules.filter_tags.output,
            ref_gz = GPO['check_VCF']['ref_gz']
    output: log_out = GPO['check_VCF']['log_out'],
            list = GPO['check_VCF']['list'],
            ref = temp(GPO['check_VCF']['ref']),
            fai = temp(GPO['check_VCF']['fai'])
    params: prefix = GPO['check_VCF']['prefix_out'],
            script = GPO['check_VCF']['script']
    envmodules: config['envmodules']['samtools']
    message: "Validate the filt VCF against ref genome and ID problematic SNPs"
    benchmark: GPO['check_VCF']['benchmark']
    log:    GPO['check_VCF']['log']
    shell:
        """
        exec > {log} 2>&1
        gunzip -c {input.ref_gz} > {output.ref}
        samtools faidx {output.ref}
        python {params.script} -r {output.ref} -v {input.vcf} -o {params.prefix} --exclude {output.list}
        """

rule exclude_SNPs:
    input:  vcf = rules.filter_tags.output,
            list = rules.check_VCF.output.list
    output: GPO['exclude_SNPs']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Exclude SNPs IDed in CheckVCF"
    benchmark: GPO['exclude_SNPs']['benchmark']
    log:    GPO['exclude_SNPs']['log']
    shell:
        """
        bcftools view -e 'ID=@{input.list}' {input.vcf} -O z -o {output} > {log} 2>&1
        """

rule idx_vcf:
    input:  rules.exclude_SNPs.output
    output: GPO['idx_vcf']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Index the final filtered VCF"
    benchmark: GPO['idx_vcf']['benchmark']
    log:    GPO['idx_vcf']['log']
    shell:
        """
        tabix -p vcf {input} > {log} 2>&1
        """

rule get_sample_list:
    input:  vcf = rules.exclude_SNPs.output,
            idx = rules.idx_vcf.output
    output: GPO['get_sample_list']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Extract final sample list"
    benchmark: GPO['get_sample_list']['benchmark']
    log:    GPO['get_sample_list']['log']
    shell:  "bcftools query -l {input.vcf} > {output} 2> {log}"

rule vcf_to_plink:
    input:  vcf = rules.exclude_SNPs.output,
            idx = rules.idx_vcf.output
    output: GPO['vcf_to_plink']['output']
    params: prefix = GPO['vcf_to_plink']['params']
    envmodules: config['envmodules']['plink']
    message: "Convert genotypes VCF to plink format"
    benchmark: GPO['vcf_to_plink']['benchmark']
    log:    GPO['vcf_to_plink']['log']
    shell:  "plink --vcf {input.vcf} --double-id --make-bed --out {params.prefix} > {log} 2>&1"

rule get_ld_pruned_snps:
    input:  rules.vcf_to_plink.output
    output: GPO['get_ld_pruned_snps']['output']
    params: input_prefix = GPO['vcf_to_plink']['params'],
            output_prefix = GPO['get_ld_pruned_snps']['params']
    envmodules: config['envmodules']['plink']
    message: "LD prune SNPs before running PCA on genotypes"
    benchmark: GPO['get_ld_pruned_snps']['benchmark']
    log:    GPO['get_ld_pruned_snps']['log']
    shell:
        """
        plink --bfile {params.input_prefix} \
              --indep-pairwise 250 5 0.2 \
              --out {params.output_prefix} > {log} 2>&1
        """

rule prune_genotypes:
    input:  bfile = rules.vcf_to_plink.output,
            included = rules.get_ld_pruned_snps.output
    output: GPO['prune_genotypes']['output']
    params: input_prefix = GPO['vcf_to_plink']['params'],
            output_prefix = GPO['prune_genotypes']['params']
    envmodules: config['envmodules']['plink']
    message: "Create genotypes plink-format file containing only pruned SNPs"
    benchmark: GPO['prune_genotypes']['benchmark']
    log:    GPO['prune_genotypes']['log']
    shell:
        """
        plink --bfile {params.input_prefix} \
              --extract {input.included} \
              --make-bed \
              --out {params.output_prefix} > {log} 2>&1
        """

rule calc_genotype_pcs:
    input:  rules.prune_genotypes.output
    output: GPO['calc_genotype_pcs']['output']
    params: input_prefix = GPO['prune_genotypes']['params'],
            output_prefix = GPO['calc_genotype_pcs']['params'],
            pcs = GPO['calc_genotype_pcs']['pcs']
    envmodules: config['envmodules']['plink']
    message: "Run PCA on pruned genotypes to get genotype-specific PC covariates"
    benchmark: GPO['calc_genotype_pcs']['benchmark']
    log:    GPO['calc_genotype_pcs']['log']
    shell:
        """
        plink --bfile {params.input_prefix} \
              --pca {params.pcs} \
              --out {params.output_prefix} > {log} 2>&1
        """

rule geno_post_report:
    input:  impute_summary = rules.impute_check_cat.output,
            eigenvec = rules.calc_genotype_pcs.output,
            rmd_script = GPO["geno_post_report"]["script"]
    output: GPO["geno_post_report"]["html"]
    params: eigenval = GPO["calc_genotype_pcs"]["params"] + ".eigenval",
            chromosomes = ",".join(str(c) for c in CHROMOSOMES),
            hwe = GPO["filter_tags"]["hwe"],
            maf = GPO["filter_tags"]["maf"],
            rsq = GPO["filter_tags"]["rsq"],
            pcs = GPO["calc_genotype_pcs"]["pcs"],
            bmark_dir = GPO["geno_post_report"]["bmark_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate genotyping post-imputation report"
    benchmark: GPO["geno_post_report"]["benchmark"]
    log:    GPO["geno_post_report"]["log"]
    shell:
        """
        Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd_script}', \
            output_file = file.path(wd, '{output}'), \
            knit_root_dir = wd, \
            params = list(impute_summary = '{input.impute_summary}', \
                          eigenval = '{params.eigenval}', \
                          chromosomes = '{params.chromosomes}', \
                          hwe = '{params.hwe}', \
                          maf = '{params.maf}', \
                          rsq = '{params.rsq}', \
                          pcs = '{params.pcs}', \
                          bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["04genotypes_post_impute"] = [GPO["geno_post_report"]["html"]]
