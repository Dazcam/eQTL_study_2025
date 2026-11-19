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
    log: GPO['impute_check']['log']
    shell:
        """
        # Create directory and unzip
        if [ ! -f "{output.vcf}.tbi" ]; then
            pwd=$(cat {params.pwd})
            unzip -P "$pwd" -n {input} -d {params.dir} > {log} 2>&1
        fi

        # Count total SNPs from dose.vcf.gz
        snp_count=$(bcftools view -H {params.dir}chr{wildcards.chr}.dose.vcf.gz | wc -l)
        
        # Create index for vcf file
        tabix -p vcf {params.dir}chr{wildcards.chr}.dose.vcf.gz >> {log} 2>&1              

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
    input:  expand(rules.impute_check.output.vcf, chr=CHROMOSOMES)
    output: GPO['impute_check_cat']['output']
    message: "Cat impute check counts for each chr into single file"
    benchmark: GPO['impute_check_cat']['benchmark']
    log:      GPO['impute_check_cat']['log']
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
    params: web_link = GPO['dwnld_dbsnp_ref']['web_link'],
            prefix = GPO['dwnld_dbsnp_ref']['prefix']
    message: "Download dbSNP reference"
    benchmark: GPO['dwnld_dbsnp_ref']['benchmark']
    log:    GPO['dwnld_dbsnp_ref']['log']
    shell:
        """
        wget {params.web_link} -O {params.prefix} > {log} 2>&1
        wget {params.web_link}.tbi -O {params.prefix}.tbi >> {log} 2>&1
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
        bcftools +fill-tags {input} -O z -- -t HWE |\
        bcftools view -e 'HWE<{params.hwe} |\
                          INFO/MAF<{params.maf} |\
                          INFO/R2<{params.rsq}' -O z -o {output} > {log} 2>&1
        """

rule check_VCF:
    input:  vcf = rules.filter_tags.output,
            ref_gz = GPO['check_VCF']['ref_gz']
    output: log_out = GPO['check_VCF']['log'],
            list = GPO['check_VCF']['list'],
            ref = temp(GPO['check_VCF']['ref']),
            fai = temp(GPO['check_VCF']['fai'])    
    params: prefix = GPO['check_VCF']['prefix_out'],
            script = GPO['check_VCF']['script']
    envmodules: config['envmodules']['samtools']
    message: "Validate the filt VCF against ref genome and ID problematic SNPs"
    benchmark: GPO['check_VCF']['benchmark']
    log:    GPO['check_VCF']['log']
    shell:  """
            gunzip -c {input.ref_gz} > {output.ref} 2>> {log}
            samtools faidx {output.ref} >> {log} 2>&1
            python {params.script} -r {output.ref} -v {input.vcf} -o {params.prefix} --exclude {output.list} > {log} 2>&1
            """

rule exclude_SNPs:
    input:  vcf = rules.filter_tags.output,
            list = rules.check_VCF.output.list
    output: GPO['exclude_SNPs']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Exclude SNPs IDed in CheckVCF"
    benchmark: GPO['exclude_SNPs']['benchmark']
    log:    GPO['exclude_SNPs']['log']
    shell:  """
            bcftools view -e 'ID=@{input.list}' {input.vcf} -O z -o {output} > {log} 2>&1
            """

rule idx_vcf:
    input:  rules.exclude_SNPs.output
    output: GPO['idx_vcf']['output']
    envmodules: config['envmodules']['bcftools']
    benchmark: GPO['idx_vcf']['benchmark']
    log:    GPO['idx_vcf']['log']
    shell:  """
            tabix -p vcf {input} > {log} 2>&1
            """

rule create_combined_log:
    input:  impute_check_cat = rules.impute_check_cat.output,
            vcf_cat = rules.vcf_cat.output,
            filter_tags = rules.filter_tags.output,
            check_vcf_log = rules.check_VCF.output.log_out,
            exclude_snps = rules.exclude_SNPs.output
    output: combined_log = GPO['create_combined_log']['output']
    envmodules: config['envmodules']['bcftools']
    message: "Summarise logs from previous rules"
    benchmark: GPO['create_combined_log']['benchmark']
    log:      GPO['create_combined_log']['log']
    shell:
        """
        exec > {log} 2>&1
        # Initialize the log file
        echo "Combined Genotype Processing Log - Started on $(date)" > {output.combined_log}
        echo "==================================================" >> {output.combined_log}
        echo "" >> {output.combined_log}

        # impute_check_cat: Copy contents of the output file
        echo "Step 1: impute_check_cat - Combining per-chromosome imputation summaries" >> {output.combined_log}
        echo "Description: This step concatenates per-chromosome summary statistics from imputation checks, including SNP counts and quality metrics like MAF and R²." >> {output.combined_log}
        echo "Contents of {input.impute_check_cat}:" >> {output.combined_log}
        cat {input.impute_check_cat} >> {output.combined_log}
        echo "" >> {output.combined_log}

        # vcf_cat: Count number of SNPs in the output VCF
        echo "Step 2: vcf_cat - Concatenating chromosome VCFs" >> {output.combined_log}
        echo "Description: This step combines all chromosome-specific VCF files into a single VCF file containing all imputed genotypes." >> {output.combined_log}
        vcf_cat_snp_count=$(bcftools view -H {input.vcf_cat} | wc -l)
        echo "Number of SNPs in {input.vcf_cat}: $vcf_cat_snp_count" >> {output.combined_log}
        echo "" >> {output.combined_log}

        # filter_tags: Count number of SNPs in the output VCF
        echo "Step 3: filter_tags - Filtering SNPs based on quality thresholds" >> {output.combined_log}
        echo "Description: This step applies filters for Hardy-Weinberg equilibrium (HWE), minor allele frequency (MAF), and imputation quality (R²) to retain high-quality SNPs." >> {output.combined_log}
        filter_tags_snp_count=$(bcftools view -H {input.filter_tags} | wc -l)
        echo "Number of SNPs in {input.filter_tags}: $filter_tags_snp_count" >> {output.combined_log}
        echo "" >> {output.combined_log}

        # check_VCF: Extract lines between REPORT and ACTION ITEM (exclusive)
        echo "Step 4: check_VCF - Validating VCF against reference" >> {output.combined_log}
        echo "Description: This step checks the filtered VCF against a reference genome to identify problematic SNPs, outputting a list of SNPs to exclude." >> {output.combined_log}
        echo "Selected contents from {input.check_vcf_log} (REPORT section):" >> {output.combined_log}
        sed -n '/---------------     REPORT     ---------------/,/---------------     ACTION ITEM     ---------------/{{/---------------     ACTION ITEM     ---------------/d;p}}' {input.check_vcf_log} >> {output.combined_log}
        echo "" >> {output.combined_log}

        # exclude_SNPs: Count number of SNPs in the output VCF
        echo "Step 5: exclude_SNPs - Excluding problematic SNPs" >> {output.combined_log}
        echo "Description: This final step removes SNPs identified as problematic in the check_VCF step, producing the final filtered VCF." >> {output.combined_log}
        exclude_snps_count=$(bcftools view -H {input.exclude_snps} | wc -l)
        echo "Number of SNPs in {input.exclude_snps}: $exclude_snps_count" >> {output.combined_log}
        echo "" >> {output.combined_log}

        # Final timestamp
        echo "Log completed on $(date)" >> {output.combined_log}
        """

rule get_sample_list:
    input:  rules.idx_vcf.output
    output: GPO['get_sample_list']['output']
    params: GPO['exclude_SNPs']['output']  
    envmodules: config['envmodules']['bcftools']
    message: "Extract final sample list"
    benchmark: GPO['get_sample_list']['benchmark']
    log:      GPO['get_sample_list']['log']
    shell:  "bcftools query -l {params} > {output} 2> {log}"


rule vcf_to_plink:
    input:  vcf = rules.exclude_SNPs.output,
            idx = rules.idx_vcf.output
    output: GPO['vcf_to_plink']['output']
    params: GPO['vcf_to_plink']['params'] 
    envmodules: config['envmodules']['plink']
    message: "Convert genotypes VCF to plink format"
    benchmark: GPO['vcf_to_plink']['benchmark']
    log:    GPO['vcf_to_plink']['log']
    shell:  "plink --vcf {input.vcf} --double-id --make-bed --out {params} > {log} 2>&1"

rule get_ld_pruned_snps:
    input:  rules.vcf_to_plink.output
    output: GPO['get_ld_pruned_snps']['output']
    params: input_prefix = GPO['vcf_to_plink']['params'],
            output_prefix = GPO['get_ld_pruned_snps']['params']
    envmodules: config['envmodules']['plink']
    message: "LD prune SNPs before running PCA on genotypes"
    benchmark: GPO['get_ld_pruned_snps']['benchmark']
    log:    GPO['get_ld_pruned_snps']['log']
    shell:  """
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
    message: "Create genotypes plink-format file containing only pruned SNPs"
    benchmark: GPO['prune_genotypes']['benchmark']
    envmodules: config['envmodules']['plink']
    log:    GPO['prune_genotypes']['log']
    shell:  """
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
    shell:  """
            plink --bfile {params.input_prefix} \
                  --pca {params.pcs} \
                  --out {params.output_prefix} > {log} 2>&1
            """

rule geno_post_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  stats = rules.calc_genotype_pcs.output,
            rmd_script = GPO["geno_post_report"]["script"]
    output: GPO["geno_post_report"]["html"]
    params: in_dir = GPO["geno_post_report"]["in_dir"],
            bmark_dir = GPO["geno_post_report"]["bmark_dir"],
            output_file = GPO["geno_post_report"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate genotyping pre-imputation report"
    benchmark: GPO["geno_post_report"]["benchmark"]
    log: GPO["geno_post_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["04genotypes_post_impute"] = [GPO["geno_post_report"]["html"]]
