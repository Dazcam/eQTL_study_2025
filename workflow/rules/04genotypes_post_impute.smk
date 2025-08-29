CHROMOSOMES = list(range(1, 23)) + ["X"]

rule impute_check:
    input:  config['geno_post_impute']['impute_check']['input']
    output: tsv = config['geno_post_impute']['impute_check']['tsv'],
            vcf = config['geno_post_impute']['impute_check']['vcf']
    params: dir = config['geno_post_impute']['impute_check']['in_dir'],
            pwd = config['geno_post_impute']['impute_check']['pwd']
    envmodules: "bcftools"
    message: "Evaluate imputation quality for each chromosome"
    benchmark: "reports/benchmarks/geno_post_impute.impute_check_{chr}.txt"
    log: config['geno_post_impute']['impute_check']['log']
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
    output: config['geno_post_impute']['impute_check_cat']['output']
    message: "Cat impute check counts for each chr into single file"
    benchmark: "reports/benchmarks/geno_post_impute.impute_check_cat.txt"
    shell:
        """
        # Write header
        echo -e "chr\tsnp_count\tmaf_lt_0.05\tmaf_gteq_0.05\trsq_lt_0.08\trsq_gteq_0.08\tmaf_gteq_0.05_rsq_gteq_0.08" > {output}

        # Combine per-chromosome summaries
        cat {input} | grep -v "^chr" >> {output}

        # Add overall summary
        awk 'NR>1 {{snp+=$2; mlt+=$3; mgteq+=$4; rlt+=$5; rgteq+=$6; mr+=$7}} END {{print "Overall\t" snp "\t" mlt "\t" mgteq "\t" rlt "\t" rgteq "\t" mr}}' {output} >> {output}
        """

rule dwnld_dbsnp_ref:
    output: config['geno_post_impute']['dwnld_dbsnp_ref']['output']
    params: web_link = config['geno_post_impute']['dwnld_dbsnp_ref']['web_link'],
            prefix = config['geno_post_impute']['dwnld_dbsnp_ref']['prefix']
    message: "Download dbSNP reference"
    benchmark: "reports/benchmarks/geno_post_impute.dwnld_dbsnp_ref.txt"
    log:    config['geno_post_impute']['dwnld_dbsnp_ref']['log']
    shell:
        """
        wget {params.web_link} -O {params.prefix} > {log} 2>&1
        wget {params.web_link}.tbi -O {params.prefix}.tbi >> {log} 2>&1
        """

rule add_rsID:
    input:  vcf = rules.impute_check.output.vcf,
            dbsnp = rules.dwnld_dbsnp_ref.output
    output: config['geno_post_impute']['add_rsID']['output']
    envmodules: "bcftools"
    message: "Add dbSNP rsIDs to each imputed chr-specific VCF file"
    benchmark: "reports/benchmarks/geno_post_impute.add_rsID_{chr}.txt"
    log:    config['geno_post_impute']['add_rsID']['log']
    shell:
        """
        bcftools annotate -a {input.dbsnp} -c ID {input.vcf} -O z -o {output} > {log} 2>&1  
        """

rule vcf_cat:
    input:  expand(rules.add_rsID.output, chr=CHROMOSOMES)
    output: config['geno_post_impute']['vcf_cat']['output']
    envmodules: "bcftools"
    message: "Cat imputed chr-specific VCF files into single file"
    benchmark: "reports/benchmarks/geno_post_impute.vcf_cat.txt"
    log:    config['geno_post_impute']['vcf_cat']['log']
    shell:
        """
        bcftools concat {input} -O z -o {output} > {log} 2>&1
        """

# Compute HWE and filter the concatenated VCF; Keep SNPs passing MAF, R2 and HWE thresholds
rule filter_tags:
    input:  rules.vcf_cat.output
    output: config['geno_post_impute']['filter_tags']['output']
    envmodules: "bcftools"
    message: "Compute HWE and filt concated VCF; Keep SNPs passing MAF, R^2 and HWE threshs"
    benchmark: "reports/benchmarks/geno_post_impute.filter_tags.txt"
    params: hwe = config['geno_post_impute']['filter_tags']['hwe'],
            maf = config['geno_post_impute']['filter_tags']['maf'],
            rsq = config['geno_post_impute']['filter_tags']['rsq']
    log:    config['geno_post_impute']['filter_tags']['log']
    shell:
        """
        bcftools +fill-tags {input} -O z -- -t HWE |\
        bcftools view -e 'HWE<{params.hwe} |\
                          INFO/MAF<{params.maf} |\
                          INFO/R2<{params.rsq}' -O z -o {output} > {log} 2>&1
        """

rule check_VCF:
    input:  vcf = rules.filter_tags.output,
            ref_gz = config['geno_post_impute']['check_VCF']['ref_gz']
    output: log_out = config['geno_post_impute']['check_VCF']['log'],
            list = config['geno_post_impute']['check_VCF']['list'],
            ref = temp(config['geno_post_impute']['check_VCF']['ref']),
            fai = temp(config['geno_post_impute']['check_VCF']['fai'])    
    params: config['geno_post_impute']['check_VCF']['prefix_out']
    envmodules: "samtools"
    message: "Validate the filt VCF against ref genome and ID problematic SNPs"
    benchmark: "reports/benchmarks/geno_post_impute.check_VCF.txt"
    log:    config['geno_post_impute']['check_VCF']['log']
    shell:  """
            gunzip -c {input.ref_gz} > {output.ref} 2>> {log}
            samtools faidx {output.ref} >> {log} 2>&1
            python scripts/checkVCF.py -r {output.ref} -v {input.vcf} -o {params} --exclude {output.list} > {log} 2>&1
            """

rule exclude_SNPs:
    input:  vcf = rules.filter_tags.output,
            list = rules.check_VCF.output.list
    output: config['geno_post_impute']['exclude_SNPs']['output']
    envmodules: "bcftools"
    message: "Exclude SNPs IDed in CheckVCF"
    benchmark: "reports/benchmarks/geno_post_impute.exclude_SNPs.txt"
    log:    config['geno_post_impute']['exclude_SNPs']['log']
    shell:  """
            bcftools view -e 'ID=@{input.list}' {input.vcf} -O z -o {output} > {log} 2>&1
            """

rule idx_vcf:
    input:  rules.exclude_SNPs.output
    output: config['geno_post_impute']['idx_vcf']['output']
    envmodules: "bcftools"
    log:    config['geno_post_impute']['idx_vcf']['log']
    shell:  """
            tabix -p vcf {input} > {log} 2>&1
            """

rule create_combined_log:
    input:  impute_check_cat = rules.impute_check_cat.output,
            vcf_cat = rules.vcf_cat.output,
            filter_tags = rules.filter_tags.output,
            check_vcf_log = rules.check_VCF.output.log_out,
            exclude_snps = rules.exclude_SNPs.output
    output: combined_log = config['geno_post_impute']['create_combined_log']['output']
    envmodules: "bcftools"
    message: "Summarise logs from previous rules"
    benchmark: "reports/benchmarks/geno_post_impute.create_combined_log.txt"
    shell:
        """
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
    output: config['geno_post_impute']['get_sample_list']['output']
    params: config['geno_post_impute']['exclude_SNPs']['output']  
    envmodules: "bcftools"
    message: "Extract final sample list"
    benchmark: "reports/benchmarks/geno_post_impute.get_sample_list.txt"
    shell:  "bcftools query -l {params} > {output}"


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

rule calc_genotype_pcs:
    input:  rules.prune_genotypes.output
    output: config['geno_post_impute']['calc_genotype_pcs']['output']
    params: input_prefix = config['geno_post_impute']['prune_genotypes']['params'],
            output_prefix = config['geno_post_impute']['calc_genotype_pcs']['params'],
            pcs = config['geno_post_impute']['calc_genotype_pcs']['pcs']
    envmodules: "plink/1.9"
    message: "Run PCA on pruned genotypes to get genotype-specific PC covariates"
    benchmark: "reports/benchmarks/geno_post_impute.calc_genotype_pcs.txt"
    log:    config['geno_post_impute']['calc_genotype_pcs']['log']
    shell:  """
            plink --bfile {params.input_prefix} \
                  --pca {params.pcs} \
                  --out {params.output_prefix} > {log} 2>&1
            """
