CHROMOSOMES = list(range(1, 23)) + ["X"]

rule impute_check:
    input:  config['genotypes']['impute_check']['input']
    output: tsv = config['genotypes']['impute_check']['output']['tsv'],
            vcf = config['genotypes']['impute_check']['output']['vcf']
    params: dir = config['genotypes']['impute_check']['params']['dir'],
            pwd = config['genotypes']['impute_check']['params']['pwd']
    envmodules: "bcftools"
    log: config['genotypes']['impute_check']['log']
    shell:
        """
        # Create directory and unzip
        if [ ! -f {params.dir}chr{wildcards.chr}.dose.vcf.gz ]; then
            unzip -P '{params.pwd}' {input} -d {params.dir} > {log} 2>&1
        else
            echo "VCF file already exists, skipping unzip" >> {log}
        fi

        # Count total SNPs from dose.vcf.gz
        snp_count=$(bcftools view -H {params.dir}chr{wildcards.chr}.dose.vcf.gz | wc -l)

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
    input:  expand(config['genotypes']['impute_check']['output'], chr=CHROMOSOMES)
    output: config['genotypes']['impute_check_cat']['output']
    shell:
        """
        # Write header
        echo -e "chr\tsnp_count\tmaf_lt_0.05\tmaf_gteq_0.05\trsq_lt_0.08\trsq_gteq_0.08\tmaf_gteq_0.05_rsq_gteq_0.08" > {output}

        # Combine per-chromosome summaries
        cat {input} | grep -v "^chr" >> {output}

        # Add overall summary
        awk 'NR>1 {{snp+=$2; mlt+=$3; mgteq+=$4; rlt+=$5; rgteq+=$6; mr+=$7}} END {{print "Overall\t" snp "\t" mlt "\t" mgteq "\t" rlt "\t" rgteq "\t" mr}}' {output} >> {output}
        """

# Rule to download dbSNP Build 156 for hg38
rule dwnld_dbsnp_ref:
    output: config['genotypes']['dwnld_dbsnp_ref']['output']
    params: web_link = config['genotypes']['dwnld_dbsnp_ref']['params']['web_link'],
            prefix = config['genotypes']['dwnld_dbsnp_ref']['params']['prefix']
    log:    config['genotypes']['dwnld_dbsnp_ref']['log']
    shell:
        """
        wget {params.web_link} -O {params.prefix} > {log} 2>&1
        wget {params.web_link}.tbi -O {params.prefix}.tbi >> {log} 2>&1
        """

# Rule to add rsIDs to each chromosome's VCF
rule add_rsID:
    input:  vcf = config['genotypes']['impute_check']['output']['vcf'],
            dbsnp = config['genotypes']['dwnld_dbsnp_ref']['output']
    output: config['genotypes']['add_rsID']['output']
    envmodules: "bcftools"
    log:    config['genotypes']['add_rsID']['log']
    shell:
        """
        bcftools annotate -a {input.dbsnp} -c ID {input.vcf} -O z -o {output} > {log} 2>&1  
        """

# Rule to concatenate all chromosome VCFs into one
rule vcf_cat:
    input:  expand(config['genotypes']['add_rsID']['output'], chr=CHROMOSOMES)
    output: config['genotypes']['vcf_cat']['output']
    envmodules: "bcftools"
    log:    config['genotypes']['vcf_cat']['log']
    shell:
        """
        bcftools concat {input} -O z -o {output} > {log} 2>&1
        """

# Rule to compute HWE and filter the concatenated VCF (using existing MAF and R2)
rule filter_tags:
    input:  config['genotypes']['vcf_cat']['output']
    output: config['genotypes']['filter_tags']['output']
    envmodules: "bcftools"
    log:    config['genotypes']['filter_tags']['log']
    params: hwe = config['genotypes']['filter_tags']['params']['hwe'], 
            maf = config['genotypes']['filter_tags']['params']['maf'],
            rsq = config['genotypes']['filter_tags']['params']['rsq']
            
    shell:
        """
        bcftools +fill-tags {input} -O z -- -t HWE |\
        bcftools view -e 'HWE<{params.hwe} |\
                          INFO/{params.maf} |\
                          INFO/R2<{params.rsq}' -O z -o {output} > {log} 2>&1
        """

# Rule to check for duplicates and excluded SNPs
rule check_VCF:
    input:  config['genotypes']['filter_tags']['output']
    output: config['genotypes']['check_VCF']['output']    
    params: config['genotypes']['check_VCF']['params']
    log:    config['genotypes']['check_VCF']['log']
    shell:
        """
        checkVCF.py -r {input} -o {params} --exclude {output} > {log} 2>&1
        """

# Rule to exclude problematic variants - do I also need to exclude duplicates??
rule exclude_SNPs:
    input:  vcf = config['genotypes']['filter_tags']['output'],
            list = config['genotypes']['check_VCF']['output']
    output: config['genotypes']['exclude_SNPs']['output']
    envmodules: "bcftools"
    log:    config['genotypes']['exclude_SNPs']['log']
    shell:
        """
        bcftools view -e 'ID=@{input.list}' {input.vcf} -O z -o {output} > {log} 2>&1
        """

# Rule to index the final VCF.gz
rule idx_vcf:
    input:  config['genotypes']['exclude_SNPs']['output']
    output: config['genotypes']['idx_vcf']['output']
    envmodules: "bcftools"
    log:    config['genotypes']['idx_vcf']['log']
    shell:
        """
        tabix -p vcf {input} > {log} 2>&1
        """
