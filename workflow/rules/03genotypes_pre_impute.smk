configfile: "../config/config.yaml"

rule convert_raw:
    input:   ped = config["geno_pre_impute"]["convert_raw"]["input"]["ped"],
             map = config["geno_pre_impute"]["convert_raw"]["input"]["map"]
    output:  config["geno_pre_impute"]["convert_raw"]["output"]
    resources: threads = 1, mem_mb = 50000, time = "0-1:00:00"
    params:  in_prefix = config["geno_pre_impute"]["convert_raw"]["in_prefix"],
             out_prefix = config["geno_pre_impute"]["convert_raw"]["out_prefix"]
    envmodules: 'plink/1.9'
    message: "Converting raw genotypes to plink format"
    benchmark: "reports/benchmarks/geno_pre_impute.convert_raw.benchmark.txt"
    log:     config["geno_pre_impute"]["convert_raw"]["log"]
    shell:
        """
        plink --file {params.in_prefix} --make-bed --out {params.out_prefix} > {log} 2>&1
        plink --bfile {params.out_prefix} --freq --out {params.out_prefix} >> {log} 2>&1
        """

rule rm_dup_sample:
    input:  rules.convert_raw.output
    output: config["geno_pre_impute"]["rm_dup_sample"]["output"]
    params: prefix_in = config["geno_pre_impute"]["convert_raw"]["out_prefix"],
            prefix_out = config["geno_pre_impute"]["rm_dup_sample"]["out_prefix"],
            remove_file = config["geno_pre_impute"]["rm_dup_sample"]["remove_file"]
    envmodules: "plink/1.9"
    benchmark: "reports/benchmarks/geno_pre_impute.rm_dup_sample.benchmark.txt"
    log:    config["geno_pre_impute"]["rm_dup_sample"]["log"]
    shell:
        """
        plink --bfile {params.prefix_in} \
              --remove {params.remove_file} \
              --make-bed --out {params.prefix_out} > {log} 2>&1
        """

rule make_kgp3_pgen:
    input:   config["geno_pre_impute"]["make_kgp3_pgen"]["input"]
    output:  config["geno_pre_impute"]["make_kgp3_pgen"]["output"]
    resources: threads = 1, mem_mb = 50000, time = "0-1:00:00"
    envmodules: "plink/2.0"
    params:  workdir = config["geno_pre_impute"]["make_kgp3_pgen"]["workdir"]
    message: "Create kgp3 pgen file"
    benchmark: "reports/benchmarks/geno_pre_impute.make_kgp3_pgen.benchmark.txt"
    log:      config["geno_pre_impute"]["make_kgp3_pgen"]["log"]
    shell:
        """
        (cd {params.workdir} && plink2 \
          --pfile all_phase3 \
          --keep kgp3.array_snps.id \
          --extract kgp3.array_snps.snplist \
          --make-pgen \
          --out kgp3.array_snps) > {log} 2>&1
        """

rule genotype_qc2hrc:
    input:  rules.rm_dup_sample.output,
            rules.make_kgp3_pgen.output
    output: config["geno_pre_impute"]["genotype-qc2hrc"]["output"]
    singularity: config["containers"]["genotype-qc2hrc"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    params: in_prefix = config["geno_pre_impute"]["genotype-qc2hrc"]["in_prefix"],
            out_prefix = config["geno_pre_impute"]["genotype-qc2hrc"]["out_prefix"],
            shortname = config["geno_pre_impute"]["genotype-qc2hrc"]["shortname"],
            outdir = config["geno_pre_impute"]["genotype-qc2hrc"]["outdir"],
            workdir = config["geno_pre_impute"]["genotype-qc2hrc"]["workdir"],
            report_dir = config["geno_pre_impute"]["genotype-qc2hrc"]["report_dir"]      
    message: "Run GenotypeQCtoHRC to prep genotypes for imputation"
    benchmark: "reports/benchmarks/geno_pre_impute.genotype_qc2hrc.benchmark.txt"
    log:    config["geno_pre_impute"]["genotype-qc2hrc"]["log"]
    shell:
            """
            (cd {params.workdir} && Rscript GenotypeQCtoHRC.R \
              --file {params.in_prefix} \
              --name {params.out_prefix} \
              --shortname {params.shortname} \
              --gh TRUE \
              --gh-ref TopMed \
              --lo TRUE \
              --lo-in 37 \
              --lo-out 38 \
              --clean TRUE) > {log} 2>&1

             cp {output} {params.report_dir}
             cp {params.outdir}eqtl_genotypes_hg19.qc3.sexcheck {params.report_dir}
             cp {params.outdir}eqtl_genotypes_hg19.qc4.pcrelate.kin {params.report_dir}
             cp {params.outdir}eqtl_genotypes_hg19.qc5.ancestry_inference.txt | cut -f 1-9 {params.report_dir}
             """

rule cat_genotypes:
    # Note here that input vcfs are aligned to hg38 despite what their name says
    input: rules.genotype_qc2hrc.output
    output: config["geno_pre_impute"]["cat_genotypes"]["output"]  
    envmodules: "bcftools/1.16.0"
    params: config["geno_pre_impute"]["cat_genotypes"]["in_dir"]
    benchmark: "reports/benchmarks/geno_pre_impute.cat_genotypes.benchmark.txt"
    shell:
        """ 
        bcftools concat -Oz -o {output} {params}eqtl_genotypes_hg19.gh.topmed.chr{{1..22}}.vcf.gz
        bcftools index -t {output}
        """

## Note we need a work around here to deal with chr specific files
# rm MAF A/T and G/C SNPs with MAFs > 0.4
rule rm_maf_ambig:
    input:  rules.cat_genotypes.output
    output: config["geno_pre_impute"]["rm_maf_ambig"]["output"]
    params: prefix_in = config["geno_pre_impute"]["rm_maf_ambig"]["prefix_in"],
            prefix_out = config["geno_pre_impute"]["rm_maf_ambig"]["prefix_out"],
            ambig_snp_lst = config["geno_pre_impute"]["rm_maf_ambig"]["ambig_snp_lst"]
    envmodules: "plink/1.9"
    benchmark: "reports/benchmarks/geno_pre_impute.rm_maf_ambig.benchmark.txt"
    log:    config["geno_pre_impute"]["rm_maf_ambig"]["log"]
    shell:
        """
        # Convert VCF to PLINK binary + frequency calc
        plink --vcf {input} --make-bed --freq --out {params.prefix_in}
        
        # Extract high-MAF ambiguous SNPs
        awk 'NR>1 && (($3=="A" && $4=="T") || ($3=="T" && $4=="A") || \
                      ($3=="G" && $4=="C") || ($3=="C" && $4=="G")) && $5>0.4 {{print $2}}' \
            {params.prefix_in}.frq > {params.ambig_snp_lst}
        
        # Remove them
        plink --bfile {params.prefix_in} --exclude {params.ambig_snp_lst} \
              --make-bed --out {params.prefix_out} > {log} 2>&1
        """

## Note	we need	a work around here to deal with	chr specific files
# rm MAF < 0.01
rule rm_rare_snps:
    input:  rules.rm_maf_ambig.output
    output: config["geno_pre_impute"]["rm_rare_snps"]["output"]
    params: prefix_in = config["geno_pre_impute"]["rm_maf_ambig"]["prefix_out"],
            prefix_out = config["geno_pre_impute"]["rm_rare_snps"]["prefix_out"]
    envmodules: "plink/1.9"
    benchmark: "reports/benchmarks/geno_pre_impute.rm_rm_rare_snps.benchmark.txt"
    log:    config["geno_pre_impute"]["rm_rare_snps"]["log"]    
    shell:
        """
        plink --bfile {params.prefix_in} \
              --maf 0.01 \
              --make-bed --out {params.prefix_out} > {log} 2>&1
        """

rule split_chrs:
    input:  rules.rm_rare_snps.output
    output: config["geno_pre_impute"]["split_chrs"]["output"]
    params: prefix_in = config["geno_pre_impute"]["rm_rare_snps"]["prefix_out"],
            prefix_vcf = config["geno_pre_impute"]["split_chrs"]["prefix_vcf"], 
            outdir = config["geno_pre_impute"]["split_chrs"]["outdir"]
    benchmark: "reports/benchmarks/geno_pre_impute.split_chrs.benchmark.txt"
    envmodules: "plink/1.9", "compiler/gnu/7/3.0", "bcftools/1.16.0"
    log: config["geno_pre_impute"]["split_chrs"]["log"]
    shell:
        """
        # Convert to VCF
        plink --bfile {params.prefix_in} --recode vcf bgz \
              --out {params.prefix_vcf} > {log} 2>&1
        bcftools index -t {params.prefix_vcf}.vcf.gz
        
        # Split by chromosome 1–22
        for chr in $(seq 1 22); do
            bcftools view -r $chr \
                {params.prefix_vcf}.vcf.gz \
                -Oz -o {params.outdir}eqtl_genotypes_hg38.gh.topmed.chr${{chr}}.vcf.gz
        done
        touch {output}
        """

rule gather_stats:
    input:  rules.split_chrs.output
    output: config["geno_pre_impute"]["gather_stats"]["output"]
    params: gQC2hrc_vcf_dir = config["geno_pre_impute"]["gather_stats"]["gQC2hrc_vcf_dir"],
            final_vcf_dir = config["geno_pre_impute"]["gather_stats"]["final_vcf_dir"],
            report_dir = config["geno_pre_impute"]["genotype-qc2hrc"]["report_dir"]
    log: config["geno_pre_impute"]["gather_stats"]["log"]
    shell:
        """
       	total=0
        for file in {params.gQC2hrc_vcf_dir}eqtl_genotypes_hg19.gh.topmed.chr*.vcf.gz
        do
          chr=$(basename "$file" .vcf.gz | sed 's/.*chr//')
          count=$(zcat "$file" | grep -vc '^#')
          echo -e "${{chr}}\t${{count}}"
          total=$((total + count))
        done > {params.report_dir}gQC2hrc_vcf_chr_counts.tsv

        # append total row
        echo -e "TOTAL\t${{total}}" >> {params.report_dir}gQC2hrc_vcf_chr_counts.tsv        
    
        total=0
        for file in {params.final_vcf_dir}eqtl_genotypes_hg38.gh.topmed.chr*.vcf.gz
        do
          chr=$(basename "$file" .vcf.gz | sed 's/.*chr//')
          count=$(zcat "$file" | grep -vc '^#')
          echo -e "${{chr}}\t${{count}}"
          total=$((total + count))
        done > {params.report_dir}final_vcf_chr_counts.tsv

        # append total row
        echo -e "TOTAL\t${{total}}" >> {params.report_dir}final_vcf_chr_counts.tsv
        touch {output}
        """

rule geno_pre_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  stats = rules.gather_stats.output,
            rmd_script = config["geno_pre_impute"]["geno_pre_report"]["rmd_script"]
    output: config["geno_pre_impute"]["geno_pre_report"]["html"]
    params: in_dir = config["geno_pre_impute"]["geno_pre_report"]["in_dir"],
            output_file = config["geno_pre_impute"]["geno_pre_report"]["out_file"]
    singularity: config["containers"]["R"]
    log: config["geno_pre_impute"]["geno_pre_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}'))" > {log} 2>&1
        """
