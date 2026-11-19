##  Settings  ------------------------------------------------------------------------------
GPR = config["geno_pre_impute"]


##  Rules  ---------------------------------------------------------------------------------

rule convert_raw:
    input:   ped = GPR["convert_raw"]["input"]["ped"],
             map = GPR["convert_raw"]["input"]["map"]
    output:  GPR["convert_raw"]["output"]
    resources: threads = 1, mem_mb = 50000, time = "0-1:00:00"
    params:  prefix_in = GPR["convert_raw"]["prefix_in"],
             prefix_out = GPR["convert_raw"]["prefix_out"]
    envmodules: config["envmodules"]["plink"]
    message: "Converting raw genotypes to plink format"
    benchmark: GPR["convert_raw"]["benchmark"]
    log:     GPR["convert_raw"]["log"]
    shell:
        """
        plink --file {params.prefix_in} --make-bed --freq --out {params.prefix_out} > {log} 2>&1
        """

rule add_sex_to_fam:
    input:  bed  = rules.convert_raw.output,
            sexfile = GPR["add_sex_to_fam"]["sexfile"]
    output: GPR["add_sex_to_fam"]["output"]
    params: prefix_in = GPR["convert_raw"]["prefix_out"],
            prefix_out = GPR["add_sex_to_fam"]["prefix_out"]
    envmodules: config["envmodules"]["plink"]
    message: "Standardise fam col 1 and col2 and add sex annotations to fam"
    benchmark: GPR["add_sex_to_fam"]["benchmark"]
    log:    GPR["add_sex_to_fam"]["log"]
    shell:
        """
        # Step 1: overwrite fam so FID = IID
        awk '{{ $1=$2; print }}' {params.prefix_in}.fam > {params.prefix_in}.fam.tmp
        mv {params.prefix_in}.fam.tmp {params.prefix_in}.fam

        # Step 2: add sex annotations
        plink --bfile {params.prefix_in} \
              --update-sex {input.sexfile} \
              --make-bed --freq \
              --out {params.prefix_out} > {log} 2>&1
        head {params.prefix_out}.fam >> {log}
        """

rule rm_dup_sample:
    input:  rules.add_sex_to_fam.output
    output: GPR["rm_dup_sample"]["output"]
    params: prefix_in = GPR["add_sex_to_fam"]["prefix_out"],
            prefix_out = GPR["rm_dup_sample"]["prefix_out"],
            remove_file = GPR["rm_dup_sample"]["remove_file"]
    envmodules: config["envmodules"]["plink"]
    message: "Remove duplicate sample: 14493"
    benchmark: GPR["rm_dup_sample"]["benchmark"]
    log:    GPR["rm_dup_sample"]["log"]
    shell:
        """
        plink --bfile {params.prefix_in} \
              --remove {params.remove_file} \
              --make-bed --freq --out {params.prefix_out} > {log} 2>&1
        echo "Samples after removal:  $(wc -l < {params.prefix_out}.fam)" >> {log}
        """

# LD-prune first so heterozygosity is computed on independent SNPs
rule check_het:
    input:  rules.rm_dup_sample.output
    output: GPR["check_het"]["output"]
    params: prefix = GPR["rm_dup_sample"]["prefix_out"],
            report_dir = GPR["genotype-qc2hrc"]["report_dir"]
    envmodules: config["envmodules"]["plink"]
    message: "Check for sample autosomal heterozygosity"
    benchmark: GPR["check_het"]["benchmark"]
    log:    GPR["check_het"]["log"]
    shell:
        """
        exec > {log} 2>&1
        plink --bfile {params.prefix} --indep-pairwise 50 5 0.2 --out {params.prefix}
        plink --bfile {params.prefix} --extract {params.prefix}.prune.in --het --out {params.prefix}
        cp {output} {params.report_dir}
        """

rule make_kgp3_pgen:
    input:   GPR["make_kgp3_pgen"]["input"]
    output:  GPR["make_kgp3_pgen"]["output"]
    resources: threads = 1, mem_mb = 50000, time = "0-1:00:00"
    envmodules: config["envmodules"]["plink"]
    params:  workdir = GPR["make_kgp3_pgen"]["workdir"]
    message: "Create kgp3 pgen file"
    benchmark: GPR["make_kgp3_pgen"]["benchmark"]
    log:      GPR["make_kgp3_pgen"]["log"]
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
    output: GPR["genotype-qc2hrc"]["output"]
    singularity: config["containers"]["genotype-qc2hrc"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    params: prefix_in = GPR["genotype-qc2hrc"]["prefix_in"],
            prefix_out = GPR["genotype-qc2hrc"]["prefix_out"],
            shortname = GPR["genotype-qc2hrc"]["shortname"],
            outdir = GPR["genotype-qc2hrc"]["outdir"],
            workdir = GPR["genotype-qc2hrc"]["workdir"],
            report_dir = GPR["genotype-qc2hrc"]["report_dir"]
    message: "Run GenotypeQCtoHRC to prep genotypes for imputation"
    benchmark: GPR["genotype-qc2hrc"]["benchmark"]
    log:    GPR["genotype-qc2hrc"]["log"]
    shell:
            """
            (cd {params.workdir} && Rscript GenotypeQCtoHRC.R \
              --file {params.prefix_in} \
              --name {params.prefix_out} \
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
             cp {params.outdir}eqtl_genotypes_hg19.qc5.ancestry_inference.txt {params.report_dir}
             """

rule cat_genotypes:
    # Note here that input vcfs are aligned to hg38 despite what their name says
    input:  rules.genotype_qc2hrc.output
    output: GPR["cat_genotypes"]["output"]
    envmodules: config["envmodules"]["bcftools"]
    params: in_dir = GPR["cat_genotypes"]["in_dir"]
    message: "Cat chr specific genotypes into a single file"
    benchmark: GPR["cat_genotypes"]["benchmark"]
    log:    GPR["cat_genotypes"]["log"]
    shell:
        """
        exec > {log} 2>&1
        bcftools concat -Oz -o {output} {params.in_dir}eqtl_genotypes_hg19.gh.topmed.chr{{1..22}}.vcf.gz
        bcftools index -t {output}
        """

# Rm strand ambiguous A/T and G/C SNPs with MAF > 0.4
rule rm_maf_ambig:
    input:  rules.cat_genotypes.output
    output: GPR["rm_maf_ambig"]["output"]
    params: prefix_in = GPR["rm_maf_ambig"]["prefix_in"],
            prefix_out = GPR["rm_maf_ambig"]["prefix_out"],
            ambig_snp_lst = GPR["rm_maf_ambig"]["ambig_snp_lst"]
    envmodules: config["envmodules"]["plink"]
    message: "Rm strand ambiguous SNPs with MAF > 0.4"
    benchmark: GPR["rm_maf_ambig"]["benchmark"]
    log:    GPR["rm_maf_ambig"]["log"]
    shell:
        """
        exec > {log} 2>&1

        # Convert VCF to PLINK binary + frequency calc
        plink --vcf {input} --make-bed --freq --out {params.prefix_in}

        # Extract high-MAF ambiguous SNPs
        awk 'NR>1 && (($3=="A" && $4=="T") || ($3=="T" && $4=="A") || \
                      ($3=="G" && $4=="C") || ($3=="C" && $4=="G")) && $5>0.4 {{print $2}}' \
            {params.prefix_in}.frq > {params.ambig_snp_lst}

        # Remove them
        plink --bfile {params.prefix_in} --exclude {params.ambig_snp_lst} \
              --make-bed --out {params.prefix_out}
        """

# Rm rare SNPs with MAF < 0.01
rule rm_rare_snps:
    input:  rules.rm_maf_ambig.output
    output: GPR["rm_rare_snps"]["output"]
    params: prefix_in = GPR["rm_maf_ambig"]["prefix_out"],
            prefix_out = GPR["rm_rare_snps"]["prefix_out"],
            report_dir = GPR["genotype-qc2hrc"]["report_dir"]
    envmodules: config["envmodules"]["plink"]
    message: "Rm rare SNPs with MAF < 0.01"
    benchmark: GPR["rm_rare_snps"]["benchmark"]
    log:    GPR["rm_rare_snps"]["log"]
    shell:
        """
        plink --bfile {params.prefix_in} \
              --maf 0.01 \
              --make-bed --out {params.prefix_out} > {log} 2>&1
        grep removed {log} > {params.report_dir}rm_rare_snps.txt
        """

rule split_chrs:
    input:  rules.rm_rare_snps.output
    output: GPR["split_chrs"]["output"]
    params: prefix_in = GPR["rm_rare_snps"]["prefix_out"],
            prefix_vcf = GPR["split_chrs"]["prefix_vcf"],
            outdir = GPR["split_chrs"]["outdir"]
    benchmark: GPR["split_chrs"]["benchmark"]
    message: "Split chrs for TOPMED imputation"
    envmodules: config["envmodules"]["plink1"], config["envmodules"]["bcftools"]
    log:    GPR["split_chrs"]["log"]
    shell:
        """
        exec > {log} 2>&1

        # Convert to VCF
        plink --bfile {params.prefix_in} --recode vcf bgz \
              --out {params.prefix_vcf}
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
    output: GPR["gather_stats"]["output"]
    params: gQC2hrc_vcf_dir = GPR["cat_genotypes"]["in_dir"],
            final_vcf_dir = GPR["split_chrs"]["outdir"],
            report_dir = GPR["genotype-qc2hrc"]["report_dir"]
    message: "Gather SNP stats for reporting"
    benchmark: GPR["gather_stats"]["benchmark"]
    log:    GPR["gather_stats"]["log"]
    shell:
        """
        exec > {log} 2>&1
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
    input:  stats = rules.gather_stats.output,
            rmd_script = GPR["geno_pre_report"]["script"]
    output: GPR["geno_pre_report"]["html"]
    params: in_dir = GPR["geno_pre_report"]["in_dir"],
            ambig_snp_lst = GPR["rm_maf_ambig"]["ambig_snp_lst"],
            bmark_dir = GPR["geno_pre_report"]["bmark_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate genotyping pre-imputation report"
    benchmark: GPR["geno_pre_report"]["benchmark"]
    log:    GPR["geno_pre_report"]["log"]
    shell:
        """
        Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd_script}', \
            output_file = file.path(wd, '{output}'), \
            knit_root_dir = wd, \
            params = list(in_dir = '{params.in_dir}', \
                          ambig_snp_lst = '{params.ambig_snp_lst}', \
                          bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["03genotypes_pre_impute"] = [GPR["geno_pre_report"]["html"]]
