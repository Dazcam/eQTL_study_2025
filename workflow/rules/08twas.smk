localrules: prep_exp_data

# Rule to download gemma binary for FUSION
rule get_gemma:
   output:  config["twas"]["get_gemma"]["output"] 
   shell:  
           """
           wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
           gunzip gemma-0.98.5-linux-static-AMD64.gz
           mv gemma-0.98.5-linux-static-AMD64 {output}
           chmod +x {output}
           """

rule prep_exp_data:
    input:  config["twas"]["prep_exp_data"]["input"]
    output: exp = config["twas"]["prep_exp_data"]["exp"],
            coord = config["twas"]["prep_exp_data"]["coord"]
    singularity: config["containers"]["R"]
    log:    config["twas"]["prep_exp_data"]["log"]
    script: "../scripts/prep_expr_for_FUSION.R"

# Rule to convert VCF to PLINK
rule convert_vcf:
    input: config["twas"]["convert_vcf"]["input"]
    output: bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"]
    params: config["twas"]["convert_vcf"]["prefix"]
    envmodules: "plink/2.0"
    log:    config["twas"]["convert_vcf"]["log"]
    shell:  """plink2 --vcf {input} --make-bed --out {params} > {log} 2>&1"""   

# Rule to restrict genotypes to LD reference SNPs
rule restrict_geno_to_ldref:
    input:  bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"],
            ldref_bim = expand("../resources/ldsr/ldsr_hg38_refs/plink_files/1000G.EUR.hg38.{chr}.bim", chr = range(1, 23))
    output: bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
            bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
            fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
            snp_list = config["twas"]["restrict_geno_to_ldref"]["snp_lst"]
    params: input = config["twas"]["restrict_geno_to_ldref"]["prefix_in"],
            output = config["twas"]["restrict_geno_to_ldref"]["prefix_out"]
    envmodules: "plink/2.0"
    log:
        "../results/00LOG/06TWAS/restrict_genotypes_to_ldref.log"
    shell:
        """
        # Extract SNP IDs from LD reference .bim files
        cat {input.ldref_bim} | cut -f 2 > {output.snp_list}
        # Restrict genotypes to LD reference SNPs
        plink2 --bfile {params.input} \
              --extract {output.snp_list} \
              --make-bed --out {params.output} > {log} 2>&1
        """

checkpoint get_gene_info:
    input:  config["twas"]["prep_exp_data"]["coord"]
    output: directory("../results/06TWAS/fusion_input/{cell_type}")
    shell:
        """
        mkdir {output}
        tail -n +2 {input} | while IFS=$'\t' read -r line; do
        gene=$(echo "$line" | cut -f4)
        echo "$line" > "{output}/${{gene}}_gene_info.txt"
        done
        """

rule extract_cis_snps:
    input:
        geno_bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
        geno_bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
        geno_fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
        gene_info = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_gene_info.txt"
    output:
        cis_bed = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed",
        cis_bim = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim",
        cis_fam = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam"
    params:
        prefix_in = config["twas"]["restrict_geno_to_ldref"]["prefix_out"],
        prefix_out = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis"
    envmodules: "plink/2.0"
    shell:
        """
        GENE_CHR=$(awk '{{print $1}}' {input.gene_info})
        GENE_START=$(awk '{{print $2}}' {input.gene_info})
        GENE_END=$(awk '{{print $3}}' {input.gene_info})
        CIS_START=$((GENE_START - 500000))
        CIS_END=$((GENE_END + 500000))
        
        # Set CIS_START to 0 if negative
        if (( CIS_START < 0 )); then
          CIS_START=0
        fi
        
        plink2 --bfile {params.prefix_in} \
               --chr ${{GENE_CHR}} \
               --from-bp ${{CIS_START}} \
               --to-bp ${{CIS_END}} \
               --make-bed \
               --out {params.prefix_out}
        """

rule prepare_gene_pheno:
    input:  config["twas"]["prep_exp_data"]["exp"],
    output: "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt"
    params: "{gene_id}"
    shell:
        """
        # Get the column index of the gene_id
        COL_IDX=$(head -n 1 {input} | tr '\t' '\n' | grep -n -w {params} | cut -f1 -d:)
        # Extract FID, IID, and the gene_id column
        cut -f 1,2,${{COL_IDX}} {input} > {output}
        """

rule compute_weights:
    input:
        geno_bed = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed",
        geno_bim = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim",
        geno_fam = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam",
        gene_pheno = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt",
        covar = "../results/03SCANPY/pseudobulk/{cell_type}_covariates.txt",
        gemma = config["twas"]["get_gemma"]["output"]
    output:
        "../results/06TWAS/weights/{cell_type}/{gene_id}.RDat"
    params:
        outdir = "../results/06TWAS/weights/{cell_type}"
    singularity: config["containers"]["R"]
    log: "../results/00LOG/06TWAS/compute_weights_{cell_type}_{gene_id}.log"
    shell:
        """
        Rscript ../resources/fusion/FUSION.compute_weights.R \
          --bfile {input.geno_bed} \
          --pheno {input.gene_pheno} \
          --covar {input.covar} \
          --hsq_p 0.01 \
          --crossval 5 \
          --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
          --PATH_gcta resources/fusion/gcta_nr_robust \
          --PATH_gemma {input.gemma} \
          --out {output} > {log} 2>&1
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.get_gene_info.get(**wildcards).output[0]
    return expand(
        "../results/06TWAS/weights/{cell_type}/{gene_id}.RDat",
        cell_type = wildcards.cell_type,
        gene_id = glob_wildcards(os.path.join(checkpoint_output, "{gene_id}_gene_info.txt")).gene_id,
    )

rule aggregate_per_cell_type:
    input: aggregate_input,
    output: "../results/06TWAS/{cell_type}/all_genes_done.txt"
    run: "cat {input} > {output}"

rule aggregate_all:
    input:  expand("../results/06TWAS/{cell_type}/all_genes_done.txt", cell_type = config['cell_types']),
    output: "../results/06TWAS/all_genes_in_all_cell_types_done.txt",
    shell:  "cat {input} > {output}"
