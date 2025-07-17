localrules: prep_exp_data

def get_all_genes(cell_types):
    genes = []
    for cell_type in cell_types:
        coord_file = f"../results/06TWAS/fusion_input/{cell_type}_gene_coord.txt"
        with open(coord_file, 'r') as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                gene_id = parts[0]
                genes.append(gene_id)
    return list(set(genes))

# Rule to download gemma binary for FUSION
get_gemma:
   ouput:  config["twas"]["get_gemma"]["output"] 
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

rule get_gene_info:
    input:  config["twas"]["prep_exp_data"]["coord"]
    output: "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_gene_info.txt"
    params: "{gene_id}"
    shell:
        """
        grep -w {params} {input} > {output}
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
        GENE_CHR=$(awk '{{print $2}}' {input.gene_info})
        GENE_START=$(awk '{{print $3}}' {input.gene_info})
        GENE_END=$(awk '{{print $4}}' {input.gene_info})
        CIS_START=$((GENE_START - 500000))
        CIS_END=$((GENE_END + 500000))
        
        plink2 --bfile {params.prefix_in} \
               --chr ${GENE_CHR} \
               --from-bp ${CIS_START} \
               --to-bp ${CIS_END} \
               --make-bed \
               --out {params.prefix_out}
        """


# Rule to compute FUSION weights for each gene
#rule compute_weights:
#    input:  geno_bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
#            geno_bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
#            geno_fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
#            exp = config["twas"]["prep_exp_data"]["exp"],
#            covar = "../results/03SCANPY/pseudobulk/{cell_type}_covariates.txt",
#            ldref = lambda wildcards: f"../resources/ldsr/ldsr_hg38_refs/plink_files/1000G.EUR.hg38.{wildcards.chr}.bed",
#            coord = config["twas"]["prep_exp_data"]["coord"],
#            gemma = config["twas"]["get_gemma"]["output"]
#    output: config["twas"]["compute_weights"]["output"]
#    params: chr = "{chr}",
#            gene = "{gene_id}",
#            outdir = config["twas"]["compute_weights"]["outdir"]
#    singularity: config["containers"]["R"]
#    log:    config["twas"]["compute_weights"]["log"]
#    shell:
#            """
#            Rscript ../resources/fusion/FUSION.compute_weights.R \
#              --bfile {input.geno_bed} \
#              --pheno {input.exp} \
#              --covar {input.covar} \
#              --hsq_p 0.01 \
#              --crossval 5 \
#              --chr {params.chr} \
#              --gene {params.gene} \
#              --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
#              --PATH_gcta resources/fusion/gcta_nr_robust \
#              --PATH_gemma ../resources/gemma/gemma \
#              --out {output} > {log} 2>&1
#            """

#rule aggregate_weights:
#    input:  lambda wildcards: expand(config["twas"]["compute_weights"]["output"],
#                                     cell_type = wildcards.cell_type,
#                                     gene_id = [g["gene_id"] for g in get_genes(wildcards)],
#                                     chr = [g["chr"] for g in get_genes(wildcards)])
#    output: config["twas"]["aggregate_weights"]["output"]
#    log:    config["twas"]["aggregate_weights"]["log"]
#    shell:
#            """
#            echo "Summary of FUSION weights for {wildcards.cell_type}" > {output}
#            ls {input} | wc -l >> {output}
#            echo "Done" >> {output}
#            """


