localrules: prep_exp_data

def get_genes(wildcards):
    coord_file = f"../results/06TWAS/fusion_input/{wildcards.cell_type}_gene_coord.txt"
    df = pd.read_csv(coord_file, sep="\t")
    return [{"gene_id": row["gene_id"], "chr": row["chr"]} for _, row in df.iterrows()]

# Rule to download gemma binary for FUSION
get_gemma:
   ouput:  config["twas"]["get_gemma"]["output"] 
   shell:  """
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

# Rule to compute FUSION weights for each gene
rule compute_weights:
    input:  geno_bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
            geno_bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
            geno_fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
            exp = config["twas"]["prep_exp_data"]["exp"],
            covar = "../results/03SCANPY/pseudobulk/{cell_type}_covariates.txt",
            ldref = lambda wildcards: f"../resources/ldsr/ldsr_hg38_refs/plink_files/1000G.EUR.hg38.{wildcards.chr}.bed",
            coord = config["twas"]["prep_exp_data"]["coord"],
            gemma = config["twas"]["get_gemma"]["output"]
    output: config["twas"]["prep_exp_data"]["output"]
    params: chr = "{chr}",
            gene = "{gene_id}",
            outdir = config["twas"]["compute_weights"]["outdir"]
    singularity: config["containers"]["R"]
    log:    config["twas"]["compute_weights"]["log"]
    shell:
            """
            Rscript ../resources/fusion/FUSION.compute_weights.R \
              --bfile {input.geno_bed} \
              --pheno {input.exp} \
              --covar {input.covar} \
              --hsq_p 0.01 \
              --crossval 5 \
              --chr {params.chr} \
              --gene {params.gene} \
              --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
              --PATH_gcta resources/fusion/gcta_nr_robust \
              --PATH_gemma ../resources/gemma/gemma \
              --out {output} > {log} 2>&1
            """

rule aggregate_weights:
    input:  lambda wildcards: expand("../results/06TWAS/weights/{cell_type}/{gene_id}.wgt.RDat",
                                     cell_type = wildcards.cell_type,
                                     gene_id = [g["gene_id"] for g in get_genes(wildcards)],
                                     chr = [g["chr"] for g in get_genes(wildcards)])
    output: "../results/06TWAS/weights/{cell_type}/summary.txt"
    log:    "../results/00LOG/06TWAS/aggregate_{cell_type}.log"
    shell:
            """
            echo "Summary of FUSION weights for {wildcards.cell_type}" > {output.summary}
            ls {input} | wc -l >> {output.summary}
            echo "Done" >> {output.summary}
            """


