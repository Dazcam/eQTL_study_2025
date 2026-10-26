configfile: "../config/config.yaml"

rule get_gemma:
   output:  config["twas"]["get_gemma"]["output"] 
   message: "Download gemma binary for FUSION"
   benchmark: "reports/benchmarks/11twas.get_gemma.benchmark.txt"   
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
    params: run_test = config['twas_run_test']
    singularity: config["containers"]["r_eqtl"]
    message: "Creating plink ready expression data and gene coordinate file for FUSION weight calculation"
    benchmark: "reports/benchmarks/11twas.prep_exp_data_{cell_type}.benchmark.txt"    
    log:    config["twas"]["prep_exp_data"]["log"]
    script: "../scripts/twas_prep_expr_for_FUSION.R"

rule convert_vcf:
    input: config["twas"]["convert_vcf"]["input"]
    output: bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"]
    params: config["twas"]["convert_vcf"]["prefix"]
    envmodules: "plink/2.0"
    message: "Converting genotypes from VCF to PLINK2 format for FUSION weight calculation"
    benchmark: "reports/benchmarks/11twas.convert_vcf.benchmark.txt"
    log:    config["twas"]["convert_vcf"]["log"]
    shell:  """plink2 --vcf {input} --make-bed --out {params} > {log} 2>&1"""   

rule get_ldref_snplist:
    input:       expand(config["twas"]["get_ldref_snplist"]["input"], chr = range(1, 23))
    output:      config["twas"]["get_ldref_snplist"]["snp_lst"]
    envmodules: "plink/2.0"
    message:    "Restrict genotyped SNPs to LD reference SNPs"
    benchmark:  "reports/benchmarks/11twas.restrict_geno_to_ldref.benchmark.txt"
    log:        config["twas"]["get_ldref_snplist"]["log"]
    shell:      """
                # Extract SNP IDs from LD reference .bim files
                cat {input} | cut -f 2 > {output} 2>&1 | tee {log}                
                """

rule prepare_covar:
    input:  
        lambda w: config["tensorQTL"]["split_covariates"]["output"].format(
            cell_type=w.cell_type, 
            norm_method="quantile",
            geno_pc=config['tensorQTL']['geno_pcs'], 
            exp_pc=config['exp_pc_map'][w.cell_type]
        )
    output: config["twas"]["prepare_covar"]["output"]
    message: "Prepare covariate gene expression data for FUSION"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/11twas.prepare_covariates.benchmark_{cell_type}.txt"
    log:  config["twas"]["prepare_covar"]["log"]
    script:   "../scripts/twas_prep_covar_for_FUSION.R"    

rule compute_weights:
    input:  bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"],
            expr_file = rules.prep_exp_data.output.exp,
            coord_file = rules.prep_exp_data.output.coord,
            ldref_snps = rules.get_ldref_snplist.output,
            covar = rules.prepare_covar.output,
            gemma = config["twas"]["get_gemma"]["output"]
    output: touch(config["twas"]["compute_weights"]["output"])
    params: prefix_in = config["twas"]["convert_vcf"]["prefix"],
            outdir = config["twas"]["compute_weights"]["outdir"],
            tmp_dir = config["twas"]["compute_weights"]["tmp_dir"],
            fusion_script = config["twas"]["compute_weights"]["fusion_script"],
            gcta = config["twas"]["compute_weights"]["gcta"],
            plink = config["twas"]["compute_weights"]["plink"],
            cis_window = config["tensorQTL"]["window"],
            hsq_p = config["twas"]["compute_weights"]["hsq_p"],
            crossval = config["twas"]["compute_weights"]["crossval"],
            models = config["twas"]["compute_weights"]["models"]
    singularity: config["containers"]["twas"]
    message: "Running FUSION weights for ALL genes in {wildcards.cell_type}"
    benchmark: "reports/benchmarks/11twas.compute_all_weights_{cell_type}.benchmark.txt"
    log: config["twas"]["compute_weights"]["log"]
    shell:
        r"""
        echo "Starting weight computation for {wildcards.cell_type}" > {log}

        COORD={input.coord_file}
        EXPR={input.expr_file}
        LDREF={input.ldref_snps}
        COVAR={input.covar}
        GEMMA={input.gemma}
        GENOPREFIX={params.prefix_in}
        OUTDIR={params.outdir}

        mkdir -p $OUTDIR
        
        while read CHR START END GENE; do
            echo "Processing $GENE" >> {log}

            WORKDIR="{params.tmp_dir}/{wildcards.cell_type}/${{GENE}}"
            mkdir -p $WORKDIR

            CIS_START=$((START - {params.cis_window}))
            CIS_END=$((END + {params.cis_window}))
            [ $CIS_START -lt 0 ] && CIS_START=0

            # Extract cis SNPs
            plink2 --bfile $GENOPREFIX \
                   --chr $CHR \
                   --from-bp $CIS_START \
                   --to-bp $CIS_END \
                   --extract $LDREF \
                   --force-intersect \
                   --make-bed \
                   --out $WORKDIR/cis || true

            for ext in bed bim fam; do
                [[ -s "$WORKDIR/cis.$ext" ]] || touch "$WORKDIR/cis.$ext"
            done

            # Extract expression column
            COL=$(head -1 $EXPR | tr '\t' '\n' | grep -n -w $GENE | cut -d: -f1)
            cut -f1,2,$COL $EXPR > $WORKDIR/pheno.txt

            # Run FUSION if SNPs exist
            OUTPREFIX="$OUTDIR/${{GENE}}"
            if [ -s "$WORKDIR/cis.bed" ]; then
                Rscript {params.fusion_script} \
                    --bfile $WORKDIR/cis \
                    --pheno $WORKDIR/pheno.txt \
                    --covar $COVAR \
                    --hsq_p {params.hsq_p} \
                    --crossval {params.crossval} \
                    --tmp $OUTPREFIX.tmp \
                    --PATH_gemma $GEMMA \
                    --PATH_gcta {params.gcta} \
                    --PATH_plink {params.plink} \
                    --verbose 2 \
                    --models {params.models} \
                    --out $OUTPREFIX >> {log} 2>&1 || true
            fi

            [ -f "$OUTPREFIX.wgt.RDat" ] || touch "$OUTPREFIX.wgt.RDat"
            
            rm -rf $WORKDIR 
            if [ ! -d "$WORKDIR" ]; then
              echo "$WORKDIR removed" >> {log}
            else
              echo "WARNING: Failed to remove $WORKDIR" >> {log}
            fi 
            
        done < <(tail -n +2 $COORD)

        echo "All genes completed" >> {log}
        """

rule make_pos_file:
    input:  weights = rules.compute_weights.output,
            coord = rules.prep_exp_data.output.coord,
            expr = rules.prep_exp_data.output.exp
    output: config["twas"]["make_pos_file"]["output"]
    params: cis_window = config["tensorQTL"]["window"]
    message:  "Generate FUSION .pos file for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/11twas.make_pos_file_{cell_type}.benchmark.txt"
    log:      config["twas"]["make_pos_file"]["log"]
    script:   "../scripts/twas_make_fusion_pos_file.py"

rule twas_weights_report:
    input:  pos = expand(rules.make_pos_file.output, cell_type = config['cell_types']),
            rmd_script = config["twas"]["twas_weights_report"]["script"]
    output: config["twas"]["twas_weights_report"]["output"]
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            log_dir = config["twas"]["twas_weights_report"]["log_dir"],
            bmark_dir = config["twas"]["twas_weights_report"]["bmark_dir"],
            output_file = config["twas"]["twas_weights_report"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message:  "Generate TWAS weights report data for report"
    benchmark: "reports/benchmarks/11twas.twas_weights_summary.benchmark.txt"
    log: config["twas"]["twas_weights_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), log_dir = '{params.log_dir}'))" > {log} 2>&1
        """
