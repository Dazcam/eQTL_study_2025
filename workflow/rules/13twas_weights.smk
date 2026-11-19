##  Settings  ------------------------------------------------------------------------------
TW = config["twas_weights"]


##  Rules  ---------------------------------------------------------------------------------

rule get_gemma:
    output:  TW["get_gemma"]["output"]
    message: "Download gemma binary for FUSION"
    benchmark: TW["get_gemma"]["benchmark"]
    log:      TW["get_gemma"]["log"]
    shell:
            """
            exec > {log} 2>&1
            wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
            gunzip gemma-0.98.5-linux-static-AMD64.gz
            mv gemma-0.98.5-linux-static-AMD64 {output}
            chmod +x {output}
            """

rule prep_exp_data:
    input:  TW["prep_exp_data"]["input"]
    output: exp = TW["prep_exp_data"]["exp"],
            coord = TW["prep_exp_data"]["coord"]
    params: run_test = config['twas_run_test']
    singularity: config["containers"]["r_eqtl"]
    message: "Creating plink ready expression data and gene coordinate file for FUSION weight calculation"
    benchmark: TW["prep_exp_data"]["benchmark"]
    log:    TW["prep_exp_data"]["log"]
    script: TW["prep_exp_data"]["script"]

rule convert_vcf:
    input: TW["convert_vcf"]["input"]
    output: bed = TW["convert_vcf"]["bed"],
            bim = TW["convert_vcf"]["bim"],
            fam = TW["convert_vcf"]["fam"]
    params: prefix = TW["convert_vcf"]["prefix"]
    envmodules: config["envmodules"]["plink"]
    message: "Converting genotypes from VCF to PLINK binary (bed/bim/fam) format for FUSION weight calculation"
    benchmark: TW["convert_vcf"]["benchmark"]
    log:    TW["convert_vcf"]["log"]
    shell:  """plink2 --vcf {input} --make-bed --out {params.prefix} > {log} 2>&1"""

rule get_ldref_snplist:
    input:       expand(TW["get_ldref_snplist"]["input"], chr = range(1, 23))
    output:      TW["get_ldref_snplist"]["snp_lst"]
    message:    "Extract SNP IDs from LD reference .bim files"
    benchmark:  TW["get_ldref_snplist"]["benchmark"]
    log:        TW["get_ldref_snplist"]["log"]
    shell:      """
                # Extract SNP IDs from LD reference .bim files
                cat {input} | cut -f 2 > {output} 2> {log}
                """

rule prepare_covar:
    input:
        lambda w: config["tensorQTL"]["split_covariates"]["output"].format(
            cell_type=w.cell_type,
            norm_method="quantile",
            geno_pc=config['tensorQTL']['geno_pcs'],
            exp_pc=config['exp_pc_map'][w.cell_type]
        )
    output: TW["prepare_covar"]["output"]
    message: "Prepare covariate gene expression data for FUSION"
    singularity: config["containers"]["r_eqtl"]
    benchmark: TW["prepare_covar"]["benchmark"]
    log:  TW["prepare_covar"]["log"]
    script:   TW["prepare_covar"]["script"]

rule compute_weights:
    input:  bed = rules.convert_vcf.output.bed,
            bim = rules.convert_vcf.output.bim,
            fam = rules.convert_vcf.output.fam,
            expr_file = rules.prep_exp_data.output.exp,
            coord_file = rules.prep_exp_data.output.coord,
            ldref_snps = rules.get_ldref_snplist.output,
            covar = rules.prepare_covar.output,
            gemma = rules.get_gemma.output
    output: touch(TW["compute_weights"]["output"])
    params: prefix_in = TW["convert_vcf"]["prefix"],
            outdir = TW["compute_weights"]["outdir"],
            tmp_dir = TW["compute_weights"]["tmp_dir"],
            fusion_script = TW["compute_weights"]["fusion_script"],
            gcta = TW["compute_weights"]["gcta"],
            plink = TW["compute_weights"]["plink"],
            cis_window = config["tensorQTL"]["window"],
            hsq_p = TW["compute_weights"]["hsq_p"],
            crossval = TW["compute_weights"]["crossval"],
            models = TW["compute_weights"]["models"]
    singularity: config["containers"]["twas"]
    message: "Running FUSION weights for ALL genes in {wildcards.cell_type}"
    benchmark: TW["compute_weights"]["benchmark"]
    log: TW["compute_weights"]["log"]
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
    output: TW["make_pos_file"]["output"]
    params: cis_window = config["tensorQTL"]["window"]
    message:  "Generate FUSION .pos file for {wildcards.cell_type}"
    benchmark: TW["make_pos_file"]["benchmark"]
    log:      TW["make_pos_file"]["log"]
    script:   TW["make_pos_file"]["script"]

rule twas_weights_report:
    input:  pos = expand(rules.make_pos_file.output, cell_type = config['cell_types']),
            coord = expand(rules.prep_exp_data.output.coord, cell_type = config['cell_types']),
            rmd_script = TW["twas_weights_report"]["script"],
            functions = TW["twas_weights_report"]["functions"]
    output: TW["twas_weights_report"]["output"]
    params: cell_types = ",".join(config["cell_types"]),
            weights_dir = lambda w: TW["compute_weights"]["outdir"],
            pos_file = lambda w: TW["make_pos_file"]["output"],
            coord_file = lambda w: TW["prep_exp_data"]["coord"],
            weights_log = lambda w: TW["compute_weights"]["log"],
            bmark_dir = TW["twas_weights_report"]["bmark_dir"],
            cis_window = config["tensorQTL"]["window"],
            hsq_p = TW["compute_weights"]["hsq_p"],
            crossval = TW["compute_weights"]["crossval"],
            models = TW["compute_weights"]["models"]
    singularity: config["containers"]["r_eqtl"]
    message:  "Generate TWAS weights report"
    benchmark: TW["twas_weights_report"]["benchmark"]
    log: TW["twas_weights_report"]["log"]
    shell:
        """
        Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd_script}', \
            output_file = file.path(wd, '{output}'), \
            knit_root_dir = wd, \
            params = list(cell_types = '{params.cell_types}', \
                          weights_dir = '{params.weights_dir}', \
                          pos_file = '{params.pos_file}', \
                          coord_file = '{params.coord_file}', \
                          weights_log = '{params.weights_log}', \
                          bmark_dir = '{params.bmark_dir}', \
                          cis_window = '{params.cis_window}', \
                          hsq_p = '{params.hsq_p}', \
                          crossval = '{params.crossval}', \
                          models = '{params.models}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["13twas_weights"] = [TW["twas_weights_report"]["output"]]
