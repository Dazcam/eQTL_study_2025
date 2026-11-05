configfile: "../config/config.yaml"

DS          = config["dev_specificity"]
NORM_METHOD = config["tensorQTL"]["norm_methods"][0]
GENO_PC     = config["tensorQTL"]["geno_pcs"]

# Pseudotime bins (NPC-to-*) excluded from the eQTL sharing pipeline
DEVSPEC_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]
DEVSPEC_FETAL_L1   = config["fetal_l1_types"]

# Jang et al. 2026 association files (downloaded by setup get_jang_singlebrain)
JANG_FULL = config["setup"]["get_jang_singlebrain"]["full_assoc"]
JANG_TOP  = config["setup"]["get_jang_singlebrain"]["top_assoc"]


def perm_files(cell_types):
    """tensorQTL permutation output per cell type, at the final expPC (exp_pc_map)."""
    return [config["tensorQTL"]["tensorqtl_perm"]["output"].format(
                cell_type=ct, norm_method=NORM_METHOD, geno_pc=GENO_PC,
                exp_pc=config["exp_pc_map"][ct])
            for ct in cell_types]


def snp_lookup_files(cell_types):
    """06 snp_lookup output (hg38 SNP coordinates) per cell type."""
    return [config["qtl_rep"]["snp_lookup"]["out_file"].format(
                cell_type=ct, norm_method=NORM_METHOD, geno_pc=GENO_PC,
                exp_pc=config["exp_pc_map"][ct])
            for ct in cell_types]


def smr_nom_files(cell_types):
    """10SMR nominal export (cat_tensorqtl_nom_snps) per cell type."""
    return [config["smr"]["cat_tensorqtl_nom_snps"]["cat_snps"].format(cell_type=ct)
            for ct in cell_types]


# Rules  -----------------------------------------------------------------------------

rule extract_eqtl_universe:
    input:  qtl_perm_files = perm_files(DEVSPEC_CELL_TYPES)
    output: DS["extract_universe"]["out_file"]
    params: cell_types       = DEVSPEC_CELL_TYPES,
            cell_types_L1    = DEVSPEC_FETAL_L1,
            jang_label_match = DS["extract_universe"]["jang_label_match"],
            eqtl_fdr         = config["eqtl_fdr"],
            geno_pc          = GENO_PC,
            norm_method      = NORM_METHOD
    singularity: config["containers"]["r_eqtl"]
    message: "Extract unique fetal eGene-eSNP universe (L1+L2)"
    benchmark: DS["extract_universe"]["benchmark"]
    log:    DS["extract_universe"]["log"]
    script: DS["extract_universe"]["script"]

rule classify_eqtl_specificity:
    input:  eqtl_universe   = rules.extract_eqtl_universe.output[0],
            jang_full_files = expand(JANG_FULL, jang_cell_type=config["jang_cell_types"]),
            jang_top_files  = expand(JANG_TOP, jang_cell_type=config["jang_cell_types"])
    output: DS["classify_specificity"]["out_file"]
    params: jang_cell_types = config["jang_cell_types"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 128000, time = "04:00:00"
    threads: 4
    message: "Classify fetal eGene-eSNP pairs as testable/untestable and shared/fetal-specific vs Jang et al. 2026 adult eQTL (all L1+L2 types)"
    benchmark: DS["classify_specificity"]["benchmark"]
    log:    DS["classify_specificity"]["log"]
    script: DS["classify_specificity"]["script"]

rule heterogeneity_test:
    input:  eqtl_universe  = rules.extract_eqtl_universe.output[0],
            classification = rules.classify_eqtl_specificity.output[0],
            pvar_file      = config["tensorQTL"]["convert_genotypes"]["pvar"],
            gene_lookup    = config["gene_lookup"]
    output: DS["heterogeneity_test"]["out_file"]
    params: fdr_thresh = DS["heterogeneity_test"]["fdr_thresh"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 2, mem_mb = 8000, time = "00:20:00"
    threads: 2
    message: "Heterogeneity test (Cochran's Q) on fetal-specific eQTL vs Jang et al. 2026 adult cell types"
    benchmark: DS["heterogeneity_test"]["benchmark"]
    log:    DS["heterogeneity_test"]["log"]
    script: DS["heterogeneity_test"]["script"]

rule prep_ld_clump_input:
    input:  eqtl_universe = rules.extract_eqtl_universe.output[0]
    output: assoc_dir      = directory(DS["prep_ld_clump_input"]["assoc_dir"]),
            singleton_snps = DS["prep_ld_clump_input"]["singleton_snps"],
            snp_list       = DS["prep_ld_clump_input"]["snp_list"]
    singularity: config["containers"]["r_eqtl"]
    message: "Prepare per-gene plink association files for LD-independence clumping"
    benchmark: DS["prep_ld_clump_input"]["benchmark"]
    log:    DS["prep_ld_clump_input"]["log"]
    script: DS["prep_ld_clump_input"]["script"]

rule run_ld_clump:
    input:  assoc_dir = rules.prep_ld_clump_input.output.assoc_dir,
            snp_list  = rules.prep_ld_clump_input.output.snp_list
    output: clump_dir = directory(DS["run_ld_clump"]["clump_dir"]),
            sentinel  = DS["run_ld_clump"]["sentinel"]
    params: bfile_prefix  = config["geno_post_impute"]["vcf_to_plink"]["params"],
            subset_prefix = DS["run_ld_clump"]["subset_prefix"],
            clump_kb = DS["run_ld_clump"]["clump_kb"],
            clump_r2 = DS["run_ld_clump"]["clump_r2"]
    envmodules: "PLINK"
    message: "LD-clump multi-SNP genes' eSNPs against cohort genotypes (independence check)"
    benchmark: DS["run_ld_clump"]["benchmark"]
    log:    DS["run_ld_clump"]["log"]
    shell:
        """
        mkdir -p {output.clump_dir}
        mkdir -p $(dirname {params.subset_prefix})

        echo "=== Subsetting genome-wide bfile to candidate SNPs ===" > {log}
        plink --bfile {params.bfile_prefix} \
              --extract {input.snp_list} \
              --make-bed \
              --out {params.subset_prefix} >> {log} 2>&1

        n_requested=$(wc -l < {input.snp_list})
        n_extracted=$(wc -l < {params.subset_prefix}.bim)
        echo "Candidate SNPs requested: $n_requested" >> {log}
        echo "Candidate SNPs found in genotype data: $n_extracted" >> {log}
        if [ "$n_requested" -ne "$n_extracted" ]; then
            echo "WARNING: $((n_requested - n_extracted)) candidate SNPs missing from genotype data" >> {log}
        fi

        echo "" >> {log}
        echo "=== Running per-gene clump (kb={params.clump_kb}, r2={params.clump_r2}) ===" >> {log}

        n_genes=0
        for assoc_file in {input.assoc_dir}/*.assoc; do
            gene_id=$(basename "$assoc_file" .assoc)
            plink --bfile {params.subset_prefix} \
                  --clump "$assoc_file" \
                  --clump-p1 1 --clump-p2 1 \
                  --clump-r2 {params.clump_r2} \
                  --clump-kb {params.clump_kb} \
                  --clump-field P \
                  --clump-snp-field SNP \
                  --out {output.clump_dir}/$gene_id >> {log} 2>&1
            n_genes=$((n_genes + 1))
        done

        n_clumped=$(find {output.clump_dir} -maxdepth 1 -name '*.clumped' 2>/dev/null | wc -l)
        echo "" >> {log}
        echo "Genes processed: $n_genes" >> {log}
        echo "Clump output files produced: $n_clumped" >> {log}
        if [ "$n_genes" -ne "$n_clumped" ]; then
            echo "WARNING: $((n_genes - n_clumped)) genes produced no .clumped output (check log above for plink errors/warnings on those genes)" >> {log}
            echo "Not writing sentinel -- treat this run as incomplete." >> {log}
            exit 1
        fi

        touch {output.sentinel}
        """

rule parse_ld_clump:
    input:  clump_dir             = rules.run_ld_clump.output.clump_dir,
            run_ld_clump_sentinel = rules.run_ld_clump.output.sentinel,
            singleton_snps        = rules.prep_ld_clump_input.output.singleton_snps
    output: DS["parse_ld_clump"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Parse LD clump results into final pruned, LD-independent eGene-eSNP universe"
    benchmark: DS["parse_ld_clump"]["benchmark"]
    log:    DS["parse_ld_clump"]["log"]
    script: DS["parse_ld_clump"]["script"]

rule ocr_enrich:
    input:  classification     = rules.classify_eqtl_specificity.output[0],
            heterogeneity_test = rules.heterogeneity_test.output[0],
            independent_snps   = rules.parse_ld_clump.output[0],
            snp_lookups        = snp_lookup_files(DEVSPEC_CELL_TYPES),
            li_s6              = DS["ocr_enrich"]["li_s6"],
            li_s7              = DS["ocr_enrich"]["li_s7"]
    output: DS["ocr_enrich"]["out_file"]
    params: seed = DS["ocr_enrich"]["seed"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 16000, time = "01:00:00"
    threads: 4
    message: "Test fetal-specific/shared/heterogeneous eQTL enrichment in adult cortical union OCRs"
    benchmark: DS["ocr_enrich"]["benchmark"]
    log:    DS["ocr_enrich"]["log"]
    script: DS["ocr_enrich"]["script"]

rule glm_enrich:
    input:  classification   = rules.classify_eqtl_specificity.output[0],
            independent_snps = rules.parse_ld_clump.output[0],
            snp_lookups      = snp_lookup_files(DEVSPEC_CELL_TYPES),
            li_s6            = DS["ocr_enrich"]["li_s6"],
            li_s7            = DS["ocr_enrich"]["li_s7"]
    output: DS["glm_enrich"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 2, mem_mb = 8000, time = "00:20:00"
    threads: 2
    message: "GLM comparison: are fetal-specific eQTL less concentrated than shared eQTL in Li et al. 2023 adult cortical OCRs?"
    benchmark: DS["glm_enrich"]["benchmark"]
    log:    DS["glm_enrich"]["log"]
    script: DS["glm_enrich"]["script"]

rule eqtl_effect_sizes:
    input:  eqtl_universe = rules.extract_eqtl_universe.output[0],
            nom_files     = smr_nom_files(DEVSPEC_CELL_TYPES),
            jang_l1_files = expand(JANG_FULL, jang_cell_type=config["jang_l1_types"]),
            gene_lookup   = config["gene_lookup"]
    output: DS["eqtl_effects"]["out_file"]
    params: cell_types    = DEVSPEC_CELL_TYPES,
            cell_types_L1 = DEVSPEC_FETAL_L1,
            jang_l1_types = config["jang_l1_types"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 64000, time = "02:00:00"
    threads: 4
    message: "Extract effect sizes for unique eGene-eSNP pairs across fetal and adult cell types"
    benchmark: DS["eqtl_effects"]["benchmark"]
    log:    DS["eqtl_effects"]["log"]
    script: DS["eqtl_effects"]["script"]

rule eqtl_effect_replication:
    input:  qtl_perm_files = perm_files(DEVSPEC_CELL_TYPES),
            nom_files      = smr_nom_files(DEVSPEC_CELL_TYPES),
            jang_files     = expand(JANG_FULL, jang_cell_type=config["jang_l1_types"]),
            pvar_file      = config["tensorQTL"]["convert_genotypes"]["pvar"]
    output: DS["eqtl_effect_replication"]["out_file"]
    params: cell_types = DEVSPEC_CELL_TYPES,
            fetal_l1   = DEVSPEC_FETAL_L1,
            jang_order = config["jang_l1_order"],
            fdr_thresh = config["eqtl_fdr"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 64000, time = "02:00:00"
    threads: 4
    message: "Discovery -> replication eQTL effect-size correlations (fetal-fetal and fetal L1-adult), no significance cutoff on the replication side"
    benchmark: DS["eqtl_effect_replication"]["benchmark"]
    log:    DS["eqtl_effect_replication"]["log"]
    script: DS["eqtl_effect_replication"]["script"]

# Benchmark files of the current rules (read by the report if present)
DEVSPEC_BENCHMARKS = [(rule, DS[key]["benchmark"]) for rule, key in [
    ("extract_eqtl_universe", "extract_universe"),
    ("classify_eqtl_specificity", "classify_specificity"),
    ("heterogeneity_test", "heterogeneity_test"),
    ("prep_ld_clump_input", "prep_ld_clump_input"),
    ("run_ld_clump", "run_ld_clump"),
    ("parse_ld_clump", "parse_ld_clump"),
    ("ocr_enrich", "ocr_enrich"),
    ("glm_enrich", "glm_enrich"),
    ("eqtl_effect_sizes", "eqtl_effects"),
    ("eqtl_effect_replication", "eqtl_effect_replication"),
    ("dev_specificity_report", "report"),
]]


rule dev_specificity_report:
    input:  eqtl_universe           = rules.extract_eqtl_universe.output[0],
            classification          = rules.classify_eqtl_specificity.output[0],
            heterogeneity_test      = rules.heterogeneity_test.output[0],
            ocr_enrich              = rules.ocr_enrich.output[0],
            glm_enrich              = rules.glm_enrich.output[0],
            eqtl_effects            = rules.eqtl_effect_sizes.output[0],
            eqtl_effect_replication = rules.eqtl_effect_replication.output[0],
            gene_lookup             = config["gene_lookup"],
            rmd_script              = DS["report"]["script"]
    output: DS["report"]["html"]
    params: eqtl_universe           = lambda wc, input: "../" + input.eqtl_universe,
            classification          = lambda wc, input: "../" + input.classification,
            heterogeneity_test      = lambda wc, input: "../" + input.heterogeneity_test,
            ocr_enrich              = lambda wc, input: "../" + input.ocr_enrich,
            glm_enrich              = lambda wc, input: "../" + input.glm_enrich,
            eqtl_effects            = lambda wc, input: "../" + input.eqtl_effects,
            eqtl_effect_replication = lambda wc, input: "../" + input.eqtl_effect_replication,
            gene_lookup             = lambda wc, input: "../" + input.gene_lookup,
            fetal_l1    = ",".join(DEVSPEC_FETAL_L1),
            jang_l1     = ",".join(config["jang_l1_order"]),
            eqtl_fdr    = config["eqtl_fdr"],
            fdr_thresh  = DS["heterogeneity_test"]["fdr_thresh"],
            clump_kb    = DS["run_ld_clump"]["clump_kb"],
            clump_r2    = DS["run_ld_clump"]["clump_r2"],
            seed        = DS["ocr_enrich"]["seed"],
            benchmarks  = ",".join(f"{rule}=../{path}" for rule, path in DEVSPEC_BENCHMARKS),
            output_file = DS["report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 64000, time = "01:00:00"
    message: "Generate eQTL sharing report"
    benchmark: DS["report"]["benchmark"]
    log:    DS["report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(eqtl_universe = '{params.eqtl_universe}', \
                          classification = '{params.classification}', \
                          heterogeneity_test = '{params.heterogeneity_test}', \
                          ocr_enrich = '{params.ocr_enrich}', \
                          glm_enrich = '{params.glm_enrich}', \
                          eqtl_effects = '{params.eqtl_effects}', \
                          eqtl_effect_replication = '{params.eqtl_effect_replication}', \
                          gene_lookup = '{params.gene_lookup}', \
                          fetal_l1 = '{params.fetal_l1}', \
                          jang_l1 = '{params.jang_l1}', \
                          eqtl_fdr = '{params.eqtl_fdr}', \
                          fdr_thresh = '{params.fdr_thresh}', \
                          clump_kb = '{params.clump_kb}', \
                          clump_r2 = '{params.clump_r2}', \
                          seed = '{params.seed}', \
                          benchmarks = '{params.benchmarks}'))" > {log} 2>&1
        """
