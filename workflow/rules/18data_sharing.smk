import os


##  Settings  ------------------------------------------------------------------------------
DSH = config['data_sharing']

# Manuscript label -> pipeline label, in config cell_types order
CT_NEW2OLD = {manuscript_label(ct): ct for ct in config['cell_types']}
assert len(CT_NEW2OLD) == len(config['cell_types']), "Relabelling produced duplicate cell type names"
MS_LABELS = list(CT_NEW2OLD)


##  Functions  -----------------------------------------------------------------------------
def nominal_files(cell_types):
    """12SMR nominal export (cat_tensorqtl_nom_snps) per cell type."""
    return [config['smr']['cat_tensorqtl_nom_snps']['cat_snps'].format(cell_type=ct)
            for ct in cell_types]


def weights_files(cell_types):
    """14CTWAS copy_fusion_weights done-file per cell type (weights copied to dest_dir)."""
    return [config['ctwas']['copy_fusion_weights']['output'].format(cell_type=ct)
            for ct in cell_types]


def susie_file(wc):
    """10SUSIE sort_susie output for the pipeline label of a manuscript-label wildcard."""
    return config['susie']['sort_susie']['output'].format(cell_type=CT_NEW2OLD[wc.cell_type])


def interleave(a, b):
    """[a1, b1, a2, b2, ...]: perm then nominal file for each cell type."""
    return [f for pair in zip(a, b) for f in pair]


##  Rules  ---------------------------------------------------------------------------------

rule data_perm:
    input:  eqtl        = perm_files(config['cell_types']),
            pvar        = config['tensorQTL']['convert_genotypes']['pvar'],
            gene_lookup = config['gene_lookup']
    output: tsv      = expand(DSH['data_perm']['out_file'], cell_type=MS_LABELS),
            summary  = DSH['data_perm']['summary'],
            dup_snps = DSH['data_perm']['dup_snps']
    params: cell_types = config['cell_types'],
            labels     = config['manuscript_labels']
    singularity: config['containers']['r_eqtl']
    resources: time = "1:00:00"
    message:   "eQTL permutation files for data sharing"
    benchmark: DSH['data_perm']['benchmark']
    log:       DSH['data_perm']['log']
    script:    DSH['data_perm']['script']

rule data_nominal:
    input:  eqtl        = nominal_files(config['cell_types']),
            pvar        = config['tensorQTL']['convert_genotypes']['pvar'],
            gene_lookup = config['gene_lookup']
    output: tsv      = expand(DSH['data_nominal']['out_file'], cell_type=MS_LABELS),
            summary  = DSH['data_nominal']['summary'],
            dup_snps = DSH['data_nominal']['dup_snps']
    params: cell_types = config['cell_types'],
            labels     = config['manuscript_labels']
    singularity: config['containers']['r_eqtl']
    resources: time = "4:00:00", mem_mb = 80000
    message:   "eQTL nominal files for data sharing"
    benchmark: DSH['data_nominal']['benchmark']
    log:       DSH['data_nominal']['log']
    script:    DSH['data_nominal']['script']

rule mk_eqtl_tar:
    input:  interleave(rules.data_perm.output.tsv, rules.data_nominal.output.tsv)
    output: DSH['data_mk_eqtl_tar']['out_file']
    params: script = DSH['data_mk_eqtl_tar']['script']
    envmodules: "pigz"
    threads:   8
    resources: threads = 8, mem_mb = 16000, time = "5:00:00"
    message:   "Tarball of eQTL permutation and nominal files"
    benchmark: DSH['data_mk_eqtl_tar']['benchmark']
    log:       DSH['data_mk_eqtl_tar']['log']
    shell:  "{params.script} {output} {threads} {input} > {log} 2>&1"

rule data_weights:
    input:  weights_files(config['cell_types'])
    output: DSH['data_weights']['out_file']
    params: script       = DSH['data_weights']['script'],
            weights_dir  = os.path.dirname(config['ctwas']['copy_fusion_weights']['dest_dir']),
            archive_root = DSH['data_weights']['archive_root'],
            pairs        = " ".join(f"{old}={new}" for new, old in CT_NEW2OLD.items())
    envmodules: "pigz"
    threads:   10
    resources: threads = 10, mem_mb = 80000, time = "3:00:00"
    message:   "Tarball of FUSION TWAS weights"
    benchmark: DSH['data_weights']['benchmark']
    log:       DSH['data_weights']['log']
    shell:  "{params.script} {output} {threads} {params.weights_dir} "
            "{params.archive_root} {params.pairs} > {log} 2>&1"

rule final_genotypes:
    input:  vcf     = config['geno_post_impute']['exclude_SNPs']['output'],
            samples = DSH['final_genotypes']['samples']
    output: DSH['final_genotypes']['out_file']
    envmodules: "BCFtools"
    resources: time = "30:00"
    message:   "Genotype VCF for the final samples"
    benchmark: DSH['final_genotypes']['benchmark']
    log:       DSH['final_genotypes']['log']
    shell:  """
            bcftools view \
              --samples-file {input.samples} \
              --output-type z \
              --output {output} \
              {input.vcf} 2> {log}
            bcftools index --tbi {output} 2>> {log}
            """

rule data_susie:
    # {cell_type} is the manuscript label; the input uses the pipeline label
    input:  susie       = susie_file,
            gene_lookup = config['gene_lookup']
    output: tsv = DSH['data_susie']['out_file']
    params: cell_type_old = lambda wc: CT_NEW2OLD[wc.cell_type],
            labels        = config['manuscript_labels']
    singularity: config['containers']['r_eqtl']
    resources: time = "0:20:00", mem_mb = 16000
    message:   "SuSiE credible sets for data sharing: {wildcards.cell_type}"
    benchmark: DSH['data_susie']['benchmark']
    log:       DSH['data_susie']['log']
    script:    DSH['data_susie']['script']

rule data_susie_mk_tar:
    input:  expand(rules.data_susie.output.tsv, cell_type=MS_LABELS)
    output: DSH['data_susie_mk_tar']['out_file']
    params: in_dir  = lambda wc, input: os.path.dirname(input[0]),
            top_dir = DSH['data_susie_mk_tar']['top_dir'],
            n_files = lambda wc, input: len(input)
    resources: time = "0:30:00", mem_mb = 4000
    message:   "Tarball of SuSiE credible sets"
    benchmark: DSH['data_susie_mk_tar']['benchmark']
    log:       DSH['data_susie_mk_tar']['log']
    shell:  """
            set -euo pipefail
            echo "Creating {output} from {params.n_files} files" > {log}
            printf '%s\\n' {input} | xargs -n1 basename > {output}.filelist
            tar --create --file {output} \\
                --directory {params.in_dir} \\
                --transform 's,^,{params.top_dir}/,' \\
                --files-from {output}.filelist \\
                >> {log} 2>&1
            rm {output}.filelist
            tar --list --file {output} | wc -l >> {log}
            echo "Done" >> {log}
            """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["18data_sharing"] = [rules.mk_eqtl_tar.output,
                             rules.data_susie_mk_tar.output,
                             rules.data_weights.output,
                             rules.final_genotypes.output]
