configfile: '../config/config.yaml'

#---  Need to handle cell type name changes for V2
def relabel_cell_type(x):
    return x.replace("Glu-UL", "Glu-A").replace("Glu-DL", "Glu-B").replace("-Q4-", "-")

# manuscript label -> pipeline label (used by any rule that needs to go from the
# manuscript/figshare cell type name back to the name used on disk everywhere upstream)
CT_NEW2OLD = {relabel_cell_type(ct): ct for ct in config["cell_types"]}
assert len(CT_NEW2OLD) == len(config["cell_types"]), "Relabelling produced duplicate cell type names"
assert not any(k for k in CT_NEW2OLD if any(s in k for s in ("UL", "DL", "Q4"))), "Old names remain"

susie_data_files = expand(config['data_sharing']['data_susie']['out_file'],
                          cell_type = list(CT_NEW2OLD))


# ----


rule data_perm:
    output: config['data_sharing']['data_perm']['out_file']
    params: in_dir           = config['data_sharing']['data_perm']['in_dir'],
            allele_file      = config['data_sharing']['data_perm']['allele_file'],
            cell_types       = config['cell_types'],
            exp_pc_map       = config['exp_pc_map'],
            norm_method      = config['tensorQTL']['norm_methods'][0],
            genPC            = config['tensorQTL']['geno_pcs'],
            gene_lookup_file = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"
    log:  config['data_sharing']['data_perm']['log']
    script: "../scripts/data_sharing_eQTL_perm.R"

rule data_nominal:
    output: config['data_sharing']['data_nominal']['out_file']
    params: in_dir           = config['data_sharing']['data_nominal']['in_dir'],
            allele_file      = config['data_sharing']['data_nominal']['allele_file'],
            cell_types       = config['cell_types'],
            gene_lookup_file = config['gene_lookup']
    singularity: config["containers"]["r_eqtl"]
    resources: time="4:00:00", threads = 10, mem_mb = 80000
    log:  config['data_sharing']['data_nominal']['log']
    script: "../scripts/data_sharing_eQTL_nominal.R"

rule mk_eqtl_tar:
    input:  perm = config['data_sharing']['data_perm']['out_file'],
            nom = config['data_sharing']['data_nominal']['out_file']
    output: config['data_sharing']['data_mk_eqtl_tar']['out_file']
    params: out_dir    = config['data_sharing']['data_mk_eqtl_tar']['out_dir'],
            cell_types  = list(CT_NEW2OLD)     # manuscript labels: Glu-A, Glu-B, NPC-to-Glu-A-Bin1, ...
    envmodules: "pigz"
    resources: time="5:00:00", threads = 8, mem_mb = 16000
    log:  config['data_sharing']['data_mk_eqtl_tar']['log']
    shell: """
           scripts/data_sharing_mk_eqtl_tar.sh \
           {params.out_dir} \
           {threads} \
           "{params.cell_types}" \
           >> {log} 2>&1
           """

rule data_weights:
    output: config['data_sharing']['data_weights']['out_file']
    params: out_dir = config['data_sharing']['data_weights']['out_dir'],
            file_name = config['data_sharing']['data_weights']['file_name'],
            cell_types = config['cell_types']
    envmodules: "pigz"
    resources: time="3:00:00", threads = 10, mem_mb = 80000
    log:  config['data_sharing']['data_weights']['log']
    shell:
           """
           scripts/data_sharing_twas_weights.sh \
             {params.out_dir} \
             {params.file_name} \
             "{params.cell_types}" \
             {threads} \
            >> {log} 2>&1
           """

rule final_genotypes:
    input:  vcf     = config["geno_post_impute"]["exclude_SNPs"]["output"],
            samples = config["data_sharing"]["final_genotypes"]["samples"]
    output: config["data_sharing"]["final_genotypes"]["out_file"]
    envmodules: "BCFtools"
    resources: time="30:00"
    log:    config["data_sharing"]["final_genotypes"]["log"]
    shell: """
           bcftools view \
             --samples-file {input.samples} \
             --output-type z \
             --output {output} \
             {input.vcf} 2>> {log}
           bcftools index --tbi {output} 2>> {log}
           """

rule data_susie:
    # {cell_type} is the manuscript label; the input uses the pipeline label
    input:  
        susie = lambda wc: config["susie"]["sort_susie"]["output"].replace(
                        "{cell_type}", CT_NEW2OLD[wc.cell_type]),
        gene_lookup = config['gene_lookup']
    output: tsv = config['data_sharing']['data_susie']['out_file']
    params: cell_type_old = lambda wc: CT_NEW2OLD[wc.cell_type]
    singularity: config["containers"]["r_eqtl"]
    resources: time="0:20:00", threads=2, mem_mb=16000
    log:    config['data_sharing']['data_susie']['log']
    script: "../scripts/data_sharing_susie_cred_sets.R"

rule data_susie_mk_tar:
    input:  susie_data_files
    output: config['data_sharing']['data_susie_mk_tar']['out_file']
    params: in_dir  = lambda wc, input: os.path.dirname(input[0]),
            top_dir = config['data_sharing']['data_susie_mk_tar']['top_dir'],
            n_files = lambda wc, input: len(input)
    resources: time="0:30:00", threads=1, mem_mb=4000
    log:    config['data_sharing']['data_susie_mk_tar']['log']
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
