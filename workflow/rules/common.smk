# Shared settings and helpers for every rule file (included first by the Snakefile)


##  Settings  ------------------------------------------------------------------------------
NORM_METHOD = config["tensorQTL"]["norm_methods"][0]
GENO_PC     = config["tensorQTL"]["geno_pcs"]

# Jang et al. 2026 association files (downloaded by setup get_jang_singlebrain)
JANG_FULL = config["setup"]["get_jang_singlebrain"]["full_assoc"]
JANG_TOP  = config["setup"]["get_jang_singlebrain"]["top_assoc"]


##  Functions  -----------------------------------------------------------------------------
def perm_files(cell_types):
    """tensorQTL permutation output per cell type, at the final expPC (exp_pc_map)."""
    return [config["tensorQTL"]["tensorqtl_perm"]["output"].format(
                cell_type=ct, norm_method=NORM_METHOD, geno_pc=GENO_PC,
                exp_pc=config["exp_pc_map"][ct])
            for ct in cell_types]


def manuscript_label(ct):
    """Pipeline -> manuscript cell type label (config manuscript_labels, in order)."""
    for old, new in config["manuscript_labels"].items():
        ct = ct.replace(old, new)
    return ct
