# Retired from 18data_sharing.smk (Nov 2026); config block data_final_scanpy_obj commented out.

rule final_scanpy_obj:
    # L1 object + the four L2 subclustering objects + one pseudotime object per trajectory
    input:
        l1 = config['scanpy']['clusters_h5ad'],
        l2 = expand(config['data_sharing']['final_scanpy_obj']['l2_h5ad'],
                    l1_pop = config['data_sharing']['final_scanpy_obj']['l2_populations']),
        pt = expand(config['pseudotime']['palantir']['h5ad'],
                    trajectory = config['trajectories'])
    output:
        h5ad    = config['data_sharing']['final_scanpy_obj']['out_h5ad'],
        mapping = config['data_sharing']['final_scanpy_obj']['out_mapping']
    params:
        l1_cluster_col        = config['data_sharing']['final_scanpy_obj']['l1_cluster_col'],
        l2_populations        = config['data_sharing']['final_scanpy_obj']['l2_populations'],   # same order as input.l2
        l1_without_l2         = config['data_sharing']['final_scanpy_obj']['l1_without_l2'],
        pt_trajectories       = config['trajectories'],   # same order as input.pt
        pt_n_bins             = config['data_sharing']['final_scanpy_obj']['pt_n_bins'],
        keep_obs              = config['data_sharing']['final_scanpy_obj']['keep_obs'],
        keep_obsm             = config['data_sharing']['final_scanpy_obj']['keep_obsm'],
        x_from_counts         = config['data_sharing']['final_scanpy_obj']['x_from_counts'],
        manuscript_cell_types = list(CT_NEW2OLD)
    conda: config["scanpy"]["env"]
    resources: time="1:00:00", mem_mb=64000
    log: config['data_sharing']['final_scanpy_obj']['log']
    script: "../scripts/data_sharing_final_scanpy_obj.py"
