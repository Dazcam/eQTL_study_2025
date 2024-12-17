import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import gc


def load_and_dwnsmple_anndata(ann_dir, n_cells):
    """
    Load an AnnData object from the specified directory and downsample it to a fixed number of cells.

    Parameters:
    ann_dir (str): The directory containing the AnnData file (assumed to be named 'anndata.h5ad').
    n_cells (int): The number of cells to retain in the downsampled AnnData object.

    Returns:
    AnnData: A downsampled AnnData object.

    Raises:
    FileNotFoundError: If the specified AnnData file does not exist in the directory.
    ValueError: If n_cells is greater than the total number of cells in the AnnData object.
    """
    try:
        # Attempt to load the AnnData object
        adata_obj = sc.read(ann_dir + 'anndata.h5ad')
    except FileNotFoundError:
        raise FileNotFoundError(f"The file 'anndata.h5ad' was not found in the directory: {ann_dir}")

    if n_cells > adata_obj.n_obs:
        raise ValueError(f"n_cells ({n_cells}) cannot be greater than the total number of cells in the AnnData object ({adata_obj.n_obs}).")

    # Perform downsampling
    random_indices = np.random.choice(adata_obj.n_obs, n_cells, replace=False)
    adata_obj = adata_obj[random_indices, :].copy()

    # Trigger garbage collection to free unused memory
    gc.collect()

    return adata_obj

def create_counts_per_sample_boxplt(ann_obj):
    """
    Create a two-panelled boxplot visualizing total counts and feature counts per sample.

    Parameters
    ----------
    ann_obj : anndata.AnnData
        An AnnData object containing the following columns in `obs`:
        - 'tscp_count': Total counts per cell.
        - 'gene_count': Number of genes per cell.
        - 'sample': Sample labels for each cell.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plots.
    axes : numpy.ndarray
        The array of axes objects for the plots.

    Raises
    ------
    ValueError
        If required columns ('tscp_count', 'gene_count', 'sample') are not present in `obs`.
    """
    required_columns = ['tscp_count', 'gene_count', 'sample']
    
    # Error check: Ensure required columns are in `obs`
    if not all(col in ann_obj.obs.columns for col in required_columns):
        missing = [col for col in required_columns if col not in ann_obj.obs.columns]
        raise ValueError(f"Missing required columns in `obs`: {missing}")

    # Create figure and axes for subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 12), gridspec_kw={'wspace': 0.2})

    # Prepare data for first plot (total counts)
    tscp_data = ann_obj.obs[['tscp_count', 'gene_count', 'sample']].copy()
    tscp_data['sample'] = tscp_data['sample'].astype('category')

    # Plot total counts
    sns.boxplot(
        x=np.log10(tscp_data['tscp_count']), y='sample', data=tscp_data,
        linewidth=0.7, palette="Set2", ax=axes[0], fliersize=2, hue='sample',
        dodge=False
    )
    axes[0].axvline(x=np.log10(500), color='red', linestyle='dotted')
    axes[0].axvline(x=np.log10(1000), color='red', linestyle='dotted')
    axes[0].set_xlabel('Log10 Total Counts')
    axes[0].set_ylabel('')
    axes[0].set_title('Tscp Counts per Sample')

    # Prepare data for second plot (feature counts)
    gene_data = ann_obj.obs[['gene_count', 'sample']].copy()

    # Order sample labels correctly
    gene_data['sample_order'] = gene_data['sample'].str.extract(r'(\d+)')[0].astype(int)
    gene_data['fc_flag'] = gene_data['sample'].str.contains('_FC$')
    gene_data = gene_data.sort_values(by=['fc_flag', 'sample_order'], ascending=[False, True])
    gene_data['sample'] = pd.Categorical(gene_data['sample'], categories=gene_data['sample'].unique(), ordered=True)

    # Plot feature counts
    sns.boxplot(
        x=np.log10(gene_data['gene_count']), y='sample', data=gene_data,
        linewidth=0.7, palette="Set2", ax=axes[1], fliersize=2, hue='sample',
        dodge=False
    )
    axes[1].axvline(x=np.log10(500), color='red', linestyle='dotted')
    axes[1].axvline(x=np.log10(1000), color='red', linestyle='dotted')
    axes[1].set_xlabel('Log10 Feature Counts')
    axes[1].set_ylabel('')
    axes[1].set_title('Gene Counts per Sample')

    # Return figure and axes objects
    return fig, axes

