import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import gc
import re
import ipywidgets as widgets
from IPython import get_ipython
from IPython.display import display, Markdown
from ipywidgets import Tab, Output
from memory_profiler import profile
import warnings
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
import seaborn as sns
from matplotlib import MatplotlibDeprecationWarning
import logging
from sklearn.decomposition import PCA
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import math

##########  LOAD FUNCTIONS  ############  
#@profile
def load_and_dwnsmpl_data(downsample_cells=None, is_anndata=False, *adata_dirs):
    """
    Load, process, and merge AnnData objects with optional downsampling.
    
    Parameters:
        downsample_cells (int, optional): Number of cells to randomly downsample to. If None, no downsampling is performed.
        *adata_dirs (str): Full paths to AnnData files. Each path should point to an .h5ad file.
        
    Returns:
        AnnData: A merged and processed AnnData object, or the single input AnnData if only one is provided.
    """
    # If no adata directories are provided, raise an error
    if len(adata_dirs) == 0:
        raise ValueError("At least one adata directory must be provided.")
    
    adata_list = []   
    
    # Extract plate numbers from the file paths
    plate_numbers = [re.search(r'plate(\d+)', adata_dir).group(1) for adata_dir in adata_dirs]

    # Report how many plates will be processed and downsampling details
    if downsample_cells is not None:
        print(f"Processing {len(adata_dirs)} plate(s) with downsampling to {downsample_cells} cells per plate.")
    else:
        print(f"Processing {len(adata_dirs)} plate(s), no downsampling applied.")

    for i, adata_dir in enumerate(adata_dirs):
        # Load the AnnData object directly from the full file path
        print(f"Loading plate {plate_numbers[i]} from {adata_dir} ...")

        if is_anndata:
            if adata_dir.endswith('.h5ad'):
                # Case 1: Full path to .h5ad file
                adata = sc.read(adata_dir)
            else:
                # Case 2: Directory path, assume 'anndata.h5ad' inside it
                h5ad_path = os.path.join(adata_dir, 'anndata.h5ad')
                print(f"Loading {h5ad_path} ...")
                adata = sc.read(h5ad_path)
        else:
            # Case 3: Parse format, load .mtx and metadata
            print(f"Loading Parse mtx, genes and metadata separately ...")
            adata = sc.read_mtx(os.path.join(adata_dir, 'count_matrix.mtx'))
            
            # Read gene and cell metadata
            gene_data = pd.read_csv(os.path.join(adata_dir, 'all_genes.csv'))
            cell_meta = pd.read_csv(os.path.join(adata_dir, 'cell_metadata.csv'))
            
            # Filter out NaN genes
            gene_data = gene_data[gene_data.gene_name.notnull()]
            notNa = gene_data.index
            notNa = notNa.to_list()
            
            # Set gene names
            adata = adata[:, notNa]
            adata.var = gene_data
            adata.var.set_index('gene_name', inplace=True)
            adata.var.index.name = None
            adata.var_names_make_unique()
            
            # Add cell metadata
            adata.obs = cell_meta
            adata.obs.set_index('bc_wells', inplace=True)
            adata.obs.index.name = None
            adata.obs_names_make_unique()
            
        # Check if 'plate' column exists and add it only if not present
        if 'plate' not in adata.obs.columns and len(adata_dirs) > 1:
            adata.obs['plate'] = 'plate' + plate_numbers[i]
        
        # Avoid duplicating plate prefix in obs_names
        if not adata.obs_names.str.startswith(f"plate{plate_numbers[i]}_").any():
            adata.obs_names = f"plate{plate_numbers[i]}_" + adata.obs_names

        # Log dimensions
        print(f"Plate {plate_numbers[i]} dimensions: {adata.shape}")
        print(f"Plate {plate_numbers[i]} matrix dimensions: {adata.X.shape}")
 
        # Optional downsampling
        if downsample_cells is not None and downsample_cells < adata.n_obs:
            random_indices = np.random.choice(adata.n_obs, downsample_cells, replace=False)
            adata = adata[random_indices, :].copy()
        
        adata_list.append(adata)

    # If more than one adata object, find common genes across all datasets
    if len(adata_list) > 1:
        common_genes = adata_list[0].var_names
        for adata in adata_list[1:]:
            common_genes = common_genes.intersection(adata.var_names)
        
        print(f"Common genes across all plates: {common_genes.shape[0]}")
        
        # Subset all adata objects for common genes
        for i, adata in enumerate(adata_list):
            adata_list[i] = adata[:, common_genes]
    else:
        print(f"Only one plate processed; no need to find common genes.")
        
    # Concatenate the AnnData objects if more than one
    if len(adata_list) > 1:
        print(f"Merging plates ...")
        adata_mrg = sc.concat(adata_list, join='inner')
    else:
        adata_mrg = adata_list[0]
    
    # Remove 'hg38' from gene names
    adata_mrg.var.index = adata_mrg.var.index.str.replace('_hg38', '')

    # Check if adata.X is sparse and if all values are integers
    # Slice a small portion of the sparse matrix (e.g., first 1000 rows and 1000 columns)
    subset = adata_mrg.X[:1000, :1000]
    is_all_integers = np.all(np.mod(subset.data, 1) == 0)
    print(f'Are all values of adata.X integers? {is_all_integers}')

    # Cleanup to release memory
    del adata_list
    gc.collect()

    print(adata_mrg)

    return adata_mrg

##########  QC FUNCTIONS  ############  
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

def get_raw_counts_all_genes(adata):
    """
    Returns the raw count matrix (all genes) from `adata.raw` or `adata.layers["counts"]`, ensuring
    that it contains integer values. If neither contains integers, raises an error.
    
    Args:
    - adata: AnnData object containing the raw counts or counts in layers.
    
    Returns:
    - A sparse matrix with raw counts in integer format, for all genes.
    """
    # Check if adata.raw exists and contains integer values
    if adata.raw is not None:
        # Check a small subset to see if data is integer type (avoiding unnecessary copying)
        raw_is_int = np.issubdtype(adata.raw.X[:10, :10].dtype, np.integer)
    else:
        raw_is_int = False

    # Check if adata.layers["counts"] exists and contains integer values
    counts_is_int = False
    if "counts" in adata.layers:
        counts_is_int = np.issubdtype(adata.layers["counts"][:10, :10].dtype, np.integer)

    # Return the raw counts matrix from the appropriate source
    if raw_is_int:
        print("adata.raw contains integers")
        return adata.raw.X  # No copy, return as-is
    elif counts_is_int:
        print("adata.layers['counts'] contains integers")
        return adata.layers["counts"]  # No copy, return as-is
    else:
        raise ValueError("Neither adata.raw nor adata.layers['counts'] contain integer values!")


##########  Filtering functions  ##########
def filter_anndata(
    adata, 
    min_genes_per_cell=300, 
    max_genes_per_cell=6000,
    min_counts_per_cell=500, 
    min_cells_per_gene=5, 
    min_cells_per_sample=10
):
    """
    Filter an AnnData object based on thresholds for cells, genes, and samples.
    
    Parameters:
    ----------
    adata : AnnData
        The AnnData object containing the single-cell data to be filtered. 
        This object will be modified in place.
    min_genes_per_cell : int, optional
        Minimum number of genes required to keep a cell (default: 300).
    max_genes_per_cell : int, optional
        Maximum number of genes allowed in a cell to keep that cell (default: 6000).
    min_counts_per_cell : int, optional
        Minimum number of total counts (reads) required to keep a cell (default: 500).
    min_cells_per_gene : int, optional
        Minimum number of cells expressing a gene to keep that gene (default: 5).
    min_cells_per_sample : int, optional
        Minimum number of cells required for a sample to be retained (default: 10).
    
    Returns:
    -------
    None
        The `adata` object is modified in place.
    """
    print(f"Dimensions before applying filter: {adata.shape}")

    if min_genes_per_cell is not None:
        # Step 1: Filter cells based on minimum number of genes
        num_cells_before = adata.shape[0]
        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        num_cells_after = adata.shape[0]

    if max_genes_per_cell is not None:
        # Step 2: Filter cells based on maximum number of genes
        num_cells_before = adata.shape[0]
        sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)
        num_cells_after = adata.shape[0]

    if min_counts_per_cell is not None:
        # Step 3: Filter cells based on minimum total counts (reads)
        num_cells_before = adata.shape[0]
        sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
        num_cells_after = adata.shape[0]

    if min_cells_per_gene is not None:
        # Step 4: Filter genes based on minimum number of cells expressing them
        num_genes_before = adata.shape[1]
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        num_genes_after = adata.shape[1]
    
    if min_cells_per_sample is not None:
        # Step 5: Filter samples based on minimum number of cells
        num_samples_before = adata.obs['sample'].nunique()
        valid_samples = adata.obs['sample'].value_counts() >= min_cells_per_sample
        valid_sample_ids = valid_samples[valid_samples].index
        valid_cells = adata.obs['sample'].isin(valid_sample_ids)  # Boolean mask for valid cells
        adata._inplace_subset_obs(valid_cells)  # In-place filtering of cells
        num_samples_after = adata.obs['sample'].nunique()
        print(f"filtered out {num_samples_before - num_samples_after} samples with < {min_cells_per_sample} cells.")

    print(f"Dimensions after applying filter: {adata.shape}")


def detect_mad_outliers_per_sample(adata, group_column, target_column, threshold=3, log=False, use_median=True, outlier_column="mad_outlier"):
    """
    Detect MAD-based outliers per sample group in an AnnData object.
    
    Parameters:
    - adata: AnnData object containing observation data.
    - group_column: Column in `adata.obs` to group data by (e.g., 'sample').
    - target_column: Column in `adata.obs` to calculate outliers on (e.g., 'total_counts').
    - threshold: Number of MADs or standard deviations to define outliers.
    - log: Whether to log-transform the target data before analysis.
    - use_median: If True, use median and MAD; otherwise, use mean and standard deviation.
    - outlier_column: Name of the output column in `adata.obs` to store outlier flags.
    
    Returns:
    - None: Modifies `adata.obs` in place by adding a new column with outlier flags.
    """
    def is_outlier(data):
        if log:
            data = np.log1p(data)
        
        if use_median:
            center = np.median(data)
            deviation = np.median(np.abs(data - center))  # MAD
        else:
            center = np.mean(data)
            deviation = np.std(data)  # Standard deviation
        
        lower_bound = center - threshold * deviation
        upper_bound = center + threshold * deviation
        
        return (data < lower_bound) | (data > upper_bound)
    
    adata.obs[outlier_column] = adata.obs.groupby(group_column)[target_column].transform(
        lambda x: is_outlier(x)
    )


def filter_cells_and_genes(
    adata,
    mito_threshold=5,
    ribo_threshold=5,
    remove_doublets=True,
    remove_mad_outliers=True,
    genes_to_remove=None,
):
    """
    Apply filters to cells and genes in the AnnData object, modifying in place.
    
    Cell filters:
        - Mitochondrial gene percentage > mito_threshold (default: 5%)
        - Ribosomal gene percentage > ribo_threshold (default: 5%)
        - MAD outliers (optional)
        - Predicted doublets (optional)
    
    Gene filters:
        - Remove specified genes
        - Default to removing mitochondrial genes and MALAT1 if not specified

    Parameters:
        adata (AnnData): The AnnData object to filter (modified in place).
        mito_threshold (float, optional): Threshold percentage for mitochondrial genes. Default is 5.
        ribo_threshold (float, optional): Threshold percentage for ribosomal genes. Default is 5.
        remove_doublets (bool, optional): Whether to remove predicted doublets. Default is True.
        remove_mad_outliers (bool, optional): Whether to remove MAD outliers. Default is True.
        genes_to_remove (list, optional): List of gene names (or lists of genes) to remove. Default is None.
    """
    # Error checking for input parameters
    if not isinstance(mito_threshold, (int, float)) or mito_threshold < 0:
        raise ValueError("mito_threshold must be a non-negative number.")
    if not isinstance(ribo_threshold, (int, float)) or ribo_threshold < 0:
        raise ValueError("ribo_threshold must be a non-negative number.")
    if not isinstance(remove_doublets, bool):
        raise ValueError("remove_doublets must be a boolean value.")
    if not isinstance(remove_mad_outliers, bool):
        raise ValueError("remove_mad_outliers must be a boolean value.")
    if genes_to_remove is not None and not isinstance(genes_to_remove, list):
        raise ValueError("genes_to_remove must be a list or None.")
    
    print(f"Applying cell filters with thresholds: mito > {mito_threshold}%, ribo > {ribo_threshold}%")
    
    # Convert 'predicted_doublet' to boolean if needed
    if remove_doublets and 'predicted_doublet' in adata.obs.columns:
        adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)

    # Apply mitochondrial and ribosomal filters
    adata.obs['mito_gt_threshold'] = adata.obs['pct_counts_mt'] > mito_threshold
    adata.obs['ribo_gt_threshold'] = adata.obs['pct_counts_ribo'] > ribo_threshold

    # Determine which filters to apply
    columns_to_check = ['mito_gt_threshold', 'ribo_gt_threshold']
    if remove_doublets:
        columns_to_check.append('predicted_doublet')
    if remove_mad_outliers:
        columns_to_check.append('mad_outlier')

    # Create a boolean column indicating if the cell is an outlier
    adata.obs['is_outlier'] = adata.obs[columns_to_check].astype(bool).any(axis=1)

    # Create a summary table for all filtering criteria
    filter_counts = {
        "mito_gt_threshold": adata.obs['mito_gt_threshold'].sum(),
        "ribo_gt_threshold": adata.obs['ribo_gt_threshold'].sum(),
    }
    if remove_doublets:
        filter_counts["predicted_doublet"] = adata.obs['predicted_doublet'].sum()
    if remove_mad_outliers:
        filter_counts["mad_outlier"] = adata.obs['mad_outlier'].sum()

    filter_counts["is_outlier"] = adata.obs['is_outlier'].sum()

    # Print the summary table
    print("Counts of outliers for removal:")
    for key, value in filter_counts.items():
        print(f"  {key}: {value}")

    # Filter cells in place
    print(f"Dimensions before cell filter: {adata.shape}")
    adata._inplace_subset_obs(adata.obs['is_outlier'] == False)  # Apply cell filter in place
    print(f"Dimensions after cell filter: {adata.shape}")

    # Applying gene filters
    print("Applying gene filters based on specified genes to remove")
    
    # Flatten the genes_to_remove into a single list, if provided
    if genes_to_remove is not None:
        flat_genes_to_remove = set(
            gene for sublist in genes_to_remove for gene in (sublist if isinstance(sublist, list) else [sublist])
        )
    else:
        flat_genes_to_remove = set()

    # Subset genes in place
    genes_to_keep = ~adata.var_names.isin(flat_genes_to_remove)
    print(f"Number of genes to remove: {len(flat_genes_to_remove - set(adata.var_names))}")
    print(f"Number of genes remaining after filter: {genes_to_keep.sum()}")
    adata._inplace_subset_var(genes_to_keep)  # Apply gene filter in place

    # Confirm changes
    print(f"Dimensions after gene filter: {adata.shape}")

def filter_genes_by_read_count(adata, min_reads=10, min_samples=100, inplace=True, verbose=False):
    """
    Optimized function to filter genes by aggregating read counts at the sample level.
    Designed for large datasets with sparse matrix support.
    
    Parameters:
    - adata: AnnData object containing gene expression data.
    - min_reads: Minimum number of reads required per gene per sample.
    - min_samples: Minimum number of samples that must meet the `min_reads` threshold.
    - inplace: If True, filters genes directly in `adata`. If False, returns a mask.
    - verbose: If True, prints detailed debugging information.

    Returns:
    - If `inplace=True`: Modifies `adata` by retaining only the filtered genes.
    - If `inplace=False`: Returns a boolean mask of the genes that meet the criteria.
    """
    if verbose:
        print(f"Dataset dimensions: {adata.shape}")
        print(f"Filtering genes with min_reads={min_reads} and min_samples={min_samples}.")

    # Ensure 'sample' exists in `obs`
    if 'sample' not in adata.obs:
        raise ValueError("The 'sample' column is missing in adata.obs. Please define the sample groupings.")

    # Map sample group indices
    unique_samples, sample_indices = np.unique(adata.obs['sample'], return_inverse=True)
    if verbose:
        print(f"Found {len(unique_samples)} unique samples.")

    # Create a sparse matrix to map samples to rows
    sample_mapping = csr_matrix((np.ones_like(sample_indices), (sample_indices, np.arange(len(sample_indices)))),
                                 shape=(len(unique_samples), adata.shape[0]))

    if verbose:
        print("Sample mapping matrix constructed.")

    # Aggregate read counts by sample using sparse matrix multiplication
    X = adata.X
    if not isinstance(X, csr_matrix):
        X = csr_matrix(X)  # Ensure the matrix is in sparse format
    sample_aggregation = sample_mapping @ X

    if verbose:
        print("Sparse matrix aggregation completed.")

    # Apply filtering criteria
    if verbose:
        print("Applying filtering criteria: minimum reads per sample and minimum samples per gene.")
    gene_mask = (sample_aggregation >= min_reads).sum(axis=0).A1 >= min_samples

    if verbose:
        print(f"Number of genes passing the filter: {gene_mask.sum()} out of {adata.shape[1]} total genes.")

    # Plot gene retention across varying sample thresholds
    if verbose or True:  # Always plot as requested
        thresholds = range(1, len(unique_samples) + 1)
        gene_counts = [(sample_aggregation >= min_reads).sum(axis=0) >= t for t in thresholds]
        gene_counts = [mask.sum() for mask in gene_counts]

        plt.figure(figsize=(8, 6))
        plt.plot(thresholds, gene_counts, marker='o')
        plt.xlabel("Minimum Samples")
        plt.ylabel("Number of Genes Passing Filter")
        plt.title("Gene Retention Across Sample Thresholds")
        plt.grid(True)
        plt.show()

    # If in-place, subset the AnnData object
    if inplace:
        if verbose:
            print("Filtering genes in place.")
        adata._inplace_subset_var(gene_mask)
        if verbose:
            print(f"New dataset dimensions after filtering: {adata.shape}")
    else:
        if verbose:
            print("Returning the gene mask.")
        return gene_mask


########## General functions ##########
def is_running_in_jupyter():
    """
    Determine if the current environment is a Jupyter Notebook.

    Returns:
        bool: True if running in Jupyter, False otherwise.
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or Jupyter Lab
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other types of shells
    except NameError:
        return False  # Standard Python interpreter

##########  VISUALISATION FUNCTIONS  ############  
def plot_celltype_and_gene_features(
    adata,
    final_genes,
    cell_type_column="cell_type",
    figsize_width=9.0,
    base_height_per_row=0.95,
    extra_height=3.0,
    ncols=3,
    cmap_name="light_dark_blue",
    umap_palette=None,
    vmin = 0,
    vmax = 'p99',
):
    """
    Create a publication-style figure with:
      - Left: large UMAP colored by cell type (panel A)
      - Right: grid of feature plots for selected genes (panel B)
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object containing UMAP coordinates and gene expression
    final_genes : list of str
        List of gene names to plot as feature plots
    cell_type_column : str, default: "cell_type"
        Column in adata.obs to use for coloring the main cell type UMAP
    figsize_width : float, default: 9.0
        Width of figure in inches (Nature Genetics double column ~174 mm ≈ 9 inches)
    base_height_per_row : float, default: 0.95
        Height per row of feature plots (tune for your preference)
    extra_height : float, default: 3.0
        Additional height for the top part (title space + main UMAP)
    ncols : int, default: 3
        Number of columns for the gene feature plot grid
    cmap_name : str, default: "light_dark_blue"
        Name of the custom colormap to create (currently only light→dark blue supported)
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure object (you can call plt.show() or save it)
    """
    # ── Custom fetaure plot colormap 
    if cmap_name == "light_dark_blue":
        cmap_name = LinearSegmentedColormap.from_list(
            "light_dark_blue",
            ["#f7fbff", "#08306b"]
        )

    # ── Layout calculations 
    n_genes = len(final_genes)
    nrows = math.ceil(n_genes / ncols)
    
    fig_height = extra_height + (nrows * base_height_per_row)
    
    # Nature Genetics style defaults
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.titlesize": 12,
    })

    # ── Create figure & grid 
    fig = plt.figure(figsize=(figsize_width, fig_height))
    
    gs = gridspec.GridSpec(
        nrows=nrows,
        ncols=ncols + 2,
        figure=fig,
        width_ratios=[1.4, 1.4] + [1] * ncols,
        wspace=0.35,
        hspace=0.35
    )

    # ── Panel A: Cell type UMAP 
    ax_main = fig.add_subplot(gs[:, :2])
    sc.pl.umap(
        adata,
        color=cell_type_column,
        ax=ax_main,
        palette=umap_palette,
        legend_loc="on data",
        legend_fontsize=11,
        legend_fontoutline=3,
        frameon=False,
        title="",
        show=False
    )
    ax_main.set_title("A", loc="left", fontweight="bold", fontsize=14)

    # ── Panel B: Gene feature plots 
    feature_axes = []
    for i, gene in enumerate(final_genes):
        row = i // ncols
        col = i % ncols
        
        ax = fig.add_subplot(gs[row, col + 2])
        feature_axes.append(ax)
        
        sc.pl.umap(
            adata,
            color=gene,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap_name,
            frameon=False,
            title=gene,
            show=False
        )
        ax.title.set_fontsize(12)

    # Add "B" label on the first feature plot
    if feature_axes:  # only if there are genes
        feature_axes[0].set_title("B", loc="left", fontweight="bold", fontsize=14)

    plt.tight_layout()
    
    return fig



def plot_gene_read_distribution(adata):
    """
    Plots the distribution of total reads per gene in an AnnData object.

    Parameters:
    - adata: AnnData object containing gene expression data.

    Outputs:
    - Prints basic statistics (total, min, max, median, mean reads per gene).
    - Displays a histogram of the distribution of total reads per gene.
    """
    # Calculate total reads per gene
    gene_total_reads = np.array(adata.X.sum(axis=0)).flatten()

    # Print basic statistics
    print(f"Total genes: {len(gene_total_reads)}")
    print(f"Minimum reads per gene: {gene_total_reads.min()}")
    print(f"Maximum reads per gene: {gene_total_reads.max()}")
    print(f"Median reads per gene: {np.median(gene_total_reads)}")
    print(f"Mean reads per gene: {gene_total_reads.mean()}")

    # Plot the distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(gene_total_reads, bins=50, kde=True, log_scale=True)
    plt.xlabel("Total Reads per Gene")
    plt.ylabel("Frequency")
    plt.title("Distribution of Total Reads per Gene")
    plt.show()

def plot_doublet_umaps(ann_obj):

    ann_obj.obs["predicted_doublet"] = ann_obj.obs["predicted_doublet"].astype("category")
    sc.pp.normalize_total(ann_obj) # Norm to median total count
    sc.pp.log1p(ann_obj)
    sc.pp.highly_variable_genes(ann_obj, n_top_genes=2000, flavor="seurat_v3")
    sc.tl.pca(ann_obj, svd_solver='arpack')
    sc.pp.neighbors(ann_obj)
    sc.tl.leiden(ann_obj)
    sc.tl.umap(ann_obj)
    sc.pl.umap(ann_obj, color=['leiden'])

    return sc.pl.umap(ann_obj, color = ["leiden", "predicted_doublet", "doublet_score"], wspace = 0.1)
        
def create_umap_visualisations(adata, resolutions, leiden_prefix="leiden", clustering_algorithm="Leiden"):
    """
    Creates a multi-panel UMAP visualization for multiple Leiden clustering resolutions.
    
    Parameters:
    - adata: AnnData object with UMAP coordinates and Leiden clusters computed.
    - resolutions: List of resolution values used for Leiden clustering.
    - leiden_prefix: Prefix used for Leiden cluster keys in adata.obs.
    - clustering_algorithm: The clustering algorithm used (default is "Leiden").

    Returns:
    - A matplotlib figure with multiple panels, one for each resolution.
    """
    # Calculate number of rows and columns for the subplots grid
    num_plots = len(resolutions)
    num_cols = 3  # Adjust based on your layout preference
    num_rows = (num_plots + num_cols - 1) // num_cols  # Ceiling division
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 5, num_rows * 5))
    
    # Flatten axes array for easy indexing
    axes = axes.flatten()

    # Loop through resolutions and plot each UMAP
    for i, res in enumerate(resolutions):
        ax = axes[i]
        sc.pl.umap(adata, color=f'{leiden_prefix}_{res}', ax=ax, show=False)
        
        # Set title with clustering algorithm and resolution
        ax.set_title(f'{clustering_algorithm} - Resolution {res}', fontsize=10)
    
    # Hide unused subplots (if any)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    
    # Adjust layout to avoid overlap
    plt.tight_layout()
    
    return fig

def plot_filtered_violin(
    adata, 
    gene_sets,  # Gene sets are now a list of tuples (name, gene_set)
    groupby_base, 
    resolutions=None, 
    row_palette=None, 
    swap_axes=True, 
    clustering_algorithm="Leiden",  # Optional: Add clustering algorithm to title
    **kwargs
):
    """
    Plot stacked violin plots for multiple gene sets, ignoring missing genes and reporting missing ones,
    with support for multiple resolutions in a multi-panel format.
    
    Parameters:
        adata (AnnData): The AnnData object containing the data.
        gene_sets (list of tuples): List of tuples where each tuple is (name, gene_set).
        groupby_base (str): Base name of the groupby key (e.g., 'leiden').
        resolutions (list, optional): List of resolutions for which to plot violin plots. Default is None.
        row_palette (list or dict, optional): Color palette for rows. Default is None.
        swap_axes (bool, optional): Whether to swap axes. Default is True.
        clustering_algorithm (str, optional): Clustering algorithm used (e.g., "Leiden").
        **kwargs: Additional keyword arguments passed to `sc.pl.stacked_violin`.

    Returns:
        - A tuple (fig, axes) containing the figure and axes objects.
    """

    # Suppress warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

    # Suppress matplotlib/seaborn logs
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

    # Separate valid and missing genes for each gene set
    valid_gene_sets = []
    missing_gene_sets = []

    for name, gene_set in gene_sets:
        valid_genes = [gene for gene in gene_set if gene in adata.var_names or gene in getattr(adata.raw, 'var_names', [])]
        missing_genes = [gene for gene in gene_set if gene not in valid_genes]
        
        valid_gene_sets.append(valid_genes)
        missing_gene_sets.append(missing_genes)

        # Report missing genes for each gene set
        if missing_genes:
            print(f"Genes not found in dataset for gene set {name}: {', '.join(missing_genes)}")
    
    # Check if any valid gene set exists
    if not any(valid_gene_sets):
        print("None of the genes in the provided gene sets are present in the dataset.")
        return None

    # Default to a single resolution if none are provided
    if resolutions is None:
        resolutions = [""]
    
    # Calculate number of rows and columns for the subplots grid
    num_plots = len(gene_sets) * len(resolutions)
    num_cols = 3  # Adjust based on your layout preference
    num_rows = (num_plots + num_cols - 1) // num_cols  # Ceiling division
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 5, num_rows * 5))
    
    # Flatten axes array for easy indexing
    axes = axes.flatten()

    # Temporarily suppress the specific warning inside the function
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # Suppress all UserWarnings
        
        # Loop through resolutions and plot each gene set
        plot_idx = 0
        for res in resolutions:
            groupby = f"{groupby_base}_{res}" if res else groupby_base
            print(f"Generating plots for resolution: {res if res else 'default'}")
        
            for i, (gene_set_tuple, valid_genes) in enumerate(zip(gene_sets, valid_gene_sets)):
                name, gene_set = gene_set_tuple
                print(f"Plotting gene set {name} with {len(valid_genes)} genes out of {len(gene_set)} provided.")
                
                ax = axes[plot_idx]
                sc.pl.stacked_violin(
                    adata, 
                    valid_genes, 
                    groupby=groupby, 
                    row_palette=row_palette, 
                    swap_axes=swap_axes, 
                    ax=ax, 
                    show=False,  # Suppress inline plotting
                #    use_raw=False,  # Use raw data if available
                    **kwargs
                )
                
                # Set title with gene set name and resolution
                ax.set_title(f'{name} - {clustering_algorithm} - Resolution {res}' if res else f'{name} - {clustering_algorithm}')
                plot_idx += 1

        # Hide unused subplots (if any)
        for j in range(plot_idx, len(axes)):
            axes[j].axis('off')

        # Rotate the legend and position it above the plot
        for ax in axes[:plot_idx]:  # Only apply to axes that have a plot
            legend = ax.get_legend()
            if legend:
                legend.set_bbox_to_anchor((0.5, 1.1))  # Position legend above the plot area
                legend.set_loc("center")
                legend.set_frame_on(False)
                for label in legend.get_texts():
                    label.set_rotation(0)  # No rotation, keeping legend horizontally aligned

    # Adjust layout to avoid overlap
    plt.tight_layout()

    return fig, axes  # Return both figure and axes for further customization

def plot_umap_grid(ann_obj, obs_columns, grid_size=(2, 2), figsize=(10, 8), save_path=None):
    """
    Generate a grid of UMAP plots for specified obs columns in a Scanpy AnnData object.

    Parameters:
        ann_obj (AnnData): The Scanpy AnnData object containing UMAP coordinates and metadata.
        obs_columns (list): List of obs columns to color UMAP plots by.
        grid_size (tuple): Grid size for the plots (rows, columns). Defaults to (2, 2).
        figsize (tuple): Figure size for the entire grid. Defaults to (10, 8).
        save_path (str): Optional path to save the figure. Defaults to None.

    Returns:
        None
    """
    num_plots = len(obs_columns)
    rows, cols = grid_size

    # Ensure enough rows and columns to fit all plots
    if rows * cols < num_plots:
        raise ValueError(f"Grid size {grid_size} is too small for {num_plots} plots.")
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    axes = np.atleast_1d(axes).flatten()

    for i, column in enumerate(obs_columns):
        sc.pl.umap(ann_obj, color=column, ax=axes[i], show=False, title=column)
    
    # Remove unused subplots
    for j in range(num_plots, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    # Save the figure if a path is specified
    if save_path:
        plt.savefig(save_path)
        print(f"Figure saved to {save_path}")

    plt.show()

def plot_stacked_figure(adata, sample_column, color_column=None, barplot=False, recalculate_columns=True):
    """
    Create a stacked figure with boxplots and an optional stacked percentage barplot.
    
    Parameters:
    - adata: AnnData object.
    - sample_column: Column in `adata.obs` to use as sample IDs for x-axis.
    - color_column: Column in `adata.obs` to determine cell types or classes for the barplot.
    - barplot: Whether to include the stacked barplot (default: False).
    - recalculate_columns: Whether to check and potentially recalculate `n_counts` and `n_genes` columns in `adata.obs` (default: True).
    """
    # Ensure sample_column is a string
    adata.obs[sample_column] = adata.obs[sample_column].astype(str)

    def detect_data_type(matrix):
        """
        Detect the type of data in adata.X.
        Returns: 'raw_counts', 'log_normalized', or 'unknown'
        """
        # Extract a small subset of the matrix
        subset = matrix[:5, :5].toarray() if issparse(matrix) else matrix[:5, :5]
        
        # Check if all values are non-negative integers (raw counts)
        if np.all(subset >= 0) and np.all(np.isclose(subset, subset.astype(int))):
            return 'raw_counts'
        
        # Check if values appear log-normalized (e.g., <= 1, often log-scaled)
        if np.all(subset <= 10):  # Arbitrary threshold for log-normalized data
            return 'log_normalized'
        
        # If it doesn't match the above patterns, return 'unknown'
        return 'unknown'

    if recalculate_columns:
        # Detect data type
        data_type = detect_data_type(adata.X)
        if data_type == 'raw_counts':
            print("Recalculating 'n_counts' and 'n_genes'...")
            recalculated_n_counts = adata.X.sum(axis=1).A1 if issparse(adata.X) else adata.X.sum(axis=1)
            recalculated_n_genes = np.array((adata.X > 0).sum(axis=1)).flatten() if issparse(adata.X) else (adata.X > 0).sum(axis=1)

            # Compare with existing columns
            if ('n_counts' in adata.obs and np.allclose(recalculated_n_counts, adata.obs['n_counts'])) and \
               ('n_genes' in adata.obs and np.all(recalculated_n_genes == adata.obs['n_genes'])):
                print("No changes detected in 'n_counts' and 'n_genes'. Columns remain unchanged.")
            else:
                print("Differences detected. Updating 'n_counts' and 'n_genes' columns in the adata object...")
                adata.obs['n_counts'] = recalculated_n_counts
                adata.obs['n_genes'] = recalculated_n_genes
        else:
            print("Skipping recalculation of 'n_counts' and 'n_genes' due to non-raw data in `adata.X`.")

    # Get sample order and group data
    sample_order = adata.obs[sample_column].unique()
    umi_counts = adata.obs.groupby(sample_column)["n_counts"].apply(list)
    gene_counts = adata.obs.groupby(sample_column)["n_genes"].apply(list)

    # Extract Scanpy's color palette for the color_column (if applicable)
    cluster_colors = None
    if color_column:
        cluster_colors = adata.uns.get(f'{color_column}_colors', None)

    # Calculate proportions for barplot
    proportions = None
    if color_column and barplot:
        proportions = adata.obs.groupby([sample_column, color_column]).size().unstack(fill_value=0)
        proportions = proportions.div(proportions.sum(axis=1), axis=0)  # Normalize to percentages

    # Initialize figure
    fig, axs = plt.subplots(3 if barplot else 2, 1, figsize=(15, 10), sharex=True, gridspec_kw={'height_ratios': [1, 1, 2] if barplot else [1, 1]})

    # Boxplot for UMIs per cell
    umi_data = pd.DataFrame({'Sample': np.repeat(sample_order, [len(umi_counts[s]) for s in sample_order]),
                             'UMIs': np.concatenate([umi_counts[s] for s in sample_order])})
    sns.boxplot(data=umi_data, x="Sample", y="UMIs", ax=axs[0], color="white", fliersize=1, linewidth=1)
    axs[0].set_ylabel("Number of UMIs per cell")
    axs[0].set_xlabel("")
    axs[0].tick_params(axis='x')
    axs[0].set_xticklabels(umi_data['Sample'].unique(), rotation=90, fontsize=7)

    # Add median annotation for UMIs
    umi_median = np.median(np.concatenate([umi_counts[s] for s in sample_order]))
    axs[0].text(0.02, 0.95, f'Median: {int(umi_median):,} UMIs', transform=axs[0].transAxes, 
                ha="left", va="top", fontsize=10, bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    # Boxplot for genes per cell
    gene_data = pd.DataFrame({'Sample': np.repeat(sample_order, [len(gene_counts[s]) for s in sample_order]),
                              'Genes': np.concatenate([gene_counts[s] for s in sample_order])})
    sns.boxplot(data=gene_data, x="Sample", y="Genes", ax=axs[1], color="white", fliersize=1, linewidth=1)
    axs[1].set_ylabel("Number of genes per cell")
    axs[1].set_xlabel("")
    axs[1].tick_params(axis='x')
    axs[1].set_xticklabels(gene_data['Sample'].unique(), rotation=90, fontsize=7)

    # Add median annotation for genes
    gene_median = np.median(np.concatenate([gene_counts[s] for s in sample_order]))
    axs[1].text(0.02, 0.95, f'Median: {int(gene_median):,} genes', transform=axs[1].transAxes, 
                ha="left", va="top", fontsize=10, bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    # Stacked barplot for cell class proportions
    if barplot and color_column and proportions is not None:
        proportions.plot(kind="bar", stacked=True, ax=axs[2], width=1, legend=False, color=cluster_colors)
        axs[2].set_ylabel("Cell class proportions")
        axs[2].set_xlabel("Sample")
        axs[2].tick_params(axis='x')
        axs[2].set_xticklabels(proportions.index, rotation=90, fontsize=7)
        if cluster_colors:
            axs[2].legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=color_column)

    plt.tight_layout()
    plt.show()

# Obsolete: Need to add excel generation to above function
def plot_and_save_cluster_percentages(adata, output_dir, clustering_param='leiden'):
    """
    Generate a stacked bar plot showing the percentage of cells in each cluster per sample,
    and save the pivot table as an Excel file.

    Parameters:
        adata: AnnData
            The annotated data matrix.
        output_dir: str
            Directory where the Excel file will be saved.
        clustering_param: str
            The clustering parameter to base the clusters on (e.g., 'leiden', 'leiden_0.5').

    Returns:
        fig: matplotlib.figure.Figure
            The generated figure.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract sample and cluster information
    df = adata.obs[['sample', clustering_param]]

    # Count the number of cells per sample and cluster
    counts = df.groupby(['sample', clustering_param]).size().reset_index(name='cell_count')

    # Pivot the table to get clusters as columns
    pivot_table = counts.pivot(index='sample', columns=clustering_param, values='cell_count').fillna(0)

    # Calculate the total number of cells per sample
    pivot_table['total_cells'] = pivot_table.sum(axis=1)

    # Convert counts to percentages
    percentage_table = pivot_table.div(pivot_table['total_cells'], axis=0).drop(columns=['total_cells']) * 100

    # Define the cluster labels and colors
    clusters = percentage_table.columns
    cluster_colors = adata.uns[f'{clustering_param}_colors']  # Use Scanpy's cluster colors

    # Plot the stacked bar plot
    fig, ax = plt.subplots(figsize=(30, 6))
    percentage_table.plot(kind='bar', stacked=True, ax=ax, color=cluster_colors[:len(clusters)], width=0.8)

    # Customize the plot
    ax.set_title('Percentage of Cells from Each Cluster per Sample', fontsize=16)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('Percentage of Cells', fontsize=14)
    ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticklabels(percentage_table.index, rotation=45, ha='right')
    plt.tight_layout()

    # Save the pivot table as an Excel file
    excel_file_path = os.path.join(output_dir, "cell_cnts_per_smpl_per_clust.xlsx")
    pivot_table.to_excel(excel_file_path)

    # Return the plot and the pivot table
    return fig


def adjust_colorbar(cb, xlabel):
    """Safely adjust colorbar placement and labels"""
    if cb:
        cb.ax.xaxis.set_ticks_position("top")  # Move ticks to top
        cb.ax.xaxis.set_label_position("top")  # Move label to top
        cb.ax.set_xlabel(xlabel)  # Set label
        cb.ax.tick_params(axis="x", rotation=0)  # Ensure horizontal labels

def plot_rank_genes_groups(adata, n_genes=5, key="t-test_ov", groupby="leiden_harmony_0.1"):
    fig, axes = plt.subplots(1, 3, figsize=(15, 10))  # 3-column layout

    # Dotplot
    sc.pl.rank_genes_groups_dotplot(
        adata, n_genes=n_genes, key=key, groupby=groupby, ax=axes[0], show=False, swap_axes=True
    )
    if axes[0].collections:
        cb = axes[0].collections[-1].colorbar
        adjust_colorbar(cb, "Mean expression\nin group")

    # Stacked Violin
    sc.pl.rank_genes_groups_stacked_violin(
        adata, n_genes=n_genes, key=key, groupby=groupby, ax=axes[1], show=False, swap_axes=True
    )
    if axes[1].collections:
        cb = axes[1].collections[-1].colorbar
        adjust_colorbar(cb, "Median expression\nin group")

    # Matrixplot
    sc.pl.rank_genes_groups_matrixplot(
        adata, n_genes=n_genes, key=key, groupby=groupby, ax=axes[2], show=False, swap_axes=True
    )
    if axes[2].collections:
        cb = axes[2].collections[-1].colorbar
        adjust_colorbar(cb, "Mean expression\nin group")

    plt.tight_layout()
    plt.show()


def classify_sample_sex(adata, output_csv="sex_assignments.csv"):
    """
    Classify samples by sex using expression of XIST, RPS4Y1, and DDX3Y from adata.raw.
    Generates a visualization with three subplots: scatter plot, boxplot, and PCA.
    
    Parameters:
        adata (AnnData): Annotated data matrix with raw counts in `adata.raw`.
        output_csv (str): Path to save the resulting CSV file.
    
    Returns:
        tuple: (pd.DataFrame, matplotlib.figure.Figure) containing the sex classification DataFrame
               and the figure with three subplots.
    """
    
    # Step 1: Access raw counts
    raw_adata = adata.raw.to_adata()
    
    # Step 2: Identify available sex genes
    sex_genes = ['XIST', 'RPS4Y1', 'DDX3Y']
    available_genes = [gene for gene in sex_genes if gene in raw_adata.var_names]
    missing_genes = [gene for gene in sex_genes if gene not in raw_adata.var_names]
    print("Available genes in raw data:", available_genes)
    print("Missing genes in raw data:", missing_genes)
    
    if not available_genes:
        raise ValueError("No sex genes found in adata.raw!")
    
    # Step 3: Normalize and log-transform
    sex_raw = raw_adata[:, available_genes]
    sc.pp.normalize_total(sex_raw, target_sum=1e4)
    sc.pp.log1p(sex_raw)
    
    # Step 4: Aggregate by sample
    sex_adata = sex_raw.to_df()
    sample_sex_exp = sex_adata.groupby(adata.obs['sample']).mean()
    
    # Step 5: Thresholds and classification
    xist_thresh = sample_sex_exp['XIST'].median() + sample_sex_exp['XIST'].std()
    y_thresh = 0.5  # adjustable if needed

    def classify_sex(row):
        if row['XIST'] > xist_thresh and row.get('RPS4Y1', 0) < y_thresh and row.get('DDX3Y', 0) < y_thresh:
            return 'Female'
        elif (row.get('RPS4Y1', 0) > y_thresh or row.get('DDX3Y', 0) > y_thresh) and row['XIST'] < xist_thresh:
            return 'Male'
        else:
            return 'Ambiguous'

    sample_sex_exp['Predicted_Sex'] = sample_sex_exp.apply(classify_sex, axis=1)
    sample_sex_exp['Sex_Numeric'] = sample_sex_exp['Predicted_Sex'].map({'Male': 1, 'Female': 2, 'Ambiguous': 0})

    # Step 6: Output DataFrame
    output_df = sample_sex_exp[['Predicted_Sex', 'Sex_Numeric']].reset_index()
    output_df.columns = ['Sample_ID', 'Predicted_Sex', 'Sex_Numeric']
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    output_df.to_csv(output_csv, index=False)
    print(f"Saved to {output_csv}")

    # Step 7: Visualization
    # Choose a Y gene to visualize (prefer RPS4Y1, fallback to DDX3Y)
    y_gene = 'RPS4Y1' if 'RPS4Y1' in available_genes else 'DDX3Y' if 'DDX3Y' in available_genes else None
    if y_gene is None:
        raise ValueError("No Y-linked gene found among available genes.")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Scatter plot: XIST vs Y gene
    sns.scatterplot(
        data=sample_sex_exp,
        x='XIST',
        y=y_gene,
        hue='Predicted_Sex',
        palette={'Male': 'blue', 'Female': 'red', 'Ambiguous': 'gray'},
        ax=axes[0]
    )
    axes[0].axvline(xist_thresh, color='red', linestyle='--', label='XIST threshold')
    axes[0].axhline(y_thresh, color='blue', linestyle='--', label=f'{y_gene} threshold')
    axes[0].set_title(f'XIST vs {y_gene}')
    axes[0].legend()

    # Boxplot: all genes by Predicted_Sex
    melted = sample_sex_exp.reset_index()[['Predicted_Sex'] + available_genes].melt(id_vars='Predicted_Sex')
    sns.boxplot(data=melted, x='variable', y='value', hue='Predicted_Sex', ax=axes[1])
    axes[1].set_title('Gene Expression by Predicted Sex')
    axes[1].set_xlabel('Gene')
    axes[1].set_ylabel('Log-normalized expression')
    axes[1].legend(loc='upper right')

    # PCA on available sex genes
    df_pca_input = sample_sex_exp[available_genes].dropna()
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(df_pca_input)
    df_pca = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=df_pca_input.index)
    df_pca['Predicted_Sex'] = sample_sex_exp.loc[df_pca.index, 'Predicted_Sex']

    sns.scatterplot(data=df_pca, x='PC1', y='PC2', hue='Predicted_Sex', palette='Set1', ax=axes[2])
    axes[2].set_title('PCA of Sex Gene Expression')

    plt.tight_layout()

    return output_df, fig