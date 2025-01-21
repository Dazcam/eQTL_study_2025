import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import gc
import re
import ipywidgets as widgets
from IPython import get_ipython
from ipywidgets import Tab, Output
from memory_profiler import profile
import warnings

##########  LOAD FUNCTIONS  ############  
#@profile
def load_and_dwnsmpl_data(downsample_cells=None, *adata_dirs):
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
        adata = sc.read(adata_dir)
        
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

def filter_genes_by_read_count(adata, min_reads=10, min_samples=100, inplace=True):
    """
    Identifies genes with at least `min_reads` in at least `min_samples` in the dataset.
    
    Parameters:
    - adata: AnnData object containing gene expression data.
    - min_reads: Minimum number of reads required per cell/sample.
    - min_samples: Minimum number of samples/cells that must meet the `min_reads` threshold.
    - inplace: If True, filters genes directly in `adata`. If False, returns a mask.

    Returns:
    - If `inplace=True`: Modifies `adata` by retaining only the filtered genes.
    - If `inplace=False`: Returns a boolean mask of the genes that meet the criteria.
    """
    # Handle both sparse and dense matrices
    if isinstance(adata.X, np.ndarray):  # Dense array
        gene_mask = (adata.X >= min_reads).sum(axis=0) >= min_samples
    else:  # Sparse matrix
        gene_mask = (adata.X >= min_reads).sum(axis=0).A1 >= min_samples

    if inplace:
        # Subset the AnnData object to retain only the filtered genes
        adata._inplace_subset_var(gene_mask)
    else:
        # Return the mask for external usage
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
def plot_doublet_umaps(ann_obj):

    ann_obj.obs["predicted_doublet"] = ann_obj.obs["predicted_doublet"].astype("category")
    sc.pp.normalize_total(ann_obj) # Norm to median total count
    sc.pp.log1p(ann_obj)
    sc.pp.highly_variable_genes(ann_obj, n_top_genes=2000, flavor="seurat_v3")
    sc.tl.pca(ann_obj, svd_solver='arpack')
    sc.pp.neighbors(ann_obj)
    sc.tl.leiden(ann_obj)
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
    axes = axes.flatten()  # Flatten axes for easy indexing

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