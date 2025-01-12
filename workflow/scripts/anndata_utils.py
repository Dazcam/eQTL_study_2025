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


##########  LOAD FUNCTIONS  ############  
@profile
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

    # Step 1: Filter cells based on minimum number of genes
    num_cells_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    num_cells_after = adata.shape[0]

    # Step 2: Filter cells based on minimum total counts (reads)
    num_cells_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    num_cells_after = adata.shape[0]

    # Step 3: Filter genes based on minimum number of cells expressing them
    num_genes_before = adata.shape[1]
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    num_genes_after = adata.shape[1]

    # Step 4: Filter samples based on minimum number of cells
    valid_samples = adata.obs['sample'].value_counts() >= min_cells_per_sample
    valid_sample_ids = valid_samples[valid_samples].index
    num_samples_before = adata.obs['sample'].nunique()
    valid_cells = adata.obs['sample'].isin(valid_sample_ids)  # Boolean mask for valid cells
    adata._inplace_subset_obs(valid_cells)  # In-place filtering of cells
    num_samples_after = adata.obs['sample'].nunique()
    print(f"filtered out {num_samples_before - num_samples_after} samples with < {min_cells_per_sample} cells.")

    print(f"Dimensions after applying filter: {adata.shape}")



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

#@profile
def load_and_process_data(downsample_cells=None, *adata_dirs):
    """
    Load, process, and merge AnnData objects with optional downsampling.
    
    Parameters:
        downsample_cells (int, optional): Number of cells to randomly downsample to. If None, no downsampling is performed.
        *adata_dirs (str): Paths to AnnData directories. Each directory should contain an 'anndata.h5ad' file.
        
    Returns:
        AnnData: A merged and processed AnnData object.
    """
    adata_list = []   
    
    # Extract plate numbers from the directory paths
    plate_numbers = [re.search(r'plate(\d+)', adata_dir).group(1) for adata_dir in adata_dirs]

    # Report how many plates will be processed and downsampling details
    if downsample_cells != None:
        print(f"Processing {len(adata_dirs)} plates with downsampling to {downsample_cells} cells per plate.")
    else:
        print(f"Processing {len(adata_dirs)} plates, no downsampling applied.")

    for i, adata_dir in enumerate(adata_dirs):
        # Load the AnnData object
        print(f"Loading plate {plate_numbers[i]} ...")
        adata = sc.read(adata_dir + 'anndata.h5ad')
        
        # Assign unique plate information to cell ids and metadata
        adata.obs['plate'] = 'plate' + plate_numbers[i]
        adata.obs_names = f"plate{plate_numbers[i]}_" + adata.obs_names

        # Log dimensions
        print(f"Plate {plate_numbers[i]} dimensions: {adata.shape}")
        print(f"Plate {plate_numbers[i]} matrix Dimensions: {adata.X.shape}")
 
        # Optional downsampling
        if downsample_cells is not None and downsample_cells < adata.n_obs:
            random_indices = np.random.choice(adata.n_obs, downsample_cells, replace=False)
            adata = adata[random_indices, :].copy()
        
        adata_list.append(adata)

    # Find common genes across all datasets
    # Optionally handle common genes only if more than one plate
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
        
    # Concatenate the AnnData objects
    print(f"Merging plates ...")
    adata_mrg = sc.concat(adata_list, join='inner')

    # Remove 'hg38' from gene names
    adata_mrg.var.index = adata_mrg.var.index.str.replace('_hg38', '')

    # Cleanup to release memory
    del adata_list
    gc.collect()

    print(adata_mrg)

    return adata_mrg


##########  VISUALISATION FUNCTIONS  ############  

def create_umap_visualisations(adata, resolutions, leiden_prefix="leiden"):
    """
    Creates UMAP visualizations for multiple Leiden clustering resolutions.
    
    Parameters:
    - adata: AnnData object with UMAP coordinates and Leiden clusters computed.
    - resolutions: List of resolution values used for Leiden clustering.
    - leiden_prefix: Prefix used for Leiden cluster keys in adata.obs.

    Returns:
    - If in a Jupyter environment, returns a tab widget with UMAP plots for each resolution.
    - Otherwise, returns a list of matplotlib figures for each resolution.
    """
    # Check if in a Jupyter environment
    in_jupyter = is_running_in_jupyter()
    
    plots = []
    tabs = []
    
    for res in resolutions:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color=f'{leiden_prefix}_{res}', ax=ax, show=False)
        
        if in_jupyter:
            # Save the figure as a widget for the tab
            output = widgets.Output()
            with output:
                plt.show()
            tabs.append((f"Resolution {res}", output))
        else:
            # Append figure for non-Jupyter environments
            plots.append(fig)
    
    if in_jupyter:
        # Create and return a tab widget
        tab_widget = widgets.Tab()
        tab_widget.children = [content for _, content in tabs]
        for i, (title, _) in enumerate(tabs):
            tab_widget.set_title(i, title)
        return tab_widget
    else:
        # Return a list of figures
        return plots


def plot_filtered_violin(
    adata, 
    gene_sets,  # Gene sets are now a list of tuples (name, gene_set)
    groupby_base, 
    resolutions=None, 
    row_palette=None, 
    swap_axes=True, 
    **kwargs
):
    """
    Plot stacked violin plots for multiple gene sets, ignoring missing genes and reporting missing ones,
    with support for multiple resolutions.

    Parameters:
        adata (AnnData): The AnnData object containing the data.
        gene_sets (list of tuples): List of tuples where each tuple is (name, gene_set).
        groupby_base (str): Base name of the groupby key (e.g., 'leiden').
        resolutions (list, optional): List of resolutions for which to plot violin plots. Default is None.
        row_palette (list or dict, optional): Color palette for rows. Default is None.
        swap_axes (bool, optional): Whether to swap axes. Default is True.
        **kwargs: Additional keyword arguments passed to `sc.pl.stacked_violin`.

    Returns:
        - A tab widget if running in a Jupyter notebook.
        - A list of matplotlib figures if not using a Jupyter notebook.
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
    
    # Determine execution environment
    in_jupyter = is_running_in_jupyter()

    # Prepare plots
    plots = []
    tabs = []

    for res in resolutions:
        groupby = f"{groupby_base}_{res}" if res else groupby_base
        print(f"Generating plots for resolution: {res if res else 'default'}")
    
        for i, (gene_set_tuple, valid_genes) in enumerate(zip(gene_sets, valid_gene_sets)):
            name, gene_set = gene_set_tuple
            print(f"Plotting gene set {name} with {len(valid_genes)} genes out of {len(gene_set)} provided.")
            
            if in_jupyter:
                # Generate plot and store as a widget for tab
                output = Output()
                with output:
                    sc.pl.stacked_violin(
                        adata, 
                        valid_genes, 
                        groupby=groupby, 
                        row_palette=row_palette, 
                        swap_axes=swap_axes, 
                        **kwargs
                    )
                tabs.append((f"{name}: {res}", output))  # Include gene set name in the tab title
            else:
                # Generate plot and store figure
                fig = sc.pl.stacked_violin(
                    adata, 
                    valid_genes, 
                    groupby=groupby, 
                    row_palette=row_palette, 
                    swap_axes=swap_axes, 
                    show=False,  # Suppress inline plotting
                    **kwargs
                )
                plots.append(fig)


    # Return appropriate output
    if in_jupyter:
        # Create a tab widget
        tab_widget = Tab()
        tab_widget.children = [content for _, content in tabs]

        # Set tab titles
        for i, (title, _) in enumerate(tabs):
            tab_widget.set_title(i, title)

        return tab_widget
    else:
        return plots
