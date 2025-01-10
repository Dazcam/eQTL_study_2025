import os
import sys
import re
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.io as sio
import scanpy.external as sce
import matplotlib.pyplot as plt
import gc
import warnings
from itables import init_notebook_mode
import numba
from memory_profiler import profile
import ipywidgets as widgets
import logging  # Import the logging module
from pathlib import Path

# Function to configure logging dynamically
def setup_logging(plate):
    log_filename = f"../../results/00LOG/03SCANPY/scanpy_qc_{plate}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_filename)
        ]
    )
    logger = logging.getLogger()
    return logger
from pathlib import Path

# Main function to encapsulate script behavior
def main():
    # Detect the environment and set plate
    if os.path.exists('/scratch/'):
        # Extract plate identifier from snakemake.output['nb']
        try:
            plate=param1
            print(f"Running on Hawk. Plate set to: {plate}")       
        except Exception as e:
            raise ValueError(f"Failed to extract plate from snakemake.output['nb']: {e}")
    else:
        plate = "plate1"  # Placeholder for local debugging or development
        print(f"Running locally. Plate set to: {plate}")

    # Setup logging
    logger = setup_logging(plate)
    logger.info("Plate variable detected: %s", plate)
    logger.info("Processing plate: %s", plate)

    # Suppress warnings
    warnings.simplefilter("ignore", FutureWarning)
    warnings.simplefilter("ignore", UserWarning)
    warnings.simplefilter("ignore", RuntimeWarning)

    # Initialize interactive table mode
    init_notebook_mode(all_interactive=True)

    # Adjust Scanpy figure defaults
    sc.settings.set_figure_params(
        dpi=100, fontsize=10, dpi_save=400,
        facecolor='white', figsize=(12, 6), format='png'
    )

    # Set root directory
    sc.settings.verbosity = 4
    if os.path.exists('/scratch/'):
        root_dir = '/scratch/c.c1477909/eQTL_study_2025/workflow/'
        os.environ['OMP_NUM_THREADS'] = '16'
    else:
        root_dir = '/Users/darren/Desktop/eQTL_study_2025/workflow/'

    script_dir = root_dir + 'scripts/'
    results_dir = root_dir + 'results/'
    data_dir = results_dir + '02PARSE/'
    plate_path = data_dir + f'combine_{plate}/all-sample/DGE_filtered/anndata.h5ad'
    scanpy_dir = results_dir + '03SCANPY/'
    sc.settings.figdir = results_dir + '/figs/'
    sys.path.append(script_dir)

    # Log important directories
    logger.info("Script initialized. Root directory: %s", root_dir)
    logger.info("Data directory: %s", data_dir)
    logger.info("Saving to directory: %s", scanpy_dir)
    logger.info("Directory exists: %s", os.path.exists(scanpy_dir))


# Ensure script runs only when executed directly
if __name__ == "__main__":
    main()

