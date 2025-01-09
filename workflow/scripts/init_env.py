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

# Configure logging - dynamically change log filename based on 'plate' value
def setup_logging(plate):
    log_filename = f"../../results/00LOG/03SCANPY/scanpy_qc_{plate}.log" # Need to be rel to this script
    logging.basicConfig(
        level=logging.INFO,  # Set the log level (INFO to capture informational messages)
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),  # Print to console
            logging.FileHandler(log_filename)  # Log to a file
        ]
    )
    logger = logging.getLogger()
    return logger

# Check if plate is defined, if not set default plate and initialize logging
if 'plate' not in locals() and 'plate' not in globals():
    plate = 'plate1'  # Set a default plate in case it's not passed in
    logger = setup_logging(plate)
    logger.info("Using default plate: %s", plate)
else:
    # Initialize logging after plate is assigned
    logger = setup_logging(plate)
    logger.info("Plate argument passed: %s", plate)
    logger.info("Processing plate: %s", plate)

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

init_notebook_mode(all_interactive=True)

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor='white', figsize=(12,6), format='png')

# Set root dir
sc.settings.verbosity = 4
if os.path.exists('/scratch/c.c1477909/eQTL_study_2025/'):
    root_dir = '/scratch/c.c1477909/eQTL_study_2025/'
    import argparse

    def parse_arguments():
        parser = argparse.ArgumentParser(description="Setup environment with plate argument.")
        parser.add_argument('--plate', type=str, required=True, help='Plate identifier (e.g., plate1, plate2)')
        return parser.parse_args()

    # Check for arguments if called as part of the script
    if __name__ == '__main__':
        args = parse_arguments()
        plate = args.plate
        logger.info("Plate argument passed: %s", plate)
    
    logger.info("Plate argument passed: %s", plate)
    logger.info("Processing plate: %s", plate)
    os.environ['OMP_NUM_THREADS'] = '16'
else:
    root_dir = '/Users/darren/Desktop/eQTL_study_2025/'
    plate = 'plate1'
    logger.info("Using default plate: %s", plate)

script_dir = root_dir + 'scripts/'
results_dir = root_dir + 'results/'
data_dir = results_dir + '02PARSE/'
plate_path = data_dir + f'combine_{plate}/all-sample/DGE_filtered/anndata.h5ad'
scanpy_dir = results_dir + '03SCANPY/'
sc.settings.figdir = results_dir + '/figs/'
sys.path.append(script_dir)  # Add the directory containing custom scripts to sys.path

logger.info("Script initialized. Root directory: %s", root_dir)
logger.info("Data directory: %s", data_dir)
logger.info(f"Saving to directory: {scanpy_dir}")
logger.info(f"Directory exists: {os.path.exists(scanpy_dir)}")
