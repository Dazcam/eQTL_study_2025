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
     os.environ['OMP_NUM_THREADS'] = '16'

script_dir = workflow_dir + 'workflow/scripts/'
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
logger.info(f"plate_path set to: {plate_path}")
