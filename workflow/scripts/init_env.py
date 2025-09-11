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
import matplotlib
import matplotlib.pyplot as plt
import gc
import warnings
from itables import init_notebook_mode
import numba
from memory_profiler import profile
import ipywidgets as widgets
import logging  # Import the logging module
from pathlib import Path

def initialize_env(plate):
    """
    Initialize the environment for the script, including logging, settings, and paths.
    Accepts `plate` as an argument.
    """
    # Function to configure logging dynamically

    
    def setup_logging(plate):
        if os.path.exists('/scratch/'):
            root_dir = '/scratch/c.c1477909/eQTL_study_2025/'
        else:
            root_dir = '/Users/darren/Desktop/eQTL_study_2025/'

        log_dir = os.path.join(root_dir, "results/00LOG/02SCANPY")
        os.makedirs(log_dir, exist_ok=True)

        if 'plate' in plate: 
            log_filename = os.path.join(log_dir, f"scanpy_qc_{plate}.log")
        else:
            log_filename = os.path.join(log_dir, f"scanpy_{plate}.log")
            
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

    # Set thread environment variable if on a specific system
    if os.path.exists('/scratch/'):
        os.environ['OMP_NUM_THREADS'] = '16'
        root_dir = '/scratch/c.c1477909/eQTL_study_2025/'
    else:
        root_dir = '/Users/darren/Desktop/eQTL_study_2025/'

    # Define and set paths
    resources_dir = root_dir + 'resources/'
    sheets_dir = resources_dir + 'sheets/'
    cameron_dir = resources_dir + 'public_datasets/cameron_2023/'
    script_dir = root_dir + 'workflow/scripts/'
    results_dir = root_dir + 'results/'
    parse_dir = results_dir + '01PARSE/'
    plate_dir = parse_dir + f'combine_{plate}/all-sample/DGE_filtered/'
    scanpy_dir = results_dir + '02SCANPY/'
    report_dir = root_dir + 'workflow/reports/02SCANPY/'
    sc.settings.figdir = root_dir + 'reports/figs/'
    sys.path.append(script_dir)

    # Log important directories
    logger.info("Script initialized. Root directory: %s", root_dir)
    logger.info("Parse directory: %s", parse_dir)
    logger.info("Scanpy directory: %s", scanpy_dir)
    logger.info("Scanpy directory exists: %s", os.path.exists(scanpy_dir))
    logger.info("Sheets directory: %s", sheets_dir)
    logger.info(f"plate_dir set to: {plate_dir}")

    # Return logger and paths as needed
    return logger, root_dir, sheets_dir, plate_dir, scanpy_dir, report_dir

