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

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

init_notebook_mode(all_interactive=True)

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(12,6), format='png')

# Set root dir

sc.settings.verbosity = 4
if 'snakemake' in locals() or 'snakemake' in globals():
    root_dir = '/scratch/c.c1477909/eQTL_study_2025/'
    os.environ['OMP_NUM_THREADS'] = '16'
else:
    root_dir = '/Users/darren/Desktop/eQTL_study_2025/'

script_dir = root_dir + 'scripts/'
results_dir = root_dir + 'results/'
data_dir = results_dir + '02PARSE/'
plate1_dir = data_dir + 'combine_plate1/all-sample/DGE_filtered/'
plate2_dir = data_dir + 'combine_plate2/all-sample/DGE_filtered/'
scanpy_dir = results_dir + '03SCANPY/'
sc.settings.figdir = results_dir + '/figs/'
sys.path.append(script_dir)  # Add the directory containing custom scripts to sys.path
