# This file imports all necessary packages and libraries used in the project.
# It includes mathematical, file handling, parallel processing, and plotting libraries.

# ------------------------------------------------------
# Math Packages
# ------------------------------------------------------
import numpy as np
from math import *  # careful: imports everything into global namespace
import dedalus.public as d3
from scipy import special
from scipy.interpolate import CubicSpline
from decimal import *

# Set decimal precision
getcontext().prec = 50

# ------------------------------------------------------
# File Handling and System Utilities
# ------------------------------------------------------
import sys
import os
import argparse
import importlib
import h5py
import csv
import gc
import re
import subprocess

# ------------------------------------------------------
# Parallel Processing
# ------------------------------------------------------
import concurrent.futures

# ------------------------------------------------------
# Timing, Warnings, Logging
# ------------------------------------------------------
import timeit
import warnings
warnings.filterwarnings("ignore")

import logging
logging.getLogger().setLevel(logging.WARNING)

# ------------------------------------------------------
# Plotting
# ------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker

# Global matplotlib settings
plt.rcParams['figure.constrained_layout.use'] = True
mpl.rcParams['axes.linewidth'] = 2  # Set axes line width globally
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Verdana', 'Liberation Sans']

# Default font for plots
font = {
    'color': 'black',
    'weight': 'normal',
    'size': 20,
}
