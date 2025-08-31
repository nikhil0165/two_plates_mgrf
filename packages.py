## -------------------------------
## Math packages
## -------------------------------

import numpy as np               # Standard numerical Python library
from math import *               # Standard Python math functions (e.g., sin, cos, sqrt, log)
import dedalus.public as d3      # Dedalus spectral PDE solver library
from scipy import special        # Special functions (e.g., Bessel, gamma functions)
from scipy.interpolate import CubicSpline  # Cubic spline interpolation

from decimal import *            # High-precision arithmetic
getcontext().prec=50             # Set decimal precision to 50 digits

## -------------------------------
## File reading and writing
## -------------------------------

import sys                        # System-specific parameters and functions
import os                         # Operating system interface (e.g., file paths)
import argparse                   # Command-line argument parsing
import importlib                  # Dynamic module import
import h5py                       # Reading/writing HDF5 files
import csv                        # Reading/writing CSV files
import gc                         # Garbage collection (manual memory management)

## -------------------------------
## For parallelising self-energy calculations
## -------------------------------

import concurrent.futures         # Multi-core parallelization using ProcessPoolExecutor/ThreadPoolExecutor
import re                         # Regular expressions
import subprocess                 # Running external shell commands

## -------------------------------
## Timing, warnings, and logging
## -------------------------------

import timeit                     # Measure execution time of code snippets
import warnings
warnings.filterwarnings("ignore") # Suppress warnings
import logging
logging.getLogger().setLevel(logging.WARNING)  # Set logging level to WARNING

## -------------------------------
## For plotting
## -------------------------------

import matplotlib.pyplot as plt   # Standard plotting library
import matplotlib as mpl
plt.rcParams['figure.constrained_layout.use'] = True  # Automatically adjust figure layout
mpl.rcParams['axes.linewidth'] = 2                    # Set global axes line width
plt.rcParams['font.family'] = 'sans-serif'           # Set global font family
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica' , 'Verdana', 'Liberation Sans']  # Preferred fonts

# Default font style for plots
font = {'color':  'black',
        'weight': 'normal',
        'size': 20,
        }
