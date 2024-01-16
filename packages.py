## Math packages

import numpy as np
from math import*
import dedalus.public as d3
from scipy import special
from scipy.interpolate import CubicSpline

from decimal import*
getcontext().prec=50

## File reading and writing
import sys
import os
import argparse
import importlib
import h5py
import csv
import gc
## For parallelising self-energy calculations

import concurrent.futures

import re
import subprocess
## timing, warnings, logging

import timeit
import warnings
warnings.filterwarnings("ignore")
import logging
logging.getLogger().setLevel(logging.WARNING)

## For plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['figure.constrained_layout.use'] = True
mpl.rcParams['axes.linewidth'] = 2 #set the value globally
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica' , 'Verdana', 'Liberation Sans']
font = {'color':  'black',
        'weight': 'normal',
        'size': 20,
        }
