# This file defines numerical parameters used across various simulations.
# These parameters include grid sizes, tolerances, quadrature points,
# and computational settings for MGRF, PB, and Green's function calculations.

from packages import *
import os

# ------------------------------------------------------
# Fourier Inversion Parameters for Green's function
# ------------------------------------------------------
s_conv = 32           # Approximation of infinity for Fourier inverse of Green's function
V_conv = log(s_conv + 1)  # Integration performed in log-space
quads = 16            # Number of Legendre-Gauss quadrature points for Fourier inverse

# ------------------------------------------------------
# Grid and Dealiasing Parameters
# ------------------------------------------------------
N_grid = 128          # Total number of grid points (should be even, 3/2 rule often used for dealiasing)
N_exc = 10            # Grid points for the exclusion zone near plates
dealias = 2           # Dealiasing factor used in Dedalus spectral computations

# ------------------------------------------------------
# Non-Constant Coefficient Cutoffs for NLBVP
# ------------------------------------------------------
ncc_cutoff_mgrf = 1e-2   # Cutoff for MGRF solver
ncc_cutoff_pb = 1e-1     # Cutoff for Poisson-Boltzmann solver
ncc_cutoff_greens = 1e-1 # Cutoff for Green's function solver

# ------------------------------------------------------
# Mixing and Quadrature Parameters
# ------------------------------------------------------
num_ratio = 0.1          # Mixing ratio of new to old in nconc_mgrf iteration
grandfe_quads = 25       # Number of Legendre-Gauss points for grand free energy integration

# ------------------------------------------------------
# Parallelization
# ------------------------------------------------------
cores = min(8, os.cpu_count())  # Number of parallel processes for Fourier inverse calculation

# ------------------------------------------------------
# Tolerances
# ------------------------------------------------------
tolerance = 1e-5         # Outer loop convergence tolerance for mgrf_2plate
tolerance_pb = 1e-7      # Tolerance for inner loop in PB solver
tolerance_num = 1e-4     # Convergence tolerance for nconc_mgrf iteration
tolerance_greens = 1e-7  # Tolerance for nonlinear Green's function problem

# ------------------------------------------------------
# Iteration Limits
# ------------------------------------------------------
iter_max = 1e7           # Maximum allowed iterations for any iterative loop

# ------------------------------------------------------
# Test output when run as a script
# ------------------------------------------------------
if __name__ == "__main__":
    for var in ['N_grid', 'N_exc', 'quads', 'cores', 'tolerance', 'tolerance_pb', 
                'tolerance_num', 'tolerance_greens', 'num_ratio', 's_conv', 'V_conv',
                'dealias', 'ncc_cutoff_mgrf', 'ncc_cutoff_pb', 'ncc_cutoff_greens', 
                'grandfe_quads', 'iter_max']:
        print(f'{var} = {globals()[var]}')
