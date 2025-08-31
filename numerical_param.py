from packages import *

# -------------------------------
# Numerical Parameters
# -------------------------------

# Convergence parameters for Fourier inverse of Green's function
s_conv = 32  # Approximation for infinity in Fourier inverse integration of Green's function
V_conv = log(s_conv + 1)  # Integration performed in log-space for Fourier inverse

# Quadrature settings
quads = 16  # Number of Legendre–Gauss quadrature points for Fourier inverse of Green's function

# Grid and spectral method parameters
N_grid = 128  # Number of Chebyshev grid points (must be even for 2 dealiasing)
dealias = 2  # Dealiasing factor used in Dedalus spectral methods

# Non-constant coefficient (NCC) cutoffs for solving nonlinear boundary value problems
ncc_cutoff_mgrf = 1e-2  # Cutoff for MGRF NLBVP
ncc_cutoff_pb = 1e-1    # Cutoff for PB NLBVP
ncc_cutoff_greens = 1e-1  # Cutoff for Green's function NLBVP

# Mixing ratios for iterative updates
num_ratio = 0.1     # Mixing ratio for new/old ion concentration in nconc_mgrf and mgrf_2plate

# Quadrature for free energy calculations
grandfe_quads = 10  # Number of Legendre–Gauss quadrature points for free energy integration

# Parallelization
cores = min(8, os.cpu_count())  # Number of parallel processes for Fourier inverse calculation

# Tolerances for iterative solvers
tolerance = pow(10, -5)       # Convergence tolerance for outermost PB–MGRF loop
tolerance_pb = pow(10, -7)    # Tolerance for inner PB/MGRF loop
tolerance_num = pow(10, -4)   # Tolerance for nconc_mgrf iteration loop
tolerance_greens = pow(10, -7)  # Tolerance for nonlinear Green's function solver

# Iteration limits
iter_max = pow(10, 7)  # Maximum number of iterations for any iterative loop

if __name__ == "__main__":
    print(f'N_grid = {N_grid}')
    print(f'quads = {quads}')
    print(f'cores = {cores}')
    print(f'tolerance = {tolerance}')
    print(f'tolerance_pb = {tolerance_pb}')
    print(f'tolerance_num = {tolerance_num}')
    print(f'tolerance_greens = {tolerance_greens}')
    print(f'num_ratio = {num_ratio}')
    print(f's_conv = {s_conv}')
    print(f'V_conv = {V_conv}')
    print(f'dealias = {dealias}')
    print(f'ncc_cutoff_mgrf = {ncc_cutoff_mgrf}')
    print(f'ncc_cutoff_pb = {ncc_cutoff_pb}')
    print(f'ncc_cutoff_greens = {ncc_cutoff_greens}')
    print(f'grandfe_quads = {grandfe_quads}')
    print(f'iter_max = {iter_max}')