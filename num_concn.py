from packages import *
from numerical_param import *
import calculate
import selfe_2plate

# -------------------------------
# Concentration Profile Functions
# -------------------------------

# Function: nconc_pb
# Purpose:
#   Compute ionic concentration profiles using the mean-field Poisson–Boltzmann (PB) approximation.
# Inputs:
#   psi_profile      : Electrostatic potential profile (array of length = # grid points)
#   valency  : Array of ion valencies
#   n_bulk   : Bulk ion concentrations
# Outputs:
#   Concentration profile (2D array, shape: [len(psi_profile), num_species])
def nconc_pb(psi_profile, valency, n_bulk):
    return n_bulk * np.exp(-np.array(valency) * psi_profile[:,np.newaxis] )


# Function: nconc_mgrf
# Purpose:
#   Compute ionic concentration profiles and coefficients for the
#   Modified Gaussian Reference Fluid (MGRF) approach.
# Inputs:
#   psi_profile          : Electrostatic potential profile
#   uself        : Self-energy profile
#   eta_profile  : Local packing fraction profile
#   uself_bulk   : Bulk self-energy
#   n_bulk       : Bulk ion concentration
#   valency      : Array of ion valencies
#   vol_ions     : Ion volumes
#   eta_bulk     : Bulk packing fraction
#   equal_vols   : Boolean, True if ion/solvent volumes are equal
# Outputs:
#   n_profile : Computed ion concentration profile
#   coeffs    : Coefficients used in further calculations
def nconc_mgrf(psi_profile,uself,eta_profile,uself_bulk, n_bulk, valency, vol_ions,eta_bulk, equal_vols):
    if equal_vols:
        # Equal ion/solvent volume case
        A = n_bulk* np.exp(-np.array(valency) * psi_profile[:,np.newaxis] - (uself - uself_bulk) + vol_ions * eta_bulk)
        coeffs = valency * n_bulk* np.exp(-(uself - uself_bulk) + vol_ions * eta_bulk)
        denom = 1 + np.sum(A * vol_ions, axis=1)   # Volume correction
        n_profile= np.true_divide(A,denom[:,np.newaxis])
        coeffs = np.true_divide(coeffs,denom[:,np.newaxis])
    else:
        # General case: different ion and solvent volumes
        n_profile = n_bulk * np.exp(-np.array(valency)*psi_profile[:,np.newaxis] - (uself - uself_bulk) - vol_ions * (eta_profile[:,np.newaxis] - eta_bulk))
        coeffs = valency* n_bulk * np.exp(-(uself - uself_bulk) - vol_ions* (eta_profile[:,np.newaxis] - eta_bulk))
    return n_profile,coeffs


# Function: nconc_complete
# Purpose:
#   Iteratively solve for self-consistent ion concentration profile,
#   including mean electrostatic potential, excluded volume, and self-energy corrections.
# Inputs:
#   psi_profile         : Electrostatic potential profile
#   nconc_guess   : Initial guess for ion concentration profile
#   uself_bulk  : Bulk self-energy
#   n_bulk      : Bulk ion concentrations
#   valency     : Ion valencies
#   rad_ions    : Ion radii
#   vol_ions    : Ion volumes
#   vol_sol     : Solvent volume
#   domain      : Domain length (slab thickness)
#   epsilon     : Dielectric constant
# Outputs:
#   n_profile       : Final converged ion concentration profile
#   uself_profile   : Final converged self-energy profile
def nconc_complete(psi_profile, nconc_guess,uself_bulk,n_bulk, valency, rad_ions, vol_ions, vol_sol, domain, epsilon):  # nconc_guess is the initial guess

    # Bulk and profile excluded volume fractions
    eta_bulk = calculate.eta_loc(n_bulk, vol_ions, vol_sol)
    eta_profile = calculate.eta_profile(nconc_guess,vol_ions,vol_sol)
    nodes = len(psi_profile)

    # Initialize concentration and guess
    n_profile = np.copy(nconc_guess)
    n_guess = np.copy(nconc_guess)

    # Initial self-energy profile
    uself_profile = selfe_2plate.uself_complete(n_profile, n_bulk,rad_ions, valency,domain, epsilon)

    # Check if ion and solvent volumes are effectively equal
    equal_vols = np.all(np.abs(vol_ions - vol_sol) < vol_sol * 1e-5)

    # Iteration loop (fixed-point iteration for self-consistency)
    convergence = 1
    p = 0
    while (convergence > tolerance_num) and (p < iter_max):
        p = p + 1

        # Update concentration profile
        if equal_vols:
            A = n_bulk* np.exp(-np.array(valency) * psi_profile[:, np.newaxis] - (uself_profile - uself_bulk) + vol_ions * eta_bulk)
            denom = 1 + np.sum(A * vol_ions, axis=1)
            n_profile = np.true_divide(A, denom[:,np.newaxis])
        else:
            n_profile = n_bulk * np.exp(-np.array(valency)*psi_profile[:,np.newaxis] - (uself_profile - uself_bulk) - vol_ions * (eta_profile[:,np.newaxis] - eta_bulk))

        # Check convergence (relative L2 norm)
        convergence = np.true_divide(np.linalg.norm(n_profile - n_guess),np.linalg.norm(n_guess))

        # Relaxation update of concentration guess
        n_guess = (num_ratio) * n_profile + (1-num_ratio) * n_guess

        # Update self-energy and excluded volume profile
        uself_profile = selfe_2plate.uself_complete(n_guess,n_bulk, rad_ions, valency, domain,epsilon)
        eta_profile = calculate.eta_profile(n_guess,vol_ions,vol_sol)

        # Print debug info every 10 iterations
        if p%10==0:
            print('num='+str(convergence))

        # Fail-safe: prevent infinite loop
        if p >= iter_max:
            print("too many iterations for convergence")

    uself_profile = selfe_2plate.uself_complete(n_profile,n_bulk, rad_ions, valency, domain,epsilon)

    return n_profile, uself_profile
