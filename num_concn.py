# This file contains functions to calculate ion concentration profiles.
# It includes methods for both mean-field Poisson-Boltzmann (PB) and MGRF approaches.

from packages import *
from numerical_param import *
import calculate
import selfe_2plate

# --------------------------------------------------------------
# Mean-field Poisson-Boltzmann concentration profile
# --------------------------------------------------------------
def nconc_pb(psi_profile, valency, n_bulk):
    """
    Compute ion concentrations for a given electrostatic potential using
    the mean-field Poisson-Boltzmann approximation.

    Parameters:
    psi_profile : electrostatic potential profile (array of size grid_points)
    valency     : list/array of ion valencies
    n_bulk      : bulk ion concentrations

    Returns:
    n_profile   : ion concentration profile [grid_points x num_ions]
    """
    return n_bulk * np.exp(-np.array(valency) * psi_profile[:, np.newaxis])


# --------------------------------------------------------------
# MGRF method: returns concentration profile and coefficients for PDE
# --------------------------------------------------------------
def nconc_mgrf(psi_profile, uself_profile, eta_profile, uself_bulk,
               n_bulk, valency, vol_ions, eta_bulk, equal_vols):
    """
    Compute ion concentrations and coefficients for MGRF solver.

    Parameters:
    psi_profile   : electrostatic potential profile
    uself_profile : self-energy profile at each grid point
    eta_profile   : local packing fraction profile
    uself_bulk    : bulk self-energy
    n_bulk        : bulk concentrations
    valency       : ion valencies
    vol_ions      : ion volumes
    eta_bulk      : bulk packing fraction
    equal_vols    : flag, True if all ion volumes are equal

    Returns:
    n_profile     : ion concentrations at each grid point
    coeffs        : coefficients for PDE RHS
    """
    if equal_vols:
        # Simplified formula when all ions have equal volume
        A = n_bulk * np.exp(-np.array(valency) * psi_profile[:, np.newaxis]
                             - (uself_profile - uself_bulk) + vol_ions * eta_bulk)
        coeffs = valency * n_bulk * np.exp(-(uself_profile - uself_bulk) + vol_ions * eta_bulk)
        denom = 1 + np.sum(A * vol_ions, axis=1)
        n_profile = np.true_divide(A, denom[:, np.newaxis])
        coeffs = np.true_divide(coeffs, denom[:, np.newaxis])
    else:
        # General case for ions with unequal volumes
        n_profile = n_bulk * np.exp(-np.array(valency) * psi_profile[:, np.newaxis]
                                    - (uself_profile - uself_bulk)
                                    - vol_ions * (eta_profile[:, np.newaxis] - eta_bulk))
        coeffs = valency * n_bulk * np.exp(-(uself_profile - uself_bulk)
                                           - vol_ions * (eta_profile[:, np.newaxis] - eta_bulk))
    return n_profile, coeffs


# --------------------------------------------------------------
# Complete concentration profile solver for given psi
# --------------------------------------------------------------
def nconc_complete(psi_profile, nconc_guess, uself_bulk, n_bulk, valency,
                   rad_ions, vol_ions, vol_sol, domain, epsilon_s, epsilon_p):
    """
    Iteratively compute self-consistent ion concentrations and self-energies
    for a given electrostatic potential using the MGRF approach.

    Parameters:
    psi_profile    : electrostatic potential profile
    nconc_guess    : initial guess for ion concentrations
    uself_bulk     : bulk self-energy
    n_bulk         : bulk concentrations
    valency        : ion valencies
    rad_ions       : ion radii
    vol_ions       : ion volumes
    vol_sol        : solvent volume
    domain         : system size
    epsilon_s,p    : solvent and plate dielectric constants

    Returns:
    n_profile      : converged ion concentrations
    uself_profile  : converged self-energy profile
    """

    # Bulk and local packing fraction
    eta_bulk = calculate.eta_loc(n_bulk, vol_ions, vol_sol)
    eta_profile = calculate.eta_profile(nconc_guess, vol_ions, vol_sol)

    # Initialize profiles
    n_profile = np.copy(nconc_guess)
    n_guess = np.copy(nconc_guess)

    # Initial self-energy profile
    uself_profile = selfe_2plate.uself_complete(n_profile, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p)
    equal_vols = np.all(np.abs(vol_ions - vol_sol) < vol_sol * 1e-5)

    # Iterative solver
    convergence = 1
    p = 0
    while convergence > tolerance_num and p < iter_max:
        p += 1
        if equal_vols:
            # Equal-volume ions
            A = n_bulk * np.exp(-np.array(valency) * psi_profile[:, np.newaxis]
                                - (uself_profile - uself_bulk) + vol_ions * eta_bulk)
            denom = 1 + np.sum(A * vol_ions, axis=1)
            n_profile = np.true_divide(A, denom[:, np.newaxis])
        else:
            # General unequal-volume case
            n_profile = n_bulk * np.exp(-np.array(valency) * psi_profile[:, np.newaxis]
                                        - (uself_profile - uself_bulk)
                                        - vol_ions * (eta_profile[:, np.newaxis] - eta_bulk))

        # Check convergence
        convergence = np.true_divide(np.linalg.norm(n_profile - n_guess), np.linalg.norm(n_guess))

        # Relaxation update
        n_guess = num_ratio * n_profile + (1 - num_ratio) * n_guess

        # Update self-energy and eta profiles
        uself_profile = selfe_2plate.uself_complete(n_guess, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p)
        eta_profile = calculate.eta_profile(n_guess, vol_ions, vol_sol)

        if p % 10 == 0:
            print('num=' + str(convergence))
        if p >= iter_max:
            print("too many iterations for convergence")

    # Final self-energy update
    uself_profile = selfe_2plate.uself_complete(n_profile, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p)

    return n_profile, uself_profile
