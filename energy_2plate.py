# This file calculates the grand free energy (GFE) for a two-plate system using 
# modified Poisson-Boltzmann (MGRF) theory. It contains functions for both
# plate-specific and bulk contributions.

from packages import *
from numerical_param import *
import selfe_2plate
import selfe_bulk

def grandfe_mgrf_2plate(psi, n_profile, uself_profile,n_bulk, valency,rad_ions, vol_ions, vol_sol, sigma_1,sigma_2, domain, epsilon_s, epsilon_p):
    """
    Compute the grand free energy for two plates using MGRF theory.

    Parameters:
    psi           : electrostatic potential profile
    n_profile     : ion concentration profile
    uself_profile : self-energy profile
    n_bulk        : bulk concentrations
    valency       : ion valencies
    rad_ions      : ion radii
    vol_ions      : ion volumes
    vol_sol       : solvent volume
    sigma_1,2     : plate surface charges
    domain        : plate separation
    epsilon_s,p   : dielectric permittivities

    Returns:
    grandfe       : grand free energy relative to bulk
    """

    # Contribution from plate surfaces
    grandfe = 0.5*psi[0]*sigma_1 + 0.5*psi[-1]*sigma_2

    # Exclude extended nodes used for numerical padding
    N_exc = np.nonzero(n_profile[:,0])[0][0]
    nodes= len(psi)
    psi = psi[N_exc:nodes-N_exc]
    n_profile = n_profile[N_exc:nodes-N_exc]
    uself_profile = uself_profile[N_exc:nodes-N_exc]

    # Bulk grand free energy using uniform concentrations
    nodes = len(n_profile)-1
    n_bulk_profile = np.multiply(np.ones((nodes, len(valency))), n_bulk)
    grandfe_bulk = grandfe_mgrf_bulk(n_bulk_profile,n_bulk, valency,rad_ions, vol_ions,vol_sol, domain, epsilon_s)
    utau = np.zeros((nodes+1, len(valency)))
    taus, weights = np.polynomial.legendre.leggauss(grandfe_quads)

    print('bulk grandfe done ')

    # Spectral setup for spatial integration
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'], size = len(n_profile), bounds = (0,domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    # Midpoint averages for integration
    n_local = 0.5 * (n_profile[:-1] + n_profile[1:])
    psi_local = 0.5 * (psi[:-1] + psi[1:])
    u_local = 0.5 * (uself_profile[:-1] + uself_profile[1:])
    vol_local = np.sum(vol_ions * n_local, axis=1)

    # Contributions from electrostatics, ideal entropy, and steric effects
    grandfe = grandfe - 0.5 * np.sum(psi_local * np.dot(valency, n_local.T) * dz)
    grandfe = grandfe - np.sum(n_local*dz[:,np.newaxis])
    grandfe = grandfe - (1 / vol_sol) * np.sum((1 - vol_local) * dz)
    grandfe = grandfe + (1 / vol_sol) * np.sum(np.log(1 - vol_local) * dz)

    # Scale factors for Gaussian quadrature integration (τ ∈ [0,1])
    tau_scales = 0.5 * (taus + 1)

    # Quadrature integration of self-energy terms
    for k, tau_scale in enumerate(tau_scales):
        scaled_n_profile = tau_scale * n_profile
        scaled_n_bulk = tau_scale * n_bulk
        utau += 0.5 * weights[k] * selfe_2plate.uself_complete(
            scaled_n_profile, scaled_n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p
        )

    # Average self-energy for integration
    utau_local = 0.5 * (utau[:-1] + utau[1:])
    grandfe = grandfe + np.sum(n_local * utau_local * dz[:, np.newaxis])
    grandfe = grandfe - np.sum(n_local * u_local * dz[:, np.newaxis])

    # Return grand free energy relative to uniform bulk
    return grandfe - grandfe_bulk

def grandfe_mgrf_bulk(n_bulk_profile,n_bulk, valency,rad_ions,vol_ions, vol_sol, domain, epsilon):
    """
    Compute bulk contribution to the grand free energy using MGRF theory.
    """

    grandfe = 0
    nodes = len(n_bulk_profile)
    vol_bulk = sum([n_bulk[i] * vol_ions[i] for i in range(len(vol_ions))])
    u_bulk = selfe_bulk.uselfb_numerical(n_bulk_profile, n_bulk, rad_ions, valency, domain, epsilon)
    utau_bulk = np.zeros_like(u_bulk)
    taus, weights = np.polynomial.legendre.leggauss(grandfe_quads)

    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'], size = nodes+1, bounds = (0,domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    # Contributions from ideal entropy and solvent volume
    grandfe = grandfe - np.sum(n_bulk * dz[:, np.newaxis])
    grandfe = grandfe - np.sum(dz*(1/vol_sol)*(1 - vol_bulk))
    grandfe = grandfe + np.sum(dz*(1/vol_sol)*np.log(1-vol_bulk))

    # Quadrature integration of bulk self-energy terms
    tau_scales = 0.5 * (taus + 1)
    for k, tau_scale in enumerate(tau_scales):
        scaled_n_bulk_profile = tau_scale * n_bulk_profile
        scaled_n_bulk = tau_scale * n_bulk
        utau_bulk += 0.5 * weights[k] * selfe_bulk.uselfb_numerical(
            scaled_n_bulk_profile, scaled_n_bulk, rad_ions, valency, domain, epsilon
        )

    grandfe = grandfe + np.sum(n_bulk * utau_bulk * dz[:,np.newaxis])
    grandfe = grandfe - np.sum(n_bulk * u_bulk * dz[:,np.newaxis])

    return grandfe

def grandfe_pb_2plate(psi, n_profile,n_bulk, valency, sigma_1,sigma_2, domain):
    """
    Grand free energy using classical Poisson-Boltzmann (PB) theory for two plates.
    """

    grandfe = 0.5*psi[0]*sigma_1 + 0.5*psi[-1]*sigma_2
    nodes = len(n_profile)-1
    n_bulk_profile = np.multiply(np.ones((nodes, len(valency))), n_bulk)
    grandfe_bulk = grandfe_pb_bulk(n_bulk_profile,n_bulk, valency, domain)

    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'], size = len(n_profile), bounds = (0,domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    n_local = 0.5 * (n_profile[:-1] + n_profile[1:])
    psi_local = 0.5 * (psi[:-1] + psi[1:])

    grandfe = grandfe - 0.5 * np.sum(psi_local * np.dot(valency, n_local.T) * dz)
    grandfe = grandfe - np.sum(n_local*dz[:,np.newaxis])

    return grandfe - grandfe_bulk

def grandfe_pb_bulk(n_bulk_profile,n_bulk, valency, domain):
    """
    Bulk grand free energy using classical mean-field PB theory.
    """
    grandfe = 0
    nodes = len(n_bulk_profile)

    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'], size = nodes+1, bounds = (0,domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)
    grandfe = grandfe - np.sum(n_bulk[:, np.newaxis] * dz)

    return grandfe
