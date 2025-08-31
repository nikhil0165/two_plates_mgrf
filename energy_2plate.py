"""
Calculates the grand free energy for the two-plate MGRF system and the corresponding bulk reference.
Includes quadrature integration for self-energy and excluded volume effects.
"""
from packages import *
from numerical_param import *
import selfe_2plate
import selfe_bulk


def grandfe_mgrf_2plate(psi, n_profile, uself_profile,
                        n_bulk, valency, rad_ions, vol_ions, vol_sol,
                        sigma_1, sigma_2, domain, epsilon):
    """
    Compute the grand free energy of the two-plate system with MGRF corrections.
    Includes electrostatics, excluded volume effects, and self-energy contributions.
    """

    # Electrostatic energy contribution from surface charges
    grandfe = 0.5*psi[0]*sigma_1 + 0.5*psi[-1]*sigma_2

    # Number of spatial nodes used for discretization
    nodes = len(n_profile)-1

    # Construct a bulk density profile for reference state
    n_bulk_profile = np.multiply(np.ones((nodes, len(valency))), n_bulk)

    # Compute the bulk grand free energy reference
    grandfe_bulk = grandfe_mgrf_bulk(n_bulk_profile, n_bulk, valency,
                                     rad_ions, vol_ions, vol_sol, domain, epsilon)

    # Initialize quadrature accumulation of self-energy contributions
    utau = np.zeros((nodes+1, len(valency)))

    # Gauss–Legendre quadrature points and weights
    taus, weights = np.polynomial.legendre.leggauss(grandfe_quads)

    # Build Dedalus spectral grid
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)  
    zbasis = d3.Chebyshev(coords['z'], size=len(n_profile), bounds=(0, domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)  # cell widths

    # Midpoint (cell-centered) fields for numerical integration
    n_local = 0.5 * (n_profile[:-1] + n_profile[1:])
    psi_local = 0.5 * (psi[:-1] + psi[1:])
    u_local = 0.5 * (uself_profile[:-1] + uself_profile[1:])
    vol_local = np.sum(vol_ions * n_local, axis=1)

    # --- Free energy contributions ---
    # Electrostatic coupling between ψ and ion charge density
    grandfe = grandfe - 0.5 * np.sum(psi_local * np.dot(valency, n_local.T) * dz)

    # Ideal gas entropy term
    grandfe = grandfe - np.sum(n_local * dz[:, np.newaxis])

    # Excluded volume correction (linear + logarithmic terms)
    grandfe = grandfe - (1 / vol_sol) * np.sum((1 - vol_local) * dz)
    grandfe = grandfe + (1 / vol_sol) * np.sum(np.log(1 - vol_local) * dz)

    # Scale factors for quadrature (τ ∈ [0,1])
    tau_scales = 0.5 * (taus + 1)

    # Quadrature integration of self-energy terms
    for k, tau_scale in enumerate(tau_scales):
        scaled_n_profile = tau_scale * n_profile
        scaled_n_bulk = tau_scale * n_bulk
        utau = utau + 0.5 * weights[k] * selfe_2plate.uself_complete(
            scaled_n_profile, scaled_n_bulk, rad_ions, valency, domain, epsilon
        )

    # Midpoint self-energy fields
    utau_local = 0.5 * (utau[:-1] + utau[1:])

    # Add self-energy contributions
    grandfe = grandfe + np.sum(n_local * utau_local * dz[:, np.newaxis])
    grandfe = grandfe - np.sum(n_local * u_local * dz[:, np.newaxis])

    # Return difference relative to bulk reference
    return grandfe - grandfe_bulk


def grandfe_mgrf_bulk(n_bulk_profile, n_bulk,
                      valency, rad_ions, vol_ions, vol_sol, domain, epsilon):
    """
    Compute the bulk grand free energy reference with MGRF corrections.
    Includes excluded volume and self-energy contributions.
    """

    grandfe = 0
    nodes = len(n_bulk_profile)

    # Bulk excluded volume fraction
    vol_bulk = sum([n_bulk[i] * vol_ions[i] for i in range(len(vol_ions))])

    # Self-energy in the bulk
    u_bulk = selfe_bulk.uselfb_numerical(n_bulk_profile, n_bulk,
                                         rad_ions, valency, domain, epsilon)

    utau_bulk = np.zeros_like(u_bulk)
    taus, weights = np.polynomial.legendre.leggauss(grandfe_quads)

    # Build Dedalus spectral grid
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=nodes+1, bounds=(0, domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    # Ideal gas entropy term
    grandfe = grandfe - np.sum(n_bulk * dz[:, np.newaxis])

    # Excluded volume correction
    grandfe = grandfe - np.sum(dz * (1/vol_sol) * (1 - vol_bulk))
    grandfe = grandfe + np.sum(dz * (1/vol_sol) * np.log(1 - vol_bulk))

    # Quadrature integration of bulk self-energy
    tau_scales = 0.5 * (taus + 1)
    for k, tau_scale in enumerate(tau_scales):
        scaled_n_profile = tau_scale * n_bulk_profile
        scaled_n_bulk = tau_scale * n_bulk
        utau_bulk = utau_bulk + 0.5 * weights[k] * selfe_bulk.uselfb_numerical(
            scaled_n_profile, scaled_n_bulk, rad_ions, valency, domain, epsilon
        )

    # Add bulk self-energy contributions
    grandfe = grandfe + np.sum(n_bulk * utau_bulk * dz[:, np.newaxis])
    grandfe = grandfe - np.sum(n_bulk * u_bulk * dz[:, np.newaxis])

    return grandfe


def grandfe_pb_2plate(psi, n_profile, n_bulk,
                      valency, sigma_1, sigma_2, domain):
    """
    Compute the grand free energy for the two-plate system under 
    standard Poisson–Boltzmann (PB) theory (no self-energy or volume effects).
    """

    # Electrostatic contribution from surface charges
    grandfe = 0.5*psi[0]*sigma_1 + 0.5*psi[-1]*sigma_2

    # Number of spatial nodes
    nodes = len(n_profile)-1

    # Construct uniform bulk profile
    n_bulk_profile = np.multiply(np.ones((nodes, len(valency))), n_bulk)

    # Compute bulk PB reference
    grandfe_bulk = grandfe_pb_bulk(n_bulk_profile, n_bulk, valency, domain)

    # Build Dedalus spectral grid
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=len(n_profile), bounds=(0, domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    # Midpoint profiles for integration
    n_local = 0.5 * (n_profile[:-1] + n_profile[1:])
    psi_local = 0.5 * (psi[:-1] + psi[1:])

    # Electrostatic and ideal entropy terms
    grandfe = grandfe - 0.5 * np.sum(psi_local * np.dot(valency, n_local.T) * dz)
    grandfe = grandfe - np.sum(n_local * dz[:, np.newaxis])

    # Return difference relative to bulk
    return grandfe - grandfe_bulk


def grandfe_pb_bulk(n_bulk_profile, n_bulk, valency, domain):
    """
    Compute the bulk grand free energy reference under PB theory.
    """

    grandfe = 0
    nodes = len(n_bulk_profile)

    # Build Dedalus spectral grid
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=nodes+1, bounds=(0, domain))
    z = np.squeeze(dist.local_grids(zbasis))
    dz = np.diff(z)

    # Ideal gas entropy term
    grandfe = grandfe - np.sum(n_bulk[:, np.newaxis] * dz)

    return grandfe
