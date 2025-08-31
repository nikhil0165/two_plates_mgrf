"""
Computes self-energy corrections for ions in the two-plate MGRF system.
Includes Born, short-range, point-charge, and long-range contributions.
"""

from packages import *
import calculate
from numerical_param import *
import greens_function_2plate
import concurrent.futures

# ----------------------------
# Born solvation energy
# ----------------------------
def uself_born(rad_ions, valency, epsilon):
    """
    Born solvation energy for all ions.
    rad_ions: ionic radii
    valency: ionic charges
    epsilon: dielectric permittivity
    """
    return np.power(valency, 2) * (1 / (8 * np.pi * epsilon * rad_ions))

# ----------------------------
# Short-range self-energy
# ----------------------------
def uself_short_single(n_position, rad_ions, valency, epsilon):
    """
    Short-range self-energy for a single grid point.
    n_position: local ion concentrations
    Uses Gaussian charge spread approximation.
    """
    kappa = calculate.kappa_loc(n_position, valency, epsilon)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(kappa * rad_ions * special.erfc(kappa * rad_ions / np.sqrt(np.pi)))
    C = np.exp(Clog)
    return (np.power(valency, 2) / (8 * np.pi * epsilon * rad_ions)) * (1 - C)

def uself_short(n_profile, rad_ions, valency, epsilon):
    """Short-range self-energy for all grid points (applies uself_short_single along profile)."""
    return np.apply_along_axis(uself_short_single, 1, n_profile, rad_ions, valency, epsilon)

# ----------------------------
# Point-charge self-energy
# ----------------------------
def uself_point(n_profile, valency, epsilon):
    """Point-charge self-energy due to local screening."""
    kappas = calculate.kappa_profile(n_profile, valency, epsilon)
    return (kappas[:, np.newaxis] / (8 * np.pi * epsilon)) * np.power(valency, 2)

# ----------------------------
# Long-range self-energy
# ----------------------------
def uself_component(n_profile, n_bulk, valency, quad_point, domain, epsilon):
    """
    Long-range self-energy contribution for one quadrature point.
    quad_point: (s, weight) from Gauss-Legendre quadrature
    """
    G_full = greens_function_2plate.Gcap_full(n_profile, n_bulk, valency, quad_point[0], domain, epsilon)
    G_free = np.ones(len(n_profile)) * (1 / (2 * epsilon * quad_point[0]))
    G_component = G_full - G_free
    return quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:, np.newaxis]

def uself_long(n_profile, n_bulk, valency, domain, epsilon):
    """
    Long-range self-energy for the profile, computed in parallel.
    Uses Gauss-Legendre quadrature in log-space.
    """
    u_long = np.zeros((len(n_profile), len(valency)))
    samples, weights = np.polynomial.legendre.leggauss(quads)
    S = np.power(e, 0.5 * V_conv * samples + 0.5 * V_conv) - 1
    v1 = v2 = 0.5 * V_conv
    weights = v1 * (np.exp(v1 * samples + v2) - 1) * np.exp(v1 * samples + v2) * weights
    quad_points = np.c_[S, weights]

    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [executor.submit(uself_component, n_profile, n_bulk, valency, qp, domain, epsilon) for qp in quad_points]
        for f in concurrent.futures.as_completed(futures):
            u_long += f.result()
    return u_long

# ----------------------------
# Complete self-energy
# ----------------------------
def uself_complete(n_profile, n_bulk, rad_ions, valency, domain, epsilon):
    """
    Total self-energy combining:
    - Short-range (Gaussian spread)
    - Point-charge
    - Long-range (Green's function)
    """
    u_short = uself_short(n_profile, rad_ions, valency, epsilon)
    u_pc = uself_point(n_profile, valency, epsilon)
    u_long = uself_long(n_profile, n_bulk, valency, domain, epsilon)
    return u_short + u_pc + u_long
