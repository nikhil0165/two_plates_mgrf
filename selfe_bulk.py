"""
Computes self-energy corrections for ions in the bulk (single-plate or unconfined) MGRF system.
Includes short-range, point-charge, and long-range (Green's function) contributions.
"""

from packages import *
import calculate
from numerical_param import *
import greens_function_bulk
import concurrent.futures

# ----------------------------
# Short-range self-energy for one grid point
# ----------------------------
def uself_short_single(n_position, rad_ions, valency, epsilon):
    """
    Short-range self-energy using Gaussian charge spread approximation.
    n_position: local ion concentrations
    rad_ions: ionic radii
    valency: ionic charges
    epsilon: dielectric permittivity
    """
    kappa = calculate.kappa_loc(n_position, valency, epsilon)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(kappa * rad_ions * special.erfc(kappa * rad_ions / np.sqrt(np.pi)))
    C = np.exp(Clog)
    return (np.power(valency, 2) / (8 * np.pi * epsilon * rad_ions)) * (1 - C)

# ----------------------------
# Point-charge self-energy
# ----------------------------
def uself_point(n_bulk_profile, valency, epsilon):
    """
    Point-charge self-energy due to local ionic screening.
    n_bulk_profile: array of ion concentrations along grid
    """
    kappas = calculate.kappa_profile(n_bulk_profile, valency, epsilon)
    return (kappas[:, np.newaxis] / (8 * np.pi * epsilon)) * np.power(valency, 2)

# ----------------------------
# Short-range self-energy for bulk profile
# ----------------------------
def uself_short(n_bulk_profile, rad_ions, valency, epsilon):
    """
    Applies short-range self-energy calculation across all grid points.
    """
    return np.apply_along_axis(uself_short_single, 1, n_bulk_profile, rad_ions, valency, epsilon)

# ----------------------------
# Long-range self-energy component for one quadrature point
# ----------------------------
def uself_component_bulk(n_bulk_profile, n_bulk, valency, quad_point, domain, epsilon):
    """
    Integrand for long-range self-energy (q**2(G-G0)) using Green's function.
    quad_point: (s, weight) from Gauss-Legendre quadrature
    """
    G_full = greens_function_bulk.Gcap_full(n_bulk_profile, n_bulk, valency, quad_point[0], domain, epsilon)
    G_free = np.ones(len(n_bulk_profile)) * (1 / (2 * epsilon * quad_point[0]))
    G_component = G_full - G_free
    return quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:, np.newaxis]

# ----------------------------
# Long-range self-energy for bulk profile
# ----------------------------
def uself_long_bulk(n_bulk_profile, n_bulk, valency, domain, epsilon):
    """
    Long-range self-energy computed in parallel over quadrature points.
    Uses Gauss-Legendre quadrature in log-space.

    Globals used: quads, cores, V_conv, e
    """
    u_long = np.zeros((len(n_bulk_profile), len(valency)))
    samples, weights = np.polynomial.legendre.leggauss(quads)
    S = np.power(e, 0.5 * V_conv * samples + 0.5 * V_conv) - 1
    v1 = v2 = 0.5 * V_conv
    weights = v1 * (np.exp(v1 * samples + v2) - 1) * np.exp(v1 * samples + v2) * weights
    quad_points = np.c_[S, weights]

    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [executor.submit(uself_component_bulk, n_bulk_profile, n_bulk, valency, qp, domain, epsilon) for qp in quad_points]
        for f in concurrent.futures.as_completed(futures):
            u_long += f.result()
    return u_long

# ----------------------------
# Complete bulk self-energy
# ----------------------------
def uselfb_numerical(n_bulk_profile, n_bulk, rad_ions, valency, domain, epsilon):
    """
    Total bulk self-energy combining:
    - Short-range
    - Point-charge
    - Long-range
    """
    u_short = uself_short(n_bulk_profile, rad_ions, valency, epsilon)
    u_pc = uself_point(n_bulk_profile, valency, epsilon)
    u_long = uself_long_bulk(n_bulk_profile, n_bulk, valency, domain, epsilon)
    return u_short + u_pc + u_long
