# This file calculates the self-energy components for bulk systems.
# It includes functions for short-range, point-charge, and long-range self-energy.

from packages import *
import calculate
from numerical_param import *
import greens_function_bulk
import concurrent.futures

def uself_short_single(n_position, rad_ions, valency, epsilon):
    """Short-range self-energy for Gaussian charge spread at one position"""
    kappa = calculate.kappa_loc(n_position, valency, epsilon)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(kappa * rad_ions * special.erfc(kappa * rad_ions / np.sqrt(np.pi)))
    C = np.exp(Clog)
    return (np.power(valency, 2) / (8 * np.pi * epsilon * rad_ions)) * (1 - C)

def uself_point(n_bulk_profile, valency, epsilon):
    """Point-charge self-energy profile"""
    kappas = calculate.kappa_profile(n_bulk_profile, valency, epsilon)
    return (kappas[:, np.newaxis] / (8 * np.pi * epsilon)) * np.power(valency, 2)

def uself_short(n_bulk_profile, rad_ions, valency, epsilon):
    """Short-range self-energy profile for Gaussian charges"""
    return np.apply_along_axis(uself_short_single, 1, n_bulk_profile, rad_ions, valency, epsilon)

def uself_component_bulk(n_bulk_profile, n_bulk, valency, quad_point, domain, epsilon):
    """Integration component of long-range self-energy for bulk"""
    G_full = greens_function_bulk.Gcap_full(n_bulk_profile, n_bulk, valency, quad_point[0], domain, epsilon)
    G_free = np.ones(len(n_bulk_profile)) * (1 / (2 * epsilon * quad_point[0]))
    G_component = G_full - G_free
    return quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:, np.newaxis]

def uself_long_bulk(n_bulk_profile, n_bulk, valency, domain, epsilon):
    """Long-range component of self-energy q²(G-G0) with parallelized quadrature"""
    u_long = np.zeros((len(n_bulk_profile), len(valency)))
    samples, weights = np.polynomial.legendre.leggauss(quads)
    S = np.exp(0.5 * V_conv * samples + 0.5 * V_conv) - 1
    v1 = 0.5 * V_conv
    v2 = 0.5 * V_conv
    weights = v1 * (np.exp(v1 * samples + v2) - 1) * np.exp(v1 * samples + v2) * weights
    quad_points = np.c_[S, weights]

    chunk_size = max(1, len(quad_points) // cores)
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        futures = []
        for i in range(0, len(quad_points), chunk_size):
            chunk = quad_points[i:i + chunk_size]
            for qp in chunk:
                futures.append(executor.submit(uself_component_bulk, n_bulk_profile, n_bulk, valency, qp, domain, epsilon))
        for f in concurrent.futures.as_completed(futures):
            u_long += f.result()
    return u_long

def uselfb_numerical(n_bulk_profile, n_bulk, rad_ions, valency, domain, epsilon):
    """Total self-energy of all ions in bulk: short + point + long range"""
    u_short = uself_short(n_bulk_profile, rad_ions, valency, epsilon)
    u_pc = uself_point(n_bulk_profile, valency, epsilon)
    u_long = uself_long_bulk(n_bulk_profile, n_bulk, valency, domain, epsilon)
    return u_short + u_pc + u_long
