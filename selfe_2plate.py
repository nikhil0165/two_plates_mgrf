# This file calculates the self-energy components for two plates.
# It includes functions for Born solvation energy, short-range, point-charge, and long-range self-energy.

from packages import *
import calculate
from numerical_param import *
import greens_function_2plate
import concurrent.futures

def uself_born(rad_ions, valency, epsilon_s):
    """Born solvation energy of all the ions"""
    return np.power(valency, 2) / (8 * np.pi * epsilon_s * rad_ions)

def uself_short_single(n_position, rad_ions, valency, epsilon_s):
    """Short-range self-energy for Gaussian charge spread at one position"""
    kappa = calculate.kappa_loc(n_position, valency, epsilon_s)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(special.erfc(kappa * rad_ions / np.sqrt(np.pi)) * kappa * rad_ions)
    C = np.exp(Clog)
    u_short = (np.power(valency, 2) / (8 * np.pi * epsilon_s * rad_ions)) * (1 - C)
    return u_short

def uself_point(n_profile, valency, epsilon_s):
    """Point-charge self-energy profile"""
    kappas = calculate.kappa_profile(n_profile, valency, epsilon_s)
    return (kappas[:, np.newaxis] / (8 * np.pi * epsilon_s)) * np.power(valency, 2)

def uself_short(n_profile, rad_ions, valency, epsilon_s):
    """Short-range self-energy profile for Gaussian charges"""
    return np.apply_along_axis(uself_short_single, 1, n_profile, rad_ions, valency, epsilon_s)

def uself_component(n_profile, n_bulk, valency, quad_point, domain, epsilon_s, epsilon_p, dist_exc):
    """Integration component of long-range self-energy q²(G - G_free)"""
    G_full = greens_function_2plate.Gcap_full(n_profile, n_bulk, valency, quad_point[0], domain, epsilon_s, epsilon_p, dist_exc)
    G_free = np.ones(len(n_profile)) * (1 / (2 * epsilon_s * quad_point[0]))
    G_component = G_full - G_free
    return quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:, np.newaxis]

def uself_long(n_profile, n_bulk, valency, domain, epsilon_s, epsilon_p, dist_exc):
    """Long-range component of self-energy q²(G-G0) using parallelized quadrature"""
    u_long = np.zeros((len(n_profile), len(valency)))
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
                futures.append(executor.submit(uself_component, n_profile, n_bulk, valency, qp, domain, epsilon_s, epsilon_p, dist_exc))
        for f in concurrent.futures.as_completed(futures):
            u_long += f.result()

    return u_long

def uself_complete(n_profile, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p):
    """Total self-energy of all ions: short-range + point-charge + long-range"""
    u_short = uself_short(n_profile, rad_ions, valency, epsilon_s)
    u_pc = uself_point(n_profile, valency, epsilon_s)
    u_long = uself_long(n_profile, n_bulk, valency, domain, epsilon_s, epsilon_p, dist_exc=np.max(rad_ions))
    return u_short + u_pc + u_long
