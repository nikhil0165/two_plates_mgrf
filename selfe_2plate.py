from packages import *
import calculate
from numerical_param import*
import greens_function_2plate

def uself_born(rad_ions, valency, epsilon_s): # born solvation energy of all the ions
    born_solv_energies = np.power(valency, 2) * (1 / (8 * np.pi * epsilon_s * rad_ions))
    return born_solv_energies

def uself_short_single(n_position,rad_ions,valency,epsilon_s): #short-range self-energy for gaussian charge spread at one point
    kappa = calculate.kappa_loc(n_position,valency,epsilon_s)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(kappa * rad_ions * special.erfc(kappa * rad_ions / np.sqrt(np.pi)))
    C = np.exp(Clog)
    u_short = (np.power(valency,2) / (8 * np.pi * epsilon_s * rad_ions)) * (1 - C)
    return u_short

def uself_point(n_profile,valency,epsilon_s):
    kappas = calculate.kappa_profile(n_profile,valency,epsilon_s)
    uself_pc = (kappas[:, np.newaxis] / (8 * np.pi * epsilon_s)) * np.power(valency, 2)
    return uself_pc

def uself_short(n_profile,rad_ions,valency, epsilon_s):  #short-range self-energy profile for gaussian charge spread
    u_short = np.apply_along_axis(uself_short_single, 1, n_profile, rad_ions, valency, epsilon_s)# iterates over rows of n_profile
    return u_short

def uself_component(n_profile,n_bulk,valency,quad_point,domain,epsilon_s,epsilon_p,dist_exc): # integration component of u_long as per legendre-gauss quadrature
    G_full= greens_function_2plate.Gcap_full(n_profile,n_bulk,valency,quad_point[0],domain,epsilon_s,epsilon_p,dist_exc)
    G_free =  np.ones(len(n_profile))*(1/(2*epsilon_s*quad_point[0]))#greens_function_2plate.Gcap_free(len(n_profile),quad_point[0],domain,epsilon_s)
    G_component = (G_full-G_free)
    u_component = quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:,np.newaxis]
    return u_component

def uself_long(n_profile,n_bulk,valency,domain,epsilon_s,epsilon_p,dist_exc):# long range component of self-energy q**2(G-Go)
    u_long = np.zeros((len(n_profile),len(valency)))
    samples, weights = np.polynomial.legendre.leggauss(quads)
    S = np.power(e,0.5*V_conv*samples + 0.5*V_conv)-1 # transformation of s into v
    v1 = 0.5*V_conv
    v2 = 0.5*V_conv
    weights = v1 * (np.exp(v1 * samples + v2) - 1) * np.exp(v1 * samples + v2) * weights
    quad_points = np.c_[S,weights]
    for i in range(0,floor(quads/cores)):
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            answer = {executor.submit(uself_component,n_profile,n_bulk,valency,quad_point,domain,epsilon_s,epsilon_p,dist_exc): quad_point for quad_point in quad_points[i*cores:(i+1)*cores]}
            for future in concurrent.futures.as_completed(answer):
                u_long = u_long + future.result()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            answer = {executor.submit(uself_component,n_profile,n_bulk,valency,quad_point,domain,epsilon_s,epsilon_p,dist_exc): quad_point for quad_point in quad_points[floor(quads/cores)*cores:]}
            for future in concurrent.futures.as_completed(answer):
                u_long = u_long + future.result()
    return u_long

def uself_complete(n_profile,n_bulk,rad_ions,valency,domain,epsilon_s,epsilon_p): # self energy of all ions
    u_short = uself_short(n_profile,rad_ions,valency,epsilon_s)
    u_pc = uself_point(n_profile,valency,epsilon_s)
    u_long = uself_long(n_profile,n_bulk,valency,domain,epsilon_s,epsilon_p,dist_exc = np.max(rad_ions))
    u_self = u_short + u_pc + u_long
    return u_self




