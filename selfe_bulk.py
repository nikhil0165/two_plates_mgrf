from packages import *
import calculate
from numerical_param import*
import greens_function_bulk

def uself_short_single(n_position,rad_ions,valency,epsilon): #short-range self-energy for gaussian charge spread at one point
    kappa = calculate.kappa_loc(n_position,valency,epsilon)
    Clog = (kappa * rad_ions) ** 2 / np.pi + np.log(kappa * rad_ions * special.erfc(kappa * rad_ions / np.sqrt(np.pi)))
    C = np.exp(Clog)
    u_short = (np.power(valency,2) / (8 * np.pi * epsilon * rad_ions)) * (1 - C)
    return u_short
    
def uself_point(n_bulk_profile,valency,epsilon):
    kappas = calculate.kappa_profile(n_bulk_profile,valency,epsilon)
    uself_pc = (kappas[:, np.newaxis] / (8 * np.pi * epsilon)) * np.power(valency, 2)
    return uself_pc

def uself_short(n_bulk_profile,rad_ions,valency, epsilon): #short-range self-energy profile for gaussian charge spread
    u_short = np.apply_along_axis(uself_short_single, 1, n_bulk_profile, rad_ions, valency, epsilon)# iterates over rows of n_profile
    return u_short

def uself_component_bulk(n_bulk_profile,n_bulk,valency,quad_point,domain,epsilon): # integration component as per legendre-gauss quadrature
    G_full= greens_function_bulk.Gcap_full(n_bulk_profile,n_bulk,valency,quad_point[0],domain,epsilon)
    G_free = np.ones(len(n_bulk_profile))*(1/(2*epsilon*quad_point[0]))#greens_function_bulk.Gcap_free(len(n_bulk_profile),quad_point[0],domain,epsilon)
    G_component = (G_full-G_free)
    u_component = quad_point[1] * (np.power(valency, 2) / (4 * np.pi)) * G_component[:,np.newaxis]
    return u_component

def uself_long_bulk(n_bulk_profile,n_bulk,valency,domain,epsilon): # long range component of self-energy q**2(G-Go)
    u_long = np.zeros((len(n_bulk_profile),len(valency)))
    samples, weights = np.polynomial.legendre.leggauss(quads)
    S = np.power(e,0.5*V_conv*samples + 0.5*V_conv)-1 # transformation of s into v
    v1 = 0.5*V_conv
    v2 = 0.5*V_conv
    weights = v1 * (np.exp(v1 * samples + v2) - 1) * np.exp(v1 * samples + v2) * weights
    quad_points = np.c_[S,weights]
    for i in range(0,floor(quads/cores)):
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            answer = {executor.submit(uself_component_bulk,n_bulk_profile,n_bulk,valency,quad_point,domain,epsilon): quad_point for quad_point in quad_points[i*cores:(i+1)*cores]}
            for future in concurrent.futures.as_completed(answer):
                u_long = u_long + future.result()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            answer = {executor.submit(uself_component_bulk,n_bulk_profile,n_bulk,valency,quad_point,domain,epsilon): quad_point for quad_point in quad_points[floor(quads/cores)*cores:]}
            for future in concurrent.futures.as_completed(answer):
                u_long = u_long + future.result()
    return u_long

def uselfb_numerical(n_bulk_profile,n_bulk,rad_ions,valency,domain,epsilon): # self energy of all ions
    u_short = uself_short( n_bulk_profile,rad_ions,valency,epsilon)
    u_pc = uself_point(n_bulk_profile,valency,epsilon)
    u_long = uself_long_bulk(n_bulk_profile,n_bulk,valency,domain,epsilon) 
    u_self = u_short + u_pc + u_long
    return u_self

