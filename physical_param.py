# This file defines the physical parameters for the simulations.
# It includes global input variables, physical constants, and derived variables.

from packages import *
import numpy as np
from math import pi

## ----------------------------
## Global Input Variables (Dimensional, SI units)
## ----------------------------

cb1_d = 0.01  # primary salt bulk concentration [M]
cb2_d = 0.0   # secondary salt bulk concentration [M]

valency1 = [1, -1]  # primary salt valencies
valency2 = [1, -1]  # secondary salt valencies

born_radius1 = 2.5  # cation radius [Å]
born_radius2 = 2.5  # anion radius [Å]
rad_sol_d = max(born_radius1, born_radius2)  # solvent radius

domain_d = 100.0  # plate separation [Å]
domain_in_d = domain_d  # initial guess separation

sigma_f1_d = -0.01  # surface charge plate 1 [C/m²]
sigma_f2_d = -0.01  # surface charge plate 2
sigma_in1_d = sigma_f1_d
sigma_in2_d = sigma_f2_d

vol_sol_d = 4/3 * pi * (rad_sol_d*1e-10)**3  # solvent molecule volume [m³]

## Assign valency, radii, and volumes depending on secondary salt
if cb2_d == 0:
    valency = np.array(valency1)
    rad_ions_d = np.array([born_radius1, born_radius2])
    vol_ions_d = np.array([vol_sol_d, vol_sol_d])
else:
    valency = np.hstack((valency1, valency2))
    rad_ions_d = np.array([born_radius1, born_radius2, born_radius2, born_radius2])
    vol_ions_d = np.array([vol_sol_d]*4)

## ----------------------------
## Physical Constants
## ----------------------------
Temp = 298            # K
ec = 1.602e-19        # elementary charge [C]
k_b = 1.38064852e-23  # Boltzmann constant [J/K]
N_A = 6.02214e23      # Avogadro number
epsilon_o_d = 8.854187e-12  # vacuum permittivity [F/m]
beta = 1 / (k_b * Temp)

epsilonr_s_d = 80
epsilonr_p_d = 80

epsilon_s_d = epsilon_o_d * epsilonr_s_d
epsilon_p_d = epsilon_o_d * epsilonr_p_d

## ----------------------------
## Characteristic scales (non-dimensionalization)
## ----------------------------
l_b = ec**2 * beta / (4*pi*epsilon_s_d)  # Bjerrum length
l_c = l_b
q_c = ec
psi_c = 1 / (beta*q_c)
epsilon_c = beta*q_c**2 / l_c
sigma_c = ec / l_c**2
vol_c = l_c**3
nconc_c = 1/vol_c
conc_c = 1/vol_c

## ----------------------------
## Derived bulk concentrations
## ----------------------------
# Primary ions
cbulk_plus_d, cbulk_neg_d = (cb1_d, cb1_d) if abs(valency[0])==abs(valency[1]) else (cb1_d*abs(valency[1]), cb1_d*abs(valency[0]))
nbulk_plus_d = cbulk_plus_d * N_A * 1e3
nbulk_neg_d = cbulk_neg_d * N_A * 1e3
n_bulk_d = [nbulk_plus_d, nbulk_neg_d]

# Secondary ions (if present)
if cb2_d != 0:
    cb2_plus_d, cb2_neg_d = (cb2_d, cb2_d) if abs(valency[2])==abs(valency[3]) else (cb2_d*abs(valency[3]), cb2_d*abs(valency[2]))
    nb2_plus_d = cb2_plus_d * N_A * 1e3
    nb2_neg_d = cb2_neg_d * N_A * 1e3
    n_bulk_d += [nb2_plus_d, nb2_neg_d]

## Debye Screening length
I = sum([valency[i]**2 * n_bulk_d[i] / len(valency) for i in range(len(valency))])
lambda_d_d = np.sqrt(epsilon_s_d / (beta * ec**2 * I))

## ----------------------------
## Scale variables for non-dimensional simulation
## ----------------------------
epsilon_s = epsilon_s_d / epsilon_c
epsilon_p = epsilon_p_d / epsilon_c
n_bulk = np.true_divide(n_bulk_d, nconc_c)
sigma_i1 = sigma_in1_d / sigma_c
sigma_f1 = sigma_f1_d / sigma_c
sigma_i2 = sigma_in2_d / sigma_c
sigma_f2 = sigma_f2_d / sigma_c
rad_ions = np.true_divide(rad_ions_d*1e-10, l_c)
rad_sol = rad_sol_d*1e-10 / l_c
vol_sol = vol_sol_d / vol_c
vol_ions = np.true_divide(vol_ions_d, vol_c)
domain = domain_d*1e-10 / l_c
domain_in = domain_in_d*1e-10 / l_c
lambda_d = lambda_d_d / l_c

## ----------------------------
## Printing variables if file is run directly
## ----------------------------
if __name__ == "__main__":
    print(f'cb1_d = {cb1_d}')
    print(f'cb2_d = {cb2_d}')
    print(f'sigma_in1_d = {sigma_in1_d}')
    print(f'sigma_in2_d = {sigma_in2_d}')
    print(f'sigma_f1_d = {sigma_f1_d}')
    print(f'sigma_f2_d = {sigma_f2_d}')
    print(f'domain_in_d = {domain_in_d}')
    print(f'domain_d = {domain_d}')
    print(f'valency = {valency}')
    print(f'rad_ions_d = {rad_ions_d}')
    print(f'rad_sol_d = {rad_sol_d}')
    print(f'epsilonr_s_d = {epsilonr_s_d}')
    print(f'epsilonr_p_d = {epsilonr_p_d}')
    print(f'electrostatic_coupling = {(2*pi*pow(valency[0],3)*(l_c**2)*np.sqrt(abs(sigma_f1_d*sigma_f2_d))/ec)}')
