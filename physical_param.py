from os import name
from packages import *

# =============================
# Global Input Variables (All quantities in SI units)
# Variables ending with _d denote dimensional (physical) units
# =============================

# -------------------------------
# Bulk salt concentrations (Molar)
# -------------------------------
cb1_d = 0.01          # Primary salt bulk concentration
cb2_d = 0.0           # Secondary salt bulk concentration (set to 0 if only one salt present)

# -------------------------------
# Ionic valencies
# -------------------------------
valency1 = [1, -1]    # Valency of primary salt (cation, anion)
valency2 = [1, -1]    # Valency of secondary salt (cation, anion)

# -------------------------------
# Ionic sizes (Born radii) in Angstroms
# -------------------------------
born_radius1 = 2.5     # Radius of cation
born_radius2 = 2.5     # Radius of anion
rad_sol_d = max(born_radius1, born_radius2)  # Solvent radius (max ion radius used as approximation)

# -------------------------------
# Plate separation and surface charges
# -------------------------------
domain_d = 100.0        # Separation between two plates in Angstroms
domain_in_d = domain_d # Initial separation for starting calculation

sigma_f1_d = -0.01     # Surface charge density of plate 1
sigma_f2_d = -0.01      # Surface charge density of plate 2
sigma_in1_d = -0.01    # Initial guess for high surface charge
sigma_in2_d = -0.01     # Initial guess for high surface charge

# -------------------------------
# Solvent volume (assumed spherical)
# -------------------------------
vol_sol_d = 4/3 * pi * pow(rad_sol_d * 1e-10, 3)  # Volume in m^3

# -------------------------------
# Assign arrays based on whether secondary salt exists
# -------------------------------
if cb2_d == 0:
    valency = np.array(valency1)
    rad_ions_d = np.array([born_radius1, born_radius2])
    vol_ions_d = np.array([vol_sol_d, vol_sol_d])
else:
    valency = np.hstack((valency1, valency2))
    rad_ions_d = np.array([born_radius1, born_radius2, born_radius2, born_radius2])
    vol_ions_d = np.array([vol_sol_d, vol_sol_d, vol_sol_d, vol_sol_d])

# =============================
# Physical constants
# =============================
Temp = 298               # Temperature in Kelvin
ec = 1.602e-19           # Elementary charge in Coulombs
k_b = 1.38064852e-23     # Boltzmann constant in J/K
N_A = 6.02214e23         # Avogadro number
epsilon_o_d = 8.854187e-12  # Vacuum permittivity (F/m)
beta = 1 / (k_b * Temp)     # 1/(k_B*T)
epsilonr_s_d = 80        # Relative permittivity of solvent
epsilon_s_d = epsilon_o_d * epsilonr_s_d  # Absolute permittivity

# -------------------------------
# Characteristic scales for non-dimensionalization
# -------------------------------
l_b = pow(ec, 2) * beta / (4 * pi * epsilon_s_d)  # Bjerrum length
l_c = l_b                                         # Characteristic length scale
q_c = ec                                         # Characteristic charge
psi_c = 1 / (beta * q_c)                         # Characteristic potential
epsilon_c = beta * q_c * q_c / l_c              # Characteristic permittivity
sigma_c = ec / l_c**2                            # Characteristic surface charge
vol_c = l_c**3                                   # Characteristic volume
nconc_c = 1 / vol_c                              # Characteristic number density
conc_c = 1 / vol_c                               # Characteristic concentration

# =============================
# Derived variables: bulk concentrations
# =============================
# Primary salt
cbulk_plus_d, cbulk_neg_d = (cb1_d, cb1_d) if abs(valency[0]) == abs(valency[1]) else (cb1_d * abs(valency[1]), cb1_d * abs(valency[0]))
nbulk_plus_d, nbulk_neg_d = cbulk_plus_d*N_A*1e3, cbulk_neg_d*N_A*1e3
n_bulk_d = [nbulk_plus_d, nbulk_neg_d]

# Secondary salt (if present)
if cb2_d != 0:
    cb2_plus_d, cb2_neg_d = (cb2_d, cb2_d) if abs(valency[2]) == abs(valency[3]) else (cb2_d * abs(valency[3]), cb2_d * abs(valency[2]))
    nb2_plus_d, nb2_neg_d = cb2_plus_d*N_A*1e3, cb2_neg_d*N_A*1e3
    n_bulk_d = [nbulk_plus_d, nbulk_neg_d, nb2_plus_d, nb2_neg_d]

# Ionic strength and Debye screening length
I = sum([(valency[i] ** 2) * n_bulk_d[i] / len(valency) for i in range(len(valency))])
lambda_d_d = np.sqrt(epsilon_s_d / (beta * ec**2 * I))  # Debye length in m

    
# =============================
# Non-dimensionalization of physical quantities
# =============================
epsilon_s = epsilon_s_d / epsilon_c
n_bulk = np.true_divide(n_bulk_d, nconc_c)
sigma_i1 = sigma_in1_d / sigma_c
sigma_f1 = sigma_f1_d / sigma_c
sigma_i2 = sigma_in2_d / sigma_c
sigma_f2 = sigma_f2_d / sigma_c
rad_ions = np.true_divide(rad_ions_d*1e-10, l_c)
rad_sol = rad_sol_d*1e-10 / l_c
cbulk_plus = cbulk_plus_d / conc_c
cbulk_neg = cbulk_neg_d / conc_c
vol_sol = vol_sol_d / vol_c
vol_ions = np.true_divide(vol_ions_d, vol_c)
domain = domain_d*1e-10 / l_c
domain_in = domain_in_d*1e-10 / l_c
lambda_d = lambda_d_d / l_c

# -------------------------------
# Print inputs if ran directly
# -------------------------------
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
    print(f'epsilonr_s_d = {epsilonr_s_d}')
