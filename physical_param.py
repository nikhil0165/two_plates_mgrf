from packages import *

## Global Input Variables, All quantities are in SI unit _d means dimensional

cb1_d = 0.1 # prinamry salt bulk concentration in M
cb2_d = 0.0 # secondary salt bulk concentration in M
valency1 = [2,-1] # valency of primary salt
valency2 = [1,-1] # valency of secondary salt
born_radius1 = 1.6# radius of cation in Angstroms
born_radius2 = 1.5 # radius of anion in Angstroms
rad_sol_d = max(born_radius1,born_radius2)

domain_d = 20.0 # separation between two plates in Angstroms
domain_in_d = domain_d # separation between two plates of the initial guess

sigma_f1_d = -0.3204 # surface charge density of plate 1
sigma_f2_d = 0.1602 # surface charge density of plate 2
sigma_in1_d = -0.3204 #initial point for starting calculation in case of high surface charge densities
sigma_in2_d = 0.1602 # initial point for starting calculation in case of high surface charge densities

print(f'cb1_d = {cb1_d}')
print(f'cb2_d = {cb2_d}')
print(f'sigma_in1_d = {sigma_in1_d}')
print(f'sigma_in2_d = {sigma_in2_d}')
print(f'sigma_f1_d = {sigma_f1_d}')
print(f'sigma_f2_d = {sigma_f2_d}')
print(f'domain_in_d = {domain_in_d}')
print(f'domain_d = {domain_d}')

vol_sol_d = 4/3*pi*pow(rad_sol_d*pow(10, -10),3)# volume of solvent molecule assuming its a sphere

if cb2_d == 0:
    valency = np.array(valency1)
    rad_ions_d = np.array([born_radius1, born_radius2])
    vol_ions_d = np.array([vol_sol_d,vol_sol_d])
else:
    valency = np.hstack((valency1,valency2))
    rad_ions_d = np.array([born_radius1, born_radius2,born_radius2,born_radius2])# rad of ions
    vol_ions_d = np.array([vol_sol_d,vol_sol_d,vol_sol_d,vol_sol_d])


## Physical constants

Temp = 298# Temperature in Kelvin
ec = 1.602 * pow(10, -19)  # electronic charge
k_b = 1.38064852 * pow(10, -23)  # boltzmann's constant
N_A = 6.02214*pow(10,23) # Avogadro number
epsilon_o_d = 8.854187 * pow(10, -12)  # permittivity in vaccuum
beta = 1 / (k_b * Temp)  # beta
epsilonr_s_d = 78.7  # relative permittivity/dielectric constant
epsilon_s_d = epsilon_o_d * epsilonr_s_d  # permittivity of the medium



## Characteristic variables - dont play with this

l_b = pow(ec,2)*(beta)*(1 / (4 * pi * epsilon_s_d)) #Bjerrum Length
l_c = l_b # characteristic length scale in the system
q_c = ec # characteristic charge
psi_c = (1/(beta*q_c))# characteristic electrostatic potential
epsilon_c = beta*q_c*q_c/l_c # characteristic dielectric permimitivity
sigma_c = (ec/pow(l_c,2)) # characteristic surface charge density
vol_c = pow(l_c,3) # characteristic volume
nconc_c = 1/vol_c # characteristic number density
conc_c = 1/vol_c # characteristic concentration


## Derived variables

# bulk concentration of primary salt ions
cbulk_plus_d, cbulk_neg_d = (cb1_d, cb1_d) if abs(valency[0]) == abs(valency[1]) else (cb1_d * abs(valency[1]), cb1_d * abs(valency[0]))
nbulk_plus_d, nbulk_neg_d = cbulk_plus_d*N_A* pow(10, 3), cbulk_neg_d*N_A* pow(10, 3) # number of concentration of ions
n_bulk_d =  [nbulk_plus_d,nbulk_neg_d]

# bulk concentration of secondary salt ions
if cb2_d != 0:
    cb2_plus_d, cb2_neg_d = (cb2_d, cb2_d) if abs(valency[2]) == abs(valency[3]) else (cb2_d * abs(valency[3]), cb2_d * abs(valency[2]))
    nb2_plus_d, nb2_neg_d = cb2_plus_d * N_A* pow(10, 3), cb2_neg_d * N_A* pow(10, 3)  # number of concentration of ions
    n_bulk_d = [nbulk_plus_d,nbulk_neg_d,nb2_plus_d,nb2_neg_d]

I = sum([(valency[i] ** 2) * n_bulk_d[i] / len(valency) for i in range(len(valency))])
lambda_d_d = np.sqrt(epsilon_s_d / (beta * pow(ec, 2) * I)) # Debye Screening length

## Scaling the variables with characteristic variables

epsilon_s = epsilon_s_d / epsilon_c
n_bulk = np.true_divide(n_bulk_d,nconc_c)
sigma_i1 = sigma_in1_d/sigma_c
sigma_f1 = sigma_f1_d/sigma_c
sigma_i2 = sigma_in2_d/sigma_c
sigma_f2 = sigma_f2_d/sigma_c
rad_ions = np.true_divide(rad_ions_d*pow(10, -10),l_c)
rad_sol  = rad_sol_d*pow(10, -10)/l_c
cbulk_plus = cbulk_plus_d/conc_c
cbulk_neg = cbulk_plus_d/conc_c
vol_sol = vol_sol_d/vol_c
vol_ions = np.true_divide(vol_ions_d,vol_c)
domain  = domain_d*pow(10,-10)/l_c
domain_in  = domain_in_d*pow(10,-10)/l_c
lambda_d = lambda_d_d/l_c 




