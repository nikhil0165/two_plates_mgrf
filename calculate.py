"""
This module provides utility functions for electrostatic and ionic calculations 
in the two-plate MGRF model.

Functions included:
- Ionic strength, charge, charge density, and Kronecker delta utilities
- Eta factors (local & profile) for steric effects
- Screening factor (kappa) calculations
- Grid interpolation for electrostatic fields
- Potential extension and residual calculations for validation
- Apparent charge distribution along the domain
"""

from packages import *
from numerical_param import*

def ionic_strength(n_position,valency):  
    # Computes ionic strength at a single spatial node
    # Formula: I = sum( (z_i^2 / N_species) * n_i )
    # where z_i = valency of ion, n_i = local ion concentration
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    return I

def charge(psi, valency, n_bulk): 
    # Returns charge density contribution from a given ion species
    # Using mean-field Poisson-Boltzmann relation:
    # q = z * n_bulk * exp(-z * psi)
    return valency * n_bulk * exp(-valency * psi)

def kron_delta(i, j): 
    # Implementation of the Kronecker delta δ_ij
    # Returns 1 if i==j, else 0
    return 1 if i == j else 0

def charge_density(n_profile,valency): 
    # Computes total charge density profile over the domain
    # q_profile = Σ_i (n_i * z_i)
    q_profile = np.dot(n_profile, valency)
    return q_profile

def eta_loc(n_position,vol_ions, vol_sol): 
    # Computes the local excluded volume (η factor) at a node
    # η = -(1/V_sol) * log(1 - Σ_i (n_i * v_i))
    vol_local = np.sum(vol_ions * n_position)
    return (-1 / vol_sol) * log(1 - vol_local)

def eta_profile(n_profile,vol_ions, vol_sol): 
    # Computes η factor for the entire concentration profile
    # Vectorized using np.apply_along_axis calling eta_loc
    eta_profile = np.apply_along_axis(eta_loc, 1, n_profile, vol_ions, vol_sol)
    return eta_profile

def kappa_loc(n_position,valency,epsilon):  
    # Computes local Debye screening factor κ at a node
    # κ = sqrt( I * N_species / ε )
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    kappa = sqrt(I * (len(valency) / epsilon))
    return kappa

def kappa_sqr(n_position,valency, epsilon):  
    # Returns κ^2 (square of screening factor) at a node
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    return (I * (len(valency) / epsilon))

def kappa_sqr_profile(n_profile,valency,epsilon): 
    # Computes κ^2 for the full concentration profile
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_profile, axis=1)
    return (I * (len(valency) / epsilon))

def kappa_profile(n_profile,valency,  epsilon): 
    # Computes κ (screening factor) for all grid points
    kappa = np.apply_along_axis(kappa_loc,1,n_profile,valency,epsilon)
    return kappa

def interpolator(psi_complete,nconc_complete,bounds,new_grid): 
    # Resamples potential (psi) and concentrations (nconc) onto a new grid
    # Uses Dedalus (d3) Chebyshev basis for interpolation

    grid_points = len(psi_complete)
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)  # serial/parallel aware distribution
    zbasis = d3.Chebyshev(coords['z'],size = grid_points,bounds = bounds)

    # Initialize fields
    n_ions = len(nconc_complete[0,:])
    nconc = np.zeros((new_grid,n_ions))
    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_complete
    psi.change_scales(new_grid/grid_points)  # rescale field to new grid size

    # First 2 ionic concentrations
    nconc0 = dist.Field(name = 'nconc0',bases = zbasis)
    nconc1 = dist.Field(name = 'nconc1',bases = zbasis)
    nconc0['g'] = nconc_complete[:,0]
    nconc1['g'] = nconc_complete[:,1]
    nconc0.change_scales(new_grid/grid_points)
    nconc1.change_scales(new_grid/grid_points)
    nconc[:,0] = nconc0['g']
    nconc[:,1] = nconc1['g']

    # For systems with 4 ionic species
    if n_ions==4:
        nconc2 = dist.Field(name = 'nconc2',bases = zbasis)
        nconc3 = dist.Field(name = 'nconc3',bases = zbasis)
        nconc2['g'] = nconc_complete[:,2]
        nconc3['g'] = nconc_complete[:,3]
        nconc2.change_scales(new_grid/grid_points)
        nconc3.change_scales(new_grid/grid_points)
        nconc[:,2] = nconc2['g']
        nconc[:,3] = nconc3['g']

    return psi['g'], nconc

def psi_extender(psi_profile,r_sol, z_lg):
    # Extends potential profile ψ(z) beyond the simulation domain
    # by linearly extrapolating slopes at boundaries
    slope1 = (psi_profile[1]-psi_profile[0])/(z_lg[1] - z_lg[0])
    slope2 = (psi_profile[-1]-psi_profile[-2])/(z_lg[-1] - z_lg[-2])
    z_ext1 = np.linspace(0,r_sol,10, endpoint=False)  # left extension
    z_ext2 = np.linspace(z_lg[-1] + r_sol,z_lg[-1] + 2*r_sol,10)[1:]  # right extension
    psi_extend1 = slope1*z_ext1 + psi_profile[0] - slope1*r_sol
    psi_extend2 = psi_profile[-1] +slope2*(z_ext2-(z_lg[-1] + r_sol))
    return np.hstack((z_ext1,(z_lg + r_sol),z_ext2)), np.hstack((psi_extend1,psi_profile,psi_extend2))

def res_2plate(psi_profile,q_profile,bounds,sigma1,sigma2,epsilon): 
    # Computes residual of Gauss' law for two charged plates
    # Boundary conditions: dψ/dz(0) = -σ1/ε, dψ/dz(L) = σ2/ε

    nodes = len(psi_profile)
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)  
    zbasis = d3.Chebyshev(coords['z'],size = nodes,bounds = bounds,dealias = dealias)

    # Define fields
    z = dist.local_grids(zbasis)
    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_profile

    grad_psi = d3.Differentiate(psi,coords['z'])
    lap_psi = d3.Laplacian(psi).evaluate()
    lap_psi.change_scales(1)

    # Boundary slopes
    slope_0 = grad_psi(z = 0).evaluate()['g'][0]
    slope_end = grad_psi(z = bounds[1]).evaluate()['g'][0]
    res = np.zeros(nodes)
    
    # Residual at boundaries
    res[0] = slope_0 + sigma1/epsilon
    res[nodes-1] = slope_end -sigma2/epsilon

    # Residual inside domain
    res[1:nodes-1] = lap_psi['g'][1:nodes-1] + q_profile[1:nodes-1]/epsilon

    # Debug print: index and value of maximum residual
    return np.max(np.abs(res))

def app_charge(sigma,q_profile,z): 
    # Computes apparent charge distribution starting from plate 1
    # Uses trapezoidal rule cumulative integration of q(z)
    cumulative_sum = np.cumsum(0.5 * (q_profile[:-1] + q_profile[1:]) * np.diff(z))
    ans = np.concatenate(([0],cumulative_sum))
    return ans + sigma