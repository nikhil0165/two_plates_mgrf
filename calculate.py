# This file contains functions for various calculations related to ionic strength, charge density, eta factor, screening factors,
# and profile extensions for electrostatic calculations in multi-ion systems.

from packages import *
from numerical_param import*

def ionic_strength(n_position,valency):  
    """
    Compute ionic strength at a given position/node.
    n_position : number of ions at this position
    valency : list/array of ion valencies
    """
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    return I

def charge(psi, valency, n_bulk): 
    """
    Compute the local charge density at a given potential using mean-field PB.
    psi : electrostatic potential
    valency : ion valencies
    n_bulk : bulk ion concentrations
    """
    return valency * n_bulk * exp(-valency * psi)

def kron_delta(i, j): 
    """
    Kronecker delta function
    Returns 1 if i == j, else 0
    """
    return 1 if i == j else 0

def charge_density(n_profile,valency): 
    """
    Compute total charge density across the spatial domain.
    n_profile : ion concentration profile
    valency : ion valencies
    """
    q_profile = np.dot(n_profile, valency)
    return q_profile

def eta_loc(n_position, vol_ions, vol_sol): 
    """
    Compute local packing fraction / excluded volume factor (eta) at a single position.
    vol_ions : array of ion volumes
    vol_sol : solvent molecular volume
    """
    vol_local = np.sum(vol_ions * n_position)
    return (-1 / vol_sol) * log(1 - vol_local)

def eta_profile(n_profile, vol_ions, vol_sol): 
    """
    Compute eta factor across the full spatial profile
    """
    eta_profile = np.apply_along_axis(eta_loc, 1, n_profile, vol_ions, vol_sol)
    return eta_profile

def kappa_loc(n_position,valency,epsilon):  
    """
    Compute local screening factor (kappa) at a single position
    epsilon : permittivity
    """
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    kappa = sqrt(I * (len(valency) / epsilon))
    return kappa

def kappa_sqr(n_position,valency, epsilon):  
    """
    Compute square of local screening factor
    """
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_position)
    return (I * (len(valency) / epsilon))

def kappa_sqr_profile(n_profile,valency,epsilon): 
    """
    Compute square of screening factor for all positions
    """
    I = np.sum((1 / len(valency)) * (np.power(valency,2)) * n_profile, axis=1)
    return (I * (len(valency) / epsilon))

def kappa_profile(n_profile,valency,  epsilon): 
    """
    Compute screening factor for all positions using vectorized local evaluation
    """
    kappa = np.apply_along_axis(kappa_loc,1,n_profile,valency,epsilon)
    return kappa

def profile_extender(psi_profile,n_profile,uself_profile, bounds,dist_exc,N_exc):
    """
    Extend profiles beyond the original domain for far-field computations.
    - Adds extra points on both sides using linear extrapolation.
    - Returns extended psi, n_profile, uself_profile, extended coordinates, and surface potentials.
    """
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)  # No mesh for serial/automatic parallelization
    zbasis = d3.Chebyshev(coords['z'],size = len(psi_profile),bounds = bounds,dealias = dealias)

    # Fields
    z = np.squeeze(dist.local_grids(zbasis))
    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_profile
    surface1_psi = psi(z = 0).evaluate()['g'][0]  # potential at first plate
    surface2_psi = psi(z = bounds[1]).evaluate()['g'][0]  # potential at second plate

    grad_psi = d3.Differentiate(psi,coords['z'])
    slope1 = grad_psi(z = 0).evaluate()['g'][0]  # gradient at first plate
    slope2 = grad_psi(z = bounds[1]).evaluate()['g'][0]  # gradient at second plate

    # Extend domain linearly outside original boundaries
    z_ext1 = np.linspace(0,dist_exc,N_exc,endpoint=False)
    z_ext2 = np.flip(np.linspace(z[-1] + 2*dist_exc,z[-1] + dist_exc,N_exc,endpoint = False))
    psi_extend1 = slope1 * z_ext1 + surface1_psi - slope1 * dist_exc
    psi_extend2 = surface2_psi + slope2*(z_ext2-(z[-1] + dist_exc))

    # Extend concentration and self-energy profiles
    n_profile = np.concatenate((np.zeros((N_exc,len(n_profile[0,:]))), n_profile), axis=0)
    n_profile = np.concatenate((n_profile,np.zeros((N_exc,len(n_profile[0,:])))), axis=0)

    uself_profile = np.concatenate((np.zeros((N_exc,len(n_profile[0,:]))),uself_profile),axis = 0)
    uself_profile = np.concatenate((uself_profile,np.zeros((N_exc,len(n_profile[0,:])))),axis = 0)
    
    # Return extended profiles and surface potentials as float
    return np.hstack((psi_extend1,psi_profile,psi_extend2)), n_profile,uself_profile,np.hstack((z_ext1,z+dist_exc,z_ext2)), [float(surface1_psi), float(surface2_psi)]

def interpolator(psi_profile,domain,points):
    """
    Interpolate psi_profile at arbitrary points using spectral Chebyshev interpolation.
    """
    grid_points = len(psi_profile)
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'],size = grid_points,bounds = (0,domain))

    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_profile

    psi_answer = np.zeros(len(points))
    for i in range(0,len(points)):
        psi_answer[i] = psi(z = points[i]).evaluate()['g'][0]

    return psi_answer

def rescaler(psi_profile,n_profile,bounds,new_grid): 
    """
    Rescale psi and n_profile onto a new grid with new_grid points.
    Works for 2 or 4 ion species.
    """
    grid_points = len(psi_profile)
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)                                                                                                      
    zbasis = d3.Chebyshev(coords['z'],size = grid_points,bounds = bounds)

    n_ions = len(n_profile[0,:])
    nconc = np.zeros((new_grid,n_ions))
    
    # Rescale psi
    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_profile
    psi.change_scales(new_grid/grid_points)

    # Rescale ion concentration profiles
    nconc0 = dist.Field(name = 'nconc0',bases = zbasis)
    nconc1 = dist.Field(name = 'nconc1',bases = zbasis)
    nconc0['g'] = n_profile[:,0]
    nconc1['g'] = n_profile[:,1]
    nconc0.change_scales(new_grid/grid_points)
    nconc1.change_scales(new_grid/grid_points)
    nconc[:,0] = nconc0.allgather_data('g')
    nconc[:,1] = nconc1.allgather_data('g')

    if n_ions==4:
        # For systems with 4 ion species
        nconc2 = dist.Field(name = 'nconc2',bases = zbasis)
        nconc3 = dist.Field(name = 'nconc3',bases = zbasis)
        nconc2['g'] = n_profile[:,2]
        nconc3['g'] = n_profile[:,3]
        nconc2.change_scales(new_grid/grid_points)
        nconc3.change_scales(new_grid/grid_points)
        nconc[:,2] = nconc2.allgather_data('g')
        nconc[:,3] = nconc3.allgather_data('g')

    return psi.allgather_data('g'), nconc,0

def res_2plate(psi_profile,q_profile,bounds,sigma1,sigma2,epsilon): 
    """
    Compute residual of Gauss's law for two-plate system.
    Returns maximum absolute residual over domain.
    """
    nodes = len(psi_profile)
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)  
    zbasis = d3.Chebyshev(coords['z'],size = nodes,bounds = bounds,dealias = dealias)

    z = dist.local_grids(zbasis)
    psi = dist.Field(name = 'psi',bases = zbasis)
    psi['g'] = psi_profile

    grad_psi = d3.Differentiate(psi,coords['z'])
    lap_psi = d3.Laplacian(psi).evaluate()
    lap_psi.change_scales(1)

    slope1 = grad_psi(z = 0).evaluate()['g'][0]
    slope2 = grad_psi(z = bounds[1]).evaluate()['g'][0]

    res = np.zeros(nodes)
    res[0] = slope1 + sigma1/epsilon
    res[nodes-1] = slope2 - sigma2/epsilon
    res[1:nodes-1] = lap_psi.allgather_data('g')[1:nodes-1] + q_profile[1:nodes-1]/epsilon

    return np.max(np.abs(res))

def app_charge(sigma,q_profile,z): 
    """
    Compute apparent charge density at plate 1.
    Integrates charge density profile cumulatively and adds plate surface charge.
    """
    cumulative_sum = np.cumsum(0.5 * (q_profile[:-1] + q_profile[1:]) * np.diff(z))
    ans = np.concatenate(([0],cumulative_sum))
    return ans + sigma
