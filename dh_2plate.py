# This file implements the Debye-Hückel (DH) theory for two parallel plates.
# It calculates the electrostatic potential (psi) and ion concentration profiles (n_profile)
# between two charged plates using linearized Poisson-Boltzmann theory.

from packages import *
import calculate
import num_concn

def dh_2plate(n_bulk,valency,sigma_1,sigma_2,grid_points,domain,epsilon):
    """
    Solve the Debye-Hückel equation for two charged plates.
    
    Parameters:
    n_bulk     : bulk ion concentrations
    valency    : ion valencies
    sigma_1    : surface charge density of plate 1
    sigma_2    : surface charge density of plate 2
    grid_points: number of discretization points along z
    domain     : plate separation distance
    epsilon    : dielectric permittivity

    Returns:
    psi_profile   : electrostatic potential profile (array)
    n_profile     : ion concentration profile (matrix)
    z             : spatial coordinates (array)
    surface_psi   : list of surface potentials [plate1, plate2]
    """
    bounds = (0,domain)
    Lz = bounds[1]
    slope1 = -sigma_1/epsilon  # boundary gradient at plate 1
    slope2 = -sigma_2/epsilon  # boundary gradient at plate 2

    # Compute squared screening factor for uniform bulk
    kappa_2 = calculate.kappa_sqr(n_bulk,valency,epsilon)

    # Spectral basis setup using Dedalus
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64) # No mesh for serial / automatic parallelization
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds)

    # Define fields for potential and auxiliary tau variables
    z = dist.local_grids(zbasis)
    psi = dist.Field(name='psi', bases=zbasis)   # electrostatic potential
    tau_1 = dist.Field(name='tau_1')            # tau for boundary lift at z=0
    tau_2 = dist.Field(name='tau_2')            # tau for boundary lift at z=Lz

    # Substitutions for derivatives and lift operators
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Set up linear boundary value problem (LBVP)
    problem = d3.LBVP([psi,tau_1, tau_2], namespace=locals())
    problem.add_equation("-lap(psi) + kappa_2*psi + lift(tau_1,-1) + lift(tau_2,-2) = 0") # DH PDE

    # Boundary conditions: slope of psi at the plates
    problem.add_equation("dz(psi)(z=0) = slope1")
    problem.add_equation("dz(psi)(z=Lz) = -slope2")

    # Build solver and solve PDE
    solver = problem.build_solver()
    solver.solve()

    # Gather the solution for psi and compute ion concentrations using PB relation
    psi_profile = psi.allgather_data('g')
    n_profile = num_concn.nconc_pb(psi_profile,valency,n_bulk)

    # Evaluate surface potentials at the plates
    surface1_psi = psi(z = 0).evaluate()['g'][0]
    surface2_psi = psi(z = bounds[1]).evaluate()['g'][0]

    # Return profiles, coordinates, and surface potentials
    return psi_profile,n_profile,np.squeeze(z), [float(surface1_psi), float(surface2_psi)]
