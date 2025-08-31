# This file implements the Poisson-Boltzmann (PB) equation for two plates.
# It calculates the electrostatic potential and ion concentration profiles iteratively.

from packages import *
from numerical_param import *
import num_concn
import calculate

def pb_2plate(psi_guess, n_bulk, valency, sigma_1, sigma_2, domain, epsilon):
    """
    Solve the mean-field Poisson-Boltzmann equation between two charged plates.

    Parameters:
        psi_guess : ndarray
            Initial guess for the electrostatic potential (from Debye-Hückel solution)
        n_bulk : list or ndarray
            Bulk concentrations of ions
        valency : list or ndarray
            Valencies of ions
        sigma_1, sigma_2 : float
            Surface charge densities on the two plates
        domain : float
            Distance between plates
        epsilon : float
            Permittivity of the medium

    Returns:
        psi_profile : ndarray
            Electrostatic potential profile
        n_profile : ndarray
            Ion concentration profiles
        z : ndarray
            Spatial grid along z
        surface_psi : list
            Potential at plate surfaces [z=0, z=Lz]
    """

    grid_points = len(psi_guess)
    bounds = (0, domain)
    Lz = bounds[1]

    # Coefficients for RHS
    coeffs = [n_bulk[i] * valency[i] / epsilon for i in range(len(valency))]
    slope1 = -sigma_1 / epsilon
    slope2 = -sigma_2 / epsilon

    # Dedalus setup
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

    # Fields
    z = dist.local_grids(zbasis)
    psi = dist.Field(name='psi', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')  # auxiliary field for BCs
    tau_2 = dist.Field(name='tau_2')

    # Operators
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Nonlinear RHS (Boltzmann term)
    boltz = lambda psi: sum(coeffs[i] * np.exp(-valency[i] * psi) for i in range(len(valency)))

    # PDE setup
    problem = d3.NLBVP([psi, tau_1, tau_2], namespace=locals())
    problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = boltz(psi)")

    # Boundary conditions
    problem.add_equation("dz(psi)(z=0) = slope1")
    problem.add_equation("dz(psi)(z=Lz) = -slope2")

    # Initial guess
    psi['g'] = psi_guess

    # Solver
    solver = problem.build_solver(ncc_cutoff=ncc_cutoff_pb)
    pert_norm = np.inf
    psi.change_scales(dealias)
    while pert_norm > tolerance_pb:
        solver.newton_iteration()
        pert_norm = sum(pert.allreduce_data_norm('c', 2) for pert in solver.perturbations)
        print(f'mean-field PB convergence = {pert_norm:.3e}')

    psi.change_scales(1)
    psi_profile = psi.allgather_data('g')

    # Ion concentrations and charge density
    n_profile = num_concn.nconc_pb(psi_profile, valency, n_bulk)
    q_profile = calculate.charge_density(n_profile, valency)

    # Surface potentials
    surface1_psi = psi(z=0).evaluate()['g'][0]
    surface2_psi = psi(z=Lz).evaluate()['g'][0]

    # Gauss's law residual check
    res = calculate.res_2plate(psi_profile, q_profile, bounds, sigma_1, sigma_2, epsilon)
    print("Gauss's law residual for mean-field PB is =", res)

    return psi_profile, n_profile, np.squeeze(z), [float(surface1_psi), float(surface2_psi)]
