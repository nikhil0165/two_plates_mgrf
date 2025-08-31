# This file computes the Green's function for bulk (homogeneous) ionic systems.
# Includes functions for free and full Green's function calculations in bulk.

from packages import *
import calculate
from numerical_param import *

def Gcap_free(grid_points, s, domain, epsilon):
    """
    Compute the free Green's function \hat{G}_0 for a bulk system.

    Parameters:
    grid_points : number of spatial grid points
    s           : wavenumber in Fourier space
    domain      : length of the domain
    epsilon     : dielectric permittivity

    Returns:
    G           : free Green's function evaluated on the grid
    """

    bounds = (0, domain)
    Lz = bounds[1]

    # Spectral method: Chebyshev basis and grid setup
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

    # Local grids and derivative/lift functions
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Fields for P(z) (log-derivative of U)
    Pz = dist.Field(name='Pz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')  # auxiliary field to enforce BC

    # Nonlinear boundary value problem for Pz
    problem = d3.NLBVP([Pz, tau_1], namespace=locals())
    problem.add_equation("-dz(Pz) + s*s + lift(tau_1,-1) = Pz**2")  # Riccati eqn for Pz

    # Boundary condition at z=0
    problem.add_equation("Pz(z=0) = s")

    # Initial guess for solver
    Pz['g'] = s

    # Solver setup and Newton iteration until convergence
    solver0 = problem.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm0 = np.inf
    Pz.change_scales(dealias)
    while pert_norm0 > tolerance_greens:
        solver0.newton_iteration()
        pert_norm0 = sum(pert0.allreduce_data_norm('c', 2) for pert0 in solver0.perturbations)

    Pz.change_scales(1)
    Pz = Pz.allgather_data('g')

    # Q(z) field (log-derivative of V) in free system
    Qz = -Pz

    # Sturm-Liouville expression for Green's function
    G = (1 / epsilon) * np.true_divide(1, Pz - Qz)

    # Clean up large arrays
    del z, Pz, Qz, tau_1, dz, lift_basis, lift, problem, solver0, pert_norm0
    gc.collect()

    return G


def Gcap_full(n_bulk_profile, n_bulk, valency, s, domain, epsilon):
    """
    Compute the full Green's function \hat{G} for a bulk system
    that includes the effect of inhomogeneous ion densities.

    Parameters:
    n_bulk_profile : ion density profile (array)
    n_bulk         : bulk ion concentrations
    valency        : ion valencies
    s              : wavenumber
    domain         : domain length
    epsilon        : dielectric permittivity

    Returns:
    G              : full Green's function on the spatial grid
    """

    grid_points = len(n_bulk_profile)
    bounds = (0, domain)
    Lz = bounds[1]

    # Spectral method setup
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

    # Local grids and derivative/lift functions
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Omega squared field: includes wavenumber + local screening factor
    omega_sqr = dist.Field(bases=zbasis)
    omega_sqr['g'] = s*s + calculate.kappa_sqr_profile(n_bulk_profile, valency, epsilon)
    omega_b = np.sqrt(s*s + calculate.kappa_sqr(n_bulk, valency, epsilon))

    # Field for P(z) (log-derivative of U)
    Pz = dist.Field(name='Pz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')

    # Nonlinear BVP for Pz
    problem = d3.NLBVP([Pz, tau_1], namespace=locals())
    problem.add_equation("-dz(Pz) + omega_sqr + lift(tau_1,-1) = Pz**2")
    problem.add_equation("Pz(z=0) = omega_b")

    # Initial guess
    Pz['g'] = omega_b

    # Solver iteration
    solver1 = problem.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm1 = np.inf
    Pz.change_scales(dealias)
    p = 0
    while pert_norm1 > tolerance_greens:
        p += 1
        solver1.newton_iteration()
        pert_norm1 = sum(pert1.allreduce_data_norm('c', 2) for pert1 in solver1.perturbations)

    Pz.change_scales(1)
    Pz = Pz.allgather_data('g')

    # Field for Q(z) (log-derivative of V)
    Qz = dist.Field(name='Qz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')

    # Nonlinear BVP for Qz
    problem1 = d3.NLBVP([Qz, tau_1], namespace=locals())
    problem1.add_equation("-dz(Qz) + omega_sqr + lift(tau_1,-1) = Qz**2")
    problem1.add_equation("Qz(z=Lz) = -omega_b")  # BC at the other end

    # Initial guess for Qz
    Qz['g'] = -omega_b

    # Solver iteration
    solver2 = problem1.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm2 = np.inf
    Qz.change_scales(dealias)
    q = 1
    while pert_norm2 > tolerance_greens:
        q += 1
        solver2.newton_iteration()
        pert_norm2 = sum(pert2.allreduce_data_norm('c', 2) for pert2 in solver2.perturbations)

    Qz.change_scales(1)
    Qz = Qz.allgather_data('g')

    # Sturm-Liouville expression for full Green's function
    G = (1 / epsilon) * np.true_divide(1, Pz - Qz)

    # Clean up
    del z, Pz, Qz, tau_1, dz, lift_basis, lift, problem, solver1, solver2, pert_norm2, pert_norm1
    gc.collect()

    return G
