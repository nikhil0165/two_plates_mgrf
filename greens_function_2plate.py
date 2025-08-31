"""
Solves for the Green's function in the two-plate geometry.
Includes both the free and full (with profile) cases using Chebyshev spectral methods.
"""
from packages import *
import calculate
from numerical_param import*


def Gcap_free(grid_points, s, domain, epsilon):  # function for \hat{Go}
    """
    Compute the free-space Green's function (no profile) in a two-plate geometry
    using Chebyshev spectral methods.
    """

    # Define domain bounds and length
    bounds = (0, domain)
    Lz = bounds[1]

    # Spectral bases and distributor for Dedalus
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)  
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

    # Grid and differentiation tools
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])   # first derivative operator
    lift_basis = zbasis.derivative_basis(2)           # lifting basis for tau terms
    lift = lambda A, n: d3.Lift(A, lift_basis, n)     # lifting operator

    # Fields for Pz (auxiliary field from Riccati transformation)
    Pz = dist.Field(name='Pz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')

    # Nonlinear BVP for Pz (comes from Riccati form of SL problem)
    problem = d3.NLBVP([Pz, tau_1], namespace=locals())
    problem.add_equation("-dz(Pz) + s*s + lift(tau_1,-1) = Pz**2")

    # Boundary condition at z=0
    problem.add_equation("Pz(z=0) = s")

    # Initial guess for Newton solver
    Pz['g'] = s

    # Build solver and iterate until convergence
    solver0 = problem.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm0 = np.inf
    Pz.change_scales(dealias)
    while pert_norm0 > tolerance_greens:
        solver0.newton_iteration()
        pert_norm0 = sum(pert0.allreduce_data_norm('c', 2) for pert0 in solver0.perturbations)

    # Extract solution
    Pz.change_scales(1)
    Pz = Pz['g']
    Qz = -Pz  # symmetry relation in free case

    # Sturm–Liouville representation of Green’s function
    G = (-1 / epsilon) * np.true_divide(1, Qz - Pz)

    # Clean up to free memory
    del z, Pz, Qz, tau_1, dz, lift_basis, lift, problem, solver0, pert_norm0
    gc.collect()

    return G


def Gcap_full(n_profile, n_bulk, valency, s, domain, epsilon):  # function for \hat{G}
    """
    Compute the full Green's function with a given charge density profile.
    Uses Chebyshev spectral methods and Riccati-type reformulation.
    """

    # Domain and grid setup
    grid_points = len(n_profile)
    bounds = (0, domain)
    Lz = bounds[1]

    # Spectral bases and distributor
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64) 
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

    # Grid and differentiation tools
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])   # derivative
    lift_basis = zbasis.derivative_basis(2)           # lifting basis
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Compute squared frequency profile (includes screening)
    omega_sqr = dist.Field(bases=zbasis)
    omega_sqr['g'] = s * s + calculate.kappa_sqr_profile(n_profile, valency, epsilon)
    omega_b = np.sqrt(s * s + calculate.kappa_sqr(n_bulk, valency, epsilon))  # bulk value

    # -------------------------
    # Solve for Pz 
    # -------------------------
    Pz = dist.Field(name='Pz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')

    problem = d3.NLBVP([Pz, tau_1], namespace=locals())
    problem.add_equation("-dz(Pz) + omega_sqr + lift(tau_1,-1) = Pz**2")

    # Boundary condition at z=0
    problem.add_equation("Pz(z=0) = s")

    # Initial guess
    Pz['g'] = s

    # Newton solver loop
    solver1 = problem.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm1 = np.inf
    Pz.change_scales(dealias)
    p = 0
    while pert_norm1 > tolerance_greens:
        p = p + 1
        solver1.newton_iteration()
        pert_norm1 = sum(pert1.allreduce_data_norm('c', 2) for pert1 in solver1.perturbations)

    # Extract converged solution for Pz
    Pz.change_scales(1)
    Pz = Pz['g']

    # -------------------------
    # Solve for Qz 
    # -------------------------
    Qz = dist.Field(name='Qz', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')

    problem1 = d3.NLBVP([Qz, tau_1], namespace=locals())
    problem1.add_equation("-dz(Qz) + omega_sqr + lift(tau_1,-1) = Qz**2")

    # Boundary condition at z=Lz
    problem1.add_equation("Qz(z=Lz) = -s")

    # Initial guess
    Qz['g'] = -s

    # Newton solver loop
    solver2 = problem1.build_solver(ncc_cutoff=ncc_cutoff_greens)
    pert_norm2 = np.inf
    Qz.change_scales(dealias)
    q = 1
    while pert_norm2 > tolerance_greens:
        q = q + 1
        solver2.newton_iteration()
        pert_norm2 = sum(pert2.allreduce_data_norm('c', 2) for pert2 in solver2.perturbations)

    # Extract converged solution for Qz
    Qz.change_scales(1)
    Qz = Qz['g']

    # -------------------------
    # Compute Green’s function
    # -------------------------
    G = (-1 / epsilon) * np.true_divide(1, Qz - Pz)

    # Clean up to free memory
    del z, Pz, Qz, tau_1, dz, lift_basis, lift, problem, solver1, solver2, pert_norm2, pert_norm1
    gc.collect()

    return G
