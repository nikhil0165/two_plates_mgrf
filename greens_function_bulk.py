from packages import *
import calculate
from numerical_param import*


def Gcap_free(grid_points,s,domain,epsilon):#\hat{Go}
    """
    Computes the free-space (bulk, no ion profile) Green's function \hat{Go} 
    using a nonlinear boundary value problem (NLBVP) in Dedalus. 
    This corresponds to solving the Riccati-transformed form of the Sturm–Liouville problem.
    """

    bounds = (0,domain)    # z-domain: [0, Lz]
    Lz = bounds[1]

    # Bases
    coords = d3.CartesianCoordinates('z')   # 1D coordinate in z
    dist = d3.Distributor(coords,dtype = np.float64)   # Distributor (for parallel / serial execution)
    zbasis = d3.Chebyshev(coords['z'],size = grid_points,bounds = bounds,dealias = dealias)

    # General fields
    z = dist.local_grids(zbasis)   # Physical grid in z
    dz = lambda A: d3.Differentiate(A,coords['z'])  # Derivative operator in z
    lift_basis = zbasis.derivative_basis(2)         # Auxiliary lifting basis
    lift = lambda A,n: d3.Lift(A,lift_basis,n)      # Lifting operator for τ fields

    # Fields for Riccati variable Pz (log-derivative of mode function U)
    Pz = dist.Field(name = 'Pz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')   # τ field to enforce boundary conditions

    # Differential equation for Pz (Riccati form):
    # -Pz' + s² + lift(τ₁) = Pz²
    problem = d3.NLBVP([Pz,tau_1],namespace = locals())
    problem.add_equation("-dz(Pz) + s*s + lift(tau_1,-1) = Pz**2")

    # Boundary condition: Pz(z=0) = s (bulk matching condition)
    problem.add_equation("Pz(z=0) = s")

    # Initial guess: start with constant Pz = s
    Pz['g'] = s

    # Solver setup
    solver0 = problem.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm0 = np.inf

    # Newton iteration until convergence
    Pz.change_scales(dealias)
    while pert_norm0 > tolerance_greens:
        solver0.newton_iteration()
        pert_norm0 = sum(pert0.allreduce_data_norm('c',2) for pert0 in solver0.perturbations)

    # Gather converged solution back to grid space
    Pz.change_scales(1)
    Pz = Pz['g']

    # Define Qz as -Pz (symmetry in Riccati formulation)
    Qz = -Pz

    # Compute free-space Green’s function via Riccati–Sturm–Liouville relation
    G = (-1 / epsilon) * np.true_divide(1,Qz - Pz)

    # Cleanup
    del z,Pz,Qz,tau_1,dz,lift_basis,lift,problem,solver0,pert_norm0
    gc.collect()

    return G


def Gcap_full(n_bulk_profile, n_bulk, valency, s, domain,epsilon): # \hat{G}
    """
    Computes the full inhomogeneous Green's function \hat{G} 
    in the presence of an ion density profile n_bulk_profile(z).
    Uses Riccati transformation for forward (Pz) and backward (Qz) solutions.
    """

    grid_points = len(n_bulk_profile)
    bounds = (0,domain)
    Lz = bounds[1]

    # Bases
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords,dtype = np.float64)
    zbasis = d3.Chebyshev(coords['z'],size = grid_points,bounds = bounds,dealias = dealias)

    # General fields
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A,coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A,n: d3.Lift(A,lift_basis,n)

    # Effective squared frequency ω²(z) = s² + κ²(z)
    omega_sqr = dist.Field(bases = zbasis)
    omega_sqr['g'] = s * s + calculate.kappa_sqr_profile(n_bulk_profile,valency,epsilon)

    # Bulk (homogeneous) ω value for boundary conditions
    omega_b = np.sqrt(s * s + calculate.kappa_sqr(n_bulk,valency,epsilon))

    # === Solve forward Riccati equation for Pz (log-derivative of U) ===
    Pz = dist.Field(name = 'Pz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')

    # Nonlinear Riccati equation:
    # -Pz' + ω²(z) + lift(τ₁) = Pz²
    problem = d3.NLBVP([Pz,tau_1],namespace = locals())
    problem.add_equation("-dz(Pz) + omega_sqr + lift(tau_1,-1) = Pz**2")

    # Boundary condition at z=0: Pz(0) = ω_b
    problem.add_equation("Pz(z=0) = omega_b")

    # Initial guess
    Pz['g'] = omega_b

    # Solve iteratively
    solver1 = problem.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm1 = np.inf
    Pz.change_scales(dealias)
    while pert_norm1 > tolerance_greens:
        solver1.newton_iteration()
        pert_norm1 = sum(pert1.allreduce_data_norm('c',2) for pert1 in solver1.perturbations)

    Pz.change_scales(1)
    Pz = Pz.allgather_data('g')[0]

    # === Solve backward Riccati equation for Qz (log-derivative of V) ===
    Qz = dist.Field(name = 'Qz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')

    # Riccati equation:
    # -Qz' + ω²(z) + lift(τ₁) = Qz²
    problem1 = d3.NLBVP([Qz,tau_1],namespace = locals())
    problem1.add_equation("-dz(Qz) + omega_sqr + lift(tau_1,-1) = Qz**2")

    # Boundary condition at z=Lz: Qz(Lz) = -ω_b
    problem1.add_equation("Qz(z=Lz) = -omega_b")

    # Initial guess
    Qz['g'] = -omega_b

    # Solve iteratively
    solver2 = problem1.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm2 = np.inf
    Qz.change_scales(dealias)
    while pert_norm2 > tolerance_greens:
        solver2.newton_iteration()
        pert_norm2 = sum(pert2.allreduce_data_norm('c',2) for pert2 in solver2.perturbations)

    Qz.change_scales(1)
    Qz = Qz.allgather_data('g')[0]

    # === Construct full Green’s function ===
    # From Riccati theory: G(z) = (-1/ε) * 1 / (Qz - Pz)
    G = (-1 / epsilon) * np.true_divide(1,Qz - Pz)

    # Cleanup unused objects
    del z,Pz,Qz,tau_1,dz,lift_basis,lift,problem,solver1,solver2,pert_norm2,pert_norm1
    gc.collect()

    return G
