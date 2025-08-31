from packages import *
import calculate
import num_concn

def dh_2plate(n_bulk, valency, sigma_1, sigma_2, grid_points, domain, epsilon):
    """
    Solves the linearized Debye–Hückel (DH) equation for electrostatic potential 
    between two parallel charged plates using Dedalus.

    Parameters
    ----------
    n_bulk : array
        Bulk ion concentrations (mol/L or number density depending on units).
    valency : array
        Valence of each ionic species.
    sigma_1 : float
        Surface charge density on plate 1 (at z=0).
    sigma_2 : float
        Surface charge density on plate 2 (at z=Lz).
    grid_points : int
        Number of Chebyshev collocation points in the domain.
    domain : float
        Separation between the plates (length of z-domain).
    epsilon : float
        Dielectric permittivity of the medium.

    Returns
    -------
    psi_profile : ndarray
        Electrostatic potential profile between the plates.
    n_profile : ndarray
        Ion concentration profile (from PB equilibrium).
    z : ndarray
        Grid points (z positions) where psi and n are evaluated.
    """

    # Domain bounds (0 to plate separation)
    bounds = (0, domain)
    Lz = bounds[1]

    # Slopes of potential at boundaries (from Gauss’ law at charged surfaces)
    slope1 = -sigma_1 / epsilon   # dψ/dz at z=0
    slope2 = -sigma_2 / epsilon   # dψ/dz at z=Lz

    # Screening parameter κ² from bulk concentrations and valencies
    kappa_2 = calculate.kappa_sqr(n_bulk, valency, epsilon)

    # Dedalus setup: 1D Chebyshev basis in z-direction
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64)  # serial/parallel distributor
    zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds)

    # Define fields
    z = dist.local_grids(zbasis)       # spatial grid (collocation points)
    psi = dist.Field(name='psi', bases=zbasis)   # electrostatic potential
    tau_1 = dist.Field(name='tau_1')   # τ fields = auxiliary variables for boundary conditions
    tau_2 = dist.Field(name='tau_2')

    # Substitution rules for differentiation and lifting (used by Dedalus)
    dz = lambda A: d3.Differentiate(A, coords['z'])  # ∂/∂z
    lift_basis = zbasis.derivative_basis(2)          # 2nd derivative basis for τ terms
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # PDE: Linearized Poisson–Boltzmann (Debye–Hückel) equation
    # -∇²ψ + κ² ψ = 0
    # with τ terms included for boundary condition enforcement
    problem = d3.LBVP([psi, tau_1, tau_2], namespace=locals())
    problem.add_equation("-lap(psi) + kappa_2*psi + lift(tau_1,-1) + lift(tau_2,-2) = 0")

    # Boundary conditions: dψ/dz fixed at each plate
    problem.add_equation("dz(psi)(z=0) = slope1")      # Plate 1
    problem.add_equation("dz(psi)(z=Lz) = -slope2")    # Plate 2

    # Build and solve the linear boundary value problem
    solver = problem.build_solver()
    solver.solve()

    # Gather potential profile across all processes
    psi_profile = psi.allgather_data('g')

    # Compute corresponding ion concentration profile via PB equilibrium
    n_profile = num_concn.nconc_pb(psi_profile, valency, n_bulk)

    # Return potential, concentrations, and spatial grid
    return psi_profile, n_profile, np.squeeze(z)
