from packages import *               # Import all standard packages, Dedalus, numpy, etc.
from numerical_param import *         # Import numerical parameters like tolerance, dealias, etc.
import num_concn                      # Module for computing concentration profiles
import calculate                      # Module for auxiliary calculations (charge density, residuals)

# -------------------------------
# Function to solve Poisson-Boltzmann equation for a two-plate system
# -------------------------------
def pb_2plate(psi_guess,n_bulk,valency, sigma_1,sigma_2, domain, epsilon):  
    # psi_guess from Debye-Hückel acts as initial guess

    grid_points = len(psi_guess)       # Number of grid points for Chebyshev spectral method
    bounds = (0,domain)                # Domain boundaries in z-direction
    Lz = bounds[1]                     

    # Coefficients for exponential terms in RHS (Boltzmann factors)
    coeffs = [n_bulk[i] * valency[i]/epsilon for i in range(len(valency))]

    # Electric field slopes at plates
    slope1 = -sigma_1/epsilon
    slope2 = -sigma_2/epsilon

    # -------------------------------
    # Spectral bases
    # -------------------------------
    coords = d3.CartesianCoordinates('z')                  # 1D Cartesian coordinate system
    dist = d3.Distributor(coords,dtype=np.float64)        # Distributor handles parallelization (or serial)
    zbasis = d3.Chebyshev(coords['z'], size= grid_points, bounds=bounds, dealias=dealias)  
                                                           # Chebyshev spectral basis for z-direction

    # -------------------------------
    # Fields
    # -------------------------------
    z = dist.local_grids(zbasis)                          # Grid points in z
    psi = dist.Field(name='psi', bases=zbasis)           # Electrostatic potential field
    tau_1 = dist.Field(name='tau_1')                     # Auxiliary field for boundary lifting
    tau_2 = dist.Field(name='tau_2')                     # Auxiliary field for boundary lifting

    # -------------------------------
    # Substitutions / helper functions
    # -------------------------------
    dz = lambda A: d3.Differentiate(A, coords['z'])      # First derivative operator
    lift_basis = zbasis.derivative_basis(2)              # Basis for lifting second derivative
    lift = lambda A,n : d3.Lift(A, lift_basis, n)       # Lifting operator for boundary conditions

    # Lambda function for RHS (Boltzmann factor), Dedalus can differentiate this automatically
    boltz = lambda psi: sum(coeffs[i] * np.exp(-valency[i] * psi) for i in range(len(valency)))

    # -------------------------------
    # PDE setup
    # -------------------------------
    problem = d3.NLBVP([psi,tau_1, tau_2], namespace=locals())  # Nonlinear BVP
    problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = boltz(psi)")  # Poisson-Boltzmann eqn

    # -------------------------------
    # Boundary Conditions
    # -------------------------------
    problem.add_equation("dz(psi)(z=0) = slope1")  # Electric field at z=0
    problem.add_equation("dz(psi)(z=Lz) = -slope2") # Electric field at z=Lz

    # -------------------------------
    # Initial guess
    # -------------------------------
    psi['g'] = psi_guess

    # -------------------------------
    # Solver
    # -------------------------------
    solver = problem.build_solver(ncc_cutoff=ncc_cutoff_pb)   # Build Dedalus solver
    pert_norm = np.inf
    psi.change_scales(dealias)                                # Apply dealiasing
    while pert_norm > tolerance_pb:                           # Newton iteration until convergence
        solver.newton_iteration()
        pert_norm = sum(pert.allreduce_data_norm('c', 2) for pert in solver.perturbations)
        print(f'mean-field PB convergence = {pert_norm:.3e}')

    psi.change_scales(1)                                     # Return to full resolution

    # -------------------------------
    # Post-processing
    # -------------------------------
    psi_profile = psi.allgather_data('g')                     # Gather solution across processors
    n_profile = num_concn.nconc_pb(psi_profile,valency,n_bulk)  # Compute ion concentrations
    q_profile = calculate.charge_density(n_profile, valency)     # Compute charge density

    res= calculate.res_2plate(psi_profile,q_profile,bounds,sigma_1,sigma_2,epsilon)  
    print("Gauss's law residual for mean-field PB is = " + str(res))

    return psi_profile, n_profile,np.squeeze(z)              # Return potential, concentrations, and grid
