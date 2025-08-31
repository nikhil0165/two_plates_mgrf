# This file implements the Modified Generalized Reaction Field (MGRF) method for two plates.
# It calculates the electrostatic potential and ion concentration profiles iteratively.

from packages import *
import num_concn
import calculate
import selfe_2plate
import selfe_bulk
from numerical_param import*

def mgrf_2plate(psi_guess, nconc_guess, n_bulk, valency, rad_ions, vol_ions, vol_sol,
                 sigma_1, sigma_2, domain, epsilon_s, epsilon_p):
    """
    Solve the electrostatics of two charged plates using the MGRF method.

    Parameters:
    psi_guess    : initial guess for electrostatic potential (from PB)
    nconc_guess  : initial guess for ion concentration profiles
    n_bulk       : bulk ion concentrations
    valency      : ion valencies
    rad_ions     : ion radii
    vol_ions     : ion volumes
    vol_sol      : solvent volume
    sigma_1,2    : plate surface charges
    domain       : system length
    epsilon_s,p  : dielectric constants (solvent, plate)

    Returns:
    psi_profile  : converged electrostatic potential profile
    n_profile    : converged ion concentration profile
    uself_profile: self-energy profile
    q_profile    : charge density profile
    Z            : spatial grid
    res          : Gauss law residual
    surface_psi  : surface potentials at plates
    """

    grid_points = len(psi_guess)
    bounds = (0, domain)
    Lz = bounds[1]
    slope1 = -sigma_1 / epsilon_s
    slope2 = -sigma_2 / epsilon_s

    # Initialize potential and self-energy profiles
    psi_profile = np.copy(psi_guess)
    uself_guess = selfe_2plate.uself_complete(nconc_guess, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p)
    eta_guess = calculate.eta_profile(nconc_guess, vol_ions, vol_sol)

    print('selfe done in the interface')

    # Bulk self-energy and eta
    n_bulk_numerical = np.multiply(np.ones((grid_points, len(valency))), n_bulk)
    uself_bulk = np.mean(selfe_bulk.uselfb_numerical(n_bulk_numerical, n_bulk, rad_ions, valency, domain, epsilon_s), axis=0)
    eta_bulk = calculate.eta_loc(n_bulk, vol_ions, vol_sol)

    print('selfe done in the bulk')

    equal_vols = np.all(np.abs(vol_ions - vol_sol) < vol_sol * 1e-5)
    print(f'equal_vols = {equal_vols}')

    n_profile = None

    # Iterative solver
    convergence_tot = np.inf
    iteration = 1
    while convergence_tot > tolerance:

        # Spectral method setup: Chebyshev basis and distributor
        coords = d3.CartesianCoordinates('z')
        dist = d3.Distributor(coords, dtype=np.float64)
        zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

        # Fields
        z = dist.local_grids(zbasis)
        psi = dist.Field(name='psi', bases=zbasis)
        tau_1 = dist.Field(name='tau_1')  # auxiliary field for BCs
        tau_2 = dist.Field(name='tau_2')  # auxiliary field for BCs

        # Substitutions
        dz = lambda A: d3.Differentiate(A, coords['z'])
        lift_basis = zbasis.derivative_basis(2)
        lift = lambda A, n: d3.Lift(A, lift_basis, n)

        # Compute coefficients for RHS of Poisson eqn
        c0 = dist.Field(bases=zbasis)
        c1 = dist.Field(bases=zbasis)
        _, coeffs = num_concn.nconc_mgrf(psi_profile, uself_guess, eta_guess,
                                         uself_bulk, n_bulk, valency, vol_ions,
                                         eta_bulk, equal_vols)
        coeffs = coeffs / epsilon_s

        # Assign coefficients to fields
        c0['g'] = np.squeeze(coeffs[:, 0])
        c1['g'] = np.squeeze(coeffs[:, 1])
        boltz0 = lambda psi: np.exp(-valency[0] * psi)
        boltz1 = lambda psi: np.exp(-valency[1] * psi)

        # For 4-component ions
        if len(valency) == 4:
            c2 = dist.Field(bases=zbasis)
            c3 = dist.Field(bases=zbasis)
            boltz2 = lambda psi: np.exp(-valency[2] * psi)
            boltz3 = lambda psi: np.exp(-valency[3] * psi)

        # PDE setup for MGRF
        problem = d3.NLBVP([psi, tau_1, tau_2], namespace=locals())
        if len(valency) == 2:
            problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = c0*boltz0(psi) + c1*boltz1(psi)")
        if len(valency) == 4:
            c2['g'] = np.squeeze(coeffs[:, 2])
            c3['g'] = np.squeeze(coeffs[:, 3])
            problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = "
                                 "c0*boltz0(psi) + c1*boltz1(psi) + c2*boltz2(psi) + c3*boltz3(psi)")

        # Boundary conditions at plates
        problem.add_equation("dz(psi)(z=0) = slope1")
        problem.add_equation("dz(psi)(z=Lz) = -slope2")

        # Initial guess for Newton iteration
        psi['g'] = psi_profile

        # Solver
        solver = problem.build_solver(ncc_cutoff=ncc_cutoff_mgrf)
        pert_norm = np.inf
        psi.change_scales(dealias)
        s = 0
        while pert_norm > tolerance_pb:
            solver.newton_iteration()
            pert_norm = sum(pert.allreduce_data_norm('c', 2) for pert in solver.perturbations)
            s += 1

        psi.change_scales(1)
        psi_profile = psi.allgather_data('g')

        if np.any(np.isnan(psi_profile)):
            print('nan in psi')

        # Update ion concentrations
        n_profile, _ = num_concn.nconc_mgrf(psi_profile, uself_guess, eta_guess,
                                            uself_bulk, n_bulk, valency, vol_ions,
                                            eta_bulk, equal_vols)

        # Check convergence
        convergence_tot = np.true_divide(np.linalg.norm(n_profile - nconc_guess), np.linalg.norm(nconc_guess))

        # Relaxation update
        nconc_guess = num_ratio * n_profile + (1 - num_ratio) * nconc_guess

        # Update self-energy and eta profiles
        uself_guess = selfe_2plate.uself_complete(nconc_guess, n_bulk, rad_ions, valency, domain, epsilon_s, epsilon_p)
        eta_guess = calculate.eta_profile(nconc_guess, vol_ions, vol_sol)

        Z = np.squeeze(z)

        # Progress print
        if iteration % 10 == 0:
            print('converg at iter = ' + str(iteration) + ' is ' + str(convergence_tot))
        iteration += 1

    # Compute full profiles after convergence
    n_profile, uself_profile = num_concn.nconc_complete(psi_profile, n_profile,
                                                        uself_bulk, n_bulk, valency,
                                                        rad_ions, vol_ions, vol_sol,
                                                        domain, epsilon_s, epsilon_p)
    eta_profile = calculate.eta_profile(n_profile, vol_ions, vol_sol)
    q_profile = calculate.charge_density(n_profile, valency)

    # Residual of Gauss's law
    res = calculate.res_2plate(psi_profile, q_profile, bounds, sigma_1, sigma_2, epsilon_s)
    print("Gauss's law residual for MGRF = " + str(res))

    # Extend profiles for edges
    psi_profile, n_profile, uself_profile, Z, surface_psi = calculate.profile_extender(
        psi_profile, n_profile, uself_profile, bounds, np.max(rad_ions), N_exc
    )

    return psi_profile, n_profile, uself_profile, q_profile, Z, res, surface_psi
