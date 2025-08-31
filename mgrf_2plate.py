from packages import *
import num_concn
import calculate
import selfe_2plate
import selfe_bulk
from numerical_param import*

def mgrf_2plate(psi_guess,nconc_guess,n_bulk,valency,rad_ions,vol_ions,vol_sol,sigma_1,sigma_2, domain, epsilon):  
    # psi_guess: initial guess for potential (from mean-field PB solution)
    # nconc_guess: initial guess for ion concentrations
    # n_bulk: bulk ion concentration
    # valency: list/array of ionic valencies
    # rad_ions: radii of ions
    # vol_ions: ionic volumes
    # vol_sol: solvent volume
    # sigma_1, sigma_2: surface charge densities at the two plates
    # domain: separation between plates
    # epsilon: dielectric permittivity

    grid_points = len(psi_guess)                # Number of collocation grid points
    bounds = (0,domain)                         # Spatial domain (z from 0 to Lz)
    Lz = bounds[1]
    slope1 = -sigma_1/epsilon                   # Boundary slope from surface charge at z=0
    slope2 = -sigma_2/epsilon                   # Boundary slope from surface charge at z=Lz

    # Initialize profiles with PB guess
    psi_profile = np.copy(psi_guess)           
    n_profile= nconc_guess
    eta_guess=calculate.eta_profile(nconc_guess,vol_ions,vol_sol)   # local volume fraction
    uself_guess = selfe_2plate.uself_complete(nconc_guess,n_bulk,rad_ions,valency,domain,epsilon) # initial self-energy profile

    print('selfe_done')

    # Bulk reference properties
    n_bulk_numerical = np.multiply(np.ones((grid_points,len(valency))),n_bulk)  
    uself_bulk = np.mean(selfe_bulk.uselfb_numerical(n_bulk_numerical, n_bulk, rad_ions, valency, domain,epsilon), axis=0) # bulk self-energy
    eta_bulk = calculate.eta_loc(n_bulk, vol_ions, vol_sol) # bulk packing fraction

    # Check if ions and solvent have equal volumes (simplifies calculation)
    equal_vols = np.all(np.abs(vol_ions - vol_sol) < vol_sol * 1e-5)

    Z = None  # Will hold the spatial grid for return
    convergence_tot = np.inf   # Overall convergence criterion
    iteration = 1                     # Iteration counter

    # ----------------- Outer self-consistency loop -----------------
    while(convergence_tot  > tolerance):

        # Dedalus bases & distributor setup
        coords = d3.CartesianCoordinates('z')
        dist = d3.Distributor(coords, dtype=np.float64)  
        zbasis = d3.Chebyshev(coords['z'], size=grid_points, bounds=bounds, dealias=dealias)

        # Field variables
        z = dist.local_grids(zbasis)        # physical grid
        psi = dist.Field(name='psi', bases=zbasis)   # electrostatic potential
        tau_1 = dist.Field(name='tau_1')    # auxiliary fields for boundary conditions
        tau_2 = dist.Field(name='tau_2')

        # Substitutions for differentiation and lifting
        dz = lambda A: d3.Differentiate(A, coords['z'])
        lift_basis = zbasis.derivative_basis(2)
        lift = lambda A, n: d3.Lift(A, lift_basis, n)

        # Fields for expansion coefficients in RHS of Poisson equation
        c0 = dist.Field(bases = zbasis)
        c1 = dist.Field(bases = zbasis)

        # Update concentration closure with current guesses
        _, coeffs = num_concn.nconc_mgrf(psi_profile, uself_guess, eta_guess, uself_bulk, n_bulk, valency, vol_ions,eta_bulk, equal_vols)
        coeffs = coeffs/epsilon   # Normalize by permittivity

        # Store coefficients into Dedalus fields
        c0['g'] = np.squeeze(coeffs[:, 0])
        c1['g'] = np.squeeze(coeffs[:, 1])

        # Boltzmann factors for ions
        boltz0 = lambda psi: np.exp(-valency[0] * psi)
        boltz1 = lambda psi: np.exp(-valency[1] * psi)

        # If symmetric 4-component electrolyte, add extra coefficients
        if len(valency) == 4:
            c2 = dist.Field(bases = zbasis)
            c3 = dist.Field(bases = zbasis)
            boltz2 = lambda psi: np.exp(-valency[2] * psi)
            boltz3 = lambda psi: np.exp(-valency[3] * psi)

        # ----------------- PDE setup -----------------
        problem = d3.NLBVP([psi, tau_1, tau_2], namespace=locals())

        # Poisson equation with MGRF-modified RHS
        if len(valency)==2:
            problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = c0*boltz0(psi) + c1*boltz1(psi)")
        if len(valency)==4:
            c2['g'] = np.squeeze(coeffs[:,2])
            c3['g'] = np.squeeze(coeffs[:,3])
            problem.add_equation("-lap(psi) + lift(tau_1,-1) + lift(tau_2,-2) = c0*boltz0(psi) + c1*boltz1(psi) + c2*boltz2(psi) + c3*boltz3(psi)")

        # Boundary conditions from surface charges
        problem.add_equation("dz(psi)(z=0) = slope1")
        problem.add_equation("dz(psi)(z=Lz) = -slope2")

        # Initial guess for psi field
        psi['g'] = psi_profile

        # Build Newton solver
        solver = problem.build_solver(ncc_cutoff=ncc_cutoff_mgrf)
        pert_norm = np.inf
        psi.change_scales(dealias)   # upscale for solving
        s = 0

        # ----------------- Newton loop for PDE -----------------
        while pert_norm > tolerance_pb:
            solver.newton_iteration()
            pert_norm = sum(pert.allreduce_data_norm('c', 2) for pert in solver.perturbations)
            s  = s +1

        # Downscale and collect updated psi
        psi.change_scales(1)
        psi_profile = psi.allgather_data('g')

        # Sanity check
        if (np.any(np.isnan(psi_profile))):
            print('nan in psi')

        # Update number concentrations with new psi
        n_profile,_ = num_concn.nconc_mgrf(psi_profile, uself_guess, eta_guess, uself_bulk, n_bulk, valency, vol_ions, eta_bulk,equal_vols)

        # Outer-loop convergence check
        convergence_tot = np.true_divide(np.linalg.norm(n_profile - nconc_guess),np.linalg.norm(n_profile))

        # Mixing update for stability
        nconc_guess = num_ratio*n_profile + (1-num_ratio)*nconc_guess

        # Update self-energy and eta profile
        uself_guess = selfe_2plate.uself_complete(nconc_guess, n_bulk,rad_ions, valency, domain,epsilon)
        eta_guess = calculate.eta_profile(nconc_guess, vol_ions, vol_sol)

        Z = np.squeeze(z)   # Save grid

        # Iteration printout
        if iteration % 10 == 0:
            print(f'Iteration {iteration}: convergence = {convergence_tot:.3e}')
        iteration += 1

    # After convergence, recompute final fields
    n_profile, uself_profile = num_concn.nconc_complete(psi_profile, n_profile, uself_bulk, n_bulk, valency, rad_ions, vol_ions, vol_sol, domain, epsilon)
    eta_profile = calculate.eta_profile(n_profile, vol_ions, vol_sol)
    q_profile = calculate.charge_density(n_profile, valency)

    # Verify Gauss’s law residual (consistency check)
    res= calculate.res_2plate(psi_profile,q_profile,bounds,sigma_1,sigma_2,epsilon)
    print("Gauss's law residual for MGRF = " + str(res))

    return psi_profile, n_profile,uself_profile,q_profile,Z, res
