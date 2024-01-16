from packages import *
import calculate
from numerical_param import*


def Gcap_free(grid_points, s, domain,epsilon): # function for \hat{Go}
    bounds = (0, domain)
    Lz = bounds[1]

    # Bases
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)  
    zbasis = d3.Chebyshev(coords['z'], size = grid_points, bounds = bounds, dealias = dealias)

    # General fields
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # Fields for G(P or log(U))
    P = dist.Field(name = 'P', bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')
    tau_2 = dist.Field(name = 'tau_2')

    # Differential equation for P (U)
    problem = d3.NLBVP([P, tau_1, tau_2], namespace = locals())
    problem.add_equation("-lap(P) + s*s + lift(tau_1,-1) + lift(tau_2,-2) = grad(P)@grad(P)")

    # Boundary conditions for P/U
    problem.add_equation("P(z=0) = 0")
    problem.add_equation("dz(P)(z=0) = s")

    # Initial guess
    P['g'] = s * np.squeeze(z)

    # Solver
    solver = problem.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm = np.inf
    P.change_scales(dealias)
    while pert_norm > tolerance_greens:
        solver.newton_iteration()
        pert_norm = sum(pert0.allreduce_data_norm('c', 2) for pert0 in solver.perturbations)

    Pz = d3.Gradient(P).evaluate()
    Pz.change_scales(1)
    Pz = Pz.allgather_data('g')[0]

    P.change_scales(1)
    Pg = P['g']

    U = [Decimal(str(value)).exp() for value in Pg]
    V = U[::-1]

    Uz = [U[i] * Decimal(str(Pz[i])) for i in range(0, grid_points)]
    Vz = Uz[::-1]
    Vz = [-x for x in Vz]

    ## Sturm-Liouville for Go

    C = [Decimal(1) / (U[i] * Vz[i] - V[i] * Uz[i]) for i in range(0, grid_points)]
    G = np.array([C[i] * U[i] * V[i] * Decimal(-1 / epsilon) for i in range(0, grid_points)], dtype = np.float64)

    del z,P,tau_1,tau_2,dz,lift_basis,lift,problem,solver,pert_norm,Pz,U,V,Uz,Vz,C
    gc.collect()

    return G


def Gcap_full(n_profile,n_bulk,valency, s, domain,epsilon):# function for \hat{G}

    grid_points = len(n_profile)
    bounds = (0, domain)
    Lz = bounds[1]

    # Bases
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype = np.float64)  
    zbasis = d3.Chebyshev(coords['z'], size = grid_points, bounds = bounds, dealias = dealias)

    # General fields
    z = dist.local_grids(zbasis)
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    omega_sqr = dist.Field(bases = zbasis)
    omega_sqr['g'] = s * s + calculate.kappa_sqr_profile(n_profile, valency, epsilon)
    omega_b = np.sqrt(s * s + calculate.kappa_sqr(n_bulk, valency, epsilon))
    omega_min = min(omega_b,np.min(np.sqrt(omega_sqr['g']))) 

    # Fields for G(P or log(U))
    P = dist.Field(name = 'P', bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')
    tau_2 = dist.Field(name = 'tau_2')

    # Differential equation for P/U
    problem1 = d3.NLBVP([P, tau_1, tau_2], namespace = locals())
    problem1.add_equation("-lap(P) + omega_sqr + lift(tau_1,-1) + lift(tau_2,-2) =grad(P)@grad(P)")

    # Boundary conditions for P
    problem1.add_equation("P(z=0) = 0")
    problem1.add_equation("dz(P)(z=0) = s")

    # Initial guess for P
    P['g'] = omega_min * np.squeeze(z)

    # Solver
    solver1 = problem1.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm1 = np.inf
    P.change_scales(dealias)
    p=0
    #print('starting P for ' + str(s))
    while pert_norm1 > tolerance_greens:
        p =p+1
        #print('P not done for ' + str(s))
        solver1.newton_iteration()
        pert_norm1 = sum(pert1.allreduce_data_norm('c', 2) for pert1 in solver1.perturbations)

    Pz = d3.Gradient(P).evaluate()
    Pz.change_scales(1)
    Pz = Pz.allgather_data('g')[0]

    P.change_scales(1)
    Pg = P['g']

    U = [Decimal(str(value)).exp() for value in Pg]
    Uz = [U[i] * Decimal(str(Pz[i])) for i in range(0, grid_points)]

    # Fields for G(Q or log(V))
    Q = dist.Field(name = 'Q', bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')
    tau_2 = dist.Field(name = 'tau_2')

    # Differential equation for Q/V
    problem2 = d3.NLBVP([Q, tau_1, tau_2], namespace = locals())
    problem2.add_equation("-lap(Q) + omega_sqr + lift(tau_1,-1) + lift(tau_2,-2) = grad(Q)@grad(Q)")

    # Boundary conditions for Q
    problem2.add_equation("Q(z=Lz) = 0")
    problem2.add_equation("dz(Q)(z=Lz) = -s")

    # Initial guess for Q
    Q['g'] = -omega_min* (np.squeeze(z) - Lz)

    # Solver
    solver2 = problem2.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm2 = np.inf
    Q.change_scales(dealias)
    q = 1
    while pert_norm2 > tolerance_greens:
        q = q+1
        solver2.newton_iteration()
        pert_norm2 = sum(pert2.allreduce_data_norm('c', 2) for pert2 in solver2.perturbations)

    Qz = d3.Gradient(Q).evaluate()
    Qz.change_scales(1)
    Qz = Qz.allgather_data('g')[0]

    Q.change_scales(1)
    Qg = Q['g']

    V = [Decimal(str(value)).exp() for value in Qg]
    Vz = [V[i] * Decimal(str(Qz[i])) for i in range(0, grid_points)]

    ## Sturm-Liouville for G

    C = [Decimal(1) / (U[i] * Vz[i] - V[i] * Uz[i]) for i in range(0, grid_points)]
    G = np.array([C[i] * U[i] * V[i] * Decimal(-1 / epsilon) for i in range(0, grid_points)], dtype = np.float64)

    del z,P,Q,tau_1,tau_2,dz,lift_basis,lift,problem1,problem2,solver1,solver2,pert_norm1,pert_norm2,Qz,Pz,Qg,Pg,U,V,Uz,Vz,C
    gc.collect()

    return G

      
