from packages import *
import calculate
from numerical_param import*


def Gcap_free(grid_points,s,domain,epsilon):#\hat{Go}
    
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

    # Fields for G(Pz or dzlog(U))
    Pz = dist.Field(name = 'Pz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')

    # Differential equation for Pz (U)
    problem = d3.NLBVP([Pz,tau_1],namespace = locals())
    problem.add_equation("-dz(Pz) + s*s + lift(tau_1,-1) = Pz**2")

    # Boundary conditions for Pz/U
    problem.add_equation("Pz(z=0) = s")

    # Initial guess
    Pz['g'] = s

    # Solver
    solver0 = problem.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm0 = np.inf
    Pz.change_scales(dealias)
    while pert_norm0 > tolerance_greens:
        solver0.newton_iteration()
        pert_norm0 = sum(pert0.allreduce_data_norm('c',2) for pert0 in solver0.perturbations)

    Pz.change_scales(1)
    Pz = Pz['g']
    Qz = -Pz

    ## Sturm-Liouville for G
    G = (1 / epsilon) * np.true_divide(1,Pz - Qz)

    del z,Pz,Qz,tau_1,dz,lift_basis,lift,problem,solver0,pert_norm0
    gc.collect()

    return G

def Gcap_full(n_bulk_profile, n_bulk, valency, s, domain,epsilon): # \hat{G}

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

    omega_sqr = dist.Field(bases = zbasis)
    omega_sqr['g'] = s * s + calculate.kappa_sqr_profile(n_bulk_profile,valency,epsilon)
    omega_b = np.sqrt(s * s + calculate.kappa_sqr(n_bulk,valency,epsilon))

    # Fields for G(Pz or log(U))
    Pz = dist.Field(name = 'Pz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')

    # Differential equation for Pz/U
    problem = d3.NLBVP([Pz,tau_1],namespace = locals())
    problem.add_equation("-dz(Pz) + omega_sqr + lift(tau_1,-1) = Pz**2")

    # Boundary conditions for Pz
    problem.add_equation("Pz(z=0) = omega_b")

    # Initial guess for Pz
    Pz['g'] = omega_b

    # Solver
    solver1 = problem.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm1 = np.inf
    Pz.change_scales(dealias)
    p = 0
    while pert_norm1 > tolerance_greens:
        p = p + 1
        solver1.newton_iteration()
        pert_norm1 = sum(pert1.allreduce_data_norm('c',2) for pert1 in solver1.perturbations)

    Pz.change_scales(1)
    Pz = Pz['g']

    # Fields for G(Qz or log(V))
    Qz = dist.Field(name = 'Qz',bases = zbasis)
    tau_1 = dist.Field(name = 'tau_1')

    # Differential equation for Qz/V
    problem1 = d3.NLBVP([Qz,tau_1],namespace = locals())
    problem1.add_equation("-dz(Qz) + omega_sqr + lift(tau_1,-1) = Qz**2")

    # Boundary conditions for Qz
    problem1.add_equation("Qz(z=Lz) = -omega_b")

    # Initial guess for Qz
    Qz['g'] = -omega_b

    # Solver
    solver2 = problem1.build_solver(ncc_cutoff = ncc_cutoff_greens)
    pert_norm2 = np.inf
    Qz.change_scales(dealias)
    q = 1
    while pert_norm2 > tolerance_greens:
        q = q + 1
        solver2.newton_iteration()
        pert_norm2 = sum(pert2.allreduce_data_norm('c',2) for pert2 in solver2.perturbations)

    Qz.change_scales(1)
    Qz = Qz['g']

    ## Sturm-Liouville for G
    G = (1 / epsilon) * np.true_divide(1,Pz - Qz)

    del z,Pz,Qz,tau_1,dz,lift_basis,lift,problem,solver1,solver2,pert_norm2,pert_norm1
    gc.collect()

    return G



