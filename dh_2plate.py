from packages import *
import calculate
import num_concn

def dh_2plate(n_bulk,valency,sigma_1,sigma_2,grid_points,domain,epsilon):

    bounds = (0,domain)
    Lz = bounds[1]
    slope1 = -sigma_1/epsilon
    slope2 = -sigma_2/epsilon

    #Making the LHS
    kappa_2 = calculate.kappa_sqr(n_bulk,valency,epsilon)

    # Bases
    coords = d3.CartesianCoordinates('z')
    dist = d3.Distributor(coords, dtype=np.float64) # No mesh for serial / automatic parallelization
    zbasis = d3.Chebyshev(coords['z'], size= grid_points,bounds =bounds)

    # Fields
    z = dist.local_grids(zbasis)
    psi = dist.Field(name='psi', bases=zbasis)
    tau_1 = dist.Field(name='tau_1')
    tau_2 = dist.Field(name='tau_2')

    # Substitutions
    dz = lambda A: d3.Differentiate(A, coords['z'])
    lift_basis = zbasis.derivative_basis(2)
    lift = lambda A, n: d3.Lift(A, lift_basis, n)

    # PDE Setup
    problem = d3.LBVP([psi,tau_1, tau_2], namespace=locals())
    problem.add_equation("-lap(psi) + kappa_2*psi + lift(tau_1,-1) + lift(tau_2,-2) = 0")

    # Boundary conditions
    problem.add_equation("dz(psi)(z=0) = slope1")
    problem.add_equation("dz(psi)(z=Lz) = -slope2")

    solver = problem.build_solver()
    solver.solve()

    psi_profile = psi.allgather_data('g')
    n_profile = num_concn.nconc_pb(psi_profile,valency,n_bulk)

    surface1_psi = psi(z = 0).evaluate()['g'][0]
    surface2_psi = psi(z = bounds[1]).evaluate()['g'][0]

    return psi_profile,n_profile,np.squeeze(z), [surface1_psi,surface2_psi]


