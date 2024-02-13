from packages import *

## Numerical Parameters

s_conv = 1e6  # approx for infinity for fourier inverse of greens function
V_conv = log(s_conv + 1)  # we do fourier inverse integration in the logspace
quads = 24  # no of legendre gauss quadrature points for fourier inverse of greens function
<<<<<<< HEAD
N_grid = 1024 # has to be even, since we often use 3/2 for dealiasing
=======
N_grid =256 # has to be even, since we often use 3/2 for dealiasing
>>>>>>> 38035e4e3cbc6229b250b5686a828f34a1c56977
dealias = 3/2  # dealiasing factor for dedalus
ncc_cutoff_mgrf = 1e-2  # some cutoff parameter for non-constant coefficients on LHS of NLBVP of MGRF
ncc_cutoff_pb = 1e-1 # some cutoff parameter for non-constant coefficients on LHS of NLBVP of PB
ncc_cutoff_greens = 1e-1 # some cutoff parameter for non-constant coefficients on LHS of NLBVP of G
num_ratio = 0.1  # mixing ratio of new to old in nconc_mgrf
selfe_ratio = 0.03 # mixing ratio of self-energy (new to old) in outermost loop of pb_mgrf
eta_ratio = 0.1  # mixing ratio of eta (new to old) in outermost loop of pb_mgrf
grandfe_quads = 25  # no of legendre gauss quadrature points for free energy calculation
cores = 24  # no of parallel processes in which you want to divide fourier inverse calculation
<<<<<<< HEAD
tolerance = pow(10,-4)  # tolerance for outermost loop for pb_mgrf
=======
tolerance = pow(10,-5)  # tolerance for outermost loop for pb_mgrf
>>>>>>> 38035e4e3cbc6229b250b5686a828f34a1c56977
tolerance_pb = pow(10,-7)  # tolerance for inner loop mgrf/outermost loop pb_
tolerance_num = pow(10,-4)  # tolerance for nconc_mgrf iteration loop
tolerance_greens = pow(10,-7)  # tolerance for nonlinear problem for greens function
iter_max = pow(10,7)  # maximum no of iterations for any  iteration loop
