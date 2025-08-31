from packages import *
from numerical_param import *
from physical_param import *
import pb_2plate
import dh_2plate
import mgrf_2plate
import energy_2plate
import selfe_2plate

if __name__ == "__main__":
    # Argument parser to accept input files for physical parameters
    parser = argparse.ArgumentParser(description='Code to calculate EDL structure using MGRF Theory with mean-field PB as an initial guess')
    parser.add_argument('input_files', nargs='+', help='Paths to the input files for physical parameters')
    args = parser.parse_args()

    # Add folder containing input files to Python path
    folder_path = os.path.dirname(args.input_files[0])
    sys.path.insert(0, folder_path)

    # Load the physical input configuration from the first file
    module_name = os.path.splitext(os.path.basename(args.input_files[0]))[0]
    input_physical = importlib.import_module(module_name)
    variables = {name: value for name, value in input_physical.__dict__.items() if not name.startswith('__')}
    (locals().update(variables))

    # Print all the physical and numerical parameters

    print(f'cb1_d = {cb1_d}')
    print(f'cb2_d = {cb2_d}')
    print(f'sigma_in1_d = {sigma_in1_d}')
    print(f'sigma_in2_d = {sigma_in2_d}')
    print(f'sigma_f1_d = {sigma_f1_d}')
    print(f'sigma_f2_d = {sigma_f2_d}')
    print(f'domain_in_d = {domain_in_d}')
    print(f'domain_d = {domain_d}')
    print(f'valency = {valency}')
    print(f'rad_ions_d = {rad_ions_d}')
    print(f'rad_sol_d = {rad_sol_d}')
    print(f'epsilonr_s_d = {epsilonr_s_d}')
    print(f'electostatic_coupling = {(2*pi*pow(valency[0],3)*(l_c**2)*sqrt(abs(sigma_f1_d*sigma_f2_d))/ec)}')

    print(f'ncc_cutoff_mgrf = {ncc_cutoff_mgrf}')
    print(f'ncc_cutoff_greens= {ncc_cutoff_greens}')
    print(f'num_ratio = {num_ratio}')
    print(f'tolerance = {tolerance}')
    print(f'tolerance_pb = {tolerance_pb}')
    print(f'tolerance_greens= {tolerance_greens}')
    print(f'N_grid = {N_grid}')


    # Start timing the computation
    start = timeit.default_timer()

    # Compute Debye-Hückel (DH) initial guess for electrostatic potential and ion concentration profiles
    psi_complete, nconc_complete, z = dh_2plate.dh_2plate(n_bulk, valency, sigma_f1, sigma_f2, N_grid, domain, epsilon_s)
    print(f'DH_done: psi_complete[0:5] = {psi_complete[0:5]}')

    # Compute mean-field Poisson-Boltzmann (PB) profiles using DH as initial guess
    psi_complete, nconc_complete, z = pb_2plate.pb_2plate(psi_complete, n_bulk, valency, sigma_f1, sigma_f2, domain, epsilon_s)
    print(f'PB_done: psi_complete[0:5] = {psi_complete[0:5]}')

    # Compute MGRF (Modified Generalized Random Field) profiles using PB as initial guess
    psi_complete, nconc_complete, uself_complete, q_complete, z, res = mgrf_2plate.mgrf_2plate(
        psi_complete, nconc_complete, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s
    )
    print(f'MGRF_done: psi_complete[0:5] = {psi_complete[0:5]}')

    # Compute grand potential (free energy) for MGRF solution
    grandfe = energy_2plate.grandfe_mgrf_2plate(
        psi_complete, nconc_complete, uself_complete, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s
    )
    print(f'Grand potential (grandfe) = {grandfe}')


    # Print elapsed time
    stop = timeit.default_timer()
    print(f'Time elapsed: {stop - start}')
    
    # Determine output directory and filename based on whether mixture of salts is present
    if cb2_d != 0:
        output_dir = os.getcwd() + f'/results-mixture{abs(valency[0])}_{abs(valency[1])}_{abs(valency[2])}_{abs(valency[3])}'
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{rad_ions_d[2]:.2f}_{rad_ions_d[3]:.2f}_{sigma_f1_d:.5f}_{sigma_f2_d:.5f}'
    else:
        output_dir = os.getcwd() + f'/results{abs(valency[0])}_{abs(valency[1])}'
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{sigma_f1_d:.5f}_{sigma_f2_d:.5f}'

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Save all results in an HDF5 file (SI units and nondimensional)
    with h5py.File(output_dir + f'/mgrf_{file_name}.h5', 'w') as file:

        # Scalar parameters as attributes
        file.attrs['ec_charge'] = ec
        file.attrs['char_length'] = l_b
        file.attrs['beta'] = beta
        file.attrs['epsilon_s'] = epsilonr_s_d
        file.attrs['epsilon_p'] = epsilonr_s_d
        file.attrs['cb1'] = cb1_d
        file.attrs['cb2'] = cb2_d
        file.attrs['domain'] = domain_d

        # Numerical parameters as attributes
        file.attrs['s_conv'] = s_conv
        file.attrs['N_grid'] = N_grid
        file.attrs['quads'] = quads
        file.attrs['grandfe_quads'] = grandfe_quads
        file.attrs['dealias'] = dealias
        file.attrs['ncc_cutoff_pb'] = ncc_cutoff_pb
        file.attrs['ncc_cutoff_mgrf'] = ncc_cutoff_mgrf
        file.attrs['ncc_cutoff_greens'] = ncc_cutoff_greens
        file.attrs['num_ratio'] = num_ratio
        file.attrs['tolerance'] = tolerance
        file.attrs['tolerance_pb'] = tolerance_pb
        file.attrs['tolerance_num'] = tolerance_num
        file.attrs['tolerance_greens'] = tolerance_greens

        # Parameter arrays
        file.create_dataset('valency', data=valency)
        file.create_dataset('radii', data=rad_ions_d)
        file.create_dataset('volumes', data=np.concatenate((vol_ions_d, [vol_sol_d])))
        file.create_dataset('surface_charges', data=np.array([sigma_f1_d, sigma_f2_d]))

        # Spatial profiles (SI units)
        file.create_dataset('z_d', data=z*l_c)
        file.create_dataset('psi_d', data=psi_complete*psi_c)
        file.create_dataset('nconc_d', data=nconc_complete*nconc_c/N_A)
        file.create_dataset('uself_d', data=uself_complete*(1/beta))
        file.create_dataset('charge_d', data=q_complete*(nconc_c*ec))

        # Spatial profiles (nondimensional)
        file.create_dataset('z', data=z)
        file.create_dataset('psi', data=psi_complete)
        file.create_dataset('nconc', data=nconc_complete)
        file.create_dataset('uself', data=uself_complete)
        file.create_dataset('charge', data=q_complete)

        # Free energy and Gauss's law residual
        file.attrs['grandfe'] = grandfe  # nondimensional
        file.attrs['grandfe_d'] = grandfe*(1/beta)  # SI units
        file.attrs['residual'] = res

