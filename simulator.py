import mgrf_2plate
from numerical_param import *
import energy_2plate
from physical_param import *
import calculate

if __name__ == "__main__":

    # Argument parser to accept input files
    parser = argparse.ArgumentParser(description='Code to calculate EDL structure using MGRF Theory with another MGRF solution as an initial guess')
    parser.add_argument('input_files', nargs='+', help='Paths to the input files for physical parameters')
    args = parser.parse_args()

    folder_path = os.path.dirname(args.input_files[0])
    sys.path.insert(0, folder_path)

    # Load the physical input configuration from the first file in the list
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

    # Determine file directory and name based on presence of secondary salt
    if cb2_d != 0:
        file_dir = os.getcwd() + f'/results-mixture{abs(valency[0])}_{abs(valency[1])}_{abs(valency[2])}_{abs(valency[3])}'
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_in_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{rad_ions_d[2]:.2f}_{rad_ions_d[3]:.2f}_{sigma_in1_d:.5f}_{sigma_in2_d:.5f}'
    else:
        file_dir = os.getcwd() + f'/results{abs(valency[0])}_{abs(valency[1])}'
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_in_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{sigma_in1_d:.5f}_{sigma_in2_d:.5f}'

    # Load previously computed MGRF solution
    with h5py.File(file_dir + f'/mgrf_{file_name}.h5', 'r') as file:
        psi_complete = np.array(file['psi'])
        nconc_complete = np.array(file['nconc'])
        uself_complete = np.array(file['uself'])
        grandfe = file.attrs['grandfe']

    # Start timer
    start = timeit.default_timer()

    # Compute MGRF profiles using previous solution as initial guess
    psi_complete, nconc_complete, uself_complete, q_complete, z, res = mgrf_2plate.mgrf_2plate(
        psi_complete, nconc_complete, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s
    )
    print('MGRF_done')
    print(f'psi_complete[0:5] = {psi_complete[0:5]}')

    # Compute grand potential (free energy) for new MGRF solution
    grandfe = energy_2plate.grandfe_mgrf_2plate(
        psi_complete, nconc_complete, uself_complete, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s
    )
    print(f'grandfe = {grandfe}')

    stop = timeit.default_timer()
    print(f'Time elapsed: {stop - start}')

    # Prepare filename for saving new results
    if cb2_d != 0:
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{rad_ions_d[2]:.2f}_{rad_ions_d[3]:.2f}_{sigma_f1_d:.5f}_{sigma_f2_d:.5f}'
    else:
        file_name = f'{cb1_d:.9f}_{cb2_d:.5f}_{float(domain_d):.2f}_{rad_ions_d[0]:.2f}_{rad_ions_d[1]:.2f}_{sigma_f1_d:.5f}_{sigma_f2_d:.5f}'

    # Writing updated MGRF solution to HDF5 file
    with h5py.File(file_dir + f'/mgrf_{file_name}.h5', 'w') as file:
        # Scalar parameters
        file.attrs['ec_charge'] = ec
        file.attrs['char_length'] = l_b
        file.attrs['beta'] = beta
        file.attrs['epsilon_s'] = epsilonr_s_d
        file.attrs['epsilon_p'] = epsilonr_s_d
        file.attrs['cb1'] = cb1_d
        file.attrs['cb2'] = cb2_d
        file.attrs['domain'] = domain_d

        # Numerical parameters
        file.attrs['s_conv'] = s_conv
        file.attrs['N_grid'] = N_grid
        file.attrs['quads'] = quads
        file.attrs['grandfe_quads'] = grandfe_quads
        file.attrs['dealias'] = dealias
        file.attrs['ncc_cutoff_pb'] = ncc_cutoff_pb
        file.attrs['ncc_cutoff_mgrf'] = ncc_cutoff_mgrf
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
        file.create_dataset('n_bulk', data=n_bulk)

        # Free energy and residual
        file.attrs['grandfe'] = grandfe
        file.attrs['grandfe_d'] = grandfe*(1/beta)
        file.attrs['residual'] = res
