# EDL simulation using DH, PB, and MGRF theories with full output
from packages import *
from numerical_param import *
import pb_2plate
import dh_2plate
import mgrf_2plate
import energy_2plate
import selfe_2plate
import calculate
from physical_param import *

if __name__ == "__main__":
    start = timeit.default_timer()

    # Argument parser to accept input files
    parser = argparse.ArgumentParser(description='Calculate EDL using MGRF Theory with PB as initial guess')
    parser.add_argument('input_files', nargs='+', help='Paths to input files for physical parameters')
    args = parser.parse_args()

    folder_path = os.path.dirname(args.input_files[0])
    sys.path.insert(0, folder_path)

    # Load physical input configuration
    module_name = os.path.splitext(os.path.basename(args.input_files[0]))[0]
    input_physical = importlib.import_module(module_name)
    variables = {name: value for name, value in input_physical.__dict__.items() if not name.startswith('__')}
    locals().update(variables)

    # Print physical & numerical parameters
    print(f'cb1_d = {cb1_d}, cb2_d = {cb2_d}')
    print(f'sigma_f1_d = {sigma_f1_d}, sigma_f2_d = {sigma_f2_d}')
    print(f'domain_d = {domain_d}, valency = {valency}')
    print(f'rad_ions_d = {rad_ions_d}, rad_sol_d = {rad_sol_d}')
    print(f'epsilonr_s_d = {epsilonr_s_d}, epsilonr_p_d = {epsilonr_p_d}')
    print(f'electrostatic_coupling = {(2*pi*pow(valency[0],3)*(l_c**2)*sqrt(abs(sigma_f1_d*sigma_f2_d))/ec)}')
    print(f'N_grid = {N_grid}, tolerance = {tolerance}, tolerance_pb = {tolerance_pb}')

    # DH initial guess
    psi_profile, n_profile, z, surface_psi = dh_2plate.dh_2plate(n_bulk, valency, sigma_f1, sigma_f2, N_grid, domain, epsilon_s)
    print('DH_done', f'surface_psi = {surface_psi}')

    # PB refinement
    psi_profile, n_profile, z, surface_psi = pb_2plate.pb_2plate(psi_profile, n_bulk, valency, sigma_f1, sigma_f2, domain, epsilon_s)
    print('PB_done', f'surface_psi = {surface_psi}')

    # MGRF iteration
    psi_profile, n_profile, uself_profile, q_profile, z, res, surface_psi = mgrf_2plate.mgrf_2plate(
        psi_profile, n_profile, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s, epsilon_p)
    print('MGRF_done', f'surface_psi = {surface_psi}')

    N_exc = np.nonzero(n_profile[:,0])[0][0]

    # Grand free energy
    grandfe = energy_2plate.grandfe_mgrf_2plate(
        psi_profile, n_profile, uself_profile, n_bulk, valency, rad_ions, vol_ions, vol_sol, sigma_f1, sigma_f2, domain, epsilon_s, epsilon_p)
    print(f'grandfe = {grandfe}')

    elapsed_time = timeit.default_timer() - start
    print(f'Time = {elapsed_time}')

    # Output directory & file name
    if cb2_d != 0:
        output_dir = os.path.join(os.getcwd(), f'results-mixture{abs(valency[0])}_{abs(valency[1])}_{abs(valency[2])}_{abs(valency[3])}')
        file_name = f"{cb1_d}_{cb2_d}_{domain_d}_{rad_ions_d[0]}_{rad_ions_d[1]}_{rad_ions_d[2]}_{rad_ions_d[3]}_{sigma_f1_d}_{sigma_f2_d}_{epsilonr_s_d}_{epsilonr_p_d}_{N_grid}"
    else:
        output_dir = os.path.join(os.getcwd(), f'results{abs(valency[0])}_{abs(valency[1])}')
        file_name = f"{cb1_d}_{cb2_d}_{domain_d}_{rad_ions_d[0]}_{rad_ions_d[1]}_{sigma_f1_d}_{sigma_f2_d}_{epsilonr_s_d}_{epsilonr_p_d}_{N_grid}"

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Write results to HDF5
    with h5py.File(os.path.join(output_dir, f'mgrf_{file_name}.h5'), 'w') as file:
        # Scalar attributes
        file.attrs.update({
            'ec_charge': ec,
            'char_length': l_b,
            'beta': beta,
            'epsilon_s': epsilonr_s_d,
            'epsilon_p': epsilonr_s_d,
            'cb1': cb1_d,
            'cb2': cb2_d,
            'domain': domain,
            'domain_d': domain * l_c,
            's_conv': s_conv,
            'N_grid': len(psi_profile)-2*N_exc,
            'N_exc': N_exc,
            'quads': quads,
            'grandfe_quads': grandfe_quads,
            'dealias': dealias,
            'ncc_cutoff_pb': ncc_cutoff_pb,
            'ncc_cutoff_mgrf': ncc_cutoff_mgrf,
            'ncc_cutoff_greens': ncc_cutoff_greens,
            'num_ratio': num_ratio,
            'tolerance': tolerance,
            'tolerance_pb': tolerance_pb,
            'tolerance_num': tolerance_num,
            'tolerance_greens': tolerance_greens,
            'time': elapsed_time
        })

        # Parameter arrays
        file.create_dataset('valency', data=valency)
        file.create_dataset('radii', data=rad_ions_d)
        file.create_dataset('volumes', data=np.concatenate((vol_ions_d, [vol_sol_d])))
        file.create_dataset('surface_charges', data=[sigma_f1_d, sigma_f2_d])

        # Spatial profiles (SI)
        file.create_dataset('z_d', data=z*l_c)
        file.create_dataset('psi_d', data=psi_profile*psi_c)
        file.create_dataset('nconc_d', data=n_profile*nconc_c/N_A)
        file.create_dataset('uself_d', data=uself_profile*(1/beta))
        file.create_dataset('charge_d', data=q_profile*(nconc_c*ec))

        # Spatial profiles (non-dimensional)
        file.create_dataset('z', data=z)
        file.create_dataset('psi', data=psi_profile)
        file.create_dataset('nconc', data=n_profile)
        file.create_dataset('uself', data=uself_profile)
        file.create_dataset('charge', data=q_profile)
        file.create_dataset('surface_psi', data=surface_psi)

        # Free energy & residual
        file.attrs.update({'grandfe': grandfe, 'grandfe_d': grandfe*(1/beta), 'residual': res})
