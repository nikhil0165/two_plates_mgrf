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

    N_exc = np.nonzero(n_profile[:,0])[0][0]

    # Grand free energy
    grandfe = energy_2plate.grandfe_pb_2plate(psi_profile, n_profile, n_bulk, valency, sigma_f1, sigma_f2, domain)
    
    print(f'grandfe = {grandfe}')

    elapsed_time = timeit.default_timer() - start
    print(f'Time = {elapsed_time}')


if cb2_d != 0:
    output_dir = os.getcwd() + '/results-pb-mixture' + str(abs(valency[0]))+ '_' + str(abs(valency[1])) + '-' + str(abs(valency[2]))+ '_' + str(abs(valency[3]))
    file_name =  str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5)) + '_' + str(round(float(domain_d), 2)) + '_' + str(round(rad_ions_d[2] / pow(10, -10), 2)) + '_' + str(round(sigma_f1_d, 5)) + '_' + str(round(sigma_f2_d, 5))
else:
    output_dir = os.getcwd() + '/results-pb' + str(abs(valency[0])) + '_' + str(abs(valency[1]))
    file_name = str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5))  + '_' + str(round(float(domain_d), 2)) + '_' + str(round(rad_ions_d[0] / pow(10, -10), 2)) + '_' + str(round(rad_ions_d[1] / pow(10, -10), 2)) + '_' + str(round(sigma_f1_d, 5)) + '_' + str(round(sigma_f2_d, 5))

# Create the output directory if it doesn't exist

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Writing everything in SI units
with h5py.File(output_dir + '/pb_' + file_name + '.h5', 'w') as file:

    # Storing scalar variables as attributes of the root group
    file.attrs['ec_charge'] = ec
    file.attrs['char_length'] = l_b
    file.attrs['beta'] = beta
    file.attrs['epsilon_s'] = epsilonr_s_d
    file.attrs['epsilon_p'] = epsilonr_s_d
    file.attrs['cb1'] = cb1_d*0.001
    file.attrs['cb2'] = cb2_d*0.001

    # Storing numerical parameters as attributes of the root group
    file.attrs['N_grid'] = N_grid
    file.attrs['dealias'] = dealias
    file.attrs['ncc_cutoff_pb'] = ncc_cutoff_pb
    file.attrs['tolerance_pb'] = tolerance_pb

    # Storing parameter arrays
    file.create_dataset('valency', data = valency)
    file.create_dataset('surface_charges', data = np.array([sigma_f1_d,sigma_f2_d]))

    # Store all spatial profiles  (SI units)
    file.create_dataset('z_d', data = z*l_c)
    file.create_dataset('psi_d', data = psi_profile*psi_c)
    file.create_dataset('nconc_d', data = n_profile*nconc_c/N_A)

    # Store all spatial profiles (non-dimensional)
    file.create_dataset('z', data = z)
    file.create_dataset('psi', data = psi_profile)
    file.create_dataset('nconc', data = n_profile)

    # Store free energy
    file.attrs['grandfe'] = grandfe # nondimensional
    file.attrs['grandfe_d'] = grandfe*(1/beta) # SI units


stop = timeit.default_timer()
print('Time: ', stop - start)



