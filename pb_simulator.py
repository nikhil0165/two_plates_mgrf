from packages import*
from physical_param import *
from numerical_param import *
import pb_2plate
import dh_2plate
import energy_2plate
import calculate
start = timeit.default_timer()

# Argument parser to accept the input files
parser = argparse.ArgumentParser(description='Code to calculate EDL structure using mean-field Theory with DH as an initial guess')
parser.add_argument('input_files', nargs='+', help='Paths to the input files for physical parameters')
args = parser.parse_args()

folder_path = os.path.dirname(args.input_files[0])
sys.path.insert(0, folder_path)

# Load the physical input configuration from the first file in the list
module_name = os.path.splitext(os.path.basename(args.input_files[0]))[0]
input_physical = importlib.import_module(module_name)
variables = {name: value for name, value in input_physical.__dict__.items() if not name.startswith('__')}
(locals().update(variables))


# The EDL structure calculations start here

psi_complete,nconc_complete,z = dh_2plate.dh_2plate(n_bulk,valency,0.25*sigma_f1,0.25*sigma_f2,N_grid,domain,epsilon_s)
print('DH_done')
print(psi_complete[0:5])

# psi_complete, nconc_complete,z = pb_2plate.pb_2plate(psi_complete,n_bulk,valency,0.25*sigma_f1,0.25*sigma_f2,domain,epsilon_s)
# print('PB_done')
# print(psi_complete[0:5])


psi_complete, nconc_complete,z = pb_2plate.pb_2plate(psi_complete,n_bulk,valency,sigma_f1,sigma_f2,domain,epsilon_s)
print('PB_done')
print(psi_complete[0:5])

q_complete = calculate.charge_density(nconc_complete,valency)

grandfe = energy_2plate.grandfe_pb_2plate(psi_complete,nconc_complete,n_bulk,valency,sigma_f1,sigma_f2,domain)
print(grandfe)

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
    file.create_dataset('psi_d', data = psi_complete*psi_c)
    file.create_dataset('nconc_d', data = nconc_complete*nconc_c/N_A)
    file.create_dataset('charge_d', data = q_complete*(nconc_c*ec))

    # Store all spatial profiles (non-dimensional)
    file.create_dataset('z', data = z)
    file.create_dataset('psi', data = psi_complete)
    file.create_dataset('nconc', data = nconc_complete)
    file.create_dataset('charge',data = q_complete)

    # Store free energy
    file.attrs['grandfe'] = grandfe # nondimensional
    file.attrs['grandfe_d'] = grandfe*(1/beta) # SI units


stop = timeit.default_timer()
print('Time: ', stop - start)



