import mgrf_2plate
from numerical_param import *
import energy_2plate
start = timeit.default_timer()

#sys.path.append('/batch_files')
# Argument parser to accept the input files
parser = argparse.ArgumentParser(description='Code to calculate EDL structure using MGRF Theory with mean-field PB as an initial guess')
parser.add_argument('input_files', nargs='+', help='Paths to the input files for physical parameters')
args = parser.parse_args()

folder_path = os.path.dirname(args.input_files[0])
sys.path.insert(0, folder_path)

# Load the physical input configuration from the first file in the list
module_name = os.path.splitext(os.path.basename(args.input_files[0]))[0]
input_physical = importlib.import_module(module_name)
variables = {name: value for name, value in input_physical.__dict__.items() if not name.startswith('__')}
(locals().update(variables))



if cb2_d != 0:
    output_dir =os.getcwd() + '/results-mixture' + str(abs(valency[0]))+ '_' + str(abs(valency[1])) + '_' + str(abs(valency[2]))+ '_' + str(abs(valency[3]))
    file_name = str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5)) + '_' + str(round(float(domain_in_d), 2)) + '_' + str(round(rad_ions_d[2] / pow(10, -10), 2)) + '_' + str(round(sigma_in1_d, 5)) + '_' + str(round(sigma_in2_d, 5))
else:
    output_dir = os.getcwd() + '/results' + str(abs(valency[0]))+ '_' + str(abs(valency[1]))
    file_name = str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5))  + '_' + str(round(float(domain_in_d), 2)) + '_' + str(round(rad_ions_d[0] / pow(10, -10), 2)) + '_' + str(round(rad_ions_d[1] / pow(10, -10), 2)) + '_' + str(round(sigma_in1_d, 5)) + '_' + str(round(sigma_in2_d, 5))

with h5py.File(output_dir + '/mgrf_' + file_name + '.h5', 'r') as file:
    # Retrieve psi and nconc
    psi_complete = np.array(file['psi'])
    nconc_complete = np.array(file['nconc'])
    uself_complete = np.array(file['uself'])

psi_complete,nconc_complete,uself_complete, q_complete, z= mgrf_2plate.mgrf_2plate(psi_complete,nconc_complete,n_bulk,valency,rad_ions,vol_ions,vol_sol,sigma_f1,sigma_f2,domain,epsilon_s)


grandfe = energy_2plate.grandfe_mgrf_2plate(psi_complete,nconc_complete,uself_complete,n_bulk,valency,rad_ions,vol_ions,vol_sol,sigma_f1,sigma_f2,domain,epsilon_s)
print(grandfe)

stop = timeit.default_timer()

print('Time: ', stop - start)
if cb2_d != 0:
    file_name = str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5)) + '_' + str(round(float(domain_d), 2)) + '_' + str(round(rad_ions_d[2] / pow(10, -10), 2)) + '_' + str(round(sigma_f1_d, 5)) + '_' + str(round(sigma_f2_d, 5))
else:
    file_name = str(round(cb1_d / pow(10, 3), 9)) + '_' + str(round(cb2_d / pow(10, 3), 5))  + '_' + str(round(float(domain_d), 2)) + '_' + str(round(rad_ions_d[0] / pow(10, -10), 2)) + '_' + str(round(rad_ions_d[1] / pow(10, -10), 2)) + '_' + str(round(sigma_f1_d, 5)) + '_' + str(round(sigma_f2_d, 5))

# # Writing everything in SI units
# with h5py.File(output_dir + '/mgrf_' + file_name + '.h5', 'w') as file:

#     # Storing scalar variables as attributes of the root group
#     file.attrs['ec_charge'] = ec
#     file.attrs['char_length'] = l_b
#     file.attrs['beta'] = beta
#     file.attrs['epsilon_s'] = epsilonr_s_d
#     file.attrs['epsilon_p'] = epsilonr_s_d
#     file.attrs['cb1'] = cb1_d * 0.001
#     file.attrs['cb2'] = cb2_d*0.001
#     file.attrs['domain'] = domain_d
    
#     # Storing numerical parameters as attributes of the root group
#     file.attrs['s_conv'] = s_conv
#     file.attrs['N_grid'] = N_grid
#     file.attrs['quads'] = quads
#     file.attrs['grandfe_quads'] = grandfe_quads
#     file.attrs['dealias'] = dealias
#     file.attrs['ncc_cutoff_pb'] = ncc_cutoff_pb
#     file.attrs['ncc_cutoff_mgrf'] = ncc_cutoff_mgrf
#     file.attrs['num_ratio'] = num_ratio
#     file.attrs['selfe_ratio'] = selfe_ratio
#     file.attrs['eta_ratio'] = eta_ratio
#     file.attrs['tolerance'] = tolerance
#     file.attrs['tolerance_pb'] = tolerance_pb
#     file.attrs['tolerance_num'] = tolerance_num
#     file.attrs['tolerance_greens'] = tolerance_greens

#     # Storing parameter arrays
#     file.create_dataset('valency', data = valency)
#     file.create_dataset('radii', data = rad_ions_d)
#     file.create_dataset('volumes', data = np.concatenate((vol_ions_d,[vol_sol_d])))
#     file.create_dataset('surface_charges', data = np.array([sigma_f1_d,sigma_f2_d]))

#     # Store all spatial profiles  (SI units)
#     file.create_dataset('z_d', data = z*l_c)
#     file.create_dataset('psi_d', data = psi_complete*psi_c)
#     file.create_dataset('nconc_d', data = nconc_complete*nconc_c/N_A)
#     file.create_dataset('uself_d', data = uself_complete*(1/beta))
#     file.create_dataset('charge_d', data = q_complete*(nconc_c*ec))

#     # Store all spatial profiles (non-dimensional)
#     file.create_dataset('z', data = z)
#     file.create_dataset('psi', data = psi_complete)
#     file.create_dataset('nconc', data = nconc_complete)
#     file.create_dataset('uself', data = uself_complete)
#     file.create_dataset('charge',data = q_complete)

#     # Store free energy
#     file.attrs['grandfe'] = grandfe # nondimensional
#     file.attrs['grandfe_d'] = grandfe*(1/beta) # SI units



