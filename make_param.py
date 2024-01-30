import numpy as np
import os
import re
import subprocess
import sys

# Define your numpy array with 'domain' values
param_pattern = r'domain_d\s*=\s*([\d.]+)'
domain_values =  np.hstack((np.arange(5.0,13,2.0),np.arange(13.0,21.0,1.0),[22.0],np.arange(23.0,41.0,2.0),[44.0,48.0,56.0,60.0]))

# Input and output file paths
slurm_param_file =str(sys.argv[1])
param_dir = str(sys.argv[2])  # Create a directory to store the duplicated files
slurm_script = 'ddl31.sh'
slurm_simulator_file = 'simulator.py'
slurm_job_name = '#SBATCH --job-name=ocr31'
slurm_output_name = '#SBATCH --output=ocr31%j.out'


# Create the output directory if it doesn't exist
if not os.path.exists(param_dir):
    os.mkdir(param_dir)

# Read the contents of the param file
with open(slurm_param_file,'r') as f:
    param_contents = f.read()

# Read the contents of the job file
with open(slurm_script, 'r') as f:
    slurm_contents = f.read()

# Find the 'domain' parameter in the contents

matches = re.findall(param_pattern,param_contents)

if len(matches) == 1:
    domain_value = float(matches[0])
    for i,new_value in enumerate(domain_values):
        modified_contents = re.sub(param_pattern,f'domain_d = {new_value}',param_contents,count = 1)

        # Create a new file with the updated 'domain' value
        p_file = f'{param_dir}/domain_{i}.py'
        with open(p_file,'w') as f:
            f.write(modified_contents)

        slurm_job_file = f'{param_dir}/job_domain_{i}.sh'
        slurm_modified_contents = slurm_contents.replace(slurm_param_file,f'{param_dir}/domain_{i}.py')
        job_name = f'#SBATCH --job-name=dom3_{i}'
        slurm_modified_contents = slurm_modified_contents.replace(slurm_job_name,f'{job_name}')
        simulator_file = slurm_simulator_file
        slurm_modified_contents = slurm_modified_contents.replace(slurm_simulator_file,f'{simulator_file}')
        output_file= f'#SBATCH --output={param_dir}/dom3_{i}_%j.out'
        slurm_modified_contents = slurm_modified_contents.replace(slurm_output_name,f'{output_file}')
        with open(slurm_job_file,'w') as f:
            f.write(slurm_modified_contents)
        # Make SLURM job script executable
        subprocess.run(['chmod', '+x', slurm_job_file])

        # Execute the SLURM job script
        subprocess.run(['sbatch', slurm_job_file])


    print(f'{len(domain_values)} parameter scan files(' + slurm_param_file +   f') and SLURM job script files have been created and submitted in the "{param_dir}" directory.')





