#!/usr/bin/env python
import os
from pathlib import PurePath, Path
from datetime import datetime
import csv
import itertools

# Helper function to make new output directories
def mkdir_p(dir):
  if not os.path.exists(dir):
    os.makedirs(dir)

# Get the project level directory
script_dir = PurePath(Path.cwd())
project_dir = str(script_dir.parent)

# Make the results directory
batch_date = datetime.now().strftime("%Y%m%d%H%M")
results_directory = f"{project_dir:s}/outputs/snapshots/{batch_date:s}"
mkdir_p(results_directory)

# Make directories to save job scripts and log sbatch output
job_directory = f"{project_dir:s}/job"
job_directory_snapshot = f"{project_dir:s}/job/{batch_date:s}"
log_directory = f"{project_dir:s}/log"
log_directory_snapshot = f"{project_dir:s}/log/{batch_date:s}"
mkdir_p(job_directory)
mkdir_p(job_directory_snapshot)
mkdir_p(log_directory)
mkdir_p(log_directory_snapshot)

# Simulation parameters
# dgp_array = [1, 2]
dgp_array = [1, 2, 3]
n_array = [500]
p_array = [10, 50]
snr_array = [0.25, 0.5, 1, 2]
ate_array = [0.5]
n_prop_array = [19]
num_simulations = 1000
on_remote = 1
job_num = 1
estimated_propensities_array = [0,1]
# residualize_xbcf_array = [0, 1]
residualize_xbcf_array = [0]
project_pi_yhat_array = [0]
n_yhat_array = [19]
use_yhat_covariate_array = [1]
grf_default_array = [0]

# Run the simulation loop - split some parameters across runs
dgp_outer_array = itertools.product(dgp_array, n_array, snr_array, ate_array, 
                                    n_prop_array, n_yhat_array, 
                                    estimated_propensities_array, 
                                    residualize_xbcf_array, project_pi_yhat_array, 
                                    use_yhat_covariate_array, grf_default_array)

for i, dgp_outer_params in enumerate(dgp_outer_array):
    # These are the "within run" parameters
    dgp_inner_array = itertools.product(p_array)
    
    # Define the file from which the SBATCH job will be executed
    job_prefix = f"sim_study_{job_num:d}"
    job_file_name = job_prefix + ".job"
    job_file = os.path.join(job_directory_snapshot, job_file_name)
    
    # Add SLURM switches to the job file
    with open(job_file, 'w+') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -N 1\n")
        fh.writelines("#SBATCH -c 20\n")
        fh.writelines("#SBATCH -t 0-04:00:00\n")
        fh.writelines("#SBATCH -p general\n")
        fh.writelines("#SBATCH -q public\n")
        fh.writelines("#SBATCH --mem=20G\n")
        fh.writelines(f"#SBATCH -o {log_directory_snapshot:s}/{job_prefix:s}.out\n")
        fh.writelines(f"#SBATCH -e {log_directory_snapshot:s}/{job_prefix:s}.err\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=asherren@asu.edu\n")
        fh.writelines("#SBATCH --export=None\n")
        fh.writelines("module purge\n")
        fh.writelines("module load r-4.2.2-gcc-11.2.0\n")
    
    for j, dgp_inner_params in enumerate(dgp_inner_array):
        # Simulation inputs unpacked from the outer array
        dgp_num = dgp_outer_params[0]
        n = dgp_outer_params[1]
        snr = dgp_outer_params[2]
        ate_true = dgp_outer_params[3]
        n_prop_submodels = dgp_outer_params[4]
        n_yhat_submodels = dgp_outer_params[5]
        estimated_propensities = dgp_outer_params[6]
        residualize_xbcf = dgp_outer_params[7]
        project_pi_yhat = dgp_outer_params[8]
        use_yhat_covariate = dgp_outer_params[9]
        grf_default = dgp_outer_params[10]
        
        # Simulation inputs unpacked from the inner array
        p = dgp_inner_params[0]
        
        # Which iteration of the inner array parameters
        script_iteration = j + 1
        
        with open(job_file, 'a+') as fh:
            # Run the simulation runner as many times as there are 
            # parameter combinations in the inner array
            fh.writelines(f"Rscript --verbose simulation_study.R {dgp_num:d} {n:d} {p:d} {snr:.2f} {num_simulations:d} {ate_true:.2f} {n_prop_submodels:d} {on_remote:d} {batch_date:s} {job_num:d} {estimated_propensities:d} {residualize_xbcf:d} {project_pi_yhat:d} {script_iteration:d} {n_yhat_submodels:d} {use_yhat_covariate:d} {grf_default:d}\n")

    # Execute the script
    os.system(f"sbatch --comment={job_file:s} {job_file:s}")

    # Increment the job number up by 1
    job_num += 1
