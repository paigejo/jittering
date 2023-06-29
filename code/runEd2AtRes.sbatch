#!/bin/bash
#SBATCH --job-name runEd2AtRes      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition CPUQ    # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 3-00:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=10000     # Memory limit for each tasks (in MB) 6850
#SBATCH -o runEd2AtRes_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e runEd2AtRes_%A_%a.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=john.paige@ntnu.no # Email to which notifications will be sent
#SBATCH --account=share-ie-imf # see https://www.hpc.ntnu.no/idun/getting-started-on-idun/accounting/ share-ie-imf


# --constraint="switch3&intel&pec6420" # use sinfo -o "%20N  %10c  %10m  %25f  %10G " to show features
# --nodelist=idun-03-0[1-9],idun-03-[10-12] # idun-03-0[1-9],idun-03-[10-12] see https://www.hpc.ntnu.no/idun/hardware/
# --account=ie-imf # see https://www.hpc.ntnu.no/idun/getting-started-on-idun/accounting/
# --exclusive=user # don't allow other users to share nodes with these jobs
# --constraint="switch3&intel&pec6420" # use sinfo -o "%20N  %10c  %10m  %25f  %10G " to show features

# for interactive job:
# https://www.hpc.ntnu.no/idun/getting-started-on-idun/running-jobs/
# srun --account=ie-imf --nodes=1 --partition=CPUQ --time=02:00:00 --pty bash 
# module purge
# module load R/4.0.3-foss-2020b
# module swap JasPer/2.0.14-GCCcore-10.2.0 JasPer/2.0.24-GCCcore-10.2.0
# module load GDAL/3.2.1-foss-2020b
module load R/4.2.1-foss-2022a
module load GDAL/3.5.0-foss-2022a

Rscript --verbose runEd2AtRes.R ${SLURM_ARRAY_TASK_ID} > runEd2AtRes_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout