#!/bin/bash
#SBATCH --job-name prepValidationFE
#SBATCH --partition CPUQ
#SBATCH --ntasks 1
#SBATCH --time 0-00:30
#SBATCH --mem-per-cpu=8000
#SBATCH -o prepValidationFE_%j.out
#SBATCH -e prepValidationFE_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=john.paige@ntnu.no
#SBATCH --account=share-ie-imf

# One-shot prep: builds fullMDMInputs_FE.RData and compiles the four TMB templates.
# Run this BEFORE submitting runValidationFE.sh, or chain them via:
#   PREP=$(sbatch --parsable prepValidationFE.sbatch)
#   sbatch --dependency=afterok:$PREP --array=1-80 runValidationFE.sbatch

module load R/4.2.1-foss-2022a
module load GDAL/3.5.0-foss-2022a

Rscript --verbose prepValidationFE.R > prepValidationFE_${SLURM_JOB_ID}.Rout
