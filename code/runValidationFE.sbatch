#!/bin/bash
#SBATCH --job-name runValidationFE
#SBATCH --partition CPUQ
#SBATCH --ntasks 1
#SBATCH --time 0-12:00
#SBATCH --mem-per-cpu=8000
#SBATCH -o runValidationFE_%A_%a.out
#SBATCH -e runValidationFE_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=john.paige@ntnu.no
#SBATCH --account=share-ie-imf

# Job array indexing:
#   i = 1..4  : models  M_D, M_M, M_DM, M_DM_comm
#   j = 1..20 : folds   1..10 = DHS held-out folds; 11..20 = MICS held-out folds
#   index    = (i-1)*20 + j   (decoded by getJobIndices in code/validationCommon.R)
# Total array size = 4 * 20 = 80; jobs whose (model, side) combination is invalid
# (M_D + MICS fold, M_M + DHS fold) exit cleanly without doing work.

module load R/4.2.1-foss-2022a
module load GDAL/3.5.0-foss-2022a

Rscript --verbose runValidationFE.R ${SLURM_ARRAY_TASK_ID} > runValidationFE_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout
