#!/bin/bash
#SBATCH --job-name runValidationFE
#SBATCH --partition CPUQ
#SBATCH --ntasks 1
#SBATCH --time 0-12:00
#SBATCH --mem-per-cpu=16000
#SBATCH -o runValidationFE_%A_%a.out
#SBATCH -e runValidationFE_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=john.paige@ntnu.no
#SBATCH --account=share-ie-imf

# Job array indexing (decodeValJob over .buildValJobList; models = VALMODELS:
#   Md, M_D, M_dm, M_DM, M_M). 20-fold CV (10 DHS + 10 MICS held-out folds) for
#   every model, collapsed for single-dataset models (the unused survey -> one
#   "predict-all" job, since its folds share a fit). Each cluster job fits once
#   and predicts ONE survey's held-out clusters (incl. cross-survey for the
#   single-dataset models). Per model: Md/M_D/M_M = 11 cluster jobs, M_dm/M_DM =
#   20 -> 73 cluster jobs; + 5 areal (leave-out folds 6:10) = 78 total.
# Concurrency capped at 16 (see .sh). Areal jobs run predGrid -> higher mem.
# Per-survey + denominator-weighted combined scores are assembled at collation
# time (scoreValidationFE.R -> .scoreModelCluster).

module load R/4.2.1-foss-2022a
module load GDAL/3.5.0-foss-2022a

Rscript --verbose runValidationFE.R ${SLURM_ARRAY_TASK_ID} > runValidationFE_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout
