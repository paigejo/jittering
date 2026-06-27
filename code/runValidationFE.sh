#!/bin/bash
# 73 cluster jobs (20-fold CV; collapsed to 11 jobs for single-dataset models) +
# 5 areal = 78 array tasks; cap concurrency at 16 cores (shared cluster).
sbatch --array=1-78%16 runValidationFE.sbatch
