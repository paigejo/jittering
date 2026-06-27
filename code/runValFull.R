# Driver: run the validation in parallel via callr (no SLURM), mirroring the sim
# study's runScoresFull_*.R -> scoreSimStudyFull.R pattern. From the repo root:
#   Rscript code/runValFull.R          # FE family (cluster + areal)
#   Rscript code/runValFull.R BYM2     # BYM2 family (areal; cluster pending)
# Worker count auto-detects: local 4 (BYM2 2), cluster 16.
source("code/runValidationFull.R")
.args  <- commandArgs(trailingOnly = TRUE)
family <- if(length(.args) >= 1) .args[1] else "FE"
runValidationFull(family = family, nsim = 10000, arealNsim = 2000)
