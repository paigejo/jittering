# Driver: run the validation in parallel via callr (no SLURM), mirroring the sim
# study's runScoresFull_*.R -> scoreSimStudyFull.R pattern. From the repo root:
#   Rscript code/runValFull.R                 # FE family (cluster + areal)
#   Rscript code/runValFull.R BYM2            # BYM2 family
#   Rscript code/runValFull.R BYM2 regenerate # force re-fit/predict (ignore cached
#                                             # predsJob/arealPred); needed to apply
#                                             # PREDICTION-side fixes (e.g. Md cluster
#                                             # w_bym2Star). Scoring-only fixes apply
#                                             # on a normal (cached) re-run.
# Worker count auto-detects: local 4 (BYM2 2), cluster 16.
source("code/runValidationFull.R")
.args  <- commandArgs(trailingOnly = TRUE)
family <- if(length(.args) >= 1) .args[1] else "FE"
regen  <- any(tolower(.args) %in% c("regenerate", "regen", "true", "regenerate=true"))
cat(sprintf("[runValFull] family=%s  regenerate=%s\n", family, regen))
runValidationFull(family = family, nsim = 10000, arealNsim = 2000, regenerate = regen)
