# Parallel-pool scoring on the first 10 sims for both generative models.
#   - subarea scoring removed (speed)
#   - NDRAWS = 1000 (was 5000)
#   - LPT load-balanced callr workers (default 6) — M_DM_BYM2 fits dominate
# regenerate = TRUE so all 10 × 8 × 2 = 160 score files are rewritten uniformly
# with the new fields (scoresHyper, fitTime, scoreTime) and the smaller NDRAWS.
# BYM2 and SPDE generative models run SEQUENTIALLY (one at a time) so the
# combined memory load stays ≤ nWorkers R processes.

source("setup.R")
options(error = traceback)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/modM_DSep.R")
source("code/modM_MSep.R")
source("code/modM_DMSep.R")
source("code/modMdSep.R")
source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/scoreSimStudy.R")

scoreSimStudyParallel(model = "bym2", nsim = 10, nWorkers = 6,
                      NDRAWS = 1000, regenerate = TRUE)
scoreSimStudyParallel(model = "spde", nsim = 10, nWorkers = 6,
                      NDRAWS = 1000, regenerate = TRUE)

cat("\nDone.\n")
