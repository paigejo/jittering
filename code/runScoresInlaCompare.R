# Drive scoreInlaCompare for both generative models (BYM2 then SPDE),
# 4 workers each, regenerate=TRUE (so results are written from scratch
# with the fixed Pt-step Gaussian path AND the new INLA-style path).
source("setup.R")
options(error = traceback)
source("code/modFED.R");   source("code/modFEM.R")
source("code/modFEMD.R");  source("code/modM_DSep.R")
source("code/modM_MSep.R");source("code/modM_DMSep.R")
source("code/modMdSep.R"); source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/scoreSimStudy.R")
source("code/inlaStyleDraws.R")
source("code/scoreInlaCompare.R")

# nWorkers = 1 (serial). BYM2 fits + INLA-style CCD walks have per-process
# memory peaks > 12 GB; even 2 parallel workers crashed the workstation.
# BYM2 generative is already complete on disk; only SPDE remains.
# regenerate = FALSE skips the 2 SPDE files already saved.
scoreInlaCompare(model = "spde", nsim = 10, nWorkers = 1,
                 NDRAWS = 1000, regenerate = FALSE)

cat("\nDone.\n")
