# Cluster driver: score 100 BYM2-simulated sims for all 8 models.
# Platform-aware worker counts (32 FE / 8 BYM2 on Linux; 8 / 1 on Mac/Win).
# Per-fit posterior draws default to useInla = "auto" (Gauss with INLA-style
# fallback for non-PSD jointPrecision). Set regenerate = TRUE to overwrite.
source("setup.R")
options(error = traceback)
source("code/modFED.R");   source("code/modFEM.R")
source("code/modFEMD.R");  source("code/modM_DSep.R")
source("code/modM_MSep.R");source("code/modM_DMSep.R")
source("code/modMdSep.R"); source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/inlaStyleDraws.R")
source("code/scoreSimStudy.R")
source("code/scoreSimStudyFull.R")

scoreSimStudyFull(model      = "bym2",
                  nsim       = 100,
                  regenerate = FALSE,
                  useInla    = "auto")

cat("\nDone.\n")
