# Cluster driver: score 100 SPDE-simulated sims for all 8 models.
# Same structure as runScoresFull_BYM2.R but for the SPDE generative model.
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

scoreSimStudyFull(model      = "spde",
                  nsim       = 100,
                  regenerate = FALSE,
                  useInla    = "auto")

cat("\nDone.\n")
