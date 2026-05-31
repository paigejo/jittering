# Cluster driver: score 100 BYM2 then 100 SPDE simulations for all 8 models.
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

scoreSimStudyFullBoth(nsim = 100, regenerate = FALSE, useInla = "auto")

cat("\nDone.\n")
