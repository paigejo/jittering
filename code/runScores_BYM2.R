# Score 8 sim-study models on 10 BYM2-simulated populations.
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

scoreSimStudy(model = "bym2", nsim = 10, regenerate = FALSE)

cat("\nDone.\n")
