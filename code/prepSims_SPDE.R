# Generate 10 SPDE-simulated populations + K=16/21 DHS integration points.
source("setup.R")
options(error = traceback)
source("code/simData.R")
source("code/makeIntegrationPoints.R")
source("code/testInfrastructure.R")

prepareSims(model = "spde", nsim = 10, regenerate = TRUE)

cat("\nDone.\n")
