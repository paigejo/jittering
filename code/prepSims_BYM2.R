# Generate 10 BYM2-simulated populations + K=16/21 DHS integration points.
source("setup.R")
options(error = recover)
source("code/simData.R")
source("code/makeIntegrationPoints.R")
source("code/testInfrastructure.R")

prepareSims(model = "bym2", nsim = 10, regenerate = TRUE)

cat("\nDone.\n")
