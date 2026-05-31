# Generate 100 sims + K=16/21 DHS integration points for both generative models.
#
# Default: BYM2 and SPDE run SEQUENTIALLY in this single R process (safe for
# memory-constrained machines — running them in parallel crashed the workstation).
# To run them concurrently as two child R processes instead, set parallel = TRUE
# below (requires `callr`).
#
# simData1 / simData1BYM2 checkpoint after every iteration, so this script is
# idempotent and crash-safe: re-running picks up wherever the previous attempt
# left off (matching seed=123).

source("setup.R")
options(error = traceback)
source("code/simData.R")
source("code/makeIntegrationPoints.R")
source("code/testInfrastructure.R")

prepareSimsForModels(models = c("bym2", "spde"),
                     nsim = 100, seed = 123,
                     KDHSu = 16, KDHSr = 21,
                     regenerate = FALSE,
                     parallel = TRUE)

cat("\nDone.\n")
