# BYM2-from-scratch reproducibility test:
#   - Generate nsim BYM2-simulated populations + surveys (regenerate=TRUE)
#   - Build K=16/21 DHS integration points for each sim
#   - Fit fitMM (BYM2 MICS-only) on sim 1, report MLE vs truth + NLL profile in phi
# All heavy lifting lives in testInfrastructure.R::testFitMM.
source("setup.R")
options(error = recover)

res <- testFitMM(model = "bym2", nsim = 1, regenerate = TRUE)

cat("\nDone.\n")
