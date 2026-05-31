# Drive scoreInlaCompare for BYM2 then SPDE generative model.
# nWorkers = 1 (serial) — BYM2 fits + INLA-style CCD walks have per-process
# memory peaks > 12 GB; even 2 parallel workers crashed the workstation.
source("setup.R")
options(error = recover)

scoreInlaCompare(model = "bym2", nsim = 10, nWorkers = 1,
                 NDRAWS = 1000, regenerate = FALSE)
scoreInlaCompare(model = "spde", nsim = 10, nWorkers = 1,
                 NDRAWS = 1000, regenerate = FALSE)

cat("\nDone.\n")
