# Cluster driver: score 100 SPDE-simulated sims for all 8 models.
# Same structure as runScoresFull_BYM2.R but for the SPDE generative model.
source("setup.R")
options(error = recover)

scoreSimStudyFull(model = "spde", nsim = 100, regenerate = FALSE, useInla = "auto")

cat("\nDone.\n")
