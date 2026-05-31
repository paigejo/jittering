# Score 8 sim-study models on 10 SPDE-simulated populations.
source("setup.R")
options(error = recover)

scoreSimStudy(model = "spde", nsim = 10, regenerate = FALSE)

cat("\nDone.\n")
