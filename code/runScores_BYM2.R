# Score 8 sim-study models on 10 BYM2-simulated populations.
source("setup.R")
options(error = recover)

scoreSimStudy(model = "bym2", nsim = 10, regenerate = FALSE)

cat("\nDone.\n")
