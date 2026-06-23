# Cluster driver: score 100 BYM2 then 100 SPDE simulations for all 8 models.
source("setup.R")
options(error = recover)

scoreSimStudyFullBoth(nsim = 100, regenerate = TRUE, useInla = "auto")

cat("\nDone.\n")
