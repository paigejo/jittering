# Cluster driver: score 100 BYM2-simulated sims for all 8 models.
# Platform-aware worker counts (32 FE / 8 BYM2 on Linux; 8 / 1 on Mac/Win).
# Per-fit posterior draws default to useInla = "auto" (Gauss with INLA-style
# fallback for non-PSD jointPrecision). Set regenerate = TRUE to overwrite.
source("setup.R")
options(error = recover)

scoreSimStudyFull(model = "bym2", nsim = 100, regenerate = FALSE, useInla = "auto")

cat("\nDone.\n")
