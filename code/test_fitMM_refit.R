# Refit fitMM on cached BYM2 sim 1 (reuses simulated populations + DHS int pts).
# Useful for re-running just the fit after R-side patches (e.g. random=...).
source("setup.R")
options(error = recover)

res <- testFitMM(model = "bym2", nsim = 1, regenerate = FALSE)

cat("\nDone.\n")
