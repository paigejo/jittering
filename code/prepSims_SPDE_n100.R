# Extend SPDE sims from nsim=10 to nsim=100:
#   1. regenerate simPopsSurveys_SPDE.RData with nsim=100 (seed=123, same as before
#      so the first 10 sims are byte-identical to the existing file)
#   2. build K=16/21 DHS int pts for any missing sim — sims 1-10 stay cached.
source("setup.R")
options(error = recover)

env <- new.env()
simulateSurveys("spde", nsim = 100, seed = 123, regenerate = TRUE, envir = env)
surveysDHS <- env$surveysDHS

for(i in 1:100)
    buildSimDhsIntPts(i, surveysDHS, "spde", KDHSu = 16, KDHSr = 21,
                      regenerate = FALSE)

cat("\nDone.\n")
