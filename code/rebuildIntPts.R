# Rebuild ALL 200 DHS integration points after the 2026-06-13 regeneration.
# The targetPopMat alignment fix changed the simulator's RNG path, so every
# survey in both worlds differs from the previous files (including BYM2 sims
# 1-9, which were byte-identical under the earlier doSmoothRisk-only change
# but are NOT anymore). Rebuild all 100 BYM2 + 100 SPDE = 200.
# Each build ~3 min, so ~10 hr total.

source("code/setup.R")

rebuildAll <- function(model) {
    cat(sprintf("\n--- %s | %s surveys (all 1-100) ---\n",
                format(Sys.time()), toupper(model)))
    env <- new.env()
    simulateSurveys(model, nsim = 100, seed = 123, regenerate = FALSE, envir = env)
    surveysDHS <- env$surveysDHS
    t0 <- proc.time()[3]
    for(i in 1:100) {
        buildSimDhsIntPts(i, surveysDHS, model,
                          KDHSu = 16, KDHSr = 21, regenerate = TRUE)
        if(i %% 10 == 0) {
            elapsed <- (proc.time()[3] - t0) / 60
            eta     <- elapsed * (100 - i) / i
            cat(sprintf("  %s [%d/100]  elapsed %.1f min  eta %.1f min\n",
                        toupper(model), i, elapsed, eta))
        }
    }
}

rebuildAll("bym2")
rebuildAll("spde")
cat(sprintf("\n--- %s | done ---\n", format(Sys.time())))
