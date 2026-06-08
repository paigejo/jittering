# Rebuild stale DHS integration points after regenerating simPopsSurveys with
# doSmoothRisk=TRUE:
#   - BYM2 sims 1-9   : unchanged, skip
#   - BYM2 sims 10-100: surveys drifted (sim 10 hit a rejection-resample), rebuild
#   - SPDE sims 1-100 : all surveys changed, rebuild all
#
# Total: 91 BYM2 + 100 SPDE = 191 files.  Each build is ~30s, so ~95 min.

source("code/setup.R")

# ---- BYM2 ----
cat(sprintf("\n--- %s | BYM2 surveys (sims 10-100 only) ---\n", format(Sys.time())))
simEnv <- new.env()
simulateSurveys("bym2", nsim = 100, seed = 123, regenerate = FALSE, envir = simEnv)
surveysDHSb <- simEnv$surveysDHS

t0 <- proc.time()[3]
for(i in 10:100) {
    buildSimDhsIntPts(i, surveysDHSb, "bym2",
                      KDHSu = 16, KDHSr = 21,
                      regenerate = TRUE)
    if((i - 9) %% 10 == 0) {
        done    <- i - 9
        elapsed <- (proc.time()[3] - t0) / 60
        eta     <- elapsed * (91 - done) / done
        cat(sprintf("  BYM2 [%d/91 (sim %d)]  elapsed %.1f min  eta %.1f min\n",
                    done, i, elapsed, eta))
    }
}

# ---- SPDE ----
cat(sprintf("\n--- %s | SPDE surveys (all 1-100) ---\n", format(Sys.time())))
spdeEnv <- new.env()
simulateSurveys("spde", nsim = 100, seed = 123, regenerate = FALSE, envir = spdeEnv)
surveysDHSs <- spdeEnv$surveysDHS

t0 <- proc.time()[3]
for(i in 1:100) {
    buildSimDhsIntPts(i, surveysDHSs, "spde",
                      KDHSu = 16, KDHSr = 21,
                      regenerate = TRUE)
    if(i %% 10 == 0) {
        elapsed <- (proc.time()[3] - t0) / 60
        eta     <- elapsed * (100 - i) / i
        cat(sprintf("  SPDE [%d/100 (sim %d)]  elapsed %.1f min  eta %.1f min\n",
                    i, i, elapsed, eta))
    }
}

cat(sprintf("\n--- %s | done ---\n", format(Sys.time())))
