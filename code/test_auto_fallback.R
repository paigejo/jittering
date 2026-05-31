# Quick smoke test of the new posteriorDraws auto-fallback.
# Fit M_D_BYM2 sim 1 (known degenerate) and verify:
#   1. useInla = "auto"  --> Gauss attempted, falls back to INLA on Cholesky fail
#   2. useInla = FALSE   --> Gauss attempted, error if Cholesky fails
#   3. useInla = TRUE    --> INLA forced regardless
source("setup.R")
options(error = traceback)
source("code/modM_DSep.R")
source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/inlaStyleDraws.R")

simIdx <- 1
simEnv <- new.env()
simulateSurveys("bym2", nsim = simIdx, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
inp <- buildInputsForSim(simIdx, simEnv$surveysDHS, simEnv$surveysMICS,
                         micsEnv$intPtsMICS, "bym2",
                         KMICS = 100, KDHSu = 16, KDHSr = 21)

cat("\n=== Fitting M_D_BYM2 sim 1 ===\n")
t0 <- proc.time()[3]
res <- fitMD(datDHS = inp$datDHS, datMICS = inp$datMICS,
             inputsMDM = inp$inputsMDM,
             KMICS = 100, KDHSurb = 16, KDHSrur = 21,
             Qgh = 10, getSDs = TRUE, verbose = FALSE)
cat(sprintf("fitMD time: %.1f min  pdHess=%s\n",
            (proc.time()[3]-t0)/60,
            as.character(res$TMBsd$pdHess)))

# --- Mode 1: auto ---------------------------------------------------------
cat("\n=== posteriorDraws(useInla = \"auto\") ===\n")
t0 <- proc.time()[3]
d_auto <- tryCatch(posteriorDraws(res, NDRAWS = 200, useInla = "auto"),
                   error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  time: %.1f min  rows=%s  cols=%s\n",
            (proc.time()[3]-t0)/60,
            if(is.null(d_auto)) "FAIL" else nrow(d_auto),
            if(is.null(d_auto)) "FAIL" else ncol(d_auto)))

# --- Mode 2: explicit FALSE (force Gauss) --------------------------------
cat("\n=== posteriorDraws(useInla = FALSE) ===\n")
t0 <- proc.time()[3]
d_g <- tryCatch(posteriorDraws(res, NDRAWS = 200, useInla = FALSE),
                error = function(e) { cat("EXPECTED ERROR (Gauss can't handle non-PSD):",
                                          conditionMessage(e), "\n"); NULL })
cat(sprintf("  time: %.1f min  result: %s\n",
            (proc.time()[3]-t0)/60,
            if(is.null(d_g)) "errored as expected" else "succeeded"))

# --- Mode 3: explicit TRUE (force INLA) ----------------------------------
cat("\n=== posteriorDraws(useInla = TRUE) ===\n")
t0 <- proc.time()[3]
d_i <- tryCatch(posteriorDraws(res, NDRAWS = 200, useInla = TRUE),
                error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  time: %.1f min  rows=%s  cols=%s\n",
            (proc.time()[3]-t0)/60,
            if(is.null(d_i)) "FAIL" else nrow(d_i),
            if(is.null(d_i)) "FAIL" else ncol(d_i)))

cat("\nDone.\n")
