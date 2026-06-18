# Localize Md's "non-conformable arguments" inside predGrid.
# No tryCatch: options(error=traceback) prints the full call stack on death.

source("code/setup.R")
options(error = quote({ traceback(2); q(status = 1) }))

SIMIDX <- 17
simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)

cat(sprintf("\n=== %s | fitting Md ===\n", format(Sys.time())))
res <- .fitOne("Md", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

cat("\n--- structure of the fit ---\n")
cat("par.fixed names: ", paste(names(res$TMBsd$par.fixed), collapse=", "), "\n")
cat("par.random name counts:\n")
print(table(names(res$TMBsd$par.random)))
cat("jointPrecision: ", ifelse(is.null(res$TMBsd$jointPrecision), "NULL",
                               paste(dim(res$TMBsd$jointPrecision), collapse=" x ")), "\n")

cat("\n--- posterior draw rownames (via posteriorDraws, what predGrid consumes) ---\n")
draws <- posteriorDraws(res, NDRAWS = 50, useInla = "auto")
print(table(rownames(draws)))

cat(sprintf("\n=== %s | predGrid (no tryCatch; traceback on error) ===\n",
            format(Sys.time())))
grid <- predGrid(res$TMBsd, popMat = popMatNGAThresh, nsim = 1000,
                 obj = res$TMBobj, admLevel = "stratMICS", res = res)
cat("predGrid SUCCEEDED. names(grid): ", paste(names(grid), collapse=", "), "\n")
cat(sprintf("pixel preds mean = %.4f\n", mean(grid$preds, na.rm = TRUE)))

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
