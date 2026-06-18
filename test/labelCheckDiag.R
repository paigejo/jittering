# Cheap diagnostic: fit Md_FE (no random effects → fast) on sim 1 BYM2, then
# run .checkDrawLabels() to see whether sd(draws[name,]) agrees with
# summary(SD)[name, "Std. Error"] per parameter name. If they do, both TMB
# code paths agree on what variable lives at what name — labels canonical.
# If not, the warning() in .checkDrawLabels prints the offending names.

source("code/setup.R")
options(warn = 1)  # print warnings as they happen, not at end

cat(sprintf("\n=== %s | label-check on Md_FE sim1 ===\n", format(Sys.time())))

simEnv <- new.env()
simulateSurveys("bym2", nsim = 1, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

inp <- buildInputsForSim(simIdx = 1,
                         surveysDHS  = simEnv$surveysDHS,
                         surveysMICS = simEnv$surveysMICS,
                         intPtsMICS  = micsEnv$intPtsMICS,
                         model       = "bym2",
                         KMICS = 100, KDHSu = 16, KDHSr = 21)

cat(sprintf("\n=== %s | fitting Md_FE ===\n", format(Sys.time())))
t0 <- proc.time()[3]
res <- .fitOne("Md_FE", inp,
               KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS  = c("urban","access","elev","distRiversLakes","normPop"))
cat(sprintf("[fit] %.1f sec\n", proc.time()[3] - t0))

SD <- res$TMBsd
cat("\n=== sdreport overview ===\n")
cat("class(SD): ", class(SD), "\n")
cat("names(par.fixed):  ", paste(names(SD$par.fixed),  collapse=", "), "\n")
cat("names(par.random): ", paste(names(SD$par.random), collapse=", "), "\n")
cat("dim(jointPrecision): ",
    if(!is.null(SD$jointPrecision)) paste(dim(SD$jointPrecision), collapse=" x ") else "NULL",
    "\n")
cat("has rownames(jointPrecision)? ",
    !is.null(SD$jointPrecision) && !is.null(rownames(SD$jointPrecision)), "\n")
if(!is.null(SD$jointPrecision) && !is.null(rownames(SD$jointPrecision))) {
    jpN <- rownames(SD$jointPrecision)
    cnN <- names(c(SD$par.fixed, SD$par.random))
    cat("len(jointPrecision rownames): ", length(jpN), "\n")
    cat("len(c(par.fixed, par.random)): ", length(cnN), "\n")
    if(length(jpN) == length(cnN)) {
        cat("identical(jpRownames, c(par.fixed, par.random)-names)?  ",
            identical(jpN, cnN), "\n")
        firstDiff <- which(jpN != cnN)[1]
        if(!is.na(firstDiff))
            cat("first mismatch at idx", firstDiff,
                ":  jp=", jpN[firstDiff], "  canonical=", cnN[firstDiff], "\n")
    }
}

cat(sprintf("\n=== %s | sampling draws ===\n", format(Sys.time())))
draws <- posteriorDraws(res, NDRAWS = 1000, useInla = "auto")

cat("\n=== per-name SD comparison ===\n")
tab <- .checkDrawLabels(draws, res, relTol = 8)
if(!is.null(tab)) print(round(tab[order(abs(tab$relDiff), decreasing = TRUE), ], 5))

cat(sprintf("\n=== %s | DONE ===\n", format(Sys.time())))
