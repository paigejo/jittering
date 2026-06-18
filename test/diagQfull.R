# ONE capture of everything needed to pin the Q divergence. At a CENTRE theta
# and an AXIAL theta, for BOTH object constructions, dump: last.par outer
# block (is it actually set to theta?), inner mode mu, and the full Q matrix.
# Also test whether forcing the hyper block of the par passed to spHess fixes
# the reuse Q.
source("code/setup.R"); options(warn = 1)
SIMIDX <- 17
simEnv <- new.env(); simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env(); load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2", KMICS = 100, KDHSu = 16, KDHSr = 21)
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))
SD <- res$TMBsd; obj <- res$TMBobj
dataList <- obj$env$data; paramsTemplate <- obj$env$parameters; DLL <- obj$env$DLL
hyperNames <- c("log_tau","logit_phi","log_tauEps")
thetaC <- SD$par.fixed[hyperNames]

objB <- MakeADFun(data = dataList, parameters = paramsTemplate,
                  random = c("alpha","beta","w_bym2Free","u_bym2Free"),
                  DLL = DLL, silent = TRUE)
randIdx <- objB$env$random
outerNames <- names(objB$par)
cat("objB$par (outer) names:\n"); print(outerNames)

check <- function(label, theta) {
    cat(sprintf("\n############ %s ############\n", label))
    # mapped/rebuild
    oA <- .makeFixedHyperObj(dataList, paramsTemplate, theta, DLL, innerWarm = NULL)
    invisible(oA$fn(numeric())); muA <- oA$env$last.par.best
    QA <- as.matrix(oA$env$spHess(muA, random = TRUE))
    # free/reuse
    thB <- as.numeric(theta[outerNames]); invisible(objB$fn(thB))
    lp <- objB$env$last.par
    cat("last.par outer block (should equal theta):\n")
    print(round(lp[setdiff(seq_along(lp), randIdx)], 6))
    cat("theta:\n"); print(round(as.numeric(theta), 6))
    muB <- lp[randIdx]
    QB <- as.matrix(objB$env$spHess(lp, random = TRUE))
    # free, but FORCE hypers in the par passed to spHess
    lpForce <- lp; lpForce[setdiff(seq_along(lp), randIdx)] <- thB
    QBf <- as.matrix(objB$env$spHess(lpForce, random = TRUE))
    cat(sprintf("max|muA - muB|           = %.3e\n", max(abs(muA - muB))))
    cat(sprintf("max|QA - QB(last.par)|    = %.3e\n", max(abs(QA - QB))))
    cat(sprintf("max|QA - QB(forced hyp)|  = %.3e\n", max(abs(QA - QBf))))
}
check("CENTRE", thetaC)
thetaX <- thetaC; thetaX["logit_phi"] <- thetaC["logit_phi"] + 3.0
check("AXIAL (logit_phi +3)", thetaX)
cat("\n=== done ===\n")
