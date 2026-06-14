# Test whether save/restore of last.par + last.par.best makes a REUSED TMB
# object give the correct Q at an axial theta after it's been "dirtied" by
# prior evaluations (the walk context). Ground truth = fresh rebuild object.
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
thetaX <- thetaC; thetaX["logit_phi"] <- thetaC["logit_phi"] + 3.0

# Ground truth: fresh rebuild obj at thetaX
oA <- .makeFixedHyperObj(dataList, paramsTemplate, thetaX, DLL, innerWarm = NULL)
invisible(oA$fn(numeric())); QA <- as.matrix(oA$env$spHess(oA$env$last.par.best, random = TRUE))

# Reused obj
objB <- MakeADFun(data = dataList, parameters = paramsTemplate,
                  random = c("alpha","beta","w_bym2Free","u_bym2Free"),
                  DLL = DLL, silent = TRUE)
on <- names(objB$par)
LP0  <- objB$env$last.par
LPB0 <- objB$env$last.par.best

Qax <- function() as.matrix(objB$env$spHess(objB$env$last.par, random = TRUE))

# Dirty the object with several evaluations across the parameter space
for(dz in c(0, 2, -2, 1.5, -3))
    invisible(objB$fn(as.numeric((thetaC + c(0, dz, 0))[on])))

# (1) NO reset: evaluate at axial, get Q
invisible(objB$fn(as.numeric(thetaX[on]))); Q_noreset <- Qax()

# (2) WITH reset of last.par + last.par.best to pre-walk values, then axial
objB$env$last.par      <- LP0
objB$env$last.par.best <- LPB0
invisible(objB$fn(as.numeric(thetaX[on]))); Q_reset <- Qax()

cat(sprintf("\nmax|Q_noreset - QA(fresh)| = %.3e\n", max(abs(Q_noreset - QA))))
cat(sprintf("max|Q_reset   - QA(fresh)| = %.3e\n", max(abs(Q_reset   - QA))))
cat("\n=== done ===\n")
