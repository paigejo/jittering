# Localize the reuse-vs-rebuild difference to mu_k and/or Q_k at a single
# CCD point. At the SAME theta, compute the inner mode and precision two ways:
#   (A) rebuild: .makeFixedHyperObj (hypers mapped to constants), fn(numeric())
#   (B) reuse:   shared free-hyper walkObj, fn(theta)
# Compare mu (per group) and Q. If these match, the divergence is downstream;
# if they differ, this is the structural cause and we see exactly where.

source("code/setup.R")
options(warn = 1)

SIMIDX <- 17
simEnv <- new.env(); simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env(); load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2", KMICS = 100, KDHSu = 16, KDHSr = 21)
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

SD <- res$TMBsd
hyperNames <- c("log_tau","logit_phi","log_tauEps")
theta <- SD$par.fixed[hyperNames]      # centre theta = hyper modes
cat("theta (centre):\n"); print(theta)

obj  <- res$TMBobj
dataList <- obj$env$data; paramsTemplate <- obj$env$parameters; DLL <- obj$env$DLL

# (A) rebuild mapped obj
objA <- .makeFixedHyperObj(dataList, paramsTemplate, theta, DLL, innerWarm = NULL)
nllA <- as.numeric(objA$fn(numeric()))
muA  <- objA$env$last.par.best
QA   <- objA$env$spHess(muA, random = TRUE)

# (B) reuse free-hyper obj
objB <- MakeADFun(data = dataList, parameters = paramsTemplate,
                  random = c("alpha","beta","w_bym2Free","u_bym2Free"),
                  DLL = DLL, silent = TRUE)
nllB <- as.numeric(objB$fn(as.numeric(theta[names(objB$par)])))
fullB <- objB$env$last.par
muB  <- fullB[objB$env$random]
QB   <- objB$env$spHess(fullB, random = TRUE)

cat(sprintf("\nNLL: rebuild=%.6f  reuse=%.6f  diff=%.3e\n", nllA, nllB, nllA - nllB))
cat(sprintf("length(muA)=%d  length(muB)=%d\n", length(muA), length(muB)))
cat(sprintf("names match: %s\n", identical(names(muA), names(muB))))

cat("\n-- mu per-group max abs diff (rebuild vs reuse) --\n")
for(nm in unique(names(muA))) {
    a <- muA[names(muA) == nm]; b <- muB[names(muB) == nm]
    cat(sprintf("  %-12s n=%-3d  max|diff|=%.3e\n", nm, length(a), max(abs(a - b))))
}

QA <- as.matrix(QA); QB <- as.matrix(QB)
cat(sprintf("\nQ dims: A %s  B %s\n", paste(dim(QA), collapse="x"), paste(dim(QB), collapse="x")))
if(all(dim(QA) == dim(QB)))
    cat(sprintf("max|QA - QB| = %.3e   (relative to max|QA|=%.3e: %.3e)\n",
                max(abs(QA - QB)), max(abs(QA)), max(abs(QA - QB))/max(abs(QA))))

cat(sprintf("\n=== done ===\n"))
