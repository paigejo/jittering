# Confirm the spHess aliasing hypothesis and verify a copy fixes it.
# If spHess reuses one internal value buffer, then a reference stored from an
# early call will CHANGE after a later call (aliasing), while an explicit copy
# stays put.
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

objB <- MakeADFun(data = dataList, parameters = paramsTemplate,
                  random = c("alpha","beta","w_bym2Free","u_bym2Free"),
                  DLL = DLL, silent = TRUE)
on <- names(objB$par)

# Point 1 (centre): store a REFERENCE and a COPY (two copy strategies)
invisible(objB$fn(as.numeric(thetaC[on])))
Q1_ref  <- objB$env$spHess(objB$env$last.par, random = TRUE)   # no copy
Q1_cpA  <- Q1_ref + 0                                          # arithmetic copy
Q1_cpB  <- as.matrix(objB$env$spHess(objB$env$last.par, random = TRUE))  # dense copy
ref_sum_before <- sum(abs(Q1_ref@x)); cpA_before <- sum(abs(Q1_cpA@x)); cpB_before <- sum(abs(Q1_cpB))

# Point 2 (axial): another spHess call on the SAME object
invisible(objB$fn(as.numeric(thetaX[on])))
Q2 <- objB$env$spHess(objB$env$last.par, random = TRUE)
q2_sum <- sum(abs(Q2@x))

ref_sum_after <- sum(abs(Q1_ref@x)); cpA_after <- sum(abs(Q1_cpA@x)); cpB_after <- sum(abs(Q1_cpB))

cat(sprintf("\nQ1_ref  sum |x|: before=%.4f  after point2=%.4f   %s\n",
            ref_sum_before, ref_sum_after,
            ifelse(abs(ref_sum_before-ref_sum_after) > 1e-6, "CHANGED -> ALIASED", "stable")))
cat(sprintf("Q1_cpA  sum |x|: before=%.4f  after point2=%.4f   %s   (arithmetic copy Q+0)\n",
            cpA_before, cpA_after,
            ifelse(abs(cpA_before-cpA_after) > 1e-6, "CHANGED", "stable -> copy OK")))
cat(sprintf("Q1_cpB  sum |x|: before=%.4f  after point2=%.4f   %s   (as.matrix copy)\n",
            cpB_before, cpB_after,
            ifelse(abs(cpB_before-cpB_after) > 1e-6, "CHANGED", "stable -> copy OK")))
cat(sprintf("point2 Q sum |x| = %.4f\n", q2_sum))
cat(sprintf("\nIf ALIASED: Q1_ref after (%.4f) == point2 (%.4f)? %s\n",
            ref_sum_after, q2_sum, abs(ref_sum_after - q2_sum) < 1e-6))
cat("\n=== done ===\n")
