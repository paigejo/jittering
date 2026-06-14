# Diff per-CCD-point Q diagonals AND the full closest vector between reuse
# and rebuild. mu_k already proven identical (5e-14); this checks the two
# remaining deterministic inputs to the inner draws.
source("code/setup.R"); options(warn = 1)
SEED <- 20260614; SIMIDX <- 17
simEnv <- new.env(); simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env(); load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2", KMICS = 100, KDHSu = 16, KDHSr = 21)
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))
set.seed(SEED); dn <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)
source("code/inlaStyleDraws_OLD.R")
set.seed(SEED); do <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)

qN <- attr(dn, "ccd_Qdiag"); qO <- attr(do, "ccd_Qdiag")
cat(sprintf("\nccd_Qdiag dims: new %s old %s\n",
            paste(dim(qN), collapse="x"), paste(dim(qO), collapse="x")))
if(all(dim(qN) == dim(qO))) {
    d <- abs(qN - qO)
    cat(sprintf("max|Qdiag diff| overall: %.3e\n", max(d, na.rm=TRUE)))
    cat("per-CCD-point (col) max|Qdiag diff|:\n"); print(round(apply(d, 2, max, na.rm=TRUE), 8))
}
cN <- attr(dn, "closest_vec"); cO <- attr(do, "closest_vec")
cat(sprintf("\nclosest vector identical: %s   # differing draws: %d / %d\n",
            identical(cN, cO), sum(cN != cO), length(cN)))

# And the final inner draws themselves
inN <- dn[!(rownames(dn) %in% c("log_tau","logit_phi","log_tauEps")), ]
inO <- do[!(rownames(do) %in% c("log_tau","logit_phi","log_tauEps")), ]
cat(sprintf("\ninner draws max|diff|: %.3e   # draw-columns differing >1e-8: %d / %d\n",
            max(abs(inN - inO)),
            sum(apply(abs(inN - inO), 2, max) > 1e-8), ncol(inN)))
cat("\n=== done ===\n")
