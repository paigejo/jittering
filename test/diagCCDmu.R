# Diff the per-CCD-point inner modes (ccd_mu) across the WHOLE walk between
# reuse and rebuild. If these are identical, mu_k is not the cause and the
# divergence is in Cholesky/solve; if they differ, we see which point/param.
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

muN <- attr(dn, "ccd_mu"); muO <- attr(do, "ccd_mu")
cat(sprintf("\nccd_mu dims: new %s  old %s\n",
            paste(dim(muN), collapse="x"), paste(dim(muO), collapse="x")))
if(all(dim(muN) == dim(muO))) {
    d <- abs(muN - muO)
    cat(sprintf("max|ccd_mu diff| overall: %.3e\n", max(d)))
    cat("per-CCD-point (column) max|diff|:\n")
    print(round(apply(d, 2, max), 6))
    # which rows (params) drive the worst column
    wc <- which.max(apply(d, 2, max))
    cat(sprintf("\nworst CCD point = col %d; per-param max there:\n", wc))
    rn <- rownames(muN)
    for(nm in unique(rn))
        cat(sprintf("  %-12s %.3e\n", nm, max(d[rn == nm, wc])))
}
cat("\n=== done ===\n")
