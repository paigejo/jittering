# Compare the CCD structure + RNG-sensitive assignment between reuse and
# rebuild draws (same seed). inlaStyleDraws attaches: n_ccd, ccd_nll,
# closest_counts. If n_ccd or the dNLL set differ, max.col(ties="random")
# consumes different RNG and desyncs the inner sampling -> the observed
# "hypers match / inner differ" signature.

source("code/setup.R")
options(warn = 1)

SEED <- 20260614; SIMIDX <- 17
simEnv <- new.env(); simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env(); load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2", KMICS = 100, KDHSu = 16, KDHSr = 21)
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

set.seed(SEED); dr_new <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)
source("code/inlaStyleDraws_OLD.R")
set.seed(SEED); dr_old <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)

cat("\n================ CCD structure: reuse vs rebuild ================\n")
cat(sprintf("n_ccd:   new=%s  old=%s\n", attr(dr_new,"n_ccd"), attr(dr_old,"n_ccd")))
nllN <- attr(dr_new,"ccd_nll"); nllO <- attr(dr_old,"ccd_nll")
cat(sprintf("length(ccd_nll): new=%d old=%d\n", length(nllN), length(nllO)))
if(length(nllN) == length(nllO))
    cat(sprintf("max|ccd_nll diff| = %.3e\n", max(abs(nllN - nllO))))
ccN <- attr(dr_new,"closest_counts"); ccO <- attr(dr_old,"closest_counts")
cat(sprintf("length(closest_counts): new=%d old=%d\n", length(ccN), length(ccO)))
if(length(ccN) == length(ccO)) {
    cat(sprintf("closest_counts identical: %s   max|diff|=%d\n",
                identical(ccN, ccO), as.integer(max(abs(ccN - ccO)))))
    cat("new counts: ", paste(ccN, collapse=","), "\n")
    cat("old counts: ", paste(ccO, collapse=","), "\n")
}
cat("\n=== done ===\n")
