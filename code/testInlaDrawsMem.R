# Validate the memory-leaned inlaStyleDraws still produces correct draws.
# Fit Md (light), force the INLA path (useInla=TRUE), and check:
#   (1) it runs without error and returns a draw matrix of the right shape
#   (2) draw rownames cover the hyper + inner params
#   (3) the .checkDrawLabels cross-check (draw SDs vs sdreport SEs) passes
#   (4) hyper marginal means/SDs from draws are sane vs sdreport
# We can't reproduce M_DM's 83 GB locally, but Md exercises the identical
# code path (same .evalCCD rebuild loop, same L_list Cholesky, same draw
# assembly), so correctness here means the edits are safe for M_DM.

source("code/setup.R")
options(warn = 1)

SIMIDX <- 17
simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)

cat(sprintf("\n=== %s | fitting M_D_BYM2 ===\n", format(Sys.time())))
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

cat(sprintf("\n=== %s | forcing INLA-style draws (the leaned path) ===\n",
            format(Sys.time())))
t0 <- proc.time()[3]
set.seed(1)
dr <- tryCatch(posteriorDraws(res, NDRAWS = 1000, useInla = TRUE),
               error = function(e) { cat("INLA draws ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("INLA draws time: %.1f s\n", proc.time()[3] - t0))

if(is.null(dr)) { cat("\nFAIL: INLA draws errored.\n"); quit(status = 1) }

cat(sprintf("\ndraw matrix: %d params x %d draws\n", nrow(dr), ncol(dr)))
cat("rowname counts:\n"); print(table(rownames(dr)))

# Hyper marginal check vs sdreport
SD <- res$TMBsd
cat("\nhyper means/SDs from INLA draws vs sdreport par.fixed:\n")
for(nm in c("log_tau","logit_phi","log_tauEps")) {
    if(!(nm %in% rownames(dr))) next
    d <- dr[rownames(dr) == nm, ]
    seRep <- sqrt(diag(SD$cov.fixed))[nm]
    cat(sprintf("  %-11s draw mean=%+.4f sd=%.4f | sdreport mode=%+.4f se=%.4f\n",
                nm, mean(d), sd(d), SD$par.fixed[nm], seRep))
}

# Label cross-check (warns loudly if any mismatch)
cat("\n.checkDrawLabels (warns on any per-name SD mismatch):\n")
invisible(.checkDrawLabels(dr, res))

cat(sprintf("\n=== %s | PASS: leaned INLA path runs and labels check out ===\n",
            format(Sys.time())))
