# Validate the tape-reuse INLA walk against the rebuild-each-step version
# (git HEAD, = the leaned version bitwise per earlier test). Same fit, same
# seed; compare draws with a TOLERANCE (not bitwise): reuse warm-starts the
# inner optimization sequentially instead of rebuilding+centre-warm-starting,
# so inner modes can differ at optimizer tolerance (~1e-6..1e-8). The inner
# effects are well-identified given theta, so modes should still agree
# closely; theta draws share the same RNG so differ only through the
# slightly-different skew scales.
#
# PASS criteria: per-parameter marginal means/SDs agree to a small relative
# tolerance, and the bulk of element-wise diffs are tiny. (Exact bitwise is
# NOT expected and not required.)

source("code/setup.R")              # NEW reuse version
options(warn = 1)

SEED   <- 20260614
SIMIDX <- 17
simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)

cat(sprintf("\n=== %s | fitting M_D_BYM2 (shared by both draw versions) ===\n",
            format(Sys.time())))
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

cat(sprintf("\n=== %s | NEW reuse-tape INLA draws ===\n", format(Sys.time())))
set.seed(SEED); t <- proc.time()[3]
dr_new <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)
t_new <- proc.time()[3] - t
cat(sprintf("  NEW (reuse) walk time: %.1f s\n", t_new))

cat(sprintf("\n=== %s | OLD rebuild-each-step INLA draws (git HEAD) ===\n",
            format(Sys.time())))
source("code/inlaStyleDraws_OLD.R")
set.seed(SEED); t <- proc.time()[3]
dr_old <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)
t_old <- proc.time()[3] - t
cat(sprintf("  OLD (rebuild) walk time: %.1f s\n", t_old))

cat("\n================ REUSE vs REBUILD ================\n")
cat(sprintf("walk time: NEW %.1f s   OLD %.1f s   speedup %.1fx\n",
            t_new, t_old, t_old / t_new))
cat(sprintf("dims: new %s  old %s   rownames identical: %s\n",
            paste(dim(dr_new), collapse="x"), paste(dim(dr_old), collapse="x"),
            identical(rownames(dr_new), rownames(dr_old))))

if(all(dim(dr_new) == dim(dr_old))) {
    d <- abs(dr_new - dr_old)
    cat(sprintf("element-wise abs diff: max=%.3e  median=%.3e  mean=%.3e\n",
                max(d), median(d), mean(d)))
    cat("\nper-group marginal mean / sd (new vs old) and rel diff:\n")
    worst <- 0
    for(nm in unique(rownames(dr_new))) {
        a <- dr_new[rownames(dr_new) == nm, ]; b <- dr_old[rownames(dr_old) == nm, ]
        rmean <- abs(mean(a) - mean(b)) / (abs(mean(b)) + 1e-8)
        rsd   <- abs(sd(a)   - sd(b))   / (abs(sd(b))   + 1e-8)
        worst <- max(worst, rmean, rsd)
        cat(sprintf("  %-12s mean %+.5f/%+.5f (rel %.1e)   sd %.5f/%.5f (rel %.1e)\n",
                    nm, mean(a), mean(b), rmean, sd(a), sd(b), rsd))
    }
    cat(sprintf("\nworst relative marginal diff: %.2e   %s\n", worst,
                ifelse(worst < 1e-2, "PASS (< 1%)",
                       ifelse(worst < 5e-2, "OK (< 5%)", "INVESTIGATE"))))
} else cat("DIM MISMATCH\n")

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
