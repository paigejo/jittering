# Confirm the memory-leaned INLA draws are statistically identical to the
# previous (HEAD) version. Same fit, same RNG seed, draws computed by NEW
# code then by OLD code (sourced from code/inlaStyleDraws_OLD.R = git HEAD).
#
# Expectation: bitwise identical. The CCD walk (.evalCCD x ~25) consumes no
# R RNG (TMB inner Newton is deterministic), so the first rnorm is Step 4's
# split-Gaussian sampling, identical in both; factorization is supernodal in
# both on the same Q; only memory management differs. If max abs diff is ~0
# (or < 1e-10 fp), the leaning is provably transparent to results.

source("code/setup.R")          # loads NEW inlaStyleDraws.R
options(warn = 1)

SEED   <- 4321
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

cat(sprintf("\n=== %s | NEW INLA draws (leaned, one-factor-at-a-time) ===\n",
            format(Sys.time())))
set.seed(SEED)
dr_new <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)

cat(sprintf("\n=== %s | sourcing OLD inlaStyleDraws (git HEAD) ===\n",
            format(Sys.time())))
source("code/inlaStyleDraws_OLD.R")   # overwrites posteriorDraws/.evalCCD/inlaStyleDraws
set.seed(SEED)
dr_old <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)

cat("\n================ OLD vs NEW INLA draws ================\n")
cat(sprintf("dims: new %s   old %s\n",
            paste(dim(dr_new), collapse="x"), paste(dim(dr_old), collapse="x")))
cat(sprintf("rownames identical: %s\n", identical(rownames(dr_new), rownames(dr_old))))

if(all(dim(dr_new) == dim(dr_old))) {
    d <- max(abs(dr_new - dr_old))
    cat(sprintf("max abs diff (all params x draws): %.3e   %s\n", d,
                ifelse(d == 0, "BITWISE IDENTICAL",
                       ifelse(d < 1e-10, "identical (fp)",
                              ifelse(d < 1e-4, "tiny", "DIFFERS")))))
    # Per-parameter-group marginal comparison (robust even if not bitwise)
    cat("\nper-group marginal mean / sd (new vs old):\n")
    for(nm in unique(rownames(dr_new))) {
        a <- dr_new[rownames(dr_new) == nm, ]
        b <- dr_old[rownames(dr_old) == nm, ]
        cat(sprintf("  %-12s mean %+.5f / %+.5f   sd %.5f / %.5f\n",
                    nm, mean(a), mean(b), sd(a), sd(b)))
    }
} else {
    cat("DIM MISMATCH — cannot diff directly\n")
}

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
