# Timing under the final intended config (leaned INLA draws; bias.correct
# dropped). On one M_D_BYM2 fit, measure:
#   (1) sdreport WITH bias.correct+sd  vs  WITHOUT  (the fit-side cost of
#       dropping bias.correct)
#   (2) Gaussian posterior draws  vs  (3) INLA-style draws (leaned)
# Draw timing is independent of bias.correct (that's a fit/sdreport step).

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
obj <- res$TMBobj

# (1) sdreport timing: with vs without bias.correct
t <- proc.time()[3]
SD_bc <- TMB::sdreport(obj, getJointPrecision = TRUE,
                       bias.correct = TRUE, bias.correct.control = list(sd = TRUE))
t_bc <- proc.time()[3] - t
t <- proc.time()[3]
SD_nobc <- TMB::sdreport(obj, getJointPrecision = TRUE)
t_nobc <- proc.time()[3] - t

# (2) Gaussian draws (may fail if jointPrecision non-PD -> that's why INLA exists)
res$TMBsd <- SD_nobc
t <- proc.time()[3]
gaussOK <- TRUE
dr_g <- tryCatch(posteriorDraws(res, NDRAWS = 1000, useInla = FALSE),
                 error = function(e) { cat("Gaussian FAILED:", conditionMessage(e), "\n"); gaussOK <<- FALSE; NULL })
t_gauss <- proc.time()[3] - t

# (3) INLA draws (leaned)
t <- proc.time()[3]
dr_i <- posteriorDraws(res, NDRAWS = 1000, useInla = TRUE)
t_inla <- proc.time()[3] - t

cat("\n================ TIMING (M_D_BYM2 sim 17) ================\n")
cat(sprintf("sdreport WITH bias.correct+sd : %7.1f s\n", t_bc))
cat(sprintf("sdreport WITHOUT bias.correct  : %7.1f s   (%.1fx faster)\n",
            t_nobc, t_bc / t_nobc))
cat(sprintf("Gaussian draws (useInla=FALSE) : %7.1f s   %s\n",
            t_gauss, if(gaussOK) "" else "[FAILED - non-PD jointPrecision]"))
cat(sprintf("INLA draws (useInla=TRUE, leaned): %7.1f s\n", t_inla))
if(gaussOK)
    cat(sprintf("\nINLA / Gaussian draw-time ratio : %.0fx\n", t_inla / t_gauss))
cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
