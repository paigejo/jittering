# Does sdreport's bias.correct affect ANY scored output? Fit a light BYM2
# model once, run sdreport TWICE on the same fitted obj (with the current
# bias.correct=TRUE + sd=TRUE, and without bias.correct), then score each
# with an identical RNG seed and compare scoresFE / scoresHyper / scoresArea.
#
# Rationale: bias.correct only adds bias-corrected ADREPORT values
# (log_tau, logit_phi) to the sdreport object; it does NOT change
# par.fixed / par.random / jointPrecision / cov.fixed, which is all the
# scoring path consumes. So with a fixed seed the scores should be IDENTICAL
# (max abs diff 0), proving the option can be turned off with no effect on
# results — only on memory/time.
#
# Md is light enough to fit locally; the bias.correct code path is the same
# for M_DM_BYM2, so the conclusion transfers.

source("code/setup.R")
options(warn = 1)

SIMIDX <- 17
MOD    <- "Md"

simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)
truths   <- modelTruths("bym2")
areaPops <- simEnv$areaPops_smoothRisk

cat(sprintf("\n=== %s | fitting %s sim %d (gives us the fitted obj) ===\n",
            format(Sys.time()), MOD, SIMIDX))
res <- .fitOne(MOD, ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))
obj <- res$TMBobj

# --- two sdreports on the SAME obj ---
cat(sprintf("\n=== %s | sdreport WITH bias.correct + sd=TRUE (current) ===\n",
            format(Sys.time())))
t0 <- proc.time()[3]
SD_bc <- TMB::sdreport(obj, getJointPrecision = TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
cat(sprintf("  sdreport(bias.correct) time: %.1f s\n", proc.time()[3] - t0))

cat(sprintf("\n=== %s | sdreport WITHOUT bias.correct ===\n", format(Sys.time())))
t1 <- proc.time()[3]
SD_nobc <- TMB::sdreport(obj, getJointPrecision = TRUE)
cat(sprintf("  sdreport(no bias.correct) time: %.1f s\n", proc.time()[3] - t1))

# --- the core objects the scoring path uses must be identical ---
cat("\n=== core sdreport objects: bias.correct vs not ===\n")
cat(sprintf("  max|par.fixed  diff| : %.3e\n",
            max(abs(SD_bc$par.fixed  - SD_nobc$par.fixed))))
cat(sprintf("  max|par.random diff| : %.3e\n",
            max(abs(SD_bc$par.random - SD_nobc$par.random))))
cat(sprintf("  max|cov.fixed  diff| : %.3e\n",
            max(abs(SD_bc$cov.fixed  - SD_nobc$cov.fixed))))
cat(sprintf("  max|jointPrec  diff| : %.3e\n",
            max(abs(as.matrix(SD_bc$jointPrecision) - as.matrix(SD_nobc$jointPrecision)))))

# --- score each, identical seed before each scoring call ---
scoreWithSeed <- function(SD, seed = 99) {
    r <- res; r$TMBsd <- SD
    set.seed(seed); fe    <- .scoreFE(r, 1000)
    set.seed(seed); hyper <- .scoreHyper(r, truths, 1000)
    set.seed(seed); area  <- .scoreArea(r, areaPops, SIMIDX, 1000)
    list(FE = fe, Hyper = hyper, Area = area)
}
cat(sprintf("\n=== %s | scoring with bias.correct SD ===\n", format(Sys.time())))
sc_bc   <- scoreWithSeed(SD_bc)
cat(sprintf("\n=== %s | scoring with no-bias.correct SD ===\n", format(Sys.time())))
sc_nobc <- scoreWithSeed(SD_nobc)

cmp <- function(a, b, nm) {
    if(is.null(a) && is.null(b)) { cat(sprintf("  %-10s both NULL\n", nm)); return(invisible()) }
    if(is.null(a) || is.null(b)) { cat(sprintf("  %-10s ONE NULL (bc=%s nobc=%s)\n",
                                                nm, !is.null(a), !is.null(b))); return(invisible()) }
    d <- max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
    cat(sprintf("  %-10s max abs diff = %.3e   %s\n", nm, d,
                ifelse(d == 0, "IDENTICAL", ifelse(d < 1e-8, "identical (fp)", "DIFFERS"))))
}
cat("\n================ SCORE COMPARISON (bias.correct vs not) ================\n")
cmp(sc_bc$FE,    sc_nobc$FE,    "scoresFE")
cmp(sc_bc$Hyper, sc_nobc$Hyper, "scoresHyper")
cmp(sc_bc$Area,  sc_nobc$Area,  "scoresArea")

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
