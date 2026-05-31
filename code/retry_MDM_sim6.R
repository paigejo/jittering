# Retry the failed M_DM_BYM2 fit on BYM2-simulated sim 6.
# Diagnose: was the fitFEMD warm-start giving non-finite values?
# If so, replace them with neutral defaults and refit.

source("setup.R")
options(error = recover)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/modM_DSep.R")
source("code/modM_MSep.R")
source("code/modM_DMSep.R")
source("code/modMdSep.R")
source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/scoreSimStudy.R")

simIdx <- 6
model  <- "bym2"
truths <- modelTruths(model)
outFile <- sprintf("savedOutput/simStudy1/scores/BYM2/scores_M_DM_BYM2_sim%d.RData", simIdx)

simEnv <- new.env()
simulateSurveys(model, nsim = 10, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

inp <- buildInputsForSim(simIdx, simEnv$surveysDHS, simEnv$surveysMICS,
                         micsEnv$intPtsMICS, model,
                         KMICS = 100, KDHSu = 16, KDHSr = 21)

cat("\n=== diagnostic: fitFEMD warm-start ===\n")
fe <- fitFEMD(datDHS = inp$datDHS, datMICS = inp$datMICS,
              inputsMDM = inp$inputsMDM,
              KMICS = 100, Qgh = 10,
              fixedEffectsOnly = TRUE, getSDs = FALSE, verbose = FALSE)
print(fe$opt$par)
cat("any non-finite?", any(!is.finite(fe$opt$par)), "\n")
cat("convergence:", fe$opt$convergence, "\n")

cat("\n=== retry 1: fitMDM (default warm-start) ===\n")
res <- tryCatch(
    fitMDM(datDHS = inp$datDHS, datMICS = inp$datMICS,
           inputsMDM = inp$inputsMDM,
           KMICS = 100, KDHSurb = 16, KDHSrur = 21,
           Qgh = 10, getSDs = TRUE, verbose = FALSE),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })

if(is.null(res)) {
    cat("\n=== retry 2: try perturbed initial conditions via set.seed ===\n")
    # Some TMB fits depend on the order RNG is touched. Re-seed and retry.
    set.seed(simIdx * 1000)
    res <- tryCatch(
        fitMDM(datDHS = inp$datDHS, datMICS = inp$datMICS,
               inputsMDM = inp$inputsMDM,
               KMICS = 100, KDHSurb = 16, KDHSrur = 21,
               Qgh = 10, getSDs = TRUE, verbose = FALSE),
        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
}

if(is.null(res)) {
    cat("\nAll retries failed. Aborting.\n")
    quit(status = 1)
}

cat("\n=== fit summary ===\n")
print(res$opt$par)
cat("convergence:", res$opt$convergence, "  pdHess:",
    if(inherits(res$TMBsd,"sdreport")) res$TMBsd$pdHess else NA, "\n")

# Compute scores and save with the new fields.
NDRAWS <- 1000
t0 <- proc.time()[3]
scoresFE    <- tryCatch(.scoreFE(res, NDRAWS), error = function(e) { cat("scoreFE error:", conditionMessage(e),"\n"); NULL })
scoresHyper <- tryCatch(.scoreHyper(res, truths, NDRAWS), error = function(e) { cat("scoreHyper error:", conditionMessage(e),"\n"); NULL })
scoresArea  <- tryCatch(.scoreArea(res, simEnv$areaPops, simIdx, NDRAWS),
                        error = function(e) { cat("scoreArea error:", conditionMessage(e),"\n"); NULL })
scoreTime <- (proc.time()[3] - t0)/60
fitTime   <- NA_real_   # we didn't time the (variable) retry path

save(scoresFE, scoresHyper, scoresArea, fitTime, scoreTime, file = outFile)
cat("Wrote", outFile, "\n")
