# Test the INLA-style smooth hyperparameter integration on BYM2 sim 1, models
# M_M and M_D. Compares:
#   (a) standard joint-Gaussian draws from cov.fixed / jointPrecision
#   (b) INLA-style draws (multi-step axial walk + per-axis skew-normal +
#       per-CCD-point cached Gaussian for x|θ)
# Then exercises predGrid(useInla=TRUE) for an area-level check on M_M.

source("setup.R")
options(error = traceback)
source("code/modFEM.R")
source("code/modFED.R")
source("code/modFEMD.R")
source("code/modM_DSep.R")
source("code/modM_MSep.R")
source("code/modM_DMSep.R")
source("code/modMdSep.R")
source("code/makeInputsTMB.R")
source("code/modBYM2.R")
source("code/testInfrastructure.R")
source("code/scoreSimStudy.R")
source("code/inlaStyleDraws.R")

simIdx <- 1
model  <- "bym2"
NDRAWS <- 1000
truths <- modelTruths(model)

# Load sim 1
simEnv <- new.env()
simulateSurveys(model, nsim = simIdx, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
inp <- buildInputsForSim(simIdx, simEnv$surveysDHS, simEnv$surveysMICS,
                         micsEnv$intPtsMICS, model,
                         KMICS = 100, KDHSu = 16, KDHSr = 21)

# Print marginal summaries for a draws matrix
print_marginals <- function(label, draws, truths) {
    cat(sprintf("\n----- %s -----\n", label))
    rn <- rownames(draws)
    spec <- list(
        list(name="alpha",       row="alpha",      transform=identity,
             truth=truths$alpha),
        list(name="beta_urban",  row="beta",       row_idx=1, transform=identity,
             truth=1.0),
        list(name="beta_normPop",row="beta",       row_idx=5, transform=identity,
             truth=0.5),
        list(name="sigma_BYM2",  row="log_tau",    transform=function(x) exp(-0.5*x),
             truth=truths$sigmaSpatial),
        list(name="phi",         row="logit_phi",  transform=plogis,
             truth=truths$phi),
        list(name="sigma_eps",   row="log_tauEps", transform=function(x) exp(-0.5*x),
             truth=truths$sigmaEps))
    for(s in spec) {
        idx <- which(rn == s$row)
        if(length(idx) == 0) { cat(sprintf("  %-14s  not found\n", s$name)); next }
        if(!is.null(s$row_idx)) idx <- idx[s$row_idx]
        if(length(idx) == 0) { cat(sprintf("  %-14s  index out of range\n", s$name)); next }
        d <- s$transform(draws[idx, ])
        cat(sprintf("  %-14s  mean=%.4f  sd=%.4f  q025=%.4f  q975=%.4f   truth=%s\n",
                    s$name, mean(d), sd(d), quantile(d, .025), quantile(d, .975),
                    if(is.na(s$truth)) "NA" else sprintf("%.4f", s$truth)))
    }
}

# Both Gaussian and INLA-style baselines now live in inlaStyleDraws.R's
# posteriorDraws() dispatcher. No duplicated logic here.

# ============================================================================
# M_M (MICS-only BYM2)
# ============================================================================
cat("\n##### M_M_BYM2 sim 1 #####\n")
t0 <- proc.time()[3]
res_MM <- fitMM(datDHS = inp$datDHS, datMICS = inp$datMICS,
                inputsMDM = inp$inputsMDM,
                KMICS = 100, Qgh = 10, getSDs = TRUE, verbose = FALSE)
cat(sprintf("fitMM time: %.1f min\n", (proc.time()[3]-t0)/60))
cat(sprintf("fitMM mode: log_tau=%.3f  logit_phi=%.3f  log_tauEps=%.3f  pdHess=%s\n",
            res_MM$opt$par["log_tau"], res_MM$opt$par["logit_phi"],
            res_MM$opt$par["log_tauEps"],
            as.character(res_MM$TMBsd$pdHess)))

g_MM <- posteriorDraws(res_MM, NDRAWS = NDRAWS, useInla = FALSE)
print_marginals("M_M Gaussian (posteriorDraws useInla=FALSE)", g_MM, truths)

t0 <- proc.time()[3]
i_MM <- posteriorDraws(res_MM, NDRAWS = NDRAWS, useInla = TRUE,
                       deltaZ = 1.0, deltaPi = 2.5, maxAxialSteps = 4)
cat(sprintf("inlaStyleDraws time: %.1f min\n", (proc.time()[3]-t0)/60))
print_marginals("M_M INLA-style (posteriorDraws useInla=TRUE)", i_MM, truths)

# predGrid comparison (area-level)
cat("\n--- predGrid (useInla=FALSE vs useInla=TRUE) ---\n")
t0 <- proc.time()[3]
grid_MM_g <- predGrid(res_MM$TMBsd, popMat = popMatNGAThresh, nsim = NDRAWS,
                     obj = res_MM$TMBobj, admLevel = "stratMICS",
                     useInla = FALSE)
cat(sprintf("  Gaussian predGrid: %.1f s\n", proc.time()[3]-t0))
t0 <- proc.time()[3]
grid_MM_i <- predGrid(res_MM$TMBsd, popMat = popMatNGAThresh, nsim = NDRAWS,
                     obj = res_MM$TMBobj, admLevel = "stratMICS",
                     useInla = TRUE, res = res_MM)
cat(sprintf("  INLA  predGrid: %.1f s\n", proc.time()[3]-t0))

# Aggregate to area level for both methods
agg_MM_g <- predArea(grid_MM_g, areaVarName = "area", orderedAreas = adm1@data$NAME_1)
agg_MM_i <- predArea(grid_MM_i, areaVarName = "area", orderedAreas = adm1@data$NAME_1)
truth_area <- simEnv$areaPops[, simIdx]
estG <- agg_MM_g$aggregationResults
estI <- agg_MM_i$aggregationResults
cat(sprintf("Area-level RMSE   Gaussian = %.4f   INLA = %.4f\n",
            sqrt(mean((estG$p - truth_area)^2)),
            sqrt(mean((estI$p - truth_area)^2))))
cat(sprintf("Area-level CI95 width  Gaussian = %.4f   INLA = %.4f\n",
            mean(estG[["upper_2.5%"]] - estG[["lower_97.5%"]], na.rm=TRUE),
            mean(estI[["upper_2.5%"]] - estI[["lower_97.5%"]], na.rm=TRUE)))

# ============================================================================
# M_D (DHS-only BYM2)
# ============================================================================
cat("\n\n##### M_D_BYM2 sim 1 #####\n")
t0 <- proc.time()[3]
res_MD <- fitMD(datDHS = inp$datDHS, datMICS = inp$datMICS,
                inputsMDM = inp$inputsMDM,
                KMICS = 100, KDHSurb = 16, KDHSrur = 21,
                Qgh = 10, getSDs = TRUE, verbose = FALSE)
cat(sprintf("fitMD time: %.1f min\n", (proc.time()[3]-t0)/60))
cat(sprintf("fitMD mode: log_tau=%.3f  logit_phi=%.3f  log_tauEps=%.3f  pdHess=%s\n",
            res_MD$opt$par["log_tau"], res_MD$opt$par["logit_phi"],
            res_MD$opt$par["log_tauEps"],
            as.character(res_MD$TMBsd$pdHess)))

g_MD <- tryCatch(posteriorDraws(res_MD, NDRAWS = NDRAWS, useInla = FALSE),
                 error = function(e) { cat("gaussian draws failed:", conditionMessage(e), "\n"); NULL })
if(!is.null(g_MD)) print_marginals("M_D Gaussian (posteriorDraws useInla=FALSE)", g_MD, truths)

t0 <- proc.time()[3]
i_MD <- tryCatch(posteriorDraws(res_MD, NDRAWS = NDRAWS, useInla = TRUE,
                                deltaZ = 1.0, deltaPi = 2.5, maxAxialSteps = 4),
                 error = function(e) { cat("INLA draws failed:", conditionMessage(e), "\n"); NULL })
cat(sprintf("inlaStyleDraws time: %.1f min\n", (proc.time()[3]-t0)/60))
if(!is.null(i_MD)) print_marginals("M_D INLA-style (posteriorDraws useInla=TRUE)", i_MD, truths)

cat("\n\nDone.\n")
