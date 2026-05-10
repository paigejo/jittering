#!/usr/bin/env Rscript
# FE + nugget model: TMBstan MCMC with laplace=TRUE on BYM2-simulated data
# Fixed effects for alpha, beta (urban, normPop) sampled by HMC
# Random effects (nuggets) integrated out by Laplace at each HMC step
#
# True parameters used in BYM2 simulation:
#   alpha (intercept) = -1.25
#   gamma (urban)      =  1.00
#   betaNormPop        =  0.50
#   sigmaEpsilon       = sqrt(1.5)  (nugget SD)
#   => log_tauEps      = log(1/1.5) = -0.405

library(TMB)
library(tmbstan)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# Different save file from laplace=FALSE run
saveFile = "savedOutput/testMM_BYM2sim_FEnugget_MCMC_laplaceTrue_progress.RData"

cat("===============================================\n")
cat("FE + NUGGET: TMBstan laplace=TRUE on BYM2 DATA\n")
cat("===============================================\n\n")

# Load BYM2 simulated data
cat("Loading BYM2 simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
cat("Loading DHS integration points...\n")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

cat("Data summary:\n")
cat("  DHS clusters:", nrow(ed), "\n")
cat("  MICS clusters:", nrow(edMICS), "\n\n")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)
trueLogTauEps = log(1/1.5)

# Load previous results
result_mcmc_laplace = NULL
mcmc_samples = NULL
if(file.exists(saveFile)) {
  cat("Loading saved results from", saveFile, "\n")
  load(saveFile)
}

# ---- Build TMB object via fitMM_IIDonly (optimization only, to get obj) ----
if(is.null(result_mcmc_laplace)) {
  cat("Building TMB object via fitMM_IIDonly...\n")

  result_build = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                covariates=c("urban", "normPop"),
                                fixedEffectsOnly=TRUE,
                                getSDs=FALSE, doMCMC=FALSE)
  obj = result_build$TMBobj

  cat("\nOuter parameters (sampled by HMC):", names(obj$par), "\n")
  cat("Number of outer params:", length(obj$par), "\n\n")

  # ---- TMBstan with laplace=TRUE ----
  cat("Starting TMBstan with laplace=TRUE...\n")
  cat("HMC samples only outer params; nuggets integrated out by Laplace.\n\n")

  startTime = proc.time()[3]
  tryCatch({
    fit = tmbstan(obj=obj, silent=FALSE, laplace=TRUE,
                  iter=4000, warmup=2000, chains=1)
    endTime = proc.time()[3]
    totalTime = endTime - startTime
    cat(sprintf("\nMCMC took %.1f minutes\n", totalTime/60))

    mcmc_samples = as.matrix(fit)
    result_mcmc_laplace = list(fit=fit, totalTime=totalTime)
    save(result_mcmc_laplace, mcmc_samples, file=saveFile)
    cat("Results saved to", saveFile, "\n")
  }, error = function(e) {
    cat("\nMCMC FAILED:\n", conditionMessage(e), "\n")
  })
} else {
  cat("MCMC already completed (loaded from save).\n")
}

# ---- Print results ----
if(!is.null(mcmc_samples)) {
  cat("\n\nMCMC Results (laplace=TRUE, BYM2 data)\n")
  cat("===================================================\n")

  parNames = colnames(mcmc_samples)
  cat("Sampled parameters:", paste(parNames, collapse=", "), "\n")

  feParams = c("log_tauEps", "alpha", "beta[1]", "beta[2]")
  feLabels = c("log_tauEps", "alpha", "urban", "normPop")
  truths = c(trueLogTauEps, trueAlpha, trueUrban, trueNormPop)

  cat("\nParameter          Truth     Mean      SD        2.5%      97.5%\n")
  cat("----------------------------------------------------------------\n")
  for(i in 1:length(feParams)) {
    if(feParams[i] %in% parNames) {
      samp = mcmc_samples[, feParams[i]]
      cat(sprintf("%-18s %6.3f  %8.4f  %8.4f  %8.4f  %8.4f\n",
                  feLabels[i], truths[i], mean(samp), sd(samp),
                  quantile(samp, 0.025), quantile(samp, 0.975)))
    }
  }

  # sigmaEps derived
  if("log_tauEps" %in% parNames) {
    sigmaEps_samp = exp(-0.5 * mcmc_samples[, "log_tauEps"])
    cat(sprintf("\n%-18s %6.3f  %8.4f  %8.4f  %8.4f  %8.4f\n",
                "sigmaEps", trueSigmaEps, mean(sigmaEps_samp), sd(sigmaEps_samp),
                quantile(sigmaEps_samp, 0.025), quantile(sigmaEps_samp, 0.975)))
  }

  # normPop coverage
  if("beta[2]" %in% parNames) {
    samp = mcmc_samples[, "beta[2]"]
    ci = quantile(samp, c(0.025, 0.975))
    covers = (trueNormPop >= ci[1]) && (trueNormPop <= ci[2])
    cat(sprintf("\nnormPop: 95%% CI [%.3f, %.3f] %s truth (%.2f)\n",
                ci[1], ci[2], ifelse(covers, "COVERS", "MISSES"), trueNormPop))
  }

  # ---- Trace and pair plots ----
  dir.create("Figures/testMCMC", recursive=TRUE, showWarnings=FALSE)

  pairCols = feParams[feParams %in% parNames]
  pairLabels = feLabels[feParams %in% parNames]

  pdf("Figures/testMCMC/testMM_BYM2sim_FEnugget_MCMC_laplaceTrue_pairs.pdf", width=10, height=10)
  pairs(mcmc_samples[, pairCols, drop=FALSE], labels=pairLabels,
        pch=".", col=rgb(0,0,0,0.15),
        main="FE+Nugget MCMC (laplace=TRUE) on BYM2 data: Pairs")
  dev.off()
  cat("\nPair plot saved.\n")

  nTrace = length(pairCols)
  pdf("Figures/testMCMC/testMM_BYM2sim_FEnugget_MCMC_laplaceTrue_trace.pdf", width=10, height=2.5*nTrace)
  par(mfrow=c(nTrace, 1), mar=c(2.5, 4, 1.5, 1))
  for(k in 1:nTrace) {
    samp = mcmc_samples[, pairCols[k]]
    plot(samp, type="l", col="steelblue",
         ylab=pairLabels[k], xlab="Iteration", main=pairLabels[k], cex.main=0.9)
    abline(h=truths[k], col="red", lwd=2, lty=2)
  }
  dev.off()
  cat("Trace plot saved.\n")
}

cat("\n===============================================\n")
cat("FE+NUGGET MCMC (laplace=TRUE) ON BYM2 DATA: COMPLETE\n")
cat("===============================================\n")
