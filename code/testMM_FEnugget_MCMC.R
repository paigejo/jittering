#!/usr/bin/env Rscript
# FE + nugget model: TMB optimization + TMBstan MCMC
# Fixed effects for alpha, beta (urban, normPop)
# Random effects for nuggets only (no spatial RE)
# laplace=FALSE: all nuggets sampled by HMC

library(TMB)
library(tmbstan)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

saveFile = "savedOutput/testMM_FEnugget_MCMC_progress.RData"

cat("===============================================\n")
cat("FE + NUGGET MODEL: TMB + TMBstan\n")
cat("===============================================\n\n")

# Load simulated data
cat("Loading simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

cat("  DHS survey: n =", nrow(ed), "\n")
cat("  MICS survey: n =", nrow(edMICS), "\n\n")

# Load previous results if available
result_opt_FE = NULL
result_mcmc_FE = NULL
mcmc_samples = NULL
if(file.exists(saveFile)) {
  cat("Loading saved results from", saveFile, "\n")
  load(saveFile)
}

# ---- TMB Optimization ----
if(is.null(result_opt_FE)) {
  cat("FE+Nugget: TMB Optimization (urban + normPop)\n")
  cat("-------------------------------------------\n")
  
  tryCatch({
    result_opt_FE = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                  covariates=c("urban", "normPop"),
                                  fixedEffectsOnly=TRUE,
                                  getSDs=TRUE, doMCMC=FALSE)
    save(result_opt_FE, result_mcmc_FE, mcmc_samples, file=saveFile)
    cat("\nTMB optimization completed and saved.\n")
  }, error = function(e) {
    cat("\nTMB optimization FAILED:\n", conditionMessage(e), "\n")
  })
} else {
  cat("TMB optimization already completed (loaded from save).\n")
}

# Print TMB results
if(!is.null(result_opt_FE)) {
  cat("  Time:", result_opt_FE$totalTime/60, "minutes\n")
  SD0 = result_opt_FE$TMBsd
  if(inherits(SD0, "sdreport")) {
    ss = summary(SD0, "fixed")
    cat("\nTMB Fixed effects (with SE):\n")
    print(round(ss, 4))
  }
}

# ---- TMBstan MCMC ----
if(is.null(result_mcmc_FE)) {
  cat("\n\nFE+Nugget: TMBstan MCMC (laplace=FALSE)\n")
  cat("-------------------------------------------\n")
  cat("Sampling all parameters including", 741+1442, "nuggets\n")
  cat("This will take a while...\n\n")
  
  tryCatch({
    result_mcmc_FE = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                    covariates=c("urban", "normPop"),
                                    fixedEffectsOnly=TRUE,
                                    getSDs=FALSE, doMCMC=TRUE)
    # Extract and save samples
    mcmc_samples = as.matrix(result_mcmc_FE$TMBsd)
    save(result_opt_FE, result_mcmc_FE, mcmc_samples, file=saveFile)
    cat("\nMCMC completed and saved.\n")
  }, error = function(e) {
    cat("\nMCMC FAILED:\n", conditionMessage(e), "\n")
  })
} else {
  cat("\nMCMC already completed (loaded from save).\n")
}

# ---- Print MCMC results ----
if(!is.null(mcmc_samples)) {
  cat("\n\nMCMC Results\n")
  cat("-------------------------------------------\n")
  
  parNames = colnames(mcmc_samples)
  
  # Fixed-effect-like params
  feParams = c("log_tauEps", "alpha", "beta[1]", "beta[2]")
  feLabels = c("log_tauEps", "alpha", "urban", "normPop")
  truths = c(log(1/1.5), -1.25, 1.00, 0.50)
  
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
    cat(sprintf("%-18s %6.3f  %8.4f  %8.4f  %8.4f  %8.4f\n",
                "sigmaEps", sqrt(1.5), mean(sigmaEps_samp), sd(sigmaEps_samp),
                quantile(sigmaEps_samp, 0.025), quantile(sigmaEps_samp, 0.975)))
  }
  
  # ---- Pair plots ----
  cat("\nGenerating pair plots...\n")
  pairCols = feParams[feParams %in% parNames]
  pairLabels = feLabels[feParams %in% parNames]
  
  pdf("Figures/testMCMC/testMM_FEnugget_MCMC_pairs.pdf", width=10, height=10)
  pairs(mcmc_samples[, pairCols, drop=FALSE], labels=pairLabels,
        pch=".", col=rgb(0,0,0,0.15),
        main="FE+Nugget MCMC: Posterior Pair Plot")
  # Add truth lines
  dev.off()
  cat("Pair plot saved to Figures/testMCMC/testMM_FEnugget_MCMC_pairs.pdf\n")
  
  # ---- Trace plots ----
  cat("Generating trace plots...\n")
  nTrace = length(pairCols)
  pdf("Figures/testMCMC/testMM_FEnugget_MCMC_trace.pdf", width=10, height=2.5*nTrace)
  par(mfrow=c(nTrace, 1), mar=c(2.5, 4, 1.5, 1))
  for(k in 1:nTrace) {
    samp = mcmc_samples[, pairCols[k]]
    plot(samp, type="l", col="steelblue",
         ylab=pairLabels[k], xlab="Iteration", main=pairLabels[k], cex.main=0.9)
    abline(h=truths[k], col="red", lwd=2, lty=2)
  }
  dev.off()
  cat("Trace plot saved to Figures/testMCMC/testMM_FEnugget_MCMC_trace.pdf\n")
}

cat("\n===============================================\n")
cat("FE+NUGGET TEST COMPLETE\n")
cat("===============================================\n")
