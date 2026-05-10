#!/usr/bin/env Rscript
# Simple test of IID spatial effects model
# Tests whether data simulation or BYM2 structure is the issue
#
# Optimization results (sim dataset 1):
#   log_tau:             0.0770  => sigma_u = 0.9622 (truth: sqrt(0.5) = 0.7071)
#   alpha (intercept):  -2.1736  (truth: -1.25)
#   beta[urban]:         1.5628  (truth:  1.00)
#   beta[access]:        1.6855  (truth:  0.00)
#   beta[elev]:          0.2374  (truth:  0.00)
#   beta[distRivers]:    0.5357  (truth:  0.00)
#   beta[normPop]:       2.7472  (truth:  0.50)
#   41 spatial random effects, Hessian PD, ~10.8 min

library(TMB)
library(tmbstan)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

saveFile = "savedOutput/testMM_IIDonly_nugget_progress.RData"

# Load/generate simulation data
cat("===============================================\n")
cat("IID SPATIAL EFFECTS MODEL TEST\n")
cat("===============================================\n\n")

# Load simulated data
cat("Loading simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")

# Load integration points
cat("Loading DHS integration points...\n")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

# Prepare data
ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

cat("\nData summary:\n")
cat("  DHS survey: n =", nrow(ed), "\n")
cat("  MICS survey: n =", nrow(edMICS), "\n")
cat("  Number of areas:", length(unique(edMICS$subarea)), "\n\n")

# Load previous results if available
result_opt = NULL
result_mcmc = NULL
if(file.exists(saveFile)) {
  cat("Loading saved results from", saveFile, "\n")
  load(saveFile)
}

# ---- Test 1: Optimization ----
if(is.null(result_opt)) {
  cat("IID Model with Optimization (no MCMC)\n")
  cat("-------------------------------------------\n")
  
  tryCatch({
    result_opt = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                               getSDs=TRUE, doMCMC=FALSE)
    save(result_opt, result_mcmc, file=saveFile)
    cat("\nOptimization completed and saved.\n")
  }, error = function(e) {
    cat("\nOptimization FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("Optimization already completed (loaded from save).\n")
}

# Print optimization results
if(!is.null(result_opt)) {
  cat("  Total time:", result_opt$totalTime/60, "minutes\n")
  
  SD0 = result_opt$TMBsd
  if(inherits(SD0, "sdreport")) {
    pf = SD0$par.fixed
    pfn = names(pf)
    pr = SD0$par.random
    rn = names(pr)
    
    cat("\nFixed parameters:\n")
    print(round(pf, 4))
    
    cat("\nRandom effects counts:\n")
    print(table(rn))
    
    alpha_est = pr[rn == "alpha"]
    cat("\nalpha (intercept) estimate:", round(as.numeric(alpha_est), 4), 
        " (truth: -1.25)\n")
    
    beta_est = pr[rn == "beta"]
    trueBeta = c(1.00, 0.00, 0.00, 0.00, 0.50)
    trueBetaNames = c("urban", "access", "elev", "distRiversLakes", "normPop")
    nBeta = length(beta_est)
    if(nBeta > 0) {
      cat("\nBeta estimates vs truth:\n")
      betaComp = data.frame(
        name = trueBetaNames[1:nBeta],
        truth = trueBeta[1:nBeta],
        estimate = round(as.numeric(beta_est), 4)
      )
      print(betaComp)
    }
    
    log_tau_est = pf[pfn == "log_tau"]
    sigma_u_est = exp(-0.5 * log_tau_est)
    cat("\nsigma_u:", round(as.numeric(sigma_u_est), 4), 
        " (truth: sqrt(0.5) =", round(sqrt(0.5), 4), ")\n")
    
    log_tauEps_est = pf[pfn == "log_tauEps"]
    if(length(log_tauEps_est) > 0) {
      sigmaEps_est = exp(-0.5 * log_tauEps_est)
      cat("sigmaEps:", round(as.numeric(sigmaEps_est), 4), 
          " (truth: sqrt(1.5) =", round(sqrt(1.5), 4), ")\n")
    }
  }
}

# ---- Test 2: MCMC via TMBstan ----
runMCMC = FALSE  # Set to TRUE to run MCMC (takes ~5 hours)
if(runMCMC && is.null(result_mcmc)) {
  cat("\n\nIID Model with MCMC (TMBstan)\n")
  cat("-------------------------------------------\n")
  
  tryCatch({
    result_mcmc = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                getSDs=FALSE, doMCMC=TRUE)
    # Save posterior samples separately (stanfit objects lose external pointers on save)
    mcmc_samples = as.matrix(result_mcmc$TMBsd)
    save(result_opt, result_mcmc, mcmc_samples, file=saveFile)
    cat("\nMCMC completed and saved.\n")
  }, error = function(e) {
    cat("\nMCMC FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("\n\nMCMC already completed (loaded from save).\n")
}

# Print MCMC results
if(!is.null(result_mcmc)) {
  cat("  Total time:", result_mcmc$totalTime/60, "minutes\n")
  
  # Use saved samples matrix if stanfit is corrupted (external pointers lost on reload)
  if(!exists("mcmc_samples")) {
    fit = result_mcmc$TMBsd
    mcmc_samples = as.matrix(fit)
  }
  posteriorSamples = mcmc_samples
  parNames = colnames(posteriorSamples)
  
  cat("\nMCMC parameter names:\n")
  print(parNames[!grepl("^u_spatial", parNames)])
  
  # Alpha
  if("alpha" %in% parNames) {
    alpha_samples = posteriorSamples[,"alpha"]
    cat("\nalpha (intercept): mean =", round(mean(alpha_samples), 4),
        ", sd =", round(sd(alpha_samples), 4),
        " (truth: -1.25)\n")
  }
  
  # Beta
  beta_idx = grep("^beta\\[", parNames)
  if(length(beta_idx) > 0) {
    trueBeta = c(1.00, 0.00, 0.00, 0.00, 0.50)
    trueBetaNames = c("urban", "access", "elev", "distRiversLakes", "normPop")
    cat("\nBeta estimates vs truth (posterior mean +/- sd):\n")
    nBeta = min(length(beta_idx), length(trueBeta))
    for(j in 1:nBeta) {
      samp = posteriorSamples[, beta_idx[j]]
      cat(sprintf("  %-18s truth=%5.2f  est=%7.4f +/- %6.4f\n",
                  trueBetaNames[j], trueBeta[j], mean(samp), sd(samp)))
    }
  }
  
  # log_tau / sigma_u
  if("log_tau" %in% parNames) {
    log_tau_samples = posteriorSamples[,"log_tau"]
    sigma_u_samples = exp(-0.5 * log_tau_samples)
    cat("\nsigma_u: mean =", round(mean(sigma_u_samples), 4),
        ", sd =", round(sd(sigma_u_samples), 4),
        " (truth: sqrt(0.5) =", round(sqrt(0.5), 4), ")\n")
  }
  
  # Pair plots
  cat("\nGenerating pair plots...\n")
  
  # Pick first u_spatial effect for inclusion
  u_idx = grep("^u_spatial\\[", parNames)
  u1_name = parNames[u_idx[1]]
  
  pairCols = c("log_tau", "alpha", "beta[1]", "beta[2]", "beta[3]", 
               "beta[4]", "beta[5]", u1_name)
  pairCols = pairCols[pairCols %in% parNames]
  pairLabels = c("log_tau", "alpha", "urban", "access", "elev", 
                 "distRivers", "normPop", "u_spatial[1]")
  pairLabels = pairLabels[1:length(pairCols)]
  
  pdf("Figures/testMCMCtMCMC/testMM_IIDonly_pairs.pdf", width=14, height=14)
  pairs(posteriorSamples[, pairCols], labels=pairLabels,
        pch=".", col=rgb(0,0,0,0.15),
        main="IID Model MCMC: Posterior Pair Plot")
  dev.off()
  cat("Pair plot saved to Figures/testMCMCtMCMC/testMM_IIDonly_pairs.pdf\n")
  
  # Trace plots
  cat("Generating trace plots...\n")
  nTrace = length(pairCols)
  pdf("Figures/testMCMCtMCMC/testMM_IIDonly_trace.pdf", width=10, height=2.5*nTrace)
  par(mfrow=c(nTrace, 1), mar=c(2.5, 4, 1.5, 1))
  for(k in 1:nTrace) {
    plot(posteriorSamples[, pairCols[k]], type="l", col="steelblue",
         ylab=pairLabels[k], xlab="", main=pairLabels[k], cex.main=0.9)
  }
  dev.off()
  cat("Trace plot saved to Figures/testMCMCtMCMC/testMM_IIDonly_trace.pdf\n")
}

cat("\n===============================================\n")
cat("TEST COMPLETE (full covariates)\n")
cat("===============================================\n")

# ---- Test 3: 2-covariate model (urban + normPop only), TMB only ----
saveFile2cov = "savedOutput/testMM_IIDonly_2cov_nugget_progress.RData"
result_opt_2cov = NULL
if(file.exists(saveFile2cov)) {
  cat("\nLoading saved 2-covariate results from", saveFile2cov, "\n")
  load(saveFile2cov)
}

if(is.null(result_opt_2cov)) {
  cat("\n\nIID Model (urban + normPop only) with Optimization\n")
  cat("-------------------------------------------\n")
  
  tryCatch({
    result_opt_2cov = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                    covariates=c("urban", "normPop"),
                                    getSDs=TRUE, doMCMC=FALSE)
    save(result_opt_2cov, file=saveFile2cov)
    cat("\n2-covariate optimization completed and saved.\n")
  }, error = function(e) {
    cat("\n2-covariate optimization FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("\n2-covariate optimization already completed (loaded from save).\n")
}

if(!is.null(result_opt_2cov)) {
  cat("  Total time:", result_opt_2cov$totalTime/60, "minutes\n")
  
  SD0 = result_opt_2cov$TMBsd
  if(inherits(SD0, "sdreport")) {
    pf = SD0$par.fixed
    pfn = names(pf)
    pr = SD0$par.random
    rn = names(pr)
    
    alpha_est = pr[rn == "alpha"]
    beta_est = pr[rn == "beta"]
    log_tau_est = pf[pfn == "log_tau"]
    sigma_u_est = exp(-0.5 * log_tau_est)
    
    cat("\n2-covariate IID model estimates vs truth:\n")
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "alpha", -1.25, as.numeric(alpha_est)))
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "urban", 1.00, as.numeric(beta_est[1])))
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "normPop", 0.50, as.numeric(beta_est[2])))
    cat(sprintf("  %-18s truth=%6.4f  est=%7.4f\n", "sigma_u", sqrt(0.5), as.numeric(sigma_u_est)))
    
    log_tauEps_est = pf[pfn == "log_tauEps"]
    if(length(log_tauEps_est) > 0) {
      sigmaEps_est = exp(-0.5 * log_tauEps_est)
      cat(sprintf("  %-18s truth=%6.4f  est=%7.4f\n", "sigmaEps", sqrt(1.5), as.numeric(sigmaEps_est)))
    }
  }
}

cat("\n===============================================\n")
cat("ALL TESTS COMPLETE\n")
cat("===============================================\n")

# ---- Test 4: Fixed effects only (no random effects, u_spatial=0) ----
# This is a GLM through the integration-point machinery,
# so geomasking is still handled via integration points.
saveFileFE = "savedOutput/testMM_IIDonly_FEonly_nugget_progress.RData"
result_opt_FE = NULL
if(file.exists(saveFileFE)) {
  cat("\nLoading saved fixed-effects-only results from", saveFileFE, "\n")
  load(saveFileFE)
}

if(is.null(result_opt_FE)) {
  cat("\n\nIID Model - Fixed Effects Only (urban + normPop, u_spatial=0)\n")
  cat("-------------------------------------------\n")
  
  tryCatch({
    result_opt_FE = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                  covariates=c("urban", "normPop"),
                                  fixedEffectsOnly=TRUE,
                                  getSDs=TRUE, doMCMC=FALSE)
    save(result_opt_FE, file=saveFileFE)
    cat("\nFixed-effects-only optimization completed and saved.\n")
  }, error = function(e) {
    cat("\nFixed-effects-only optimization FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("\nFixed-effects-only optimization already completed (loaded from save).\n")
}

if(!is.null(result_opt_FE)) {
  cat("  Total time:", result_opt_FE$totalTime/60, "minutes\n")
  
  SD0 = result_opt_FE$TMBsd
  if(inherits(SD0, "sdreport")) {
    pf = SD0$par.fixed
    pfn = names(pf)
    
    cat("\nFixed-effects-only model estimates vs truth:\n")
    alpha_est = pf[pfn == "alpha"]
    beta_est = pf[pfn == "beta"]
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "alpha", -1.25, as.numeric(alpha_est)))
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "urban", 1.00, as.numeric(beta_est[1])))
    cat(sprintf("  %-18s truth=%6.2f  est=%7.4f\n", "normPop", 0.50, as.numeric(beta_est[2])))
    
    log_tauEps_est = pf[pfn == "log_tauEps"]
    if(length(log_tauEps_est) > 0) {
      sigmaEps_est = exp(-0.5 * log_tauEps_est)
      cat(sprintf("  %-18s truth=%6.4f  est=%7.4f\n", "sigmaEps", sqrt(1.5), as.numeric(sigmaEps_est)))
    }
  }
}

cat("\n===============================================\n")
cat("ALL TESTS COMPLETE (including fixed-effects-only)\n")
cat("===============================================\n")
