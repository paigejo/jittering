#!/usr/bin/env Rscript
# Test IID spatial + nugget + 2 covariates (urban, normPop) + intercept model
# on data simulated from BYM2 model (simData1BYM2)
#
# True parameters used in BYM2 simulation:
#   beta0 (intercept) = -1.25
#   gamma (urban)      =  1.00
#   betaRest           = c(0, 0, 0, 0.5) => normPop = 0.5
#   sigmaBYM2          = sqrt(0.5)  (spatial SD)
#   sigmaEpsilon       = sqrt(1.5)  (nugget SD)
#   phi                = 0.8        (BYM2 mixing)

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

saveFile = "savedOutput/testMM_BYM2sim_IIDnugget_progress.RData"

cat("===============================================\n")
cat("IID+NUGGET MODEL ON BYM2-SIMULATED DATA\n")
cat("===============================================\n\n")

# Load BYM2 simulated data
cat("Loading BYM2 simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")

# Use first simulated dataset
ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

cat("Data summary:\n")
cat("  DHS clusters:", nrow(ed), "\n")
cat("  MICS clusters:", nrow(edMICS), "\n")
cat("  DHS total obs:", sum(ed$N), "\n")
cat("  MICS total obs:", sum(edMICS$N), "\n\n")

# Load DHS integration points (same geometry as SPDE sim study)
cat("Loading DHS integration points...\n")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaU = sqrt(0.5)
trueSigmaEps = sqrt(1.5)

# Load previous results if available
result_opt = NULL
if(file.exists(saveFile)) {
  cat("Loading saved results from", saveFile, "\n")
  load(saveFile)
}

# ---- Fit: IID + nugget, 2 covariates (urban + normPop) ----
if(is.null(result_opt)) {
  cat("\nFitting IID + nugget model (urban + normPop) via TMB optimization...\n")
  cat("-------------------------------------------\n")
  
  startTime = proc.time()[3]
  tryCatch({
    result_opt = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                               covariates=c("urban", "normPop"),
                               getSDs=TRUE, doMCMC=FALSE)
    save(result_opt, file=saveFile)
    cat("\nOptimization completed and saved.\n")
  }, error = function(e) {
    cat("\nOptimization FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("Optimization already completed (loaded from save).\n")
}

# ---- Print results ----
if(!is.null(result_opt)) {
  cat("\n  Total time:", round(result_opt$totalTime/60, 2), "minutes\n")
  
  SD0 = result_opt$TMBsd
  if(inherits(SD0, "sdreport")) {
    pf = SD0$par.fixed
    pfn = names(pf)
    pr = SD0$par.random
    rn = names(pr)
    se_fixed = sqrt(diag(SD0$cov.fixed))
    
    # Extract estimates
    alpha_est = as.numeric(pr[rn == "alpha"])
    beta_est = as.numeric(pr[rn == "beta"])
    
    log_tau_est = as.numeric(pf[pfn == "log_tau"])
    sigma_u_est = exp(-0.5 * log_tau_est)
    
    log_tauEps_est = as.numeric(pf[pfn == "log_tauEps"])
    sigmaEps_est = exp(-0.5 * log_tauEps_est)
    
    logit_phi_est = pf[pfn == "logit_phi"]
    
    # Standard errors for random effects
    se_rand = sqrt(diag(SD0$jointPrecision))
    # Actually, get SEs from sdreport summary
    summ = summary(SD0, select="random")
    alpha_se = summ[rownames(summ) == "alpha" | names(summ[,1]) == "alpha", 2]
    if(length(alpha_se) == 0) {
      # try alternative
      alpha_se = NA
    }
    
    cat("\n===============================================\n")
    cat("PARAMETER ESTIMATES vs TRUTH (BYM2 sim data)\n")
    cat("===============================================\n\n")
    
    cat("Note: True DGP is BYM2 (phi=0.8, 41 areas), model is IID spatial\n\n")
    
    cat(sprintf("%-18s %10s %10s %10s\n", "Parameter", "Truth", "Estimate", "Bias"))
    cat(sprintf("%-18s %10s %10s %10s\n", "---------", "-----", "--------", "----"))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f\n", "alpha (intercept)", trueAlpha, alpha_est, alpha_est - trueAlpha))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f\n", "beta[urban]", trueUrban, beta_est[1], beta_est[1] - trueUrban))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f\n", "beta[normPop]", trueNormPop, beta_est[2], beta_est[2] - trueNormPop))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f\n", "sigma_u (spatial)", trueSigmaU, sigma_u_est, sigma_u_est - trueSigmaU))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f\n", "sigmaEps (nugget)", trueSigmaEps, sigmaEps_est, sigmaEps_est - trueSigmaEps))
    
    # Log-scale parameters
    cat("\nLog-scale fixed parameters:\n")
    cat(sprintf("  log_tau    = %7.4f (SE = %7.4f)\n", log_tau_est, se_fixed[pfn == "log_tau"]))
    cat(sprintf("  log_tauEps = %7.4f (SE = %7.4f)\n", log_tauEps_est, se_fixed[pfn == "log_tauEps"]))
    if(length(logit_phi_est) > 0) {
      cat(sprintf("  logit_phi  = %7.4f (SE = %7.4f)\n", logit_phi_est, se_fixed[pfn == "logit_phi"]))
    }
    
    # Random effects summary
    cat("\nRandom effects counts:\n")
    print(table(rn))
    
    # Spatial effects summary
    u_spatial = pr[rn == "u_spatial"]
    cat(sprintf("\nu_spatial: mean=%7.4f, sd=%7.4f, range=[%7.4f, %7.4f]\n",
                mean(u_spatial), sd(u_spatial), min(u_spatial), max(u_spatial)))
    
    nuggetUrb = pr[rn == "nuggetUrbMICS"]
    nuggetRur = pr[rn == "nuggetRurMICS"]
    cat(sprintf("nuggetUrb: mean=%7.4f, sd=%7.4f (n=%d)\n", mean(nuggetUrb), sd(nuggetUrb), length(nuggetUrb)))
    cat(sprintf("nuggetRur: mean=%7.4f, sd=%7.4f (n=%d)\n", mean(nuggetRur), sd(nuggetRur), length(nuggetRur)))
    
  } else {
    cat("\nsdreport failed or not available.\n")
  }
}

cat("\n===============================================\n")
cat("TEST COMPLETE\n")
cat("===============================================\n")
