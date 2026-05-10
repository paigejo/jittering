#!/usr/bin/env Rscript
# Test covariate + nugget ONLY model (no spatial random effects)
# on data simulated from BYM2 model (simData1BYM2)
#
# Uses fixedEffectsOnly=TRUE which maps out u_spatial and log_tau,
# so alpha and beta become fixed (not random) effects optimized directly.
#
# True parameters used in BYM2 simulation:
#   beta0 (intercept) = -1.25
#   gamma (urban)      =  1.00
#   betaRest           = c(0, 0, 0, 0.5) => normPop = 0.5
#   sigmaEpsilon       = sqrt(1.5)  (nugget SD)

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

saveFile = "savedOutput/testMM_BYM2sim_nuggetOnly_progress.RData"

cat("===============================================\n")
cat("COVARIATE + NUGGET ONLY MODEL ON BYM2-SIM DATA\n")
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

# Load DHS integration points
cat("Loading DHS integration points...\n")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

# Load previous results if available
result_nugget = NULL
if(file.exists(saveFile)) {
  cat("Loading saved results from", saveFile, "\n")
  load(saveFile)
}

# ---- Fit: covariates + nugget only (no spatial RE) ----
if(is.null(result_nugget)) {
  cat("\nFitting covariate + nugget only model (urban + normPop) via TMB...\n")
  cat("  fixedEffectsOnly=TRUE => no spatial random effects\n")
  cat("-------------------------------------------\n")
  
  startTime = proc.time()[3]
  tryCatch({
    result_nugget = fitMM_IIDonly(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS,
                                  covariates=c("urban", "normPop"),
                                  fixedEffectsOnly=TRUE,
                                  getSDs=TRUE, doMCMC=FALSE)
    save(result_nugget, file=saveFile)
    cat("\nOptimization completed and saved.\n")
  }, error = function(e) {
    cat("\nOptimization FAILED with error:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("Optimization already completed (loaded from save).\n")
}

# ---- Print results ----
if(!is.null(result_nugget)) {
  cat("\n  Total time:", round(result_nugget$totalTime/60, 2), "minutes\n")
  
  SD0 = result_nugget$TMBsd
  if(inherits(SD0, "sdreport")) {
    pf = SD0$par.fixed
    pfn = names(pf)
    se_fixed = sqrt(diag(SD0$cov.fixed))
    
    # Fixed effects: alpha, beta are now in par.fixed
    alpha_est = as.numeric(pf[pfn == "alpha"])
    beta_est = as.numeric(pf[pfn == "beta"])
    alpha_se = as.numeric(se_fixed[pfn == "alpha"])
    beta_se = as.numeric(se_fixed[pfn == "beta"])
    
    log_tauEps_est = as.numeric(pf[pfn == "log_tauEps"])
    log_tauEps_se = as.numeric(se_fixed[pfn == "log_tauEps"])
    sigmaEps_est = exp(-0.5 * log_tauEps_est)
    sigmaEps_se = 0.5 * sigmaEps_est * log_tauEps_se  # delta method
    
    cat("\n===============================================\n")
    cat("PARAMETER ESTIMATES vs TRUTH (BYM2 sim data)\n")
    cat("Model: Covariate + Nugget ONLY (no spatial RE)\n")
    cat("===============================================\n\n")
    
    cat("Note: True DGP is BYM2 (phi=0.8, 41 areas), model has NO spatial effects\n\n")
    
    cat(sprintf("%-18s %10s %10s %10s %10s\n", "Parameter", "Truth", "Estimate", "SE", "Bias"))
    cat(sprintf("%-18s %10s %10s %10s %10s\n", "---------", "-----", "--------", "------", "----"))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f %10.4f\n", "alpha (intercept)", trueAlpha, alpha_est, alpha_se, alpha_est - trueAlpha))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f %10.4f\n", "beta[urban]", trueUrban, beta_est[1], beta_se[1], beta_est[1] - trueUrban))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f %10.4f\n", "beta[normPop]", trueNormPop, beta_est[2], beta_se[2], beta_est[2] - trueNormPop))
    cat(sprintf("%-18s %10.4f %10.4f %10.4f %10.4f\n", "sigmaEps (nugget)", trueSigmaEps, sigmaEps_est, sigmaEps_se, sigmaEps_est - trueSigmaEps))
    
    # Also print log-scale
    cat("\nLog-scale fixed parameters:\n")
    cat(sprintf("  log_tauEps = %7.4f (SE = %7.4f)\n", log_tauEps_est, log_tauEps_se))
    
    # Random effects summary (only nuggets)
    pr = SD0$par.random
    rn = names(pr)
    cat("\nRandom effects counts:\n")
    print(table(rn))
    
    nuggetUrb = pr[rn == "nuggetUrbMICS"]
    nuggetRur = pr[rn == "nuggetRurMICS"]
    cat(sprintf("\nnuggetUrb: mean=%7.4f, sd=%7.4f (n=%d)\n", mean(nuggetUrb), sd(nuggetUrb), length(nuggetUrb)))
    cat(sprintf("nuggetRur: mean=%7.4f, sd=%7.4f (n=%d)\n", mean(nuggetRur), sd(nuggetRur), length(nuggetRur)))
    
    # ---- Compare to IID+nugget model ----
    cat("\n===============================================\n")
    cat("COMPARISON: Nugget-only vs IID+Nugget model\n")
    cat("===============================================\n\n")
    
    iidFile = "savedOutput/testMM_BYM2sim_IIDnugget_progress.RData"
    if(file.exists(iidFile)) {
      load(iidFile)  # loads result_opt
      SD_iid = result_opt$TMBsd
      pf_iid = SD_iid$par.fixed
      pfn_iid = names(pf_iid)
      se_iid = sqrt(diag(SD_iid$cov.fixed))
      
      # IID model: alpha, beta are random
      sr_iid = summary(SD_iid, select="random")
      rn_iid = names(SD_iid$par.random)
      alpha_iid = sr_iid[rn_iid == "alpha", 1]
      alpha_iid_se = sr_iid[rn_iid == "alpha", 2]
      beta_iid = sr_iid[rn_iid == "beta", 1]
      beta_iid_se = sr_iid[rn_iid == "beta", 2]
      
      lt_iid = pf_iid[pfn_iid == "log_tau"]
      sigma_u_iid = exp(-0.5 * lt_iid)
      sigma_u_iid_se = 0.5 * sigma_u_iid * se_iid[pfn_iid == "log_tau"]
      
      le_iid = pf_iid[pfn_iid == "log_tauEps"]
      sigmaEps_iid = exp(-0.5 * le_iid)
      sigmaEps_iid_se = 0.5 * sigmaEps_iid * se_iid[pfn_iid == "log_tauEps"]
      
      cat(sprintf("%-18s %10s %22s %22s\n", "Parameter", "Truth", "Nugget-only", "IID+Nugget"))
      cat(sprintf("%-18s %10s %22s %22s\n", "---------", "-----", "-----------", "----------"))
      cat(sprintf("%-18s %10.4f %10.4f (SE %5.4f) %10.4f (SE %5.4f)\n", "alpha", trueAlpha, alpha_est, alpha_se, alpha_iid, alpha_iid_se))
      cat(sprintf("%-18s %10.4f %10.4f (SE %5.4f) %10.4f (SE %5.4f)\n", "beta[urban]", trueUrban, beta_est[1], beta_se[1], beta_iid[1], beta_iid_se[1]))
      cat(sprintf("%-18s %10.4f %10.4f (SE %5.4f) %10.4f (SE %5.4f)\n", "beta[normPop]", trueNormPop, beta_est[2], beta_se[2], beta_iid[2], beta_iid_se[2]))
      cat(sprintf("%-18s %10.4f %10s %16.4f (SE %5.4f)\n", "sigma_u", sqrt(0.5), "    NA", sigma_u_iid, sigma_u_iid_se))
      cat(sprintf("%-18s %10.4f %10.4f (SE %5.4f) %10.4f (SE %5.4f)\n", "sigmaEps", trueSigmaEps, sigmaEps_est, sigmaEps_se, sigmaEps_iid, sigmaEps_iid_se))
    } else {
      cat("IID+Nugget results not found at", iidFile, "\n")
    }
  } else {
    cat("\nsdreport failed or not available.\n")
  }
}

cat("\n===============================================\n")
cat("TEST COMPLETE\n")
cat("===============================================\n")
