#!/usr/bin/env Rscript
# Test BYM2 model on simulated BYM2 data (MICS only, no TMBstan)
# Fix log_tau and logit_phi to true values; only estimate log_tauEps (nugget)
# Run on sim 1 and sim 2.

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modBYM2.R")
source("code/modM_MSepMarg.R")

saveFile1 = "savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim1.RData"
saveFile2 = "savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim2.RData"

cat("===============================================\n")
cat("BYM2 MODEL WITH FIXED HYPERPARAMS ON BYM2 SIM\n")
cat("===============================================\n\n")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaU = sqrt(0.5)
trueSigmaEps = sqrt(1.5)
truePhi = 0.8

# Compute true hyperparameter values on TMB scale
true_log_tau = log(1 / trueSigmaU^2)      # log(2)
true_logit_phi = log(truePhi / (1-truePhi)) # log(4)

cat(sprintf("True hyperparams: log_tau=%.4f, logit_phi=%.4f\n", true_log_tau, true_logit_phi))
cat(sprintf("  (sigma_u=%.4f, phi=%.4f)\n\n", trueSigmaU, truePhi))

# Load BYM2 simulated data
cat("Loading BYM2 simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")

# ============================================================
# SIM 1
# ============================================================
result_sim1 = NULL
if(file.exists(saveFile1)) {
  cat("Loading saved sim 1 results...\n")
  load(saveFile1)
} else {
  edMICS1 = surveysMICS[[1]]
  cat(sprintf("Sim 1: %d MICS clusters, %d total obs\n", nrow(edMICS1), sum(edMICS1$N)))
  
  cat("\n--- Fitting BYM2 (fixed hyper) on sim 1 ---\n")
  tryCatch({
    result_sim1 = fitMMMarg(datMICS=edMICS1,
                            covariates=c("urban", "normPop"),
                            fixedHyperparams=list(log_tau=true_log_tau,
                                                  logit_phi=true_logit_phi),
                            getSDs=TRUE, doMCMC=FALSE)
    save(result_sim1, file=saveFile1)
    cat("Sim 1 optimization completed and saved.\n")
  }, error = function(e) {
    cat("Sim 1 FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# SIM 2
# ============================================================
result_sim2 = NULL
if(file.exists(saveFile2)) {
  cat("Loading saved sim 2 results...\n")
  load(saveFile2)
} else {
  edMICS2 = surveysMICS[[2]]
  cat(sprintf("Sim 2: %d MICS clusters, %d total obs\n", nrow(edMICS2), sum(edMICS2$N)))
  
  cat("\n--- Fitting BYM2 (fixed hyper) on sim 2 ---\n")
  tryCatch({
    result_sim2 = fitMMMarg(datMICS=edMICS2,
                            covariates=c("urban", "normPop"),
                            fixedHyperparams=list(log_tau=true_log_tau,
                                                  logit_phi=true_logit_phi),
                            getSDs=TRUE, doMCMC=FALSE)
    save(result_sim2, file=saveFile2)
    cat("Sim 2 optimization completed and saved.\n")
  }, error = function(e) {
    cat("Sim 2 FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# Extract estimates
# ============================================================
extractEst = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed
  pfn = names(pf)
  se_f = sqrt(diag(SD0$cov.fixed))
  
  # alpha and beta are random effects in this model
  sr = summary(SD0, select="random")
  rn = names(SD0$par.random)
  alpha_est = sr[rn == "alpha", 1]
  alpha_se = sr[rn == "alpha", 2]
  beta_est = sr[rn == "beta", 1]
  beta_se = sr[rn == "beta", 2]
  
  # sigmaEps from log_tauEps (the only outer parameter)
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  
  list(alpha=alpha_est, alpha_se=alpha_se,
       beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2],
       sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

# ============================================================
# Print results
# ============================================================
cat("\n===============================================\n")
cat("RESULTS: BYM2 (fixed log_tau, logit_phi)\n")
cat("===============================================\n\n")

fmtRow = function(name, truth, est1, se1, est2, se2) {
  if(is.na(est1)) { s1 = sprintf("%10s", "NA") } else { s1 = sprintf("%7.4f (SE %5.4f)", est1, se1) }
  if(is.na(est2)) { s2 = sprintf("%10s", "NA") } else { s2 = sprintf("%7.4f (SE %5.4f)", est2, se2) }
  sprintf("%-18s %7.4f  %s  %s", name, truth, s1, s2)
}

if(!is.null(result_sim1) && !is.null(result_sim2)) {
  e1 = extractEst(result_sim1)
  e2 = extractEst(result_sim2)
  
  cat(sprintf("Sim 1 time: %.1f min | Sim 2 time: %.1f min\n\n", e1$time/60, e2$time/60))
  
  cat(sprintf("%-18s %7s  %22s  %22s\n", "Parameter", "Truth", "Sim 1", "Sim 2"))
  cat(sprintf("%-18s %7s  %22s  %22s\n", "---------", "-----", "-----", "-----"))
  cat(fmtRow("alpha", trueAlpha, e1$alpha, e1$alpha_se, e2$alpha, e2$alpha_se), "\n")
  cat(fmtRow("beta[urban]", trueUrban, e1$beta1, e1$beta1_se, e2$beta1, e2$beta1_se), "\n")
  cat(fmtRow("beta[normPop]", trueNormPop, e1$beta2, e1$beta2_se, e2$beta2, e2$beta2_se), "\n")
  cat(fmtRow("sigmaEps", trueSigmaEps, e1$sigmaEps, e1$sigmaEps_se, e2$sigmaEps, e2$sigmaEps_se), "\n")
  
  cat("\nBias summary:\n")
  cat(sprintf("  %-18s %+8.4f  %+8.4f\n", "alpha", e1$alpha - trueAlpha, e2$alpha - trueAlpha))
  cat(sprintf("  %-18s %+8.4f  %+8.4f\n", "beta[urban]", e1$beta1 - trueUrban, e2$beta1 - trueUrban))
  cat(sprintf("  %-18s %+8.4f  %+8.4f\n", "beta[normPop]", e1$beta2 - trueNormPop, e2$beta2 - trueNormPop))
  cat(sprintf("  %-18s %+8.4f  %+8.4f\n", "sigmaEps", e1$sigmaEps - trueSigmaEps, e2$sigmaEps - trueSigmaEps))
} else {
  if(!is.null(result_sim1)) {
    e1 = extractEst(result_sim1)
    cat("Sim 1 results:\n")
    cat(sprintf("  alpha: %.4f (SE %.4f)\n", e1$alpha, e1$alpha_se))
    cat(sprintf("  beta[urban]: %.4f (SE %.4f)\n", e1$beta1, e1$beta1_se))
    cat(sprintf("  beta[normPop]: %.4f (SE %.4f)\n", e1$beta2, e1$beta2_se))
    cat(sprintf("  sigmaEps: %.4f (SE %.4f)\n", e1$sigmaEps, e1$sigmaEps_se))
  }
  if(!is.null(result_sim2)) {
    e2 = extractEst(result_sim2)
    cat("Sim 2 results:\n")
    cat(sprintf("  alpha: %.4f (SE %.4f)\n", e2$alpha, e2$alpha_se))
    cat(sprintf("  beta[urban]: %.4f (SE %.4f)\n", e2$beta1, e2$beta1_se))
    cat(sprintf("  beta[normPop]: %.4f (SE %.4f)\n", e2$beta2, e2$beta2_se))
    cat(sprintf("  sigmaEps: %.4f (SE %.4f)\n", e2$sigmaEps, e2$sigmaEps_se))
  }
}

cat("\n===============================================\n")
cat("TEST COMPLETE\n")
cat("===============================================\n")
