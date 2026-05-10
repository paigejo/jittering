#!/usr/bin/env Rscript
# Test both models on SECOND simulated dataset from BYM2 simulation
# 1) IID spatial + nugget + 2 covariates
# 2) Covariate + nugget only (no spatial RE)
# Then compare with sim 1 results.

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

saveFileIID = "savedOutput/testMM_BYM2sim2_IIDnugget.RData"
saveFileNug = "savedOutput/testMM_BYM2sim2_nuggetOnly.RData"

cat("===============================================\n")
cat("MODELS ON BYM2-SIMULATED DATA: SIM 2\n")
cat("===============================================\n\n")

# Load BYM2 simulated data
cat("Loading BYM2 simulated data...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")

# Use SECOND simulated dataset (MICS only)
edMICS = surveysMICS[[2]]

cat("Data summary (sim 2):\n")
cat("  MICS clusters:", nrow(edMICS), "\n")
cat("  MICS total obs:", sum(edMICS$N), "\n\n")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaU = sqrt(0.5)
trueSigmaEps = sqrt(1.5)

# ============================================================
# MODEL 1: IID spatial + nugget
# ============================================================
result_iid = NULL
if(file.exists(saveFileIID)) {
  cat("Loading saved IID+nugget results...\n")
  load(saveFileIID)
} else {
  cat("\n--- Fitting IID + nugget model (sim 2) ---\n")
  startTime = proc.time()[3]
  tryCatch({
    result_iid = fitMM_IIDonly(datMICS=edMICS,
                                covariates=c("urban", "normPop"),
                                fixedEffectsOnly=FALSE,
                                getSDs=TRUE, doMCMC=FALSE)
    save(result_iid, file=saveFileIID)
    cat("IID+nugget optimization completed and saved.\n")
  }, error = function(e) {
    cat("IID+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# MODEL 2: Covariate + nugget only
# ============================================================
result_nug = NULL
if(file.exists(saveFileNug)) {
  cat("Loading saved nugget-only results...\n")
  load(saveFileNug)
} else {
  cat("\n--- Fitting nugget-only model (sim 2) ---\n")
  startTime = proc.time()[3]
  tryCatch({
    result_nug = fitMM_IIDonly(datMICS=edMICS,
                                covariates=c("urban", "normPop"),
                                fixedEffectsOnly=TRUE,
                                getSDs=TRUE, doMCMC=FALSE)
    save(result_nug, file=saveFileNug)
    cat("Nugget-only optimization completed and saved.\n")
  }, error = function(e) {
    cat("Nugget-only FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# Extract estimates helper
# ============================================================
extractEst = function(res, isNuggetOnly=FALSE) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed
  pfn = names(pf)
  se_f = sqrt(diag(SD0$cov.fixed))
  
  if(isNuggetOnly) {
    alpha_est = as.numeric(pf[pfn == "alpha"])
    alpha_se = as.numeric(se_f[pfn == "alpha"])
    beta_est = as.numeric(pf[pfn == "beta"])
    beta_se = as.numeric(se_f[pfn == "beta"])
    sigma_u_est = NA; sigma_u_se = NA
  } else {
    sr = summary(SD0, select="random")
    rn = names(SD0$par.random)
    alpha_est = sr[rn == "alpha", 1]
    alpha_se = sr[rn == "alpha", 2]
    beta_est = sr[rn == "beta", 1]
    beta_se = sr[rn == "beta", 2]
    lt = pf[pfn == "log_tau"]; se_lt = se_f[pfn == "log_tau"]
    sigma_u_est = exp(-0.5 * lt); sigma_u_se = 0.5 * sigma_u_est * se_lt
  }
  
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  
  list(alpha=alpha_est, alpha_se=alpha_se,
       beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2],
       sigma_u=sigma_u_est, sigma_u_se=sigma_u_se,
       sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

# ============================================================
# Print sim 2 results
# ============================================================
cat("\n===============================================\n")
cat("SIM 2 RESULTS\n")
cat("===============================================\n\n")

if(!is.null(result_iid)) {
  e_iid = extractEst(result_iid, FALSE)
  cat(sprintf("IID+nugget time: %.1f min\n", e_iid$time/60))
}
if(!is.null(result_nug)) {
  e_nug = extractEst(result_nug, TRUE)
  cat(sprintf("Nugget-only time: %.1f min\n", e_nug$time/60))
}

fmtRow = function(name, truth, est, se, est2, se2) {
  if(is.na(est)) {
    s1 = sprintf("%10s", "NA")
  } else {
    s1 = sprintf("%7.4f (SE %5.4f)", est, se)
  }
  if(is.na(est2)) {
    s2 = sprintf("%10s", "NA")
  } else {
    s2 = sprintf("%7.4f (SE %5.4f)", est2, se2)
  }
  sprintf("%-18s %7.4f  %s  %s", name, truth, s1, s2)
}

if(!is.null(result_iid) && !is.null(result_nug)) {
  cat(sprintf("\n%-18s %7s  %22s  %22s\n", "Parameter", "Truth", "IID+Nugget", "Nugget-only"))
  cat(sprintf("%-18s %7s  %22s  %22s\n", "---------", "-----", "----------", "-----------"))
  cat(fmtRow("alpha", trueAlpha, e_iid$alpha, e_iid$alpha_se, e_nug$alpha, e_nug$alpha_se), "\n")
  cat(fmtRow("beta[urban]", trueUrban, e_iid$beta1, e_iid$beta1_se, e_nug$beta1, e_nug$beta1_se), "\n")
  cat(fmtRow("beta[normPop]", trueNormPop, e_iid$beta2, e_iid$beta2_se, e_nug$beta2, e_nug$beta2_se), "\n")
  cat(fmtRow("sigma_u", trueSigmaU, e_iid$sigma_u, e_iid$sigma_u_se, e_nug$sigma_u, e_nug$sigma_u_se), "\n")
  cat(fmtRow("sigmaEps", trueSigmaEps, e_iid$sigmaEps, e_iid$sigmaEps_se, e_nug$sigmaEps, e_nug$sigmaEps_se), "\n")
}

# ============================================================
# Compare with sim 1 results
# ============================================================
cat("\n===============================================\n")
cat("COMPARISON: SIM 1 vs SIM 2\n")
cat("===============================================\n\n")

sim1_iid_file = "savedOutput/testMM_BYM2sim_IIDnugget_progress.RData"
sim1_nug_file = "savedOutput/testMM_BYM2sim_nuggetOnly_progress.RData"

if(file.exists(sim1_iid_file) && file.exists(sim1_nug_file)) {
  load(sim1_iid_file)  # result_opt
  load(sim1_nug_file)  # result_nugget
  e_iid1 = extractEst(result_opt, FALSE)
  e_nug1 = extractEst(result_nugget, TRUE)
  
  fmtComp = function(name, truth, e1, se1, e2, se2) {
    b1 = e1 - truth; b2 = e2 - truth
    sprintf("%-18s %7.4f  %7.4f (%+6.3f)  %7.4f (%+6.3f)", name, truth, e1, b1, e2, b2)
  }
  
  # IID+Nugget across sims
  cat("--- IID+Nugget model ---\n")
  cat(sprintf("%-18s %7s  %18s  %18s\n", "Parameter", "Truth", "Sim1 (bias)", "Sim2 (bias)"))
  cat(sprintf("%-18s %7s  %18s  %18s\n", "---------", "-----", "-----------", "-----------"))
  cat(fmtComp("alpha", trueAlpha, e_iid1$alpha, e_iid1$alpha_se, e_iid$alpha, e_iid$alpha_se), "\n")
  cat(fmtComp("beta[urban]", trueUrban, e_iid1$beta1, e_iid1$beta1_se, e_iid$beta1, e_iid$beta1_se), "\n")
  cat(fmtComp("beta[normPop]", trueNormPop, e_iid1$beta2, e_iid1$beta2_se, e_iid$beta2, e_iid$beta2_se), "\n")
  cat(fmtComp("sigma_u", trueSigmaU, e_iid1$sigma_u, e_iid1$sigma_u_se, e_iid$sigma_u, e_iid$sigma_u_se), "\n")
  cat(fmtComp("sigmaEps", trueSigmaEps, e_iid1$sigmaEps, e_iid1$sigmaEps_se, e_iid$sigmaEps, e_iid$sigmaEps_se), "\n")
  
  cat("\n--- Nugget-only model ---\n")
  cat(sprintf("%-18s %7s  %18s  %18s\n", "Parameter", "Truth", "Sim1 (bias)", "Sim2 (bias)"))
  cat(sprintf("%-18s %7s  %18s  %18s\n", "---------", "-----", "-----------", "-----------"))
  cat(fmtComp("alpha", trueAlpha, e_nug1$alpha, e_nug1$alpha_se, e_nug$alpha, e_nug$alpha_se), "\n")
  cat(fmtComp("beta[urban]", trueUrban, e_nug1$beta1, e_nug1$beta1_se, e_nug$beta1, e_nug$beta1_se), "\n")
  cat(fmtComp("beta[normPop]", trueNormPop, e_nug1$beta2, e_nug1$beta2_se, e_nug$beta2, e_nug$beta2_se), "\n")
  cat(fmtComp("sigmaEps", trueSigmaEps, e_nug1$sigmaEps, e_nug1$sigmaEps_se, e_nug$sigmaEps, e_nug$sigmaEps_se), "\n")
} else {
  cat("Sim 1 results not found.\n")
}

cat("\n===============================================\n")
cat("TEST COMPLETE\n")
cat("===============================================\n")
