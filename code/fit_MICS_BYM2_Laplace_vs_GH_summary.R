#!/usr/bin/env Rscript
# MICS-only BYM2 comparison: Laplace nugget integration vs GH nugget integration

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MSep.R")

load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

# Ensure expected names
nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]
  toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}

truePhi = 0.8
trueSigmaBYM2 = sqrt(0.5)
trueSigmaEps = sqrt(1.5)

extract_hyper = function(outObj) {
  if(is.null(outObj$TMBsd) || !inherits(outObj$TMBsd, "sdreport")) {
    return(list(ok=FALSE, msg="No sdreport"))
  }
  fix = summary(outObj$TMBsd, "fixed")
  rn = rownames(fix)

  get_est_se = function(paramName) {
    idx = which(rn == paramName)
    if(length(idx) == 0) return(c(NA_real_, NA_real_))
    c(fix[idx[1], "Estimate"], fix[idx[1], "Std. Error"])
  }

  lt = get_est_se("log_tau")
  lp = get_est_se("logit_phi")
  le = get_est_se("log_tauEps")

  log_tau = lt[1]; se_log_tau = lt[2]
  logit_phi = lp[1]; se_logit_phi = lp[2]
  log_tauEps = le[1]; se_log_tauEps = le[2]

  sigmaBYM2 = exp(-0.5 * log_tau)
  se_sigmaBYM2 = 0.5 * sigmaBYM2 * se_log_tau

  phi = 1/(1 + exp(-logit_phi))
  se_phi = phi * (1 - phi) * se_logit_phi

  sigmaEps = exp(-0.5 * log_tauEps)
  se_sigmaEps = 0.5 * sigmaEps * se_log_tauEps

  list(
    ok=TRUE,
    log_tau=log_tau, se_log_tau=se_log_tau,
    logit_phi=logit_phi, se_logit_phi=se_logit_phi,
    log_tauEps=log_tauEps, se_log_tauEps=se_log_tauEps,
    sigmaBYM2=sigmaBYM2, se_sigmaBYM2=se_sigmaBYM2,
    phi=phi, se_phi=se_phi,
    sigmaEps=sigmaEps, se_sigmaEps=se_sigmaEps
  )
}

cat("\n============================================================\n")
cat("MICS-only BYM2: Laplace vs GH nugget integration\n")
cat("============================================================\n")

# Laplace (old)
cat("\nRunning Laplace (fitMM_old)...\n")
lap = fitMM_old(
  datDHS=ed, datMICS=edMICS,
  getSDs=TRUE,
  doMCMC=FALSE,
  maxit=1000
)

# GH (new)
cat("\nRunning GH (fitMM, Q=10)...\n")
gh = fitMM(
  datDHS=ed, datMICS=edMICS,
  Qgh=10,
  getSDs=TRUE,
  maxit=1000
)

lap_h = extract_hyper(lap)
gh_h = extract_hyper(gh)

lap_conv = if(!is.null(lap$opt$convergence)) lap$opt$convergence else NA
lap_nll = if(!is.null(lap$opt$value)) lap$opt$value else NA
lap_time = if(!is.null(lap$totalTime)) lap$totalTime else NA

gh_conv = if(!is.null(gh$opt$convergence)) gh$opt$convergence else NA
gh_nll = if(!is.null(gh$opt$value)) gh$opt$value else NA
gh_time = if(!is.null(gh$totalTime)) gh$totalTime else NA

cat("\n================================================================\n")
cat("MICS-only BYM2: Laplace vs GH (with uncertainty)\n")
cat("================================================================\n")
cat(sprintf("%-22s %10s %8s %10s %8s %10s\n", "Parameter", "Laplace", "(SE)", "GH", "(SE)", "Truth"))
cat(sprintf("%-22s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "phi", lap_h$phi, lap_h$se_phi, gh_h$phi, gh_h$se_phi, truePhi))
cat(sprintf("%-22s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaBYM2", lap_h$sigmaBYM2, lap_h$se_sigmaBYM2, gh_h$sigmaBYM2, gh_h$se_sigmaBYM2, trueSigmaBYM2))
cat(sprintf("%-22s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaEps", lap_h$sigmaEps, lap_h$se_sigmaEps, gh_h$sigmaEps, gh_h$se_sigmaEps, trueSigmaEps))
cat(sprintf("%-22s %10.4f %8.4f %10.4f %8.4f\n",
            "logit_phi", lap_h$logit_phi, lap_h$se_logit_phi, gh_h$logit_phi, gh_h$se_logit_phi))

cat(sprintf("\n%-22s %10.0f %8s %10.0f\n", "Convergence", lap_conv, "", gh_conv))
cat(sprintf("%-22s %10.4f %8s %10.4f\n", "NLL", lap_nll, "", gh_nll))
cat(sprintf("%-22s %10.1f %8s %10.1f\n", "Time (s)", lap_time, "", gh_time))

# z-style distance from truth
z_lap_phi = abs(lap_h$phi - truePhi) / lap_h$se_phi
z_gh_phi = abs(gh_h$phi - truePhi) / gh_h$se_phi
cat("\n--- Phi accuracy in SE units ---\n")
cat(sprintf("Laplace: |phi - truth| / SE = %.2f\n", z_lap_phi))
cat(sprintf("GH:      |phi - truth| / SE = %.2f\n", z_gh_phi))

save(lap, gh, lap_h, gh_h,
  file="savedOutput/fitMICS_BYM2_Laplace_vs_GH_summary.RData")
cat("\nSaved to savedOutput/fitMICS_BYM2_Laplace_vs_GH_summary.RData\n")
cat("Done.\n")
