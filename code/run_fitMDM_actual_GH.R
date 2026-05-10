#!/usr/bin/env Rscript
# Run DHS+MICS BYM2 GH fusion model on actual survey data (defaults ed/edMICS)

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_DMSep.R")

cat("\n============================================================\n")
cat("Running fitMDM() on actual survey data (defaults)\n")
cat("GH fusion, 2n-2 parameterization, default K settings\n")
cat("============================================================\n")

options(error = function() {
  cat("\n--- traceback() ---\n")
  traceback(2)
  quit(status = 1)
})

start = proc.time()[3]
out = tryCatch(
  fitMDM(
    getSDs = TRUE,
    maxit = 1000,
    Qgh = 10
  ),
  error = function(e) {
    cat("\nfitMDM() error:\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  }
)
elapsed = proc.time()[3] - start

conv = NA
nll = NA
if(!is.null(out$opt)) {
  if(!is.null(out$opt$convergence)) conv = out$opt$convergence
  if(!is.null(out$opt$value)) nll = out$opt$value
}

cat(sprintf("Convergence: %s\n", as.character(conv)))
cat(sprintf("NLL: %s\n", as.character(nll)))
cat(sprintf("Total time (s): %.1f\n", elapsed))

if(!is.null(out$TMBsd) && inherits(out$TMBsd, "sdreport")) {
  fix = summary(out$TMBsd, "fixed")
  rn = rownames(fix)

  idx_tau = which(rn == "log_tau")
  idx_phi = which(rn == "logit_phi")
  idx_eps = which(rn == "log_tauEps")

  if(length(idx_tau) > 0) {
    log_tau = fix[idx_tau[1], "Estimate"]
    se_log_tau = fix[idx_tau[1], "Std. Error"]
    sigmaBYM2 = exp(-0.5 * log_tau)
    se_sigmaBYM2 = 0.5 * sigmaBYM2 * se_log_tau
    cat(sprintf("sigmaBYM2: %.4f (SE %.4f)\n", sigmaBYM2, se_sigmaBYM2))
  }

  if(length(idx_phi) > 0) {
    logit_phi = fix[idx_phi[1], "Estimate"]
    se_logit_phi = fix[idx_phi[1], "Std. Error"]
    phi = 1/(1 + exp(-logit_phi))
    se_phi = phi * (1 - phi) * se_logit_phi
    cat(sprintf("phi: %.4f (SE %.4f)\n", phi, se_phi))
  }

  if(length(idx_eps) > 0) {
    log_tauEps = fix[idx_eps[1], "Estimate"]
    se_log_tauEps = fix[idx_eps[1], "Std. Error"]
    sigmaEps = exp(-0.5 * log_tauEps)
    se_sigmaEps = 0.5 * sigmaEps * se_log_tauEps
    cat(sprintf("sigmaEps: %.4f (SE %.4f)\n", sigmaEps, se_sigmaEps))
  }
}

save(out, elapsed, conv, nll,
     file="savedOutput/fitMDM_actual_GH.RData")
cat("Saved to savedOutput/fitMDM_actual_GH.RData\n")
