#!/usr/bin/env Rscript
# Test FE+nugget MICS-only GH model on ed/edMICS with KMICS=100, Q=10
# Uses fitMFEM from modM_FEnug.R

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_FEnug.R")

# ── True values (from simStudy1) ─────────────────────────────────
trueAlpha   = -1.25
trueUrban   = 1.00
trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

# ── Load BYM2 simulated data ─────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS_sim = surveysMICS[[1]]

# ── Run on simulated edMICS (KMICS=100, Q=10) ───────────────────
cat("\n============================================================\n")
cat("  FE+nugget GH on simulated edMICS (KMICS=100, Q=10)\n")
cat("============================================================\n\n")

fit1 = fitMFEM(datDHS=ed, datMICS=edMICS_sim,
               KMICS=100, Qgh=10,
               covariates=c("urban", "normPop"),
               fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d\n", fit1$opt$convergence))
cat(sprintf("NLL: %.4f\n", fit1$opt$objective))
cat(sprintf("Time: %.1f s (opt) + %.1f s (sd)\n", fit1$totalTime, fit1$sdTime))

if(!is.null(fit1$TMBsd)) {
  est1 = summary(fit1$TMBsd, "fixed")
  cat("\nFixed parameter estimates:\n")
  print(est1)

  pe = est1[,"Estimate"]; se = est1[,"Std. Error"]; nms = rownames(est1)
  alphaEst  = pe[nms == "alpha"]
  betaEst   = pe[nms == "beta"]
  alphaSE   = se[nms == "alpha"]
  betaSE    = se[nms == "beta"]
  logTauEps = pe[nms == "log_tauEps"]
  sigmaEps  = exp(-0.5 * logTauEps)

  cat(sprintf("\n--- Simulated edMICS Results (KMICS=100, Q=10) ---\n"))
  cat(sprintf("  alpha:      %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
              alphaEst, alphaSE, trueAlpha, abs(alphaEst - trueAlpha) / alphaSE))
  cat(sprintf("  beta[urb]:  %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
              betaEst[1], betaSE[1], trueUrban, abs(betaEst[1] - trueUrban) / betaSE[1]))
  cat(sprintf("  beta[pop]:  %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
              betaEst[2], betaSE[2], trueNormPop, abs(betaEst[2] - trueNormPop) / betaSE[2]))
  cat(sprintf("  sigmaEps:   %7.4f              truth: %7.4f\n", sigmaEps, trueSigmaEps))
} else {
  cat("Warning: sdreport failed for fit1\n")
}

# ── Run on real ed/edMICS (KMICS=100, Q=10) ─────────────────────
cat("\n\n============================================================\n")
cat("  FE+nugget GH on real edMICS (KMICS=100, Q=10)\n")
cat("============================================================\n\n")

fit2 = fitMFEM(datDHS=ed, datMICS=edMICS,
               KMICS=100, Qgh=10,
               covariates=c("urban", "normPop"),
               fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d\n", fit2$opt$convergence))
cat(sprintf("NLL: %.4f\n", fit2$opt$objective))
cat(sprintf("Time: %.1f s (opt) + %.1f s (sd)\n", fit2$totalTime, fit2$sdTime))

if(!is.null(fit2$TMBsd)) {
  est2 = summary(fit2$TMBsd, "fixed")
  cat("\nFixed parameter estimates:\n")
  print(est2)

  pe2 = est2[,"Estimate"]; se2 = est2[,"Std. Error"]; nms2 = rownames(est2)
  alphaEst2  = pe2[nms2 == "alpha"]
  betaEst2   = pe2[nms2 == "beta"]
  alphaSE2   = se2[nms2 == "alpha"]
  betaSE2    = se2[nms2 == "beta"]
  logTauEps2 = pe2[nms2 == "log_tauEps"]
  sigmaEps2  = exp(-0.5 * logTauEps2)

  cat(sprintf("\n--- Real edMICS Results (KMICS=100, Q=10) ---\n"))
  cat(sprintf("  alpha:      %7.4f (SE %.4f)\n", alphaEst2, alphaSE2))
  cat(sprintf("  beta[urb]:  %7.4f (SE %.4f)\n", betaEst2[1], betaSE2[1]))
  cat(sprintf("  beta[pop]:  %7.4f (SE %.4f)\n", betaEst2[2], betaSE2[2]))
  cat(sprintf("  sigmaEps:   %7.4f\n", sigmaEps2))
} else {
  cat("Warning: sdreport failed for fit2\n")
}

cat("\nDone.\n")
