#!/usr/bin/env Rscript
# Compare FE+nugget: GH quadrature vs Laplace approximation
# Uses fitMFEM (GH) and fitMFEM_laplace (Laplace) on simulated MICS data
# Both use KMICS=100 integration points and same covariates

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")
source("code/modM_FEnug.R")

# ── Load BYM2 simulated data ─────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

# Fix names
nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}
if(!("Stratum" %in% names(edMICS)))
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

# ── Build shared inputs once ─────────────────────────────────────
cat("Building integration point inputs (KMICS=100)...\n")
inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

# ═════════════════════════════════════════════════════════════════
# (1) GH model (nuggets integrated out analytically)
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  FE+nugget with GAUSS-HERMITE (Q=25)\n")
cat("============================================================\n")
fitGH = fitMFEM(datDHS=ed, datMICS=edMICS, inputsMDM=inputsMDM,
                covariates=c("urban", "normPop"),
                fixedEffectsOnly=TRUE, Qgh=25, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d\n", fitGH$opt$convergence))
cat(sprintf("NLL: %.4f\n", fitGH$opt$objective))
cat(sprintf("Time: %.1f s (opt) + %.1f s (sd)\n", fitGH$totalTime, fitGH$sdTime))

if(!is.null(fitGH$TMBsd)) {
  estGH = summary(fitGH$TMBsd, "fixed")
  cat("\nFixed parameter estimates (GH):\n")
  print(estGH)
} else {
  cat("Warning: sdreport failed for GH model\n")
}

# ═════════════════════════════════════════════════════════════════
# (2) Laplace model (nuggets as TMB random effects)
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  FE+nugget with LAPLACE approximation\n")
cat("============================================================\n")
fitLap = fitMFEM_laplace(datDHS=ed, datMICS=edMICS, inputsMDM=inputsMDM,
                         covariates=c("urban", "normPop"),
                         fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d\n", fitLap$opt$convergence))
cat(sprintf("NLL: %.4f\n", fitLap$opt$value))
cat(sprintf("Time: %.1f s (opt) + %.1f s (sd)\n", fitLap$totalTime, fitLap$sdTime))

if(!is.null(fitLap$TMBsd)) {
  estLap = summary(fitLap$TMBsd, "fixed")
  cat("\nFixed parameter estimates (Laplace):\n")
  print(estLap)
} else {
  cat("Warning: sdreport failed for Laplace model\n")
}

# ═════════════════════════════════════════════════════════════════
# (3) Summary comparison
# ═════════════════════════════════════════════════════════════════
cat("\n\n============================================================\n")
cat("  COMPARISON: GH vs Laplace for FE+nugget (MICS-only)\n")
cat("============================================================\n\n")

truth = c(alpha=trueAlpha, beta1=trueUrban, beta2=trueNormPop, sigmaEps=trueSigmaEps)

# Extract GH estimates
if(!is.null(fitGH$TMBsd)) {
  peGH = summary(fitGH$TMBsd, "fixed")
  ghAlpha = peGH["alpha", "Estimate"]
  ghBeta = peGH[grep("^beta$", rownames(peGH)), "Estimate"]
  ghLogTauEps = peGH["log_tauEps", "Estimate"]
  ghSigmaEps = exp(-0.5 * ghLogTauEps)
  ghSE_alpha = peGH["alpha", "Std. Error"]
  ghSE_beta = peGH[grep("^beta$", rownames(peGH)), "Std. Error"]
}

# Extract Laplace estimates
if(!is.null(fitLap$TMBsd)) {
  peLap = summary(fitLap$TMBsd, "fixed")
  lapAlpha = peLap["alpha", "Estimate"]
  lapBeta = peLap[grep("^beta$", rownames(peLap)), "Estimate"]
  lapLogTauEps = peLap["log_tauEps", "Estimate"]
  lapSigmaEps = exp(-0.5 * lapLogTauEps)
  lapSE_alpha = peLap["alpha", "Std. Error"]
  lapSE_beta = peLap[grep("^beta$", rownames(peLap)), "Std. Error"]
}

cat(sprintf("%-12s %10s %10s %10s %10s %10s\n",
            "Parameter", "Truth", "GH", "Laplace", "GH_SE", "Lap_SE"))
cat(sprintf("%-12s %10s %10s %10s %10s %10s\n",
            "---------", "-----", "------", "-------", "-----", "------"))

if(!is.null(fitGH$TMBsd) && !is.null(fitLap$TMBsd)) {
  cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              "alpha", trueAlpha, ghAlpha, lapAlpha, ghSE_alpha, lapSE_alpha))
  cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              "beta[urban]", trueUrban, ghBeta[1], lapBeta[1], ghSE_beta[1], lapSE_beta[1]))
  cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              "beta[pop]", trueNormPop, ghBeta[2], lapBeta[2], ghSE_beta[2], lapSE_beta[2]))
  cat(sprintf("%-12s %10.4f %10.4f %10.4f %10s %10s\n",
              "sigmaEps", trueSigmaEps, ghSigmaEps, lapSigmaEps, "—", "—"))

  cat(sprintf("\n%-12s %10s %10.4f %10.4f\n", "NLL", "—",
              fitGH$opt$objective, fitLap$opt$value))
  cat(sprintf("%-12s %10s %10.1f %10.1f\n", "Time (s)", "—",
              fitGH$totalTime, fitLap$totalTime))
}

cat("\nDone.\n")
