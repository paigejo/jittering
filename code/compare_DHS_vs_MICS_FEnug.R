#!/usr/bin/env Rscript
# Compare DHS-only vs MICS-only FE+nugget GH models on real data (ed / edMICS)
# All 5 covariates: urban, access, elev, distRiversLakes, normPop
# GH Q=10, DHS K=16/21, MICS KMICS=100

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modFEM.R")
source("code/modFED.R")

Qgh <- 10

# ── Ensure column names ──────────────────────────────────────────
nameTab <- rbind(c("N", "n"), c("N", "ns"), c("Z", "y"), c("Z", "ys"))
for(i in 1:nrow(nameTab)) {
  fromN <- nameTab[i,1]; toN <- nameTab[i,2]
  if(!(toN %in% names(ed)))     ed[[toN]]     <- ed[[fromN]]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] <- edMICS[[fromN]]
}
if(!("Stratum" %in% names(edMICS)))
  edMICS$Stratum <- adm2ToStratumMICS(edMICS$subarea)


# ═════════════════════════════════════════════════════════════════
# (1)  MICS-only FE+nugget GH via fitFEM (all covariates)
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  MICS-only FE+nug GH Q=10, KMICS=100 — all covariates\n")
cat("============================================================\n\n"); flush.console()

fitMICS <- fitFEM(datDHS=ed, datMICS=edMICS, KMICS=100, Qgh=Qgh,
                  fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("MICS convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs (sd)\n",
            fitMICS$opt$convergence, fitMICS$opt$objective,
            fitMICS$totalTime, fitMICS$sdTime))

# Extract MICS estimates
if(!is.null(fitMICS$TMBsd)) {
  peM <- summary(fitMICS$TMBsd, "fixed")
  mics_alpha    <- peM["alpha", "Estimate"]
  mics_alpha_se <- peM["alpha", "Std. Error"]
  betaIdxM      <- which(rownames(peM) == "beta")
  mics_beta     <- peM[betaIdxM, "Estimate"]
  mics_beta_se  <- peM[betaIdxM, "Std. Error"]
  mics_covNames <- fitMICS$covNames
  mics_sigEps   <- exp(-0.5 * peM["log_tauEps", "Estimate"])
  mics_NLL      <- fitMICS$opt$objective
  cat("\nMICS fixed parameter estimates:\n"); print(peM)
} else {
  stop("MICS sdreport failed")
}


# ═════════════════════════════════════════════════════════════════
# (2)  DHS-only FE+nugget GH via fitFED (all covariates)
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  DHS-only FE+nug GH Q=10, K=16/21 — all covariates\n")
cat("============================================================\n\n"); flush.console()

fitDHS <- fitFED(datDHS=ed, KDHSurb=16, KDHSrur=21, Qgh=Qgh,
                 fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("DHS convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs (sd)\n",
            fitDHS$opt$convergence, fitDHS$opt$objective,
            fitDHS$totalTime, fitDHS$sdTime))

# Extract DHS estimates
if(!is.null(fitDHS$TMBsd)) {
  peD <- summary(fitDHS$TMBsd, "fixed")
  dhs_alpha    <- peD["alpha", "Estimate"]
  dhs_alpha_se <- peD["alpha", "Std. Error"]
  betaIdxD     <- which(rownames(peD) == "beta")
  dhs_beta     <- peD[betaIdxD, "Estimate"]
  dhs_beta_se  <- peD[betaIdxD, "Std. Error"]
  dhs_covNames <- fitDHS$covNames
  dhs_sigEps   <- exp(-0.5 * peD["log_tauEps", "Estimate"])
  dhs_NLL      <- fitDHS$opt$objective
  cat("\nDHS fixed parameter estimates:\n"); print(peD)
} else {
  stop("DHS sdreport failed")
}


# ═════════════════════════════════════════════════════════════════
# (3)  Side-by-side comparison
# ═════════════════════════════════════════════════════════════════
cat("\n\n============================================================\n")
cat("  COMPARISON: MICS vs DHS  —  FE+nug GH Q=10, all covariates\n")
cat("============================================================\n\n")

canonNames <- c("urban", "access", "elev", "distRiversLakes", "normPop")

cat(sprintf("%-25s %16s %16s\n", "Parameter", "MICS (K=100)", "DHS (K=16/21)"))
cat(paste0(rep("-", 57), collapse=""), "\n")

# alpha
cat(sprintf("%-25s %8.4f (%5.4f) %8.4f (%5.4f)\n", "alpha",
            mics_alpha, mics_alpha_se, dhs_alpha, dhs_alpha_se))

# betas in canonical order
for(cn in canonNames) {
  mI <- match(cn, mics_covNames)
  dI <- match(cn, dhs_covNames)
  mTxt <- if(!is.na(mI)) sprintf("%8.4f (%5.4f)", mics_beta[mI], mics_beta_se[mI]) else sprintf("%16s", "—")
  dTxt <- if(!is.na(dI)) sprintf("%8.4f (%5.4f)", dhs_beta[dI], dhs_beta_se[dI])  else sprintf("%16s", "—")
  cat(sprintf("%-25s %s %s\n", paste0("beta[", cn, "]"), mTxt, dTxt))
}

# sigmaEps, NLL
cat(sprintf("%-25s %16.4f %16.4f\n", "sigmaEps", mics_sigEps, dhs_sigEps))
cat(sprintf("%-25s %16.4f %16.4f\n", "NLL", mics_NLL, dhs_NLL))
cat(sprintf("%-25s %16.1f %16.1f\n", "Time (s)",
            fitMICS$totalTime + fitMICS$sdTime, fitDHS$totalTime + fitDHS$sdTime))

cat("\nDone.\n")
