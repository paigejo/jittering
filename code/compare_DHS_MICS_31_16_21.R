#!/usr/bin/env Rscript
# Compare MICS (KMICS=100) vs DHS K=16/21 and K=31/41 for FE+nug GH models
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modFEM.R")
source("code/modFED.R")

Qgh <- 10

# ensure names
nameTab <- rbind(c("N", "n"), c("N", "ns"), c("Z", "y"), c("Z", "ys"))
for(i in 1:nrow(nameTab)) {
  fromN <- nameTab[i,1]; toN <- nameTab[i,2]
  if(!(toN %in% names(ed)))     ed[[toN]]     <- ed[[fromN]]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] <- edMICS[[fromN]]
}
if(!("Stratum" %in% names(edMICS))) edMICS$Stratum <- adm2ToStratumMICS(edMICS$subarea)

# 1) MICS-only fit
cat("\n=== MICS-only FE+nug GH Q=10, KMICS=100 ===\n")
fitMICS <- fitFEM(datDHS=ed, datMICS=edMICS, KMICS=100, Qgh=Qgh,
                  beta_pri=c(0, sqrt(1000)),
                  fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat(sprintf("MICS convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs (sd)\n",
            fitMICS$opt$convergence, fitMICS$opt$objective, fitMICS$totalTime, fitMICS$sdTime))
if(is.null(fitMICS$TMBsd)) stop("MICS sdreport failed")

# 2) DHS K=16/21
cat("\n=== DHS-only FE+nug GH Q=10, K=16/21 ===\n")
fitDHS_16_21 <- fitFED(datDHS=ed, KDHSurb=16, KDHSrur=21, Qgh=Qgh,
                       beta_pri=c(0, sqrt(1000)),
                       fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat(sprintf("DHS K=16/21 convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs (sd)\n",
            fitDHS_16_21$opt$convergence, fitDHS_16_21$opt$objective, fitDHS_16_21$totalTime, fitDHS_16_21$sdTime))

# 3) DHS K=31/41
cat("\n=== DHS-only FE+nug GH Q=10, K=31/41 ===\n")
fitDHS_31_41 <- fitFED(datDHS=ed, KDHSurb=31, KDHSrur=41, Qgh=Qgh,
                       beta_pri=c(0, sqrt(1000)),
                       fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat(sprintf("DHS K=31/41 convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs (sd)\n",
            fitDHS_31_41$opt$convergence, fitDHS_31_41$opt$objective, fitDHS_31_41$totalTime, fitDHS_31_41$sdTime))

# Extract fixed effect estimates from sdreport
extractFE <- function(fit) {
  pe <- summary(fit$TMBsd, "fixed")
  list(
    alpha = pe["alpha","Estimate"], alpha_se = pe["alpha","Std. Error"],
    beta = pe[grep("^beta", rownames(pe)), "Estimate"],
    beta_se = pe[grep("^beta", rownames(pe)), "Std. Error"],
    sigEps = exp(-0.5 * pe["log_tauEps","Estimate"]),
    NLL = fit$opt$objective,
    covNames = fit$covNames
  )
}

resM <- extractFE(fitMICS)
resD16 <- extractFE(fitDHS_16_21)
resD31 <- extractFE(fitDHS_31_41)

# print comparison table
cat("\n\n=== Comparison table ===\n")
canonNames <- c("urban", "access", "elev", "distRiversLakes", "normPop")
cat(sprintf("%-25s %20s %20s %20s\n", "Parameter", "MICS (K=100)", "DHS (K=16/21)", "DHS (K=31/41)"))
cat(paste0(rep("-", 90), collapse=""), "\n")
cat(sprintf("%-25s %8.4f (%6.4f) %8.4f (%6.4f) %8.4f (%6.4f)\n", "alpha",
            resM$alpha, resM$alpha_se, resD16$alpha, resD16$alpha_se, resD31$alpha, resD31$alpha_se))

for(cn in canonNames) {
  mI <- match(cn, resM$covNames)
  d16I <- match(cn, resD16$covNames)
  d31I <- match(cn, resD31$covNames)
  mTxt <- if(!is.na(mI)) sprintf("%8.4f (%6.4f)", resM$beta[mI], resM$beta_se[mI]) else sprintf("%16s", "—")
  dTxt16 <- if(!is.na(d16I)) sprintf("%8.4f (%6.4f)", resD16$beta[d16I], resD16$beta_se[d16I]) else sprintf("%16s", "—")
  dTxt31 <- if(!is.na(d31I)) sprintf("%8.4f (%6.4f)", resD31$beta[d31I], resD31$beta_se[d31I]) else sprintf("%16s", "—")
  cat(sprintf("%-25s %s %s %s\n", paste0("beta[", cn, "]"), mTxt, dTxt16, dTxt31))
}
cat(sprintf("%-25s %20.4f %20.4f %20.4f\n", "sigmaEps", resM$sigEps, resD16$sigEps, resD31$sigEps))
cat(sprintf("%-25s %20.4f %20.4f %20.4f\n", "NLL", resM$NLL, resD16$NLL, resD31$NLL))
cat("\nDone.\n")
