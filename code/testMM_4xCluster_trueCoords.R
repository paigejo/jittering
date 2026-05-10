#!/usr/bin/env Rscript
# Fit IID+nugget, FE+nugget, and BYM2+nugget models on the 4x cluster sim 1
# using TRUE spatial coordinates (no integration point averaging over strata).
# K=1 integration point per cluster at exact location, weight=1.

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")
source("code/modBYM2.R")
source("code/modM_MSepMarg.R")

saveFileIID  = "savedOutput/testMM_4xCluster_trueCoords_IID_sim1.RData"
saveFileFE   = "savedOutput/testMM_4xCluster_trueCoords_FE_sim1.RData"
saveFileBYM2 = "savedOutput/testMM_4xCluster_trueCoords_BYM2_sim1.RData"

cat("===============================================\n")
cat("4x CLUSTER, TRUE COORDINATES: FIT 3 MODELS\n")
cat("===============================================\n\n")

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

cat(sprintf("True: alpha=%.2f, urban=%.2f, normPop=%.2f, sigmaEps=%.4f\n\n",
            trueAlpha, trueUrban, trueNormPop, trueSigmaEps))

# Load the 4x cluster simulation
load("savedOutput/simStudy1/simPopsSurveys_4xCluster.RData")
edMICS1 = surveysMICS[[1]]

cat(sprintf("Sim 1: %d MICS clusters, %d total obs (avg %.1f per cluster)\n",
            nrow(edMICS1), sum(edMICS1$N), mean(edMICS1$N)))

# ============================================================
# Build inputsMDM with true coordinates (K=1)
# ============================================================
cat("\nBuilding inputs with TRUE coordinates (K=1, no integration point averaging)...\n")

# Add Stratum
edMICS1$Stratum = adm2ToStratumMICS(edMICS1$subarea)

# Ensure expected column names
edMICS1$ns = edMICS1$N
edMICS1$n  = edMICS1$N
edMICS1$ys = edMICS1$Z
edMICS1$y  = edMICS1$Z

# Extract covariates at true cluster locations
cat("  Extracting covariates at true cluster (lon, lat)...\n")
clusterCoords = cbind(edMICS1$lon, edMICS1$lat)
covMat = getDesignMatPopNorm(clusterCoords, useThreshPopMat=TRUE)
# covMat cols: int, pop, urb, access, elev, distRiversLakes, urbanicity

# Add normPop (same as how integration points compute it)
load("savedOutput/global/popMeanSDCal.RData")
load("savedOutput/global/covariatesNorm.RData")
popVals = terra::extract(pop, clusterCoords, method="bilinear")
normPopVals = (log1p(popVals) - popMeanCal) / popSDCal
normPopVals[is.na(normPopVals)] = min(normPopVals, na.rm=TRUE)

# Build covariate matrix matching the format: urban, access, elev, distRiversLakes, normPop
clusterCovs = cbind(
  urban = as.numeric(covMat[,"urb"]),
  access = as.numeric(covMat[,"access"]),
  elev = as.numeric(covMat[,"elev"]),
  distRiversLakes = as.numeric(covMat[,"distRiversLakes"]),
  normPop = as.numeric(normPopVals)
)

# Sort data by Stratum (same as makeInputsMDM does)
stratOrder = admFinal$NAME_FINAL
edMICS1 = sortByCol(edMICS1, "Stratum", stratOrder)

# Re-sort covariates to match
sortIdx = order(match(edMICS1$Stratum, stratOrder))
# Actually sortByCol already sorted edMICS1, but we need to track covariate rows.
# Re-extract from sorted data:
clusterCoords_sorted = cbind(edMICS1$lon, edMICS1$lat)
covMat_sorted = getDesignMatPopNorm(clusterCoords_sorted, useThreshPopMat=TRUE)
popVals_sorted = terra::extract(pop, clusterCoords_sorted, method="bilinear")
normPop_sorted = (log1p(popVals_sorted) - popMeanCal) / popSDCal
normPop_sorted[is.na(normPop_sorted)] = min(normPop_sorted, na.rm=TRUE)

clusterCovs_sorted = cbind(
  urban = as.numeric(covMat_sorted[,"urb"]),
  access = as.numeric(covMat_sorted[,"access"]),
  elev = as.numeric(covMat_sorted[,"elev"]),
  distRiversLakes = as.numeric(covMat_sorted[,"distRiversLakes"]),
  normPop = as.numeric(normPop_sorted)
)

# Split by urban/rural
isUrb = edMICS1$urban
isRur = !edMICS1$urban

XUrb = as.matrix(clusterCovs_sorted[isUrb, , drop=FALSE])  # nUrb x 5
XRur = as.matrix(clusterCovs_sorted[isRur, , drop=FALSE])  # nRur x 5

nUrb = sum(isUrb)
nRur = sum(isRur)

cat(sprintf("  %d urban clusters, %d rural clusters\n", nUrb, nRur))

# Weights: K=1, weight=1 for each cluster (no averaging)
wUrban = matrix(1, nrow=nUrb, ncol=1)
wRural = matrix(1, nrow=nRur, ncol=1)

# Area indices (0-indexed for TMB)
strataMICS = stratOrder
urbStrata = edMICS1$Stratum[isUrb]
rurStrata = edMICS1$Stratum[isRur]
areaidxlocUrbanMICS = as.integer(match(urbStrata, strataMICS) - 1)
areaidxlocRuralMICS = as.integer(match(rurStrata, strataMICS) - 1)

# Observation counts
ysUrbMICS = edMICS1$ys[isUrb]
nsUrbMICS = edMICS1$ns[isUrb]
ysRurMICS = edMICS1$ys[isRur]
nsRurMICS = edMICS1$ns[isRur]

# Build inputsMDM (only MICS parts needed)
# DHS parts are needed by makeInputsMDM but we'll provide dummy values
# since we're running M_M (MICS-only) models
inputsMDM_trueCoords = list(
  intPtsMICS = list(
    XUrb = XUrb,
    XRur = XRur,
    wUrban = wUrban,
    wRural = wRural,
    strataMICS = strataMICS
  ),
  intPtsDHS = list(
    covsUrb = matrix(0, nrow=1, ncol=5),
    covsRur = matrix(0, nrow=1, ncol=5)
  ),
  areaidxlocUrbanMICS = areaidxlocUrbanMICS,
  areaidxlocRuralMICS = areaidxlocRuralMICS,
  areaidxlocUrbanDHS = integer(0),
  areaidxlocRuralDHS = integer(0),
  AUrbMICS = NULL, ARurMICS = NULL,
  AUrbDHS = NULL, ARurDHS = NULL,
  ysUrbMICS = ysUrbMICS,
  nsUrbMICS = nsUrbMICS,
  ysRurMICS = ysRurMICS,
  nsRurMICS = nsRurMICS,
  ysUrbDHS = numeric(0),
  ysRurDHS = numeric(0),
  nsUrbDHS = numeric(0),
  nsRurDHS = numeric(0)
)

cat("  inputsMDM built with K=1 true-coordinate integration points.\n")

# ============================================================
# Fit 3 models
# ============================================================

# --- MODEL 1: IID + nugget ---
result_iid = NULL
if(file.exists(saveFileIID)) {
  cat("\nLoading saved IID+nugget results...\n")
  load(saveFileIID)
} else {
  cat("\n--- Fitting IID + nugget model (true coords) ---\n")
  tryCatch({
    result_iid = fitMM_IIDonly(datMICS=edMICS1,
                                inputsMDM=inputsMDM_trueCoords,
                                KMICS=1,
                                covariates=c("urban", "normPop"),
                                fixedEffectsOnly=FALSE,
                                getSDs=TRUE, doMCMC=FALSE)
    save(result_iid, file=saveFileIID)
    cat("IID+nugget completed and saved.\n")
  }, error = function(e) {
    cat("IID+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# --- MODEL 2: FE + nugget ---
result_fe = NULL
if(file.exists(saveFileFE)) {
  cat("\nLoading saved FE+nugget results...\n")
  load(saveFileFE)
} else {
  cat("\n--- Fitting FE + nugget model (true coords) ---\n")
  tryCatch({
    result_fe = fitMM_IIDonly(datMICS=edMICS1,
                              inputsMDM=inputsMDM_trueCoords,
                              KMICS=1,
                              covariates=c("urban", "normPop"),
                              fixedEffectsOnly=TRUE,
                              getSDs=TRUE, doMCMC=FALSE)
    save(result_fe, file=saveFileFE)
    cat("FE+nugget completed and saved.\n")
  }, error = function(e) {
    cat("FE+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# --- MODEL 3: BYM2 + nugget ---
result_bym2 = NULL
if(file.exists(saveFileBYM2)) {
  cat("\nLoading saved BYM2+nugget results...\n")
  load(saveFileBYM2)
} else {
  cat("\n--- Fitting BYM2 + nugget model (true coords) ---\n")
  tryCatch({
    result_bym2 = fitMMMarg(datMICS=edMICS1,
                             inputsMDM=inputsMDM_trueCoords,
                             KMICS=1,
                             covariates=c("urban", "normPop"),
                             getSDs=TRUE, doMCMC=FALSE)
    save(result_bym2, file=saveFileBYM2)
    cat("BYM2+nugget completed and saved.\n")
  }, error = function(e) {
    cat("BYM2+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# Extract estimates (same functions as 4xCluster)
# ============================================================
extractEst_IID = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed; pfn = names(pf); se_f = sqrt(diag(SD0$cov.fixed))
  sr = summary(SD0, select="random"); rn = names(SD0$par.random)
  alpha_est = sr[rn == "alpha", 1]; alpha_se = sr[rn == "alpha", 2]
  beta_est = sr[rn == "beta", 1]; beta_se = sr[rn == "beta", 2]
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  list(alpha=alpha_est, alpha_se=alpha_se, beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2], sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

extractEst_FE = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed; pfn = names(pf); se_f = sqrt(diag(SD0$cov.fixed))
  alpha_est = as.numeric(pf[pfn == "alpha"]); alpha_se = as.numeric(se_f[pfn == "alpha"])
  beta_est = as.numeric(pf[pfn == "beta"]); beta_se = as.numeric(se_f[pfn == "beta"])
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  list(alpha=alpha_est, alpha_se=alpha_se, beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2], sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

extractEst_BYM2 = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed; pfn = names(pf); se_f = sqrt(diag(SD0$cov.fixed))
  sr = summary(SD0, select="random"); rn = names(SD0$par.random)
  alpha_est = sr[rn == "alpha", 1]; alpha_se = sr[rn == "alpha", 2]
  beta_est = sr[rn == "beta", 1]; beta_se = sr[rn == "beta", 2]
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  list(alpha=alpha_est, alpha_se=alpha_se, beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2], sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

# ============================================================
# Print results
# ============================================================
cat("\n===============================================\n")
cat("RESULTS: 4x CLUSTER, TRUE COORDINATES (sim 1)\n")
cat("===============================================\n\n")

fmtVal = function(est, se) {
  if(is.na(est)) return(sprintf("%20s", "NA"))
  sprintf("%7.4f (SE %5.4f)", est, se)
}

e_iid = if(!is.null(result_iid)) extractEst_IID(result_iid) else NULL
e_fe  = if(!is.null(result_fe))  extractEst_FE(result_fe) else NULL
e_bym2 = if(!is.null(result_bym2)) extractEst_BYM2(result_bym2) else NULL

getV = function(e, field, sefield) {
  if(is.null(e)) return(c(NA, NA))
  c(e[[field]], e[[sefield]])
}

cat(sprintf("%-15s %7s  %-22s %-22s %-22s\n", "Parameter", "Truth", "IID+nugget", "FE+nugget", "BYM2+nugget"))
cat(sprintf("%-15s %7s  %-22s %-22s %-22s\n", "---------", "-----", "----------", "---------", "-----------"))

printRow = function(name, truth, e_iid, e_fe, e_bym2, field, sefield) {
  v1 = getV(e_iid, field, sefield)
  v2 = getV(e_fe, field, sefield)
  v3 = getV(e_bym2, field, sefield)
  cat(sprintf("%-15s %7.4f  %s %s %s\n", name, truth,
              fmtVal(v1[1], v1[2]), fmtVal(v2[1], v2[2]), fmtVal(v3[1], v3[2])))
}

printRow("alpha", trueAlpha, e_iid, e_fe, e_bym2, "alpha", "alpha_se")
printRow("beta[urban]", trueUrban, e_iid, e_fe, e_bym2, "beta1", "beta1_se")
printRow("beta[normPop]", trueNormPop, e_iid, e_fe, e_bym2, "beta2", "beta2_se")
printRow("sigmaEps", trueSigmaEps, e_iid, e_fe, e_bym2, "sigmaEps", "sigmaEps_se")

cat("\nBias summary:\n")
cat(sprintf("%-15s  %-10s %-10s %-10s\n", "Parameter", "IID+nug", "FE+nug", "BYM2+nug"))
cat(sprintf("%-15s  %-10s %-10s %-10s\n", "---------", "-------", "------", "--------"))

printBias = function(name, truth, e_iid, e_fe, e_bym2, field) {
  b1 = if(!is.null(e_iid)) sprintf("%+8.4f", e_iid[[field]] - truth) else sprintf("%8s", "NA")
  b2 = if(!is.null(e_fe))  sprintf("%+8.4f", e_fe[[field]] - truth)  else sprintf("%8s", "NA")
  b3 = if(!is.null(e_bym2)) sprintf("%+8.4f", e_bym2[[field]] - truth) else sprintf("%8s", "NA")
  cat(sprintf("%-15s  %-10s %-10s %-10s\n", name, b1, b2, b3))
}

printBias("alpha", trueAlpha, e_iid, e_fe, e_bym2, "alpha")
printBias("beta[urban]", trueUrban, e_iid, e_fe, e_bym2, "beta1")
printBias("beta[normPop]", trueNormPop, e_iid, e_fe, e_bym2, "beta2")
printBias("sigmaEps", trueSigmaEps, e_iid, e_fe, e_bym2, "sigmaEps")

cat("\nTiming:\n")
if(!is.null(e_iid))  cat(sprintf("  IID+nugget:  %.1f min\n", e_iid$time/60))
if(!is.null(e_fe))   cat(sprintf("  FE+nugget:   %.1f min\n", e_fe$time/60))
if(!is.null(e_bym2)) cat(sprintf("  BYM2+nugget: %.1f min\n", e_bym2$time/60))

cat("\n===============================================\n")
cat("DONE\n")
cat("===============================================\n")
