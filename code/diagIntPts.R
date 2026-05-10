#!/usr/bin/env Rscript
# Diagnostic: compare integration-point covariates to true-location covariates
# This checks whether the integration point pipeline produces reasonable
# covariate values for each cluster.

source("code/setup.R")
source("code/makeInputsTMB.R")

load("savedOutput/simStudy1/simPopsSurveys.RData")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

ed = surveysDHS[[1]]
edMICS = surveysMICS[[1]]

# Add Stratum column if missing (simulated data may not have it)
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}

# Load precomputed MICS integration points
load("savedOutput/global/intPtsMICS_100.RData")

# Get the integration point inputs
inputsMDM = makeInputsMDM(ed, edMICS, intPtsDHS=intPtsDHS, intPtsMICS=intPtsMICS)

# ---- Integration point covariates ----
XUrb = inputsMDM$intPtsMICS$XUrb  # [K*nObsUrb x nCov] matrix
XRur = inputsMDM$intPtsMICS$XRur
wUrb = inputsMDM$intPtsMICS$wUrban  # [nObsUrb x K]
wRur = inputsMDM$intPtsMICS$wRural

cat("XUrb dim:", dim(XUrb), "\n")
cat("XRur dim:", dim(XRur), "\n")
cat("wUrb dim:", dim(wUrb), "\n")
cat("wRur dim:", dim(wRur), "\n")
cat("XUrb colnames:", colnames(XUrb), "\n")

nUrb = nrow(wUrb)
nRur = nrow(wRur)
K = ncol(wUrb)
cat("\nnUrb:", nUrb, "nRur:", nRur, "K:", K, "\n")
cat("Expected XUrb rows:", nUrb * K, "actual:", nrow(XUrb), "\n")
cat("Expected XRur rows:", nRur * K, "actual:", nrow(XRur), "\n")

# Compute weighted-average covariates at integration points
covAvgUrb = matrix(0, nUrb, ncol(XUrb))
for(j in 1:ncol(XUrb)) {
  xmat = matrix(XUrb[,j], nrow=nUrb, ncol=K)
  covAvgUrb[,j] = rowSums(xmat * wUrb)
}
colnames(covAvgUrb) = colnames(XUrb)

covAvgRur = matrix(0, nRur, ncol(XRur))
for(j in 1:ncol(XRur)) {
  xmat = matrix(XRur[,j], nrow=nRur, ncol=K)
  covAvgRur[,j] = rowSums(xmat * wRur)
}
colnames(covAvgRur) = colnames(XRur)

cat("\n=== Integration point weighted-avg covariates (first 5 urban obs) ===\n")
print(round(covAvgUrb[1:5,], 4))

cat("\n=== Integration point weighted-avg covariates (first 5 rural obs) ===\n")
print(round(covAvgRur[1:5,], 4))

# ---- True-location covariates ----
# Sort MICS data same way as makeInputsMDM does
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}
# Re-sort to match the integration point ordering
load("savedOutput/global/intPtsMICS_100.RData")
intPtsMICS = straightenMICS(intPtsMICS)
edMICS = sortByCol(edMICS, "Stratum", intPtsMICS$strataMICS)

load("savedOutput/global/popMatNGAThresh.RData")
load("savedOutput/global/popMeanSDCal.RData")

popEN = cbind(popMatNGAThresh$east, popMatNGAThresh$north)

matchPixel = function(east, north) {
  closeE = (east > popEN[,1] - 2.5) & (east <= popEN[,1] + 2.5)
  closeN = (north > popEN[,2] - 2.5) & (north <= popEN[,2] + 2.5)
  w = which(closeE & closeN)
  if(length(w) >= 1) return(w[1])
  dists = (east - popEN[,1])^2 + (north - popEN[,2])^2
  return(which.min(dists))
}

edUrb = edMICS[edMICS$urban,]
edRur = edMICS[!edMICS$urban,]

cat("\nnUrb from data:", nrow(edUrb), "nRur from data:", nrow(edRur), "\n")

# Get pixel indices for true locations
pixIsUrb = sapply(1:nrow(edUrb), function(i) matchPixel(edUrb$east[i], edUrb$north[i]))
pixIsRur = sapply(1:nrow(edRur), function(i) matchPixel(edRur$east[i], edRur$north[i]))

# Extract true covariates
trueUrbanUrb = as.numeric(popMatNGAThresh$urban[pixIsUrb])
truePopUrb = popMatNGAThresh$pop[pixIsUrb]
trueNormPopUrb = (log1p(truePopUrb) - popMeanCalThresh) / popSDCalThresh

trueUrbanRur = as.numeric(popMatNGAThresh$urban[pixIsRur])
truePopRur = popMatNGAThresh$pop[pixIsRur]
trueNormPopRur = (log1p(truePopRur) - popMeanCalThresh) / popSDCalThresh

# ---- Compare ----
cat("\n=== URBAN clusters: true vs int-pt weighted avg ===\n")
cat("urban covariate:\n")
cat("  True mean:", round(mean(trueUrbanUrb), 4), "\n")
cat("  IntPt mean:", round(mean(covAvgUrb[,"urban"]), 4), "\n")
cat("normPop covariate:\n")
cat("  True mean:", round(mean(trueNormPopUrb), 4), "\n")
cat("  IntPt mean:", round(mean(covAvgUrb[,"normPop"]), 4), "\n")

cat("\n=== RURAL clusters: true vs int-pt weighted avg ===\n")
cat("urban covariate:\n")
cat("  True mean:", round(mean(trueUrbanRur), 4), "\n")
cat("  IntPt mean:", round(mean(covAvgRur[,"urban"]), 4), "\n")
cat("normPop covariate:\n")
cat("  True mean:", round(mean(trueNormPopRur), 4), "\n")
cat("  IntPt mean:", round(mean(covAvgRur[,"normPop"]), 4), "\n")

# Per-cluster comparison (first 10)
cat("\n=== First 10 urban clusters: true normPop vs intPt normPop ===\n")
comp = data.frame(
  true_urban = trueUrbanUrb[1:10],
  intPt_urban = round(covAvgUrb[1:10, "urban"], 4),
  true_normPop = round(trueNormPopUrb[1:10], 4),
  intPt_normPop = round(covAvgUrb[1:10, "normPop"], 4)
)
print(comp)

cat("\n=== First 10 rural clusters: true normPop vs intPt normPop ===\n")
comp = data.frame(
  true_urban = trueUrbanRur[1:10],
  intPt_urban = round(covAvgRur[1:10, "urban"], 4),
  true_normPop = round(trueNormPopRur[1:10], 4),
  intPt_normPop = round(covAvgRur[1:10, "normPop"], 4)
)
print(comp)

# Check the layout of XUrb: is it [obs1_k1, obs2_k1, ..., obs1_k2, obs2_k2, ...]
# or [obs1_k1, obs1_k2, ..., obs2_k1, obs2_k2, ...]?
cat("\n=== XUrb layout check ===\n")
cat("First integration point values for obs 1:\n")
cat("  XUrb[1,]:", round(XUrb[1,], 4), "\n")
cat("  XUrb[nUrb+1,]:", round(XUrb[nUrb+1,], 4), "\n")
cat("  XUrb[2*nUrb+1,]:", round(XUrb[2*nUrb+1,], 4), "\n")
cat("Values for obs 2:\n")
cat("  XUrb[2,]:", round(XUrb[2,], 4), "\n")
cat("  XUrb[nUrb+2,]:", round(XUrb[nUrb+2,], 4), "\n")
