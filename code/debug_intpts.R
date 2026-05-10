# Debug: check what covariates the integration-point model actually sees
setwd("c:/Users/jpaige/git/jittering")
suppressMessages(suppressWarnings(source("code/setup.R")))
source("code/simData.R")

library(TMB)

# Load the BYM2 sim data
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

# Add required columns
edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
edMICS$ns = edMICS$N
edMICS$n  = edMICS$N
edMICS$ys = edMICS$Z
edMICS$y  = edMICS$Z

# Load integration points (same as what testMM_BYM2sim.R uses)
# The BYM2 sim models use default KMICS=100
intPtFile = "savedOutput/global/intPtsMICS_100.RData"
load(intPtFile)

# Replicate what makeInputsMDM does
intPtsMICS_orig = intPtsMICS  # save original
intPtsMICS = straightenMICS(intPtsMICS)
datMICS = sortByCol(edMICS, "Stratum", intPtsMICS$strataMICS)

KMICS = 100
allNumPerStrat = aggregate(datMICS$Stratum, by=list(strat=datMICS$Stratum, urb=datMICS$urban), FUN=length, drop=FALSE)
allNumPerStrat = straightenNumPerStrat(allNumPerStrat, intPtsMICS$strataMICS)
numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
numPerStratRur[is.na(numPerStratRur)] = 0

cat("numPerStratUrb:", head(numPerStratUrb), "... total:", sum(numPerStratUrb), "\n")
cat("numPerStratRur:", head(numPerStratRur), "... total:", sum(numPerStratRur), "\n")

# Build XUrb exactly as makeInputsMDM does
XUrb = intPtsMICS$XUrb
nAreas = nrow(XUrb)/KMICS
cat("nAreas:", nAreas, "\n")

areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrb[x])}))
allAreaIs = rep(areaI, KMICS)
nUrb = length(allAreaIs)/KMICS
allIntIs = rep(1:KMICS, each=nUrb)
transformIUrb = allAreaIs + (allIntIs-1)*nAreas
XUrb = XUrb[transformIUrb,]
XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]

cat("\nXUrb after expansion:\n")
cat("  dim:", dim(XUrb), "\n")
cat("  nUrb:", nUrb, "KMICS:", KMICS, "product:", nUrb*KMICS, "\n")
cat("  cols:", names(XUrb), "\n")

# Now subset to just urban + normPop (what the model does with covariates=c("urban","normPop"))
keepIdx = which(names(XUrb) %in% c("urban", "normPop"))
cat("  keepIdx:", keepIdx, "names:", names(XUrb)[keepIdx], "\n")

XUrb_sub = as.matrix(XUrb[, keepIdx])
cat("  XUrb_sub dim:", dim(XUrb_sub), "\n")
cat("  XUrb_sub urban values: first 20:", head(XUrb_sub[,"urban"], 20), "\n")
cat("  XUrb_sub urban unique:", unique(XUrb_sub[,"urban"]), "\n")
cat("  XUrb_sub normPop range:", range(XUrb_sub[,"normPop"]), "\n")

# Key check: does the model see the CORRECT normPop for each cluster?
# In the C++ template, for urban obs i, integration point j:
#   thisIndex = num_iUrbanMICS * intI + obsI
#   eta = alpha + fe_iUrbanMICS(thisIndex) + u_spatial(area[obsI]) + nugget[obsI]
# So fe_iUrbanMICS = X_beta * beta, where X_beta is [K*nUrb x nCov]
# The indexing is: for intPt j (0-indexed), cluster i (0-indexed):
#   thisIndex = nUrb * j + i
# This means the matrix is laid out as:
#   rows 0..nUrb-1: intPt 0 for all clusters
#   rows nUrb..2*nUrb-1: intPt 1 for all clusters
#   etc.

# Let's check if the normPop for cluster i is actually varying across integration points
# (it must, because different int pts are at different locations within the stratum)
cat("\n=== CHECK: normPop variation for specific clusters ===\n")

# Get normPop for first urban cluster across all K int pts
obsI = 1
normPop_obs1 = numeric(KMICS)
for(intI in 0:(KMICS-1)) {
  thisIndex = nUrb * intI + obsI
  normPop_obs1[intI+1] = XUrb_sub[thisIndex, "normPop"]
}
cat("Cluster 1 (urban) normPop across K=100 int pts:\n")
cat("  mean:", round(mean(normPop_obs1),4), "sd:", round(sd(normPop_obs1),4), "\n")
cat("  range:", round(range(normPop_obs1),4), "\n")

# Now check: what normPop does actual cluster 1 have?
urbIdx = which(datMICS$urban)
cat("  Actual cluster 1 stratum:", datMICS$Stratum[urbIdx[2]], "\n")

# Get the actual normPop at cluster 1's true location
load("savedOutput/global/covariatesNorm.RData")
load("savedOutput/global/popMeanSDCal.RData")
cl1Coords = cbind(datMICS$lon[urbIdx[2]], datMICS$lat[urbIdx[2]])
cl1Pop = terra::extract(pop, cl1Coords, method="bilinear")
cl1NormPop = (log1p(cl1Pop) - popMeanCal) / popSDCal
cat("  Actual normPop at cluster 1:", round(cl1NormPop,4), "\n")

# Compare: do ALL clusters in the SAME stratum see the SAME int pt normPop values?
strat1 = datMICS$Stratum[urbIdx[2]]
clInSameStrat = which(datMICS$Stratum == strat1 & datMICS$urban)
cat("  Clusters in same stratum (", strat1, "):", length(clInSameStrat), "\n")

if(length(clInSameStrat) >= 2) {
  # Get obs indices for first two clusters in this stratum
  obsI_a = match(clInSameStrat[1], urbIdx) 
  obsI_b = match(clInSameStrat[2], urbIdx)
  normPop_a = XUrb_sub[nUrb * 0 + obsI_a, "normPop"]
  normPop_b = XUrb_sub[nUrb * 0 + obsI_b, "normPop"]
  cat("  Cluster a normPop (intPt 0):", round(normPop_a,4), "\n")
  cat("  Cluster b normPop (intPt 0):", round(normPop_b,4), "\n")
  cat("  Same? ", normPop_a == normPop_b, "\n")
  
  # Check across ALL int pts
  allSame = TRUE
  for(intI in 0:(KMICS-1)) {
    rowA = nUrb * intI + obsI_a
    rowB = nUrb * intI + obsI_b
    if(XUrb_sub[rowA, "normPop"] != XUrb_sub[rowB, "normPop"]) {
      allSame = FALSE
      break
    }
  }
  cat("  All K int pts identical for both clusters?", allSame, "\n")
}

# Also check the WEIGHTS
wUrban = intPtsMICS_orig$wUrban  # nStrata x K (before straightening)
# After straightenMICS and expansion in makeInputsMDM:
wUrbanFull = intPtsMICS$wUrban
stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrbanFull), each=numPerStratUrb))
wUrbanExpanded = wUrbanFull[stratIndexUrbW,]
cat("\n=== WEIGHTS ===\n")
cat("wUrbanExpanded dim:", dim(wUrbanExpanded), "\n")
cat("Row sums (should be 1): range:", round(range(rowSums(wUrbanExpanded)),4), "\n")

# Check: do clusters in the same stratum get the same weights?
if(length(clInSameStrat) >= 2) {
  obsI_a_idx = match(clInSameStrat[1], urbIdx)
  obsI_b_idx = match(clInSameStrat[2], urbIdx)
  cat("Weights for cluster a:", round(wUrbanExpanded[obsI_a_idx, 1:5],4), "...\n")
  cat("Weights for cluster b:", round(wUrbanExpanded[obsI_b_idx, 1:5],4), "...\n")
  cat("Same weights?", all(wUrbanExpanded[obsI_a_idx,] == wUrbanExpanded[obsI_b_idx,]), "\n")
}
