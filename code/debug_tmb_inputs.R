#!/usr/bin/env Rscript
# Diagnostic: check the exact data makeInputsMDM passes to TMB.
# Trace through the full code path for a BYM2 sim IID model
# and look for anything that could explain normPop ≈ 0.02.

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# Load BYM2 sim data
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

# Do the same name conversions fitMM_IIDonly does  
nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

cat("=== PART 1: Check makeInputsMDM output ===\n\n")

# Recreate exactly what fitMM_IIDonly does
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}

# Make inputs (this is what happens inside fitMM_IIDonly)
inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

cat("--- intPtsMICS$XUrb ---\n")
cat("Dimensions:", dim(inputsMDM$intPtsMICS$XUrb), "\n")
cat("Column names:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")
cat("Class:", class(inputsMDM$intPtsMICS$XUrb), "\n")

# Expected: nUrb*K rows, 5 columns (urban, access, elev, distRiversLakes, normPop)
nUrb = nrow(inputsMDM$intPtsMICS$wUrban)
nRur = nrow(inputsMDM$intPtsMICS$wRural)
KMICS = 100
cat("ysUrbMICS length:", length(inputsMDM$ysUrbMICS), "\n")
cat("ysRurMICS length:", length(inputsMDM$ysRurMICS), "\n")
cat("\nnUrb clusters:", nUrb, "\n")
cat("nRur clusters:", nRur, "\n")
cat("Expected XUrb rows:", nUrb * KMICS, "\n")
cat("Actual XUrb rows:", nrow(inputsMDM$intPtsMICS$XUrb), "\n")

cat("\n--- Urban column in XUrb (should all be 1 for urban int pts) ---\n")
urbCol = inputsMDM$intPtsMICS$XUrb[, "urban"]
cat("Min:", min(urbCol), "Max:", max(urbCol), "Mean:", mean(urbCol), "\n")
cat("Fraction == 1:", mean(urbCol == 1), "\n")
cat("Fraction == 0:", mean(urbCol == 0), "\n")
if(any(urbCol != 1)) {
  cat("WARNING: Not all urban int pts have urban=1!\n")
  cat("Values found:", sort(unique(urbCol)), "\n")
}

cat("\n--- Urban column in XRur (should all be 0 for rural int pts) ---\n")
rurUrbCol = inputsMDM$intPtsMICS$XRur[, "urban"]
cat("Min:", min(rurUrbCol), "Max:", max(rurUrbCol), "Mean:", mean(rurUrbCol), "\n")
cat("Fraction == 0:", mean(rurUrbCol == 0), "\n")
if(any(rurUrbCol != 0)) {
  cat("WARNING: Not all rural int pts have urban=0!\n")
  cat("Values found:", sort(unique(rurUrbCol)), "\n")
}

cat("\n--- normPop in XUrb ---\n")
npCol = inputsMDM$intPtsMICS$XUrb[, "normPop"]
cat("Mean:", mean(npCol), "SD:", sd(npCol), "\n")
cat("Range:", range(npCol), "\n")

cat("\n--- normPop in XRur ---\n")
npColR = inputsMDM$intPtsMICS$XRur[, "normPop"]
cat("Mean:", mean(npColR), "SD:", sd(npColR), "\n")
cat("Range:", range(npColR), "\n")

cat("\n--- areaidxlocUrbanMICS ---\n")
cat("Length:", length(inputsMDM$areaidxlocUrbanMICS), "\n")
cat("Expected (nUrb*K):", nUrb * KMICS, "\n")
cat("First nUrb values (0-indexed strata for each cluster):\n")
aidx = inputsMDM$areaidxlocUrbanMICS[1:nUrb]
cat("  Range:", range(aidx), "\n")
cat("  Unique:", length(unique(aidx)), "\n")
cat("  Table:\n")
print(table(aidx))

cat("\n--- areaidxlocRuralMICS ---\n")
cat("Length:", length(inputsMDM$areaidxlocRuralMICS), "\n")
aidxR = inputsMDM$areaidxlocRuralMICS[1:nRur]
cat("First nRur range:", range(aidxR), "\n")

cat("\n--- wUrban ---\n")
cat("Dimensions:", dim(inputsMDM$intPtsMICS$wUrban), "\n")
cat("Row sums (first 5):", rowSums(inputsMDM$intPtsMICS$wUrban)[1:5], "\n")
cat("All row sums == 1?", all(abs(rowSums(inputsMDM$intPtsMICS$wUrban) - 1) < 1e-10), "\n")

cat("\n=== PART 2: Subset covariates (like fitMM_IIDonly) ===\n\n")

# Subset to urban + normPop
allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
covariates = c("urban", "normPop")
keepIdx = which(allCovNames %in% covariates)
cat("allCovNames:", allCovNames, "\n")
cat("keepIdx:", keepIdx, "\n")
cat("Selected columns:", allCovNames[keepIdx], "\n")

XUrbSub = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
XRurSub = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]
cat("XUrbSub dims:", dim(XUrbSub), "\n")
cat("XUrbSub colnames:", colnames(XUrbSub), "\n")

cat("\n=== PART 3: Check per-cluster covariates ===\n\n")

# For each urban cluster, what normPop values does it see?
# The first nUrb rows of XUrbSub are for intPt 0
# Rows nUrb+1:2*nUrb are for intPt 1, etc.
# For cluster i, its covariate at intPt j is at row nUrb*j + i (0-indexed)

# Sort edMICS the same way as makeInputsMDM does
edMICS_sorted = sortByCol(edMICS, "Stratum", inputsMDM$intPtsMICS$strataMICS)
edMICS_urb = edMICS_sorted[edMICS_sorted$urban, ]
edMICS_rur = edMICS_sorted[!edMICS_sorted$urban, ]

cat("Checking cluster 1 (0-indexed: 0)\n")
clust0_np = sapply(0:(KMICS-1), function(j) XUrbSub[nUrb*j + 1, "normPop"])
clust0_wts = inputsMDM$intPtsMICS$wUrban[1, ]
cat("  Actual cluster stratum:", edMICS_urb$Stratum[1], "\n")
cat("  Model area index:", aidx[1], "\n")
cat("  normPop across", KMICS, "int pts: mean=", mean(clust0_np), "sd=", sd(clust0_np), "\n")
cat("  normPop range:", range(clust0_np), "\n")
cat("  Weights sum:", sum(clust0_wts), "\n")
cat("  Weighted mean normPop:", sum(clust0_wts * clust0_np), "\n")

# Check a second cluster in the SAME stratum
strat1_idx = which(aidx == aidx[1])
if(length(strat1_idx) > 1) {
  i2 = strat1_idx[2]
  clust1_np = sapply(0:(KMICS-1), function(j) XUrbSub[nUrb*j + i2, "normPop"])
  cat("\nChecking cluster", i2, "(same stratum)\n")
  cat("  normPop identical to cluster 1?", identical(clust0_np, clust1_np), "\n")
}

cat("\n=== PART 4: Check y/n ordering matches covariate ordering ===\n\n")

# The y_c/n_c from inputsMDM
cat("ysUrbMICS (first 10):", inputsMDM$ysUrbMICS[1:10], "\n")
cat("nsUrbMICS (first 10):", inputsMDM$nsUrbMICS[1:10], "\n")
cat("edMICS_urb$ys (first 10):", edMICS_urb$ys[1:10], "\n")
cat("edMICS_urb$ns (first 10):", edMICS_urb$ns[1:10], "\n")
cat("ys match?", identical(inputsMDM$ysUrbMICS, edMICS_urb$ys), "\n")
cat("ns match?", identical(inputsMDM$nsUrbMICS, edMICS_urb$ns), "\n")

# Check strata match
cat("\nedMICS_urb strata (first 10):", edMICS_urb$Stratum[1:10], "\n")
cat("areaidxloc strata (first 10):", inputsMDM$intPtsMICS$strataMICS[aidx[1:10]+1], "\n")
cat("Strata match?", all(edMICS_urb$Stratum == inputsMDM$intPtsMICS$strataMICS[aidx[1:nUrb]+1]), "\n")

cat("\n=== PART 5: Actual cluster normPop vs integration point normPop ===\n\n")

# Get normPop at true cluster locations from the raster
popMat = popsResults[[1]]$popMat
popMeanCal = 7.3259
popSDCal = 1.9238

clusterCoords = cbind(edMICS_urb$east, edMICS_urb$north)
popAtClusters = raster::extract(popMat, clusterCoords)
normPopAtClusters = (log1p(popAtClusters) - popMeanCal) / popSDCal

cat("Actual cluster normPop (first 10):", round(normPopAtClusters[1:10], 4), "\n")
cat("Weighted intPt normPop (first 10):")
for(i in 1:10) {
  np_i = sapply(0:(KMICS-1), function(j) XUrbSub[nUrb*j + i, "normPop"])
  wt_i = inputsMDM$intPtsMICS$wUrban[i, ]
  cat(" ", round(sum(wt_i * np_i), 4))
}
cat("\n")

# How different are they on average?
actual_np = normPopAtClusters
weighted_np = sapply(1:nUrb, function(i) {
  np_i = sapply(0:(KMICS-1), function(j) XUrbSub[nUrb*j + i, "normPop"])
  wt_i = inputsMDM$intPtsMICS$wUrban[i, ]
  sum(wt_i * np_i)
})
cat("Mean actual cluster normPop:", mean(actual_np), "\n")
cat("Mean weighted intPt normPop:", mean(weighted_np), "\n")
cat("Correlation:", cor(actual_np, weighted_np), "\n")
cat("Mean difference:", mean(actual_np - weighted_np), "\n")

cat("\n=== PART 6: Build TMB and evaluate NLL at different beta_normPop ===\n\n")

# Create the TMB object exactly like fitMM_IIDonly
# Then manually evaluate the marginal NLL at different beta_normPop values

# Load the saved result to get the estimated outer parameters
load("savedOutput/testMM_BYM2sim_IIDnugget_sim1.RData")
SD0 = result_opt$TMBsd

# Get estimated outer params: log_tau, log_tauEps
pf = SD0$par.fixed
pfn = names(pf)
est_log_tau = as.numeric(pf[pfn == "log_tau"])
est_log_tauEps = as.numeric(pf[pfn == "log_tauEps"])
cat("Estimated log_tau:", est_log_tau, "\n")
cat("Estimated log_tauEps:", est_log_tauEps, "\n")

# Get estimated random effects
sr = summary(SD0, select="random")
rn = names(SD0$par.random)
est_alpha = sr[rn == "alpha", 1]
est_beta = sr[rn == "beta", 1]
cat("Estimated alpha:", est_alpha, "\n")
cat("Estimated beta:", est_beta, "\n")

# Build TMB object with the subsetted data
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(inputsMDM$areaidxlocUrbanMICS, inputsMDM$areaidxlocRuralMICS)) + 1

data_full = list(
  y_iUrbanMICS=inputsMDM$ysUrbMICS,
  y_iRuralMICS=inputsMDM$ysRurMICS,
  n_iUrbanMICS=inputsMDM$nsUrbMICS,
  n_iRuralMICS=inputsMDM$nsRurMICS,
  areaidxlocUrbanMICS=inputsMDM$areaidxlocUrbanMICS,
  areaidxlocRuralMICS=inputsMDM$areaidxlocRuralMICS,
  X_betaUrbanMICS=XUrbSub,
  X_betaRuralMICS=XRurSub,
  wUrbanMICS=inputsMDM$intPtsMICS$wUrban,
  wRuralMICS=inputsMDM$intPtsMICS$wRural,
  nAreas=nAreas,
  alpha_pri=c(0, 100^2),
  beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau,
  lambdaTauEps=lambdaTauEps
)

cat("data_full dimensions check:\n")
cat("  y_iUrbanMICS:", length(data_full$y_iUrbanMICS), "\n")
cat("  y_iRuralMICS:", length(data_full$y_iRuralMICS), "\n")
cat("  X_betaUrbanMICS:", dim(data_full$X_betaUrbanMICS), "\n")
cat("  X_betaRuralMICS:", dim(data_full$X_betaRuralMICS), "\n")
cat("  wUrbanMICS:", dim(data_full$wUrbanMICS), "\n")
cat("  wRuralMICS:", dim(data_full$wRuralMICS), "\n")
cat("  areaidxlocUrbanMICS length:", length(data_full$areaidxlocUrbanMICS), "\n")
cat("  nAreas:", nAreas, "\n")

unloadDynlibs()
if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
  compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modM_MIIDonly"))

# Create TMB object with random effects
rand_effs = c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS')

tmb_params = list(
  log_tau = est_log_tau,
  log_tauEps = est_log_tauEps,
  alpha = est_alpha,
  beta = est_beta,
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
  nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
)

obj = MakeADFun(data=data_full,
                parameters=tmb_params,
                random=rand_effs,
                hessian=TRUE,
                DLL='modM_MIIDonly')

# Evaluate the marginal NLL at the estimated outer params
# (this does inner optimization over alpha, beta, u, nuggets)
cat("\nMarginal NLL at estimated (log_tau, log_tauEps):")
nll_est = obj$fn(c(est_log_tau, est_log_tauEps))
cat(nll_est, "\n")

# Get the inner optimization result (random effects at mode)
cat("\nRandom effects at mode:\n")
last_par = obj$env$last.par.best
rn_all = names(last_par)
cat("  alpha:", last_par[rn_all == "alpha"], "\n")
cat("  beta:", last_par[rn_all == "beta"], "\n")

# Now, let's try a profile: for fixed outer params, 
# evaluate the NLL at a grid of beta_normPop values
# We'll do this by creating a new TMB object where beta[2] is NOT random
cat("\n=== PART 7: Profile NLL over beta_normPop ===\n\n")

# Map beta[2] as fixed, everything else stays random
beta_normPop_grid = seq(-0.5, 1.5, by=0.1)
nll_profile = rep(NA, length(beta_normPop_grid))

for(gi in seq_along(beta_normPop_grid)) {
  this_normPop = beta_normPop_grid[gi]
  
  # Create TMB with beta[2] fixed via map
  tmb_params2 = list(
    log_tau = est_log_tau,
    log_tauEps = est_log_tauEps,
    alpha = est_alpha,
    beta = c(est_beta[1], this_normPop),
    u_spatial = rep(0, nAreas),
    nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
    nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  # Map: beta[2] fixed (NA), beta[1] free (factor 1)
  mapList2 = list(beta = factor(c(1, NA)))
  rand_effs2 = c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS')
  
  tryCatch({
    obj2 = MakeADFun(data=data_full,
                     parameters=tmb_params2,
                     random=rand_effs2,
                     map=mapList2,
                     hessian=FALSE,
                     DLL='modM_MIIDonly',
                     silent=TRUE)
    
    nll_profile[gi] = obj2$fn(c(est_log_tau, est_log_tauEps))
    cat(sprintf("beta_normPop = %5.2f: NLL = %10.4f", this_normPop, nll_profile[gi]))
    
    # Get the optimized alpha and beta[1]
    lp = obj2$env$last.par.best
    rn2 = names(lp)
    cat(sprintf("  alpha=%6.3f  beta_urban=%6.3f\n",
                lp[rn2 == "alpha"], lp[rn2 == "beta"]))
    
  }, error = function(e) {
    cat(sprintf("beta_normPop = %5.2f: FAILED - %s\n", this_normPop, conditionMessage(e)))
  })
}

cat("\n--- Profile summary ---\n")
cat("Minimum NLL at beta_normPop =", beta_normPop_grid[which.min(nll_profile)], "\n")
cat("NLL at truth (0.5):", nll_profile[beta_normPop_grid == 0.5], "\n")
cat("NLL at estimate (~0):", nll_profile[which.min(abs(beta_normPop_grid))], "\n")
cat("NLL difference (truth - est):", 
    nll_profile[beta_normPop_grid == 0.5] - nll_profile[which.min(abs(beta_normPop_grid))], "\n")

cat("\nDone.\n")
