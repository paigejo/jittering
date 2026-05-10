setwd("c:/Users/jpaige/git/jittering")
suppressMessages(suppressWarnings(source("code/setup.R")))
source("code/simData.R")
cat("setup loaded\n")

cat("\n=== PART 1: RECONSTRUCT THE EXACT SIMULATION DGP ===\n")
# Re-run the simulation offset computation exactly as simData1BYM2 does
popMat = popMatNGAThresh
popMat = popMat[order(popMat$subarea),]

LLcoords_sim = cbind(popMat$lon, popMat$lat)
tempDesMat = getDesignMat(LLcoords_sim, normalized=TRUE, useThreshPopMat=TRUE)
cat("Design matrix columns:", colnames(tempDesMat), "\n")

load("savedOutput/global/covariatesNorm.RData")
popVals_sim = extract(pop, LLcoords_sim, method="bilinear")

load("savedOutput/global/popMeanSDCal.RData")
# NOTE: simData1BYM2 hardcodes popMeanCal/popSDCal regardless of useThreshPopMat flag
normPop_sim = (log1p(popVals_sim) - popMeanCal) / popSDCal
normPop_sim[is.na(normPop_sim)] = min(normPop_sim, na.rm=TRUE)

covRestMat = tempDesMat[,-c(1:3, 7)]  # remove int, pop, urb, urbanicity
covRestMat = cbind(covRestMat, normPop=normPop_sim)
cat("covRestMat columns:", colnames(covRestMat), "\n")

betaRest = c(0, 0, 0, 0.5)
offset_sim = covRestMat %*% betaRest
cat("Sim offset: mean=", round(mean(offset_sim),4), " = 0.5 * mean(normPop)=", round(0.5*mean(normPop_sim),4), "\n")
cat("  offset is just 0.5*normPop since other betas are 0:", all.equal(as.numeric(offset_sim), 0.5*normPop_sim), "\n")

cat("\n=== PART 2: CHECK pFineScalePrevalence vs DGP ===\n")
# The DGP: logit(p_pixel) = BYM2_spatial + beta0 + gamma*urban + offset
# pFineScalePrevalence = expit(logit(p_pixel))
# We can't recover BYM2_spatial directly, but we can check partial correlations

load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS1 = surveysMICS[[1]]
cat("MICS clusters:", nrow(edMICS1), "\n")

# Get normPop at cluster locations using the SAME computation as simulation
clCoords = cbind(edMICS1$lon, edMICS1$lat)
popVals_cl = extract(pop, clCoords, method="bilinear")
normPop_cl = (log1p(popVals_cl) - popMeanCal) / popSDCal
normPop_cl[is.na(normPop_cl)] = min(normPop_cl, na.rm=TRUE)

# Get the pixel-level normPop that the simulation used for each cluster
# Clusters are at pixel grid locations, so find the matching pixel index
pixKey = paste(round(popMat$lon,8), round(popMat$lat,8))
clKey = paste(round(edMICS1$lon,8), round(edMICS1$lat,8))
pixIdx = match(clKey, pixKey)
cat("Cluster-to-pixel matching: ", sum(!is.na(pixIdx)), "of", nrow(edMICS1), "matched\n")

normPop_simAtCluster = normPop_sim[pixIdx]
cat("normPop discrepancy (raster extract vs simulation grid):\n")
cat("  max abs diff:", max(abs(normPop_cl - normPop_simAtCluster), na.rm=TRUE), "\n")

# Check logit(pFineScalePrevalence)
logitP = log(edMICS1$pFineScalePrevalence / (1 - edMICS1$pFineScalePrevalence))
# DGP: logitP = BYM2_spatial + beta0 + gamma*urban + 0.5*normPop
# Residual after removing known fixed effects:
resid_partial = logitP - (-1.25) - 1.0 * edMICS1$urban - 0.5 * normPop_simAtCluster
cat("\nResidual (logitP - beta0 - gamma*urban - 0.5*normPop):\n")
cat("  This should be ~BYM2 spatial (area-level RE). Mean:", round(mean(resid_partial, na.rm=TRUE),4), "\n")
cat("  SD:", round(sd(resid_partial, na.rm=TRUE),4), " (expect ~ sigmaBYM2=", round(sqrt(0.5),4), ")\n")
cat("  Cor with normPop:", round(cor(resid_partial, normPop_simAtCluster, use="complete"),4), "(should be ~0)\n")
cat("  Cor with urban:", round(cor(resid_partial, edMICS1$urban, use="complete"),4), "(should be ~0)\n")

# Also check: the sign of the normPop relationship in the RAW data
cat("\nDirect relationship checks:\n")
cat("  Cor(logitP, normPop):", round(cor(logitP, normPop_simAtCluster, use="complete"),4), "\n")
cat("  Cor(Z/N, normPop):", round(cor(edMICS1$Z/edMICS1$N, normPop_simAtCluster, use="complete"),4), "\n")

# GLM with all covariates from simulation
cat("\nGLM check (cluster-level, using exact simulation normPop):\n")
df = data.frame(y=edMICS1$Z, n=edMICS1$N, urban=as.numeric(edMICS1$urban), 
                normPop=normPop_simAtCluster)
fit = glm(cbind(y, n-y) ~ urban + normPop, data=df, family=binomial)
cat("  Intercept:", round(coef(fit)[1],4), "(truth: -1.25)\n")
cat("  urban:", round(coef(fit)[2],4), "(truth: 1.00)\n")
cat("  normPop:", round(coef(fit)[3],4), "(truth: 0.50)\n")

# GLM using raster-extracted normPop (what the model actually uses)
df2 = data.frame(y=edMICS1$Z, n=edMICS1$N, urban=as.numeric(edMICS1$urban), 
                 normPop=normPop_cl)
fit2 = glm(cbind(y, n-y) ~ urban + normPop, data=df2, family=binomial)
cat("\nGLM check (raster-extracted normPop):\n")
cat("  Intercept:", round(coef(fit2)[1],4), "(truth: -1.25)\n")
cat("  urban:", round(coef(fit2)[2],4), "(truth: 1.00)\n")
cat("  normPop:", round(coef(fit2)[3],4), "(truth: 0.50)\n")

cat("\n=== PART 3: CHECK WHAT getDesignMatPopNorm RETURNS ===\n")
# This is what the true-coords model uses for covariates
covMat = getDesignMatPopNorm(clCoords, useThreshPopMat=TRUE)
cat("getDesignMatPopNorm columns:", colnames(covMat), "\n")
cat("  'urb' column: mean=", round(mean(covMat[,"urb"]),4), 
    " matches survey urban?", all.equal(as.numeric(covMat[,"urb"]), as.numeric(edMICS1$urban)), "\n")

# Check: normPop from getDesignMatPopNorm vs our manual computation
# getDesignMatPopNorm may compute normPop differently!
cat("  'pop' column (if exists): mean=", round(mean(covMat[,"pop"]),4), "\n")

# Check how the true-coords model computes normPop
cat("\n=== PART 4: CHECK MODEL COVARIATE MATRIX ===\n")
# In the true-coords script, we do:
# covMat = getDesignMatPopNorm(clusterCoords, useThreshPopMat=TRUE)
# popVals = terra::extract(pop, clusterCoords, method="bilinear")
# normPop = (log1p(popVals) - popMeanCal) / popSDCal
# clusterCovs = cbind(urban=covMat[,"urb"], access=covMat[,"access"], 
#                     elev=covMat[,"elev"], distRiversLakes=covMat[,"distRiversLakes"],
#                     normPop=normPop)

# But wait, getDesignMatPopNorm might ALSO have a pop-based column
# Let me check if normPop from the function differs from manual extraction
cat("normPop manual (log1p raster):", round(mean(normPop_cl),4), "\n")
cat("normPop at matching pixels:", round(mean(normPop_simAtCluster, na.rm=TRUE),4), "\n")

# What does getDesignMat (used in the simulation) return?
tempDesMat_cl = getDesignMat(clCoords, normalized=TRUE, useThreshPopMat=TRUE)
cat("getDesignMat columns:", colnames(tempDesMat_cl), "\n")
cat("getDesignMat 'pop' values at clusters: mean=", round(mean(tempDesMat_cl[,"pop"]),4), "\n")

# CRITICAL: what does the simulation use as the offset vs what the model uses as normPop?
# Simulation covRestMat: removes int, pop, urb, urbanicity from getDesignMat, adds normPop
# Model: uses urban (from getDesignMatPopNorm "urb") and normPop (from manual raster extract)
# Are these the SAME normPop?
covRestMat_cl = tempDesMat_cl[,-c(1:3, 7)]
covRestMat_cl = cbind(covRestMat_cl, normPop=normPop_cl)
cat("\ncovRestMat columns at clusters:", colnames(covRestMat_cl), "\n")
cat("  access range:", round(range(covRestMat_cl[,"access"]),4), "\n")
cat("  elev range:", round(range(covRestMat_cl[,"elev"]),4), "\n")
cat("  distRiversLakes range:", round(range(covRestMat_cl[,"distRiversLakes"]),4), "\n")
cat("  normPop range:", round(range(covRestMat_cl[,"normPop"]),4), "\n")

cat("\n=== PART 5: CHECK INTEGRATION POINT NORMPP STRUCTURE ===\n")
library(TMB)

# Load the MICS integration points used by the models
intPtFile = "savedOutput/global/intPtsMICS_100.RData"
if(file.exists(intPtFile)) {
  load(intPtFile)
  cat("Loaded intPtsMICS_100\n")
  cat("XUrb cols:", names(intPtsMICS$XUrb), "\n")
  cat("XUrb dim:", dim(intPtsMICS$XUrb), "\n")
  cat("XRur dim:", dim(intPtsMICS$XRur), "\n")
  
  # The MICS integration points have K points per stratum (not per cluster)
  # XUrb dimensions: [K * nStrata] x nVar
  KMICS = 100
  nStrata = nrow(intPtsMICS$XUrb) / KMICS
  cat("nStrata from XUrb:", nStrata, "\n")
  
  cat("\nnormPop in integration points:\n")
  cat("  XUrb normPop: mean=", round(mean(intPtsMICS$XUrb$normPop),4), 
      " sd=", round(sd(intPtsMICS$XUrb$normPop),4), "\n")
  cat("  XRur normPop: mean=", round(mean(intPtsMICS$XRur$normPop),4),
      " sd=", round(sd(intPtsMICS$XRur$normPop),4), "\n")
  
  cat("\nnormPop at ACTUAL cluster locations:\n")
  cat("  All clusters: mean=", round(mean(normPop_simAtCluster),4),
      " sd=", round(sd(normPop_simAtCluster),4), "\n")
  cat("  Urban clusters: mean=", round(mean(normPop_simAtCluster[edMICS1$urban]),4),
      " sd=", round(sd(normPop_simAtCluster[edMICS1$urban]),4), "\n")
  cat("  Rural clusters: mean=", round(mean(normPop_simAtCluster[!edMICS1$urban]),4),
      " sd=", round(sd(normPop_simAtCluster[!edMICS1$urban]),4), "\n")
  
  # Key question: After makeInputsMDM expands the integration points,
  # ALL clusters in the same stratum share the SAME K normPop values.
  # So within a stratum, there's NO variation in normPop to estimate the coefficient.
  # The normPop coefficient can only be estimated from BETWEEN-stratum variation.
  
  # Compute between-stratum vs within-stratum normPop variation
  edMICS1$Stratum = adm2ToStratumMICS(edMICS1$subarea)
  stratMeans = tapply(normPop_simAtCluster, edMICS1$Stratum, mean)
  stratSDs = tapply(normPop_simAtCluster, edMICS1$Stratum, sd)
  cat("\nnormPop: between-stratum SD of means:", round(sd(stratMeans),4), "\n")
  cat("normPop: mean within-stratum SD:", round(mean(stratSDs, na.rm=TRUE),4), "\n")
  
  # Also split by urban/rural
  stratMeansUrb = tapply(normPop_simAtCluster[edMICS1$urban], edMICS1$Stratum[edMICS1$urban], mean)
  stratMeansRur = tapply(normPop_simAtCluster[!edMICS1$urban], edMICS1$Stratum[!edMICS1$urban], mean)
  cat("Urban normPop: between-stratum SD:", round(sd(stratMeansUrb, na.rm=TRUE),4), "\n")
  cat("Rural normPop: between-stratum SD:", round(sd(stratMeansRur, na.rm=TRUE),4), "\n")
  
  # What normPop values does the integration point model actually see?
  # It sees the STRATUM-level values from the integration points, weighted by population
  # Let's look at the weighted average normPop per stratum in the integration points
  wUrb = intPtsMICS$wUrban  # nStrata x K
  wRur = intPtsMICS$wRural
  
  # For each stratum, compute weighted mean normPop across K int pts
  intPtNormPopUrb = matrix(intPtsMICS$XUrb$normPop, nrow=nStrata, ncol=KMICS)
  intPtNormPopRur = matrix(intPtsMICS$XRur$normPop, nrow=nStrata, ncol=KMICS)
  
  wtMeanNormPopUrb = rowSums(wUrb * intPtNormPopUrb)
  wtMeanNormPopRur = rowSums(wRur * intPtNormPopRur)
  
  cat("\nWeighted mean normPop per stratum (integration points):\n")
  cat("  Urban: mean=", round(mean(wtMeanNormPopUrb),4), " sd=", round(sd(wtMeanNormPopUrb),4), "\n")
  cat("  Rural: mean=", round(mean(wtMeanNormPopRur),4), " sd=", round(sd(wtMeanNormPopRur),4), "\n")
  
  cat("\nActual cluster normPop per stratum:\n")
  cat("  Urban: mean=", round(mean(stratMeansUrb, na.rm=TRUE),4), " sd=", round(sd(stratMeansUrb, na.rm=TRUE),4), "\n")
  cat("  Rural: mean=", round(mean(stratMeansRur, na.rm=TRUE),4), " sd=", round(sd(stratMeansRur, na.rm=TRUE),4), "\n")
} else {
  cat("Integration point file not found:", intPtFile, "\n")
}

cat("\n=== PART 6: ALL SAVED BYM2 SIM RESULTS ===\n")
resultFiles = list(
  "BYM2sim IID (intpts, DHS+MICS)" = "savedOutput/testMM_BYM2sim_IIDnugget_progress.RData",
  "BYM2sim2 IID (intpts, MICS)" = "savedOutput/testMM_BYM2sim2_IIDnugget.RData",
  "BYM2sim2 nuggetOnly (MICS)" = "savedOutput/testMM_BYM2sim2_nuggetOnly.RData",
  "BYM2sim nuggetOnly" = "savedOutput/testMM_BYM2sim_nuggetOnly_progress.RData",
  "BYM2sim BYM2fixedHyper sim1" = "savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim1.RData",
  "BYM2sim BYM2fixedHyper sim2" = "savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim2.RData",
  "BYM2sim trueCoords IID" = "savedOutput/testMM_BYM2sim_trueCoords_IID_sim1.RData",
  "BYM2sim trueCoords FE" = "savedOutput/testMM_BYM2sim_trueCoords_FE_sim1.RData",
  "BYM2sim trueCoords BYM2" = "savedOutput/testMM_BYM2sim_trueCoords_BYM2_sim1.RData"
)

cat(sprintf("%-40s %10s %10s %10s\n", "Model", "alpha", "urban", "normPop"))
cat(sprintf("%-40s %10s %10s %10s\n", "-----", "-----", "-----", "-------"))
cat(sprintf("%-40s %10.4f %10.4f %10.4f\n", "TRUTH", -1.25, 1.00, 0.50))

for(nm in names(resultFiles)) {
  f = resultFiles[[nm]]
  if(!file.exists(f)) { cat(sprintf("%-40s %10s\n", nm, "FILE NOT FOUND")); next }
  e = new.env()
  load(f, envir=e)
  found = FALSE
  for(o in ls(envir=e)) {
    v = get(o, envir=e)
    if(is.list(v) && "TMBsd" %in% names(v)) {
      SD0 = v$TMBsd
      # Try random effects first
      sr = summary(SD0, select="random")
      rn = names(SD0$par.random)
      if("alpha" %in% rn) {
        aEst = sr[rn == "alpha", 1]
        bEst = sr[rn == "beta", 1]
      } else {
        # Fixed effects
        sf = summary(SD0, select="fixed")
        fn = names(SD0$par.fixed)
        aEst = sf[fn == "alpha", 1]
        bEst = sf[fn == "beta", 1]
      }
      if(length(bEst) >= 2) {
        cat(sprintf("%-40s %10.4f %10.4f %10.4f\n", nm, aEst, bEst[1], bEst[2]))
      } else if(length(bEst) == 1) {
        cat(sprintf("%-40s %10.4f %10.4f %10s\n", nm, aEst, bEst[1], "N/A"))
      }
      found = TRUE
    }
  }
  if(!found) cat(sprintf("%-40s %10s\n", nm, "NO TMBsd FOUND"))
}

cat("\nDone.\n")
