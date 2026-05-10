#!/usr/bin/env Rscript
# Simulate 2 populations/surveys with 4x cluster sizes using simData1 logic,
# then fit IID+nugget, FE+nugget, and BYM2+nugget models on sim 1.
#
# Default cluster sizes: MICS=16 HH, DHS=25 HH
# 4x cluster sizes:      MICS=64 HH, DHS=100 HH

library(TMB)

setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/simData.R")
source("code/modM_MIIDonly.R")
source("code/modBYM2.R")
source("code/modM_MSepMarg.R")

# Custom survey sampler that caps HH/cluster at available HH per EA
# (SUMMER::sampleClusterSurveys uses replace=FALSE and fails if HHperClust > nHH)
sampleClusterSurveysCapped = function(thisHHpop, HHperClust, clustpa) {
  hhDat = thisHHpop[[1]]
  uniqueEAIs = sort(unique(hhDat$eaIs))
  eaUrbs = hhDat$urban[match(uniqueEAIs, hhDat$eaIs)]
  eaAreas = hhDat$area[match(uniqueEAIs, hhDat$eaIs)]
  aggHHs = stats::aggregate(hhDat$nHH, by=list(hhDat$eaIs), FUN=sum)
  eaHHs = aggHHs$x
  
  sampledEAIs = numeric(0)
  inclusionProbs = numeric(0)
  
  for(j in 1:nrow(clustpa)) {
    thisArea = clustpa$area[j]
    nUrbEA = clustpa$EAUrb[j]
    nRurEA = clustpa$EARur[j]
    thisEAIs = uniqueEAIs[eaAreas == thisArea]
    thisEAUrbs = eaUrbs[eaAreas == thisArea]
    thisEAIsUrb = thisEAIs[thisEAUrbs]
    thisEAIsRur = thisEAIs[!thisEAUrbs]
    thisEAhhsUrb = eaHHs[eaAreas == thisArea][thisEAUrbs]
    thisEAhhsRur = eaHHs[eaAreas == thisArea][!thisEAUrbs]
    
    sampUrbEAIs = numeric(0); inclusionProbsUrb = numeric(0)
    sampRurEAIs = numeric(0); inclusionProbsRur = numeric(0)
    
    if(nUrbEA != 0) {
      inclusionProbsUrb = nUrbEA * thisEAhhsUrb / sum(thisEAhhsUrb)
      sampUrbEAIs = thisEAIsUrb[as.logical(sampling::UPmidzuno(inclusionProbsUrb))]
    }
    if(nRurEA != 0) {
      inclusionProbsRur = nRurEA * thisEAhhsRur / sum(thisEAhhsRur)
      sampRurEAIs = thisEAIsRur[as.logical(sampling::UPmidzuno(inclusionProbsRur))]
    }
    
    sampledEAIs = c(sampledEAIs, c(sampUrbEAIs, sampRurEAIs))
    inclusionProbs = c(inclusionProbs, c(inclusionProbsUrb, inclusionProbsRur))
  }
  
  # Sample HH within each EA, capped at min(HHperClust, available)
  hhSubdat = hhDat[hhDat$eaIs %in% sampledEAIs, ]
  hhIsTab = stats::aggregate(1:nrow(hhSubdat), by=list(eaIs=hhSubdat$eaIs),
    function(x) sample(x, min(HHperClust, length(x)), replace=FALSE))
  hhIs = sort(unlist(hhIsTab[,2]))
  hhDatSample = hhSubdat[hhIs, ]
  
  # Aggregate full EA-level data (for inclusion prob denominator)
  eaAgg = lapply(1:ncol(hhSubdat), function(j) {
    vn = names(hhSubdat)[j]
    stats::aggregate(hhSubdat[,j], by=list(hhSubdat$eaIs),
      FUN=function(x) if(vn %in% c("N","nHH","Z")) sum(x,na.rm=TRUE) else if(is.numeric(x)) mean(x,na.rm=TRUE) else x[1])$x
  })
  names(eaAgg) = names(hhSubdat)
  eaAgg = as.data.frame(eaAgg)
  
  # Aggregate sampled HH to cluster level
  survAgg = lapply(1:ncol(hhDatSample), function(j) {
    vn = names(hhDatSample)[j]
    stats::aggregate(hhDatSample[,j], by=list(hhDatSample$eaIs),
      FUN=function(x) if(vn %in% c("N","nHH","Z")) sum(x,na.rm=TRUE) else if(is.numeric(x)) mean(x,na.rm=TRUE) else x[1])$x
  })
  names(survAgg) = names(hhDatSample)
  surveyDat = as.data.frame(survAgg)
  surveyDat$pFineScalePrevalence = surveyDat$Z / surveyDat$N
  surveyDat$pFineScalePrevalence[surveyDat$N == 0] = 0
  surveyDat$includeProbEA = inclusionProbs[match(sampledEAIs, surveyDat$eaIs)]
  surveyDat$nHHsEA = eaAgg$nHH[match(surveyDat$eaIs, eaAgg$eaIs)]
  surveyDat$includeProbHH = surveyDat$nHH / surveyDat$nHHsEA
  surveyDat$samplingWeight = surveyDat$N / (surveyDat$includeProbHH * surveyDat$includeProbEA)
  
  list(surveyDat)
}

simSaveFile = "savedOutput/simStudy1/simPopsSurveys_4xCluster.RData"
saveFileIID = "savedOutput/testMM_4xCluster_IID_sim1.RData"
saveFileFE  = "savedOutput/testMM_4xCluster_FE_sim1.RData"
saveFileBYM2 = "savedOutput/testMM_4xCluster_BYM2_sim1.RData"

cat("===============================================\n")
cat("4x CLUSTER SIZE: SIMULATE + FIT 3 MODELS\n")
cat("===============================================\n\n")

# True parameter values (from simData1 defaults)
trueAlpha = -1.25
trueUrban = 1.00    # gamma
trueNormPop = 0.50   # betaRest[4]
trueSigmaEps = sqrt(1.5)
trueMargVar = 0.5

cat(sprintf("True: alpha=%.2f, urban=%.2f, normPop=%.2f, sigmaEps=%.4f\n\n",
            trueAlpha, trueUrban, trueNormPop, trueSigmaEps))

# ============================================================
# STEP 1: Simulate 2 populations/surveys with 4x cluster sizes
# ============================================================
if(file.exists(simSaveFile)) {
  cat("Loading existing 4x cluster simulation from", simSaveFile, "\n")
  load(simSaveFile)
} else {
  cat("Simulating 2 populations with 4x cluster sizes...\n")
  cat("  MICS: 64 HH/cluster (default 16)\n")
  cat("  DHS: 25 HH/cluster (default, unchanged)\n\n")
  
  # Replicate simData1 logic with 4x cluster sizes
  nsim = 2
  margVar = 0.5
  effRange = 200
  sigmaEpsilon = sqrt(1.5)
  beta0 = -1.25
  gamma = 1
  betaRest = c(0, 0, 0, 0.5)
  useThreshPopMat = TRUE
  fixPopPerHH = NULL
  
  # Use the same data objects that simData1 uses as defaults
  popMat = popMatNGAThresh
  targetPopMat = popMatNGAedThresh
  poppsub = poppsubNGAThresh
  easpaDat = easpaNGAed
  mesh = getSPDEMesh(doPlot=FALSE)
  
  set.seed(456)  # different seed from original (123) to get different populations
  
  # order data
  popMat = popMat[order(popMat$subarea),]
  poppsub = poppsub[order(poppsub$subarea),]
  
  # construct covariate offset
  print("Constructing offset based on covariates...")
  LLcoords = cbind(popMat$lon, popMat$lat)
  tempDesMat = getDesignMat(LLcoords, normalized=TRUE, useThreshPopMat=useThreshPopMat)
  
  load("savedOutput/global/covariatesNorm.RData")
  popVals = extract(pop, LLcoords, method="bilinear")
  
  load("savedOutput/global/popMeanSDCal.RData")
  normPop = (log1p(popVals) - popMeanCal) / popSDCal
  normPop[is.na(normPop)] = min(normPop, na.rm=TRUE)
  
  covRestMat = tempDesMat[,-c(1:3, 7)]
  covRestMat = cbind(covRestMat, normPop=normPop)
  offset = covRestMat %*% betaRest
  
  # aggregation info
  tempAreasFrom = popMat$subarea
  tempAreasTo = popMat$stratumMICS
  areasFrom = sort(unique(tempAreasFrom))
  areasToI = match(areasFrom, tempAreasFrom)
  areasTo = tempAreasTo[areasToI]
  
  # simulate
  surveysDHS = list()
  surveysMICS = list()
  subareaPops = list()
  areaPops = list()
  stratumPops = list()
  
  for(i in 1:nsim) {
    print(paste0("Simulating population ", i, "/", nsim))
    
    simPop = SUMMER::simPopSPDE(nsim=1, easpa=easpaDat, popMat=popMat, 
                                targetPopMat=targetPopMat, poppsub=poppsub, 
                                spdeMesh=mesh, margVar=margVar, 
                                sigmaEpsilon=sigmaEpsilon, effRange=effRange,
                                gamma=gamma, beta0=beta0, seed=NULL, nHHSampled=25,
                                stratifyByUrban=TRUE, subareaLevel=TRUE, offset=offset,
                                doFineScaleRisk=FALSE, doSmoothRisk=FALSE, 
                                min1PerSubarea=TRUE)
    
    # aggregate to stratum level
    stratPop = SUMMER::areaPopToArea(areaLevelPop=simPop$subareaPop,
                                     areasFrom=areasFrom, areasTo=areasTo,
                                     stratifyByUrban=TRUE, doFineScaleRisk=FALSE, 
                                     doSmoothRisk=FALSE)
    
    if(i == 1) {
      subareaPops = simPop$subareaPop$aggregationResults$pFineScalePrevalence
      areaPops = simPop$areaPop$aggregationResults$pFineScalePrevalence
      stratumPops = stratPop$aggregationResults$pFineScalePrevalence
    } else {
      subareaPops = cbind(subareaPops, simPop$subareaPop$aggregationResults$pFineScalePrevalence)
      areaPops = cbind(areaPops, simPop$areaPop$aggregationResults$pFineScalePrevalence)
      stratumPops = cbind(stratumPops, stratPop$aggregationResults$pFineScalePrevalence)
    }
    
    # generate surveys with 4x cluster sizes
    print(paste0("Generating 4x surveys for population ", i, "/", nsim))
    thisEApop = simPop$eaPop$eaDatList[1]
    thisHHpop = SUMMER::getHHpop(thisEApop, fixPopPerHH=fixPopPerHH)
    
    # DHS: keep at default 25 HH/cluster (only running MICS models)
    survDHS = SUMMER::sampleClusterSurveys(1, thisHHpop, HHperClust=25, 
                                            clustpaList=list(clustpaDHSed))
    
    # MICS: 64 HH/cluster (4 * 16), capped at available HH per EA
    tempClustpa = clustpaMICSed
    names(tempClustpa)[1] = "area"
    thisHHpop[[1]]$area = adm2ToStratumMICS(thisHHpop[[1]]$subarea)
    survMICS = sampleClusterSurveysCapped(thisHHpop, HHperClust=64, 
                                           clustpa=tempClustpa)
    
    surveysDHS = c(surveysDHS, survDHS)
    surveysMICS = c(surveysMICS, survMICS)
  }
  
  save(subareaPops, areaPops, stratumPops, surveysDHS, surveysMICS,
       file=simSaveFile)
  cat("Saved 4x cluster simulation to", simSaveFile, "\n")
}

# Quick data summary
edMICS1 = surveysMICS[[1]]
cat(sprintf("\nSim 1: %d MICS clusters, %d total obs (avg %.1f per cluster)\n",
            nrow(edMICS1), sum(edMICS1$N), mean(edMICS1$N)))

# ============================================================
# STEP 2: Fit 3 models on sim 1 (MICS only)
# ============================================================

# --- MODEL 1: IID + nugget ---
result_iid = NULL
if(file.exists(saveFileIID)) {
  cat("\nLoading saved IID+nugget results...\n")
  load(saveFileIID)
} else {
  cat("\n--- Fitting IID + nugget model ---\n")
  tryCatch({
    result_iid = fitMM_IIDonly(datMICS=edMICS1,
                                covariates=c("urban", "normPop"),
                                fixedEffectsOnly=FALSE,
                                getSDs=TRUE, doMCMC=FALSE)
    save(result_iid, file=saveFileIID)
    cat("IID+nugget completed and saved.\n")
  }, error = function(e) {
    cat("IID+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# --- MODEL 2: FE + nugget (no spatial RE) ---
result_fe = NULL
if(file.exists(saveFileFE)) {
  cat("\nLoading saved FE+nugget results...\n")
  load(saveFileFE)
} else {
  cat("\n--- Fitting FE + nugget model ---\n")
  tryCatch({
    result_fe = fitMM_IIDonly(datMICS=edMICS1,
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
  cat("\n--- Fitting BYM2 + nugget model ---\n")
  tryCatch({
    result_bym2 = fitMMMarg(datMICS=edMICS1,
                             covariates=c("urban", "normPop"),
                             getSDs=TRUE, doMCMC=FALSE)
    save(result_bym2, file=saveFileBYM2)
    cat("BYM2+nugget completed and saved.\n")
  }, error = function(e) {
    cat("BYM2+nugget FAILED:", conditionMessage(e), "\n")
  })
}

# ============================================================
# Extract estimates
# ============================================================
extractEst_IID = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed
  pfn = names(pf)
  se_f = sqrt(diag(SD0$cov.fixed))
  
  sr = summary(SD0, select="random")
  rn = names(SD0$par.random)
  alpha_est = sr[rn == "alpha", 1]
  alpha_se = sr[rn == "alpha", 2]
  beta_est = sr[rn == "beta", 1]
  beta_se = sr[rn == "beta", 2]
  
  lt = pf[pfn == "log_tau"]; se_lt = se_f[pfn == "log_tau"]
  sigma_u_est = exp(-0.5 * lt); sigma_u_se = 0.5 * sigma_u_est * se_lt
  
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  
  list(alpha=alpha_est, alpha_se=alpha_se,
       beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2],
       sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

extractEst_FE = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed
  pfn = names(pf)
  se_f = sqrt(diag(SD0$cov.fixed))
  
  alpha_est = as.numeric(pf[pfn == "alpha"])
  alpha_se = as.numeric(se_f[pfn == "alpha"])
  beta_est = as.numeric(pf[pfn == "beta"])
  beta_se = as.numeric(se_f[pfn == "beta"])
  
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  
  list(alpha=alpha_est, alpha_se=alpha_se,
       beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2],
       sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

extractEst_BYM2 = function(res) {
  SD0 = res$TMBsd
  pf = SD0$par.fixed
  pfn = names(pf)
  se_f = sqrt(diag(SD0$cov.fixed))
  
  sr = summary(SD0, select="random")
  rn = names(SD0$par.random)
  alpha_est = sr[rn == "alpha", 1]
  alpha_se = sr[rn == "alpha", 2]
  beta_est = sr[rn == "beta", 1]
  beta_se = sr[rn == "beta", 2]
  
  le = pf[pfn == "log_tauEps"]; se_le = se_f[pfn == "log_tauEps"]
  sigmaEps_est = exp(-0.5 * le); sigmaEps_se = 0.5 * sigmaEps_est * se_le
  
  list(alpha=alpha_est, alpha_se=alpha_se,
       beta1=beta_est[1], beta1_se=beta_se[1],
       beta2=beta_est[2], beta2_se=beta_se[2],
       sigmaEps=sigmaEps_est, sigmaEps_se=sigmaEps_se,
       time=res$totalTime)
}

# ============================================================
# Print combined results
# ============================================================
cat("\n===============================================\n")
cat("RESULTS: 4x CLUSTER SIZE (SPDE sim, sim 1)\n")
cat("===============================================\n\n")

fmtVal = function(est, se) {
  if(is.na(est)) return(sprintf("%20s", "NA"))
  sprintf("%7.4f (SE %5.4f)", est, se)
}

cat(sprintf("%-15s %7s  %-22s %-22s %-22s\n", "Parameter", "Truth", "IID+nugget", "FE+nugget", "BYM2+nugget"))
cat(sprintf("%-15s %7s  %-22s %-22s %-22s\n", "---------", "-----", "----------", "---------", "-----------"))

e_iid = if(!is.null(result_iid)) extractEst_IID(result_iid) else NULL
e_fe  = if(!is.null(result_fe))  extractEst_FE(result_fe) else NULL
e_bym2 = if(!is.null(result_bym2)) extractEst_BYM2(result_bym2) else NULL

getV = function(e, field, sefield) {
  if(is.null(e)) return(c(NA, NA))
  c(e[[field]], e[[sefield]])
}

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
