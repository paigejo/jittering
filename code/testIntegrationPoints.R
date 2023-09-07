# script contains methods for testing integration resolution
# See "Sensitivity analysis in Bayesian generalized linear mixed models for binary data"
# a 2011 paper by Roos and Held in 2011 for more info on using Hellinger distance

# Calculates H(Bern(.5), Bern(q)). Max is .541 at q=0 and 1, min is 0 at q=.5
HellingerBern = function(q) {
  sqrt(1 - sqrt((1-q)/2) - sqrt(q/2))
}

# given H(Bern(.5), Bern(q)), calculates q
HellingerBernGetQ = function(h) {
  (1 + sqrt(1 - 4*((1 - h^2)^2 - 1/2)^2))/2
}

plotHellingerBern = function() {
  qs = seq(0, 1, l=1000)
  
  ys = HellingerBern(qs)
  
  plot(qs, ys, main="Hellinger Distance, H(Bern(.5), Bern(q))", xlab="q", ylab="H", type="l")
}

plotHellingerBernGetQ = function() {
  hs = seq(0, HellingerBern(1), l=1000)
  
  ys = HellingerBernGetQ(hs)
  
  plot(hs, ys, main="q", xlab="Hellinger Distance, H(Bern(.5), Bern(q))", ylab="q", type="l")
}

# https://en.wikipedia.org/wiki/Hellinger_distance
HellingerMVN = function(mu1, Prec1, mu2, Prec2, eigen1=NULL, eigen2=NULL) {
  # get eigendecompositions of precision matrices
  if(is.null(eigen1)) {
    eigen1 = eigen(Prec1)
  }
  if(is.null(eigen2)) {
    eigen2 = eigen(Prec2)
  }
  
  # convert to eigendecompositions of the covariance matrices
  eigen1$values = 1/eigen1$values
  eigen2$values = 1/eigen2$values
  
  # calculate eigendecomposition of average of covariance matrices
  Sig1p1 = eigen1$vectors %*% diag(eigen1$values, nrow=length(eigen1$values)) %*% t(eigen1$vectors) + 
    eigen2$vectors %*% diag(eigen2$values, nrow=length(eigen2$values)) %*% t(eigen2$vectors)
  eigen1p2 = eigen(Sig1p1/2)
  
  # calculate log determinants
  ldet1 = sum(log(eigen1$values))
  ldet2 = sum(log(eigen2$values))
  ldet1p2 = sum(log(eigen1p2$values))
  
  # calculate log of the determinant term
  ldetTerm = (1/4) * (ldet1 + ldet2) - (1/2) * ldet1p2
  
  # calculate exponential argument
  muDiff = mu1 - mu2
  rightHalf = t(eigen1p2$vectors) %*% muDiff
  expArg = -(1/8) * t(rightHalf) %*% diag(eigen1p2$values) %*% rightHalf
  
  # Battacharyya coefficient
  BC = exp(ldetTerm + expArg)
  
  # Hellinger distance
  d = sqrt(1 - BC)
  
  d
}

# https://en.wikipedia.org/wiki/Hellinger_distance
HellingerUniveriate = function(samples1, samples2) {
  if(any(is.na(c(samples1, samples2)))) {
    warning("some samples are NA, returning NA Hellinger distance")
    return(NA)
  }
  require(kdensity)
  
  d1 = kdensity(samples1, start="gaussian", kernel="gaussian")
  d2 = kdensity(samples2, start="gaussian", kernel="gaussian")
  
  # Battacharyya coefficient
  BC = integrate(function(x) {sqrt(d1(x) * d2(x))}, lower=-Inf, upper=Inf)
  BC = min(c(max(c(BC$value, 0)), 1)) # make sure BC is in [0, 1]
  Hdist = sqrt(1 - BC)
  
  Hdist
}

# TODO: do we calculate Hellinger distance from samples, or from evaluation of density?

# nSamples: if NULL, use presaved prediction samples. Otherwise, a number 
#     specifying the number of samples to take
testResModels = function(allRes=c(50, 75, 100, 125, 150, 175, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
                         nSamples=NULL, saveGridPreds=FALSE) {
  
  # first get parameter predictions and computation times in hours
  # testMat = c()
  # for(i in 4:length(allRes)) {
  #   thisRes = allRes[i]
  # 
  # 
  #   out = load(paste0("savedOutput/ed/gridPreds2_", thisRes, "_adm2Cov.RData"))
  #   thisMeans = rowMeans(gridPreds$fixedMat)
  #   testMat = rbind(testMat, c(thisMeans))
  # }
  # row.names(testMat) = allRes[4:length(allRes)]
  # 
  # print(round(testMat, 3))
  # # 
  # # browser()
  # # 
  # # Now get means and precisions of each posterior
  # muMat = c()
  # Precs = list()
  # totalTimes = c()
  # sdTimes = c()
  # for(i in 1:length(allRes)) {
  #   thisRes = allRes[i]
  # 
  #   out = load(paste0("savedOutput/ed/fit2_", thisRes, "_adm2Cov.RData"))
  #   thisMu = summary(SD0)[,1]
  #   thisPrec = SD0$jointPrecision
  #   totalTimes[i] = totalTime/60/60
  #   sdTimes[i] = sdTime/60/60
  # 
  #   # permute mu so that it corresponds correctly to the joint precision
  #   muNames = names(thisMu)
  #   PrecNames = colnames(thisPrec)
  #   # make sure the names come in this order:
  #   permI = c(which(muNames == "log_tau"),
  #             which(muNames == "logit_phi"),
  #             which(muNames == "log_tauEps"))
  #   fixedNames = c("log_tau", "logit_phi", "log_tauEps")
  #   thisMu = thisMu[muNames %in% fixedNames][permI]
  # 
  #   thisPrec = thisPrec[PrecNames %in% fixedNames, PrecNames %in% fixedNames]
  # 
  #   muMat = cbind(muMat, thisMu)
  #   Precs = c(Precs, list(thisPrec))
  # }
  # 
  # # print parameters and time taken
  # print(round(cbind(testMat, totalHrs=totalTimes, sdHrs=sdTimes), 3))
  # 
  # # D2
  # out = load("savedOutput/ed/fit_D2.RData")
  # thisMu = summary(SD0)[,1]
  # thisPrec = SD0$jointPrecision
  # 
  # # calculate Hellinger distances
  # distFromBest = numeric(length(allRes))
  # distFromLast = c(length(allRes)-1)
  # bestMu = muMat[,length(allRes)]
  # bestPrec = Precs[[length(allRes)]]
  # bestEig = eigen(bestPrec)
  # for(i in 1:length(allRes)) {
  #   print(paste0("i: ", i, "/", length(allRes)))
  #   
  #   thisMu = muMat[,i]
  #   thisPrec = Precs[[i]]
  #   if(i != length(allRes)) {
  #     thisEig = eigen(thisPrec)
  #   } else {
  #     thisEig = bestEig
  #   }
  #   
  #   distFromBest[i] = HellingerMVN(thisMu, thisPrec, bestMu, bestPrec, eigen1=thisEig, eigen2=bestEig)
  #   
  #   if(i > 1) {
  #     lastMu = muMat[,i-1]
  #     lastPrec = Precs[[i-1]]
  #     distFromLast[i-1] = HellingerMVN(thisMu, thisPrec, lastMu, lastPrec, eigen1=thisEig, eigen2=lastEig)
  #   }
  #   lastEig = thisEig
  # }
  # 
  # # calculate Hellinger distance matrix
  # distMat = matrix(nrow=length(allRes), ncol=length(allRes))
  # for(i in 1:length(allRes)) {
  #   print(paste0("i: ", i, "/", length(allRes)))
  #   
  #   thisMui = muMat[,i]
  #   thisPreci = Precs[[i]]
  #   thisEigi = eigen(thisPreci)
  #   
  #   for(j in 1:length(allRes)) {
  #     thisMuj = muMat[,j]
  #     thisPrecj = Precs[[j]]
  #     thisEigj = eigen(thisPrecj)
  #     
  #     if(i != j) {
  #       distMat[i, j] = HellingerMVN(thisMui, thisPreci, thisMuj, thisPrecj, eigen1=thisEigi, eigen2=thisEigj)
  #     } else {
  #       distMat[i, j] = 0
  #     }
  #     
  #   }
  # }
  # 
  # colnames(distMat) = allRes
  # row.names(distMat) = allRes
  # print(round(distMat, 3))
  # 
  # # plot results
  # 
  # pdf(file=paste0("figures/integration/HellingerDistTest_", min(allRes), "_", max(allRes), ".pdf"), width=5, height=5)
  # distRange = range(c(distFromBest, distFromLast))
  # plot(allRes, distFromBest, ylim=distRange, main="Hellinger distance vs resolution", 
  #      xlab="Number integration points", ylab="Hellinger distance", pch=19, col="black")
  # points(allRes[2:length(allRes)], distFromLast, pch=19, col="blue")
  # abline(h=0, lty=2, col="grey")
  # legend("right", c("Distance from best", "Distance from last"), pch=19, col=c("black", "blue"))
  # dev.off()
  # 
  # pdf(file=paste0("figures/integration/HellingerDistMat_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  # image.plot(list(x=allRes, y=allRes, z=distMat), main="Hellinger distance vs resolution", 
  #            xlab="Number integration points", ylab="Number integration points", asp=1)
  # dev.off()
  
  # resample if necessary ----
  if(!is.null(nSamples)) {
    for(i in 1:length(allRes)) {
      print(paste0("i: ", i, "/", length(allRes)))
      
      resI = allRes[i]
      out = load(paste0("savedOutput/ed/admin1PredsM_DM2_", resI, "_adm2Cov.RData"))
      
      thisN = ncol(admin1Preds$aggregationResults$p)
      
      if(thisN != nSamples) {
        out = load(paste0("savedOutput/ed/fit2_", resI, "_adm2Cov.RData"))
        
        gridPreds = predGrid(SD0, admLevel="adm2", nsim=nSamples)
        if(saveGridPreds) {
          save(gridPreds, file=paste0("savedOutput/ed/gridPreds2_", resI, "_adm2Cov.RData"))
        }
        
        stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
        admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
        admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
        save(stratPreds, file=paste0("savedOutput/ed/stratPredsM_DM2_", resI, "_adm2Cov.RData"))
        save(admin1Preds, file=paste0("savedOutput/ed/admin1PredsM_DM2_", resI, "_adm2Cov.RData"))
        save(admin2Preds, file=paste0("savedOutput/ed/admin2PredsM_DM2_", resI, "_adm2Cov.RData"))
      }
    }
  }
  
  
  # adm1 ----
  # calculate Hellinger distance matrix (predictive admin1)
  distMatAdm1Avg = matrix(nrow=length(allRes), ncol=length(allRes))
  distMatAdm1Max = matrix(nrow=length(allRes), ncol=length(allRes))
  distMatAdm190 = matrix(nrow=length(allRes), ncol=length(allRes))
  for(i in 1:length(allRes)) {
    print(paste0("i: ", i, "/", length(allRes)))
    
    resI = allRes[i]
    out = load(paste0("savedOutput/ed/admin1PredsM_DM2_", resI, "_adm2Cov.RData"))
    
    samplesI = admin1Preds$aggregationResults$p
    
    for(j in 1:i) {
      resJ = allRes[j]
      
      if(i != j) {
        out = load(paste0("savedOutput/ed/admin1PredsM_DM2_", resJ, "_adm2Cov.RData"))
        
        samplesJ = admin1Preds$aggregationResults$p
        
        nAreas = nrow(samplesI)
        dists = numeric(nAreas)
        for(k in 1:nAreas) {
          dists[k] = HellingerUniveriate(samplesI[k,], samplesJ[k,])
        }
        distMatAdm1Avg[i, j] = mean(dists)
        distMatAdm1Max[i, j] = max(dists)
        distMatAdm190[i, j] = quantile(prob=.9, dists)
      } else {
        distMatAdm1Avg[i, j] = 0
        distMatAdm1Max[i, j] = 0
        distMatAdm190[i, j] = 0
      }
    }
  }
  for(i in 1:(length(allRes)-1)) {
    for(j in (i+1):length(allRes)) {
      distMatAdm1Avg[i, j] = distMatAdm1Avg[j,i]
      distMatAdm1Max[i, j] = distMatAdm1Max[j,i]
      distMatAdm190[i, j] = distMatAdm190[j,i]
    }
  }
  
  colnames(distMatAdm1Avg) = allRes
  row.names(distMatAdm1Avg) = allRes
  colnames(distMatAdm1Max) = allRes
  row.names(distMatAdm1Max) = allRes
  colnames(distMatAdm190) = allRes
  row.names(distMatAdm190) = allRes
  print(round(distMatAdm1Avg, 3))
  print(round(distMatAdm1Max, 3))
  print(round(distMatAdm190, 3))
  
  # plot results
  browser()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm1Avg_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm1Avg), main="Hellinger distance vs resolution (Adm1 Avg)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm1Max_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm1Max), main="Hellinger distance vs resolution (Adm1 Max)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm1Q90_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm190), main="Hellinger distance vs resolution (Adm1 Q90)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
  # adm2 ----
  # calculate Hellinger distance matrix (predictive admin2)
  distMatAdm2Avg = matrix(nrow=length(allRes), ncol=length(allRes))
  distMatAdm2Max = matrix(nrow=length(allRes), ncol=length(allRes))
  distMatAdm290 = matrix(nrow=length(allRes), ncol=length(allRes))
  for(i in 1:length(allRes)) {
    print(paste0("i: ", i, "/", length(allRes)))
    
    resI = allRes[i]
    out = load(paste0("savedOutput/ed/admin1PredsM_DM2_", resI, "_adm2Cov.RData"))
    
    samplesI = admin2Preds$aggregationResults$p
    
    for(j in 1:i) {
      resJ = allRes[j]
      
      if(i != j) {
        out = load(paste0("savedOutput/ed/admin2PredsM_DM2_", resJ, "_adm2Cov.RData"))
        
        samplesJ = admin2Preds$aggregationResults$p
        
        nAreas = nrow(samplesI)
        dists = numeric(nAreas)
        for(k in 1:nAreas) {
          
          dists[k] = HellingerUniveriate(samplesI[k,], samplesJ[k,])
        }
        distMatAdm2Avg[i, j] = mean(dists, na.rm=TRUE)
        distMatAdm2Max[i, j] = max(dists, na.rm=TRUE)
        distMatAdm290[i, j] = quantile(prob=.9, dists, na.rm=TRUE)
      } else {
        distMatAdm2Avg[i, j] = 0
        distMatAdm2Max[i, j] = 0
        distMatAdm290[i, j] = 0
      }
    }
  }
  for(i in 1:(length(allRes)-1)) {
    for(j in (i+1):length(allRes)) {
      distMatAdm2Avg[i, j] = distMatAdm2Avg[j,i]
      distMatAdm2Max[i, j] = distMatAdm2Max[j,i]
      distMatAdm290[i, j] = distMatAdm290[j,i]
    }
  }
  
  colnames(distMatAdm2Avg) = allRes
  row.names(distMatAdm2Avg) = allRes
  colnames(distMatAdm2Max) = allRes
  row.names(distMatAdm2Max) = allRes
  colnames(distMatAdm290) = allRes
  row.names(distMatAdm290) = allRes
  print(round(distMatAdm2Avg, 3))
  print(round(distMatAdm2Max, 3))
  print(round(distMatAdm290, 3))
  
  # plot results
  browser()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm2Avg_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm2Avg), main="Hellinger distance vs resolution (Adm2 Avg)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm2Max_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm2Max), main="Hellinger distance vs resolution (Adm2 Max)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
  pdf(file=paste0("figures/integration/HellingerDistMatAdm2Q90_", min(allRes), "_", max(allRes), ".pdf"), width=5.25, height=5)
  image.plot(list(x=allRes, y=allRes, z=distMatAdm290), main="Hellinger distance vs resolution (Adm2 Q90)", 
             xlab="Number integration points", ylab="Number integration points", asp=1)
  dev.off()
  
}

fitResModels = function(allRes=c(100, 125, 150, 175, 200, 225, 300)) {
  optRes = max(allRes)
  
  for(i in 1:(length(allRes)-1)) {
    thisRes = allRes[i]
    
    # fit the model at the current resolution using the optimum of the given 
    # resolution as the starting point
    fitModelAtResolution(thisRes, optRes)
  }
  
  invisible(NULL)
}

# allRes=c(100, 125, 150, 175, 200, 225, 300)
# fitModelAtResolution(125, 300)
# predGridOnly: if TRUE, loads the already fitted model to regenerate gridPreds
fitModelAtResolution = function(res, optRes=NULL, predGridOnly=FALSE) {
  # script for women's secondary education in Nigeria application
  
  # load datasets ----
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/global/edMICS.RData")
  
  # set parameters ----
  # Umut settings: 
  #   5 urban rings of 15 each (61 points total)
  #   10 rural rings of 15 each (136 points total)
  # KMICS=100
  # KDHSurb = 31 # 4 rings of 10 each
  # JInnerUrban = 4
  # KDHSrur = 71 # 4 inner + 4 outer rings of 10 each
  # JInnerRural = 4
  # JOuterRural = 4
  KMICS=res
  KDHSurb = 11 # 3 rings of 5 each
  JInnerUrban = 3
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  JInnerRural = 3
  JOuterRural = 1
  
  # do some precomputation ----
  
  # make integration points if necessary
  # intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=TRUE,
  #                                           numPtsRur=KMICS, numPtsUrb=KMICS, adm2AsCovariate = TRUE, lambda=1)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE, 
  #                                         numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
  #                                         JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
  #                                         JOuterRural=JOuterRural)
  # save(intPtsMICS, file=paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData")
  load("savedOutput/global/intPtsDHS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
  
  # AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  # ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
  ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICS$Stratum, by=list(strat=edMICS$Stratum, urb=edMICS$urban), FUN=length, drop=FALSE)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  # AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  
  # numPerStratUrb = table(edMICS$Stratum[edMICS$urban])
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  startInds = seq(1, KMICS*length(admFinal@data$NAME_FINAL), by=KMICS)
  getInds = function(intPtI = 1, numPerStrat) {
    unlist(mapply(rep, startInds+intPtI-1, each=numPerStrat))
  }
  actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrb))
  XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  # ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  # numPerStratRur = table(edMICS$Stratum[!edMICS$urban])
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRur))
  XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICS$Stratum, admFinal$NAME_FINAL)
  # edMICS = edMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICS[edMICS$urban,]$ys
  nsUrbMICS = edMICS[edMICS$urban,]$ns
  ysRurMICS = edMICS[!edMICS$urban,]$ys
  nsRurMICS = edMICS[!edMICS$urban,]$ns
  
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbMICS = t(AUrbMICS)
  ARurMICS = t(ARurMICS)
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbMICS) = "numeric"
  mode(ARurMICS) = "numeric"
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRur = XRur[,-(2:3)]
  intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  
  # convert A matrices to sparse matrices
  AUrbMICS = as(AUrbMICS, "sparseMatrix")
  ARurMICS = as(ARurMICS, "sparseMatrix")
  AUrbDHS = as(AUrbDHS, "sparseMatrix")
  ARurDHS = as(ARurDHS, "sparseMatrix")
  # AUrbMICS = as.matrix(AUrbMICS)
  # ARurMICS = as.matrix(ARurMICS)
  # AUrbDHS = as.matrix(AUrbDHS)
  # ARurDHS = as.matrix(ARurDHS)
  
  
  # compile model ----
  compile( "code/modBYM2JitterFusionNugget2.cpp", 
           framework="TMBad", safebounds=FALSE)
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/TMB/include" -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppEigen/include"  -DTMB_SAFEBOUNDS -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterFusionNugget2  -DTMB_LIB_INIT=R_init_modBYM2JitterFusionNugget2  -DTMBAD_FRAMEWORK  -I/usr/local/include   -fPIC  -Wall -g -O2  -c code/modBYM2JitterFusionNugget2.cpp -o code/modBYM2JitterFusionNugget2.o
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o code/modBYM2JitterFusionNugget2.so code/modBYM2JitterFusionNugget2.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
  
  # on Idun:
  # g++ -std=gnu++14 -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/include" -DNDEBUG -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/TMB/include" -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/RcppEigen/include"   -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterFusionNugget2  -DTMB_LIB_INIT=R_init_modBYM2JitterFusionNugget2  -DTMBAD_FRAMEWORK  -I/cluster/apps/eb/software/OpenSSL/1.1/include -I/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/include -I/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/include -I/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/include -I/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/include -I/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/include -I/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Java/11.0.2/include -I/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/include -I/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/include -I/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include/flexiblas   -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c code/modBYM2JitterFusionNugget2.cpp -o code/modBYM2JitterFusionNugget2.o
  # g++ -std=gnu++14 -shared -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -L/cluster/apps/eb/software/OpenSSL/1.1/lib64 -L/cluster/apps/eb/software/OpenSSL/1.1/lib -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib64 -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib64 -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Java/11.0.2/lib64 -L/cluster/apps/eb/software/Java/11.0.2/lib -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib64 -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib -L/cluster/apps/eb/software/GCCcore/11.3.0/lib64 -L/cluster/apps/eb/software/GCCcore/11.3.0/lib -o code/modBYM2JitterFusionNugget2.so code/modBYM2JitterFusionNugget2.o -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -lR
  
  # load in TMB function inputs
  # out = load("savedOutput/global/ed2Inputs.RData")
  
  # set priors ----
  alpha_pri = c(0, 100^2)
  beta_pri = c(0, 10^2)
  
  out = load("savedOutput/global/adm2Mat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for bym2 precision
  lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('Epsilon_bym2', 'nuggetUrbMICS', 'nuggetRurMICS', 
                 'nuggetUrbDHS', 'nuggetRurDHS', 'beta', 'alpha')
  
  # collect input data
  
  data_full = list(
    y_iUrbanMICS=ysUrbMICS, # observed binomial experiment at point i (clust)
    y_iRuralMICS=ysRurMICS, # 
    n_iUrbanMICS=nsUrbMICS, # number binomial trials
    n_iRuralMICS=nsRurMICS, # 
    AprojUrbanMICS=AUrbMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralMICS=ARurMICS, # 
    X_betaUrbanMICS=intPtsMICS$XUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralMICS=intPtsMICS$XRur, # 
    wUrbanMICS=intPtsMICS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralMICS=intPtsMICS$wRural, # 
    
    y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
    y_iRuralDHS=ysRurDHS, # 
    n_iUrbanDHS=nsUrbDHS, # number binomial trials
    n_iRuralDHS=nsRurDHS, # 
    AprojUrbanDHS=AUrbDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralDHS=ARurDHS, # 
    X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralDHS=intPtsDHS$covsRur, # 
    wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralDHS=intPtsDHS$wRural, # 
    
    Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
    V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
    alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
    beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
    tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
    gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
    lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
    lambdaTau=lambdaTau, # determines PC prior for tau
    lambdaTauEps=lambdaTauEps, 
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  anyna = function(x) {any(is.na(x))}
  myDim = function(x) {
    dims = dim(x)
    if(is.null(dims)) {
      length(x)
    }
    else {
      dims
    }
  }
  
  if(!predGridOnly) {
    if(!is.null(optRes)) {
      # set initial parameters based on simple model
      out = load(paste0("savedOutput/ed/fit2_", optRes, "_adm2Cov.RData"))
      
      tmb_params <- list(alpha = SD0$par.random[grepl("alpha", names(SD0$par.random))], # intercept
                         beta = SD0$par.random[grepl("beta", names(SD0$par.random))], 
                         log_tau = SD0$par.fixed[which(names(SD0$par.fixed) == "log_tau")], # Log tau (i.e. log spatial precision, Epsilon)
                         logit_phi = SD0$par.fixed[grepl("phi", names(SD0$par.fixed))], # SPDE parameter related to the range
                         log_tauEps = SD0$par.fixed[grepl("tauEps", names(SD0$par.fixed))], # Log tau (i.e. log spatial precision, Epsilon)
                         Epsilon_bym2 = SD0$par.random[grepl("Epsilon", names(SD0$par.random))], # RE on mesh vertices
                         nuggetUrbMICS = SD0$par.random[grepl("nuggetUrbMICS", names(SD0$par.random))], 
                         nuggetRurMICS = SD0$par.random[grepl("nuggetRurMICS", names(SD0$par.random))], 
                         nuggetUrbDHS = SD0$par.random[grepl("nuggetUrbDHS", names(SD0$par.random))], 
                         nuggetRurDHS = SD0$par.random[grepl("nuggetRurDHS", names(SD0$par.random))]
      )
    } else {
      # set initial parameters based on simple/unadjusted model
      print("No initialization supplied for optimization. Starting optimization/model fitting for the unadjusted DHS model")
      initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
      initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
      initAlpha = logit(initRurP)
      initBeta1 = logit(initUrbP) - initAlpha
      
      tmb_paramsStart <- list(alpha = initAlpha, # intercept
                              beta = c(initBeta1, rep(0, ncol(intPtsDHS$covsUrb)-1)), 
                              log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              logit_phi = 0, # SPDE parameter related to the range
                              log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                              nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                              nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
      )
      
      # specify random effects
      rand_effsStart <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS', 'beta', 'alpha')
      
      # collect input data, setting only first weights as nonzero (to 1)
      wUrbanDHStemp=intPtsDHS$wUrban
      wRuralDHStemp=intPtsDHS$wRural
      wUrbanDHStemp[,1] = 1
      wRuralDHStemp[,1] = 1
      wUrbanDHStemp[,-1] = 0
      wRuralDHStemp[,-1] = 0
      
      data_start = list(
        y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
        y_iRuralDHS=ysRurDHS, # 
        n_iUrbanDHS=nsUrbDHS, # number binomial trials
        n_iRuralDHS=nsRurDHS, # 
        AprojUrbanDHS=AUrbDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
        AprojRuralDHS=ARurDHS, # 
        X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
        X_betaRuralDHS=intPtsDHS$covsRur, # 
        wUrbanDHS=wUrbanDHStemp, # nObsUrban x nIntegrationPointsUrban weight matrix
        wRuralDHS=wRuralDHStemp, # 
        
        Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
        V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
        alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
        beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
        tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
        gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
        lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
        lambdaTau=lambdaTau, # determines PC prior for tau
        lambdaTauEps=lambdaTauEps, 
        options=0 # 1 for adreport of log tau and logit phi
      )
      
      dyn.load( dynlib("code/modBYM2JitterDHS2"))
      TMB::config(tmbad.sparse_hessian_compress = 1)
      objStart <- MakeADFun(data=data_start,
                            parameters=tmb_paramsStart,
                            random=rand_effsStart,
                            hessian=TRUE,
                            DLL='modBYM2JitterDHS2')
      
      lower = rep(-10, length(objStart[['par']]))
      upper = rep( 10, length(objStart[['par']]))
      
      # make wrapper functions that print out parameters and function values
      funWrapper = function(par, badParVal=1e10) {
        if(any(par < lower) || any(par > upper)) {
          return(badParVal)
        }
        print(par)
        objVal = testObj[['fn']](par)
        parNames = names(par)
        parVals = par
        parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
        parStr = paste(parStrs, collapse=", ")
        print(paste0("objective: ", objVal, " for parameters, ", parStr))
        objVal
      }
      
      grWrapper = function(par, badParVal=1e10) {
        if(any(par < lower) || any(par > upper)) {
          return(rep(badParVal, length(par)))
        }
        print(par)
        grVal = testObj[['gr']](par)
        parNames = names(par)
        parVals = par
        parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
        parStr = paste(parStrs, collapse=", ")
        grStr = paste(grVal, collapse=", ")
        print(paste0("gradient: ", grStr, " for parameters, ", parStr))
        grVal
      }
      
      testObj = objStart
      optPar = testObj$par
      testObj = objStart
      testObj$env$inner.control = list(maxit=1000, tol10=1e-06)
      testObj$env$tracepar = TRUE
      optStart <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                        method = c("BFGS"), hessian = FALSE, control=list(reltol=1e-06))
      optParStart = optStart$par
      print("Optimization/model fitting for the unadjusted DHS model complete")
      
      # now set the initial parameters
      tmb_params <- list(alpha = testObj$env$last.par[grepl("alpha", names(testObj$env$last.par))], # intercept
                         beta = testObj$env$last.par[grepl("beta", names(testObj$env$last.par))], 
                         log_tau = testObj$env$last.par[names(testObj$env$last.par) == "log_tau"], # Log tau (i.e. log spatial precision, Epsilon)
                         logit_phi = testObj$env$last.par[grepl("logit_phi", names(testObj$env$last.par))], # SPDE parameter related to the range
                         log_tauEps = testObj$env$last.par[grepl("log_tauEps", names(testObj$env$last.par))], # Log tau (i.e. log spatial precision, Epsilon)
                         Epsilon_bym2 = testObj$env$last.par[grepl("Epsilon_bym2", names(testObj$env$last.par))], # RE on mesh vertices
                         nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                         nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS)), 
                         nuggetUrbDHS = testObj$env$last.par[grepl("nuggetUrbDHS", names(testObj$env$last.par))], 
                         nuggetRurDHS = testObj$env$last.par[grepl("nuggetRurDHS", names(testObj$env$last.par))]
      )
      
      dyn.unload( dynlib("code/modBYM2JitterDHS2"))
    }
    
    
    # make TMB fun and grad ----
    # dyn.load( dynlib("code/modBYM2JitterFusionNugget2sparse"))
    dyn.load( dynlib("code/modBYM2JitterFusionNugget2"))
    TMB::config(tmbad.sparse_hessian_compress = 1)
    obj <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     random=rand_effs,
                     hessian=TRUE,
                     DLL='modBYM2JitterFusionNugget2')
    # objFull <- MakeADFun(data=data_full,
    #                      parameters=tmb_params,
    #                      hessian=TRUE,
    #                      DLL='modBYM2JitterFusionNugget2')
    
    lower = rep(-10, length(obj[['par']]))
    upper = rep( 10, length(obj[['par']]))
    
    # make wrapper functions that print out parameters and function values
    funWrapper = function(par, badParVal=1e10) {
      if(any(par < lower) || any(par > upper)) {
        return(badParVal)
      }
      print(par)
      objVal = testObj[['fn']](par)
      parNames = names(par)
      parVals = par
      parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
      parStr = paste(parStrs, collapse=", ")
      print(paste0("objective: ", objVal, " for parameters, ", parStr))
      objVal
    }
    
    grWrapper = function(par, badParVal=1e10) {
      if(any(par < lower) || any(par > upper)) {
        return(rep(badParVal, length(par)))
      }
      print(par)
      grVal = testObj[['gr']](par)
      parNames = names(par)
      parVals = par
      parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
      parStr = paste(parStrs, collapse=", ")
      grStr = paste(grVal, collapse=", ")
      print(paste0("gradient: ", grStr, " for parameters, ", parStr))
      grVal
    }
    
    # * Run TMB ----
    {
      # tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
      tolSeq = 1e-06
      testObj = obj
      optPar = testObj$par
      startTime = proc.time()[3]
      for(thisTol in tolSeq) {
        testObj = obj
        testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
        testObj$env$tracepar = TRUE
        print(paste0("optimizing for tol = ", thisTol, "."))
        opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                      method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
        optPar = opt1$par
        if(!is.null(opt1$message)) {
          print(paste0("error for tol = ", thisTol, ". Message:"))
          print(opt1$message)
          next
        }
        else {
          print(paste0("completed optimization for tol = ", thisTol, ""))
          
          ## Get standard errors
          print("getting standard errors...")
          sdTime = system.time({
            SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                 bias.correct = TRUE,
                                 bias.correct.control = list(sd = TRUE))
            
            if(!SD0$pdHess) {
              # try recalculating for fixed parameters numerically
              warning("Using alternative method for hessian calculation...")
              Hess = numDeriv::hessian( func=testObj$fn, x=optPar )
              SD0 <- sdreport( testObj, hessian.fixed=Hess,
                               getJointPrecision=TRUE,
                               bias.correct = TRUE,
                               bias.correct.control = list(sd = TRUE) )
              
              # system.time(HessTest <- numDeriv::jacobian(func=objFull$gr, x=testObj$env$last.par.best))
              # 15.82288 minutes, not PD
              # system.time(HessTest2 <- numDeriv::hessian(func=objFull$fn, x=testObj$env$last.par.best))
            }
          }
          )[3]
          # SD0
          print(sdTime/60)
          # 3.9447 minutes for intern=FALSE
          
          if(SD0$pdHess) {
            print("Optimization and PD hess calculation done!")
            break
          }
          else {
            print("Hessan not PD. Rerunning optimization with stricter tol...")
            
          }
        }
        
        
        
      }
      endTime = proc.time()[3]
      sdTime/60
      totalTime = endTime - startTime
      print(paste0("optimization took ", totalTime/60, " minutes"))
      # optimization took ~511.460316666667 minutes (for intern=FALSE)
    }
    
    
    # opt0 <- nlminb(start       =    obj[['par']],
    #                objective   =    obj[['fn']],
    #                gradient    =    obj[['gr']],
    #                lower = rep(-10, length(obj[['par']])),
    #                upper = rep( 10, length(obj[['par']])),
    #                control     =    list(trace=1))
    
    # * TMB Posterior Sampling ----
    
    
    ## summary(SD0, 'report')
    ## summary(SD0, 'fixed')
    
    tmbReport = obj$report()
    if(trueObjectSize(obj) < .3*2^30) {
      save(SD0, obj, tmbReport, totalTime, sdTime, file=paste0("savedOutput/ed/fit2_", res, "_adm2Cov.RData"))
    } else {
      save(SD0, tmbReport, totalTime, sdTime, file=paste0("savedOutput/ed/fit2_", res, "_adm2Cov.RData"))
    }
  } else {
    out = load(paste0("savedOutput/ed/fit2_", res, "_adm2Cov.RData"))
  }
  
  gridPreds = predGrid(SD0, admLevel="adm2")
  # \begin{table}[ht]
  # \centering
  # \begin{tabular}{rrrrrr}
  # \hline
  # & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\
  # \hline
  # (Int) & -1.79 & -1.93 & -1.88 & -1.70 & -1.64 \\
  # urb & 0.93 & 0.79 & 0.84 & 1.01 & 1.06 \\
  # access & -0.04 & -0.23 & -0.17 & 0.08 & 0.15 \\
  # elev & 0.13 & -0.11 & -0.03 & 0.30 & 0.41 \\
  # distRiversLakes & 0.04 & -0.36 & -0.20 & 0.32 & 0.44 \\
  # popValsNorm & 0.82 & 0.64 & 0.70 & 0.95 & 1.00 \\
  # sigmaSq & 0.90 & 0.75 & 0.80 & 1.00 & 1.07 \\
  # phi & 0.13 & 0.11 & 0.12 & 0.14 & 0.15 \\
  # sigmaEpsSq & 0.50 & 0.46 & 0.47 & 0.53 & 0.54 \\
  # \hline
  # \end{tabular}
  # \end{table}
  save(gridPreds, file=paste0("savedOutput/ed/gridPreds2_", res, "_adm2Cov.RData"))
  
  if(!predGridOnly) {
    out = load(paste0("savedOutput/ed/gridPreds2_", res, "_adm2Cov.RData"))
    
    stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
    admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
    admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
    save(stratPreds, file=paste0("savedOutput/ed/stratPredsM_DM2_", res, "_adm2Cov.RData"))
    save(admin1Preds, file=paste0("savedOutput/ed/admin1PredsM_DM2_", res, "_adm2Cov.RData"))
    save(admin2Preds, file=paste0("savedOutput/ed/admin2PredsM_DM2_", res, "_adm2Cov.RData"))
    out = load(paste0("savedOutput/ed/stratPredsM_DM2_", res, "_adm2Cov.RData"))
    out = load(paste0("savedOutput/ed/admin1PredsM_DM2_", res, "_adm2Cov.RData"))
    out = load(paste0("savedOutput/ed/admin2PredsM_DM2_", res, "_adm2Cov.RData"))
    
    summaryTabBYM2(SD0, popMat=popMatNGAThresh, 
                   gridPreds=gridPreds)
    # \begin{table}[ht]
    # \centering
    # \begin{tabular}{rrrr}
    # \hline
    # & Est & Q0.025 & Q0.975 \\
    # \hline
    # X.Int. & -1.79 & -1.93 & -1.64 \\
    # beta & 0.93 & 0.79 & 1.06 \\
    # beta.1 & -0.04 & -0.23 & 0.15 \\
    # beta.2 & 0.13 & -0.11 & 0.41 \\
    # beta.3 & 0.04 & -0.36 & 0.44 \\
    # beta.4 & 0.82 & 0.64 & 1.00 \\
    # sigmaSq & 0.90 & 0.75 & 1.07 \\
    # phi & 0.13 & 0.11 & 0.15 \\
    # sigmaEpsSq & 0.50 & 0.46 & 0.54 \\
    # \hline
    # \end{tabular}
    # \end{table}
    plotPreds(SD0, obj, popMat=popMatNGAThresh, 
              gridPreds=gridPreds, arealPreds=NULL, 
              plotNameRoot=paste0("edFusionM_DM2_", res, "_adm2Cov"))
    plotPreds(SD0, obj, popMat=popMatNGAThresh, 
              gridPreds=gridPreds, arealPreds=stratPreds, 
              plotNameRoot=paste0("edFusionM_DM2_", res, "_adm2Cov"), plotNameRootAreal="Strat")
    plotPreds(SD0, obj, popMat=popMatNGAThresh, 
              gridPreds=gridPreds, arealPreds=admin1Preds, 
              plotNameRoot=paste0("edFusionM_DM2_", res, "_adm2Cov"), plotNameRootAreal="Admin1")
    plotPreds(SD0, obj, popMat=popMatNGAThresh, 
              gridPreds=gridPreds, arealPreds=admin2Preds, 
              plotNameRoot=paste0("edFusionM_DM2_", res, "_adm2Cov"), plotNameRootAreal="Admin2")
  }
  
  invisible(NULL)
}


