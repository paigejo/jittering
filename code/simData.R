# simulate data from an SPDE model for simulation study 1

simData1 = function(nsim=100, margVar=.5, effRange=200, sigmaEpsilon=sqrt(1.5), 
                    beta0=-1.25, gamma=1, betaRest=c(0, 0, 0, .5), 
                    mesh=getSPDEMesh(), easpaDat=easpaNGAed, 
                    popMat=popMatNGAThresh, targetPopMat=popMatNGAedThresh, 
                    poppsub=poppsubNGAThresh, nHHMICS=16, nHHDHS=25, seed=123, 
                    useThreshPopMat=TRUE, fixPopPerHH=NULL, 
                    eaSampleStrat="pps", regenPop=FALSE) {
  set.seed(seed)
  
  # make sure everything is ordered nicely
  popMat = popMat[order(popMat$subarea),]
  poppsub = poppsub[order(poppsub$subarea),]
  
  # construct logit offset vector based on covariates in betaRest
  # first get the design matrix
  print("Constructing offset based on covariates...")
  LLcoords = cbind(popMat$lon, popMat$lat)
  tempDesMat = getDesignMat(LLcoords, normalized=TRUE, useThreshPopMat=useThreshPopMat)
  # cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, 
  #       distRiversLakes=distVals, urbanicity=urbanicityVals)
  
  # normalized population density is calculated differently
  load("savedOutput/global/covariatesNorm.RData")
  popVals = extract(pop, LLcoords, method="bilinear")
  
  load("savedOutput/global/popMeanSDCal.RData")
  popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
  popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
  normPop=(log1p(popVals)-popMeanCal)/popSDCal
  normPop[is.na(normPop)] = min(normPop, na.rm=TRUE)
  
  # get final design matrix
  covRestMat = tempDesMat[,-c(1:3, 7)] # remove int, pop, urb, urbanicity
  covRestMat = cbind(covRestMat, normPop=normPop) # add in normalized population density
  
  # calculate offset
  offset = covRestMat %*% betaRest
  
  # get aggregation info from admin2 areas to MICS strata
  tempAreasFrom = popMat$subarea
  tempAreasTo = popMat$stratumMICS
  areasFrom = sort(unique(tempAreasFrom))
  areasToI = match(areasFrom, tempAreasFrom)
  areasTo = tempAreasTo[areasToI]
  
  # simulate populations and surveys
  print("Simulating populations and surveys...")
  surveysDHS = list()
  surveysMICS = list()
  subareaPops = list()
  areaPops = list()
  stratumPops = list()
  for(i in 1:nsim) {
    # simulate population at pixel, EA levels 
    
    simPop = 
      SUMMER::simPopSPDE(nsim=1, easpa=easpaDat, popMat=popMat, targetPopMat=targetPopMat, 
                         poppsub=poppsub, spdeMesh=mesh, 
                         margVar=margVar, sigmaEpsilon=sigmaEpsilon, effRange=effRange, 
                         gamma=gamma, beta0=beta0, seed=NULL, nHHSampled=nHHSampled, 
                         stratifyByUrban=TRUE, subareaLevel=TRUE, offset=offset, 
                         doFineScaleRisk=FALSE, doSmoothRisk=FALSE, min1PerSubarea=TRUE
      )
    
    # calculate stratum level population information
    stratPop = SUMMER::areaPopToArea(areaLevelPop=simPop$subareaPop, 
                                     areasFrom=areasFrom, 
                                     areasTo=areasTo, 
                                     stratifyByUrban=TRUE, doFineScaleRisk=FALSE, doSmoothRisk=FALSE)
    
    # append population information
    if(i == 1) {
      subareaPops = simPop$subareaPop$aggregationResults$pFineScalePrevalence
      areaPops = simPop$areaPop$aggregationResults$pFineScalePrevalence
      stratumPops = stratPop$aggregationResults$pFineScalePrevalence
    } else {
      # cbind the new pop info to the full set of populations
      browser()
      subareaPops = cbind(subareaPops, 
                          simPop$subareaPop$aggregationResults$pFineScalePrevalence)
      areaPops = cbind(areaPops, 
                       simPop$areaPop$aggregationResults$pFineScalePrevalence)
      stratumPops = cbind(stratumPops, 
                          stratPop$aggregationResults$pFineScalePrevalence)
    }
    
    # generate surveys
    # get EA level population information for population i
    thisEApop = simPop$eaPop$eaDatList[i]
    
    # get associated HH level population information
    thisHHpop = getHHpop(thisEApop, fixPopPerHH=fixPopPerHH)
    
    # sample DHS survey for this population
    survDHS = sampleClusterSurveys(1, thisHHpop, HHperClust=25)
    
    # now sample the MICS survey. Do some gymnastics to make sure it works for MICS strata
    tempClustpa = clustpaMICSed
    names(tempClustpa)[1] = "area"
    
    thisHHpop[[1]]$area = adm2ToStratumMICS(thisHHpop[[1]]$subarea)
    
    survMICS = sampleClusterSurveys(1, thisHHpop, HHperClust=16, clustpaList=list(tempClustpa))
    
    # concatenate results
    surveysDHS = c(surveysDHS, survDHS)
    surveysMICS = c(surveysMICS, survMICS)
  }
  
  save(subareaPops, areaPops, stratumPops, surveysDHS, surveysMICS, 
       file="savedOutput/simStudy1/simPopsSurveys.RData")
  
  invisible(NULL)
}

# simulates household level population data given EA level population data
# 
# 
getHHpop = function(popSim, fixPopPerHH=NULL) {
  # first get ea level data from popSim
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      eaPopDat = popSim
    } else {
      stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
    }
  } else {
    stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  if(length(eaPopDat) > 1) {
    warning("length(eaPopDat) > 1, so there may be a lot of household level data, and memory issues accordingly...")
  }
  
  HHdat = list()
  for(i in 1:length(eaPopDat)) {
    eaDat = eaPopDat[[i]]
    
    # first get the number of households
    numHouseholds = eaDat$nHH
    
    # now expand the eaDat table to be in the long format, where each row is a house
    rowsLong = rep(1:nrow(eaDat), numHouseholds)
    thisHHdat = eaDat[rowsLong, ]
    thisHHdat$eaIs = rowsLong
    urbanHH = thisHHdat$urban
    
    # function for randomly spreading people among households in long form data:
    extendThisDat = function(xs, nHH) {
      revenMultinom = function(sizeK) {
        size = sizeK[1]
        k = sizeK[2]
        prob = rep(1/k, k)
        rmultinom(1, size, prob)
      }
      unlist(apply(cbind(xs, nHH), 1, revenMultinom))
    }
    
    # generate how many of the target population are in each cluster
    if(is.null(fixPopPerHH)) {
      lived = extendThisDat(eaDat$N - eaDat$Z, numHouseholds)
      died = extendThisDat(eaDat$Z, numHouseholds)
    } else if(fixPopPerHH == 1) {
      extendDatEven = function(Ns, Zs) {
        
        spreadAmongHHs = function(thisRow) {
          thisN = thisRow[1]
          thisZ = thisRow[2]
          
          # spread population evenly among households
          # nHH = thisRow$nHH
          # hhI = sample(1:nHH, nHH, replace=FALSE)
          c(rep(1, thisZ), rep(0, thisN - thisZ))
        }
        c(unlist(apply(cbind(Ns, Zs), 1, spreadAmongHHs)))
      }
      died = extendDatEven(eaDat$N, eaDat$Z)
      lived = 1 - died
    } else {
      stop("If fixPopPerHH is not NULL it must be 1. Other values not currently supported")
    }
    
    thisHHdat$N = died + lived
    thisHHdat$Z = died
    thisHHdat$nHH = 1
    thisHHdat$pFineScalePrevalence = thisHHdat$Z / thisHHdat$N
    thisHHdat$pFineScalePrevalence[thisHHdat$N == 0] = 0 # NaN otherwise
    
    HHdat = c(HHdat, list(thisHHdat))
  }
  
  HHdat
}

# function for sampling clusters from simulated population
sampleClusterSurveys = function(n=NULL, popSim=NULL, HHperClust=25, fixPopPerHH=NULL, 
                                eaSampleStrat=c("pps", "srs"), clustpaList=list(clustpaDHSed), 
                                seed=NULL) {
  
  eaSampleStrat = match.arg(eaSampleStrat)
  
  # set random seed if supplied
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # initialize EA and HH level population info to NULL
  hhPopDat = NULL
  eaPopDat = NULL
  
  # first get ea level data from popSim if need be
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      if(all(popSim[[1]]$nHH == 1)) {
        hhPopDat = popSim
      } else {
        eaPopDat = popSim
      }
    }
  } else {
    stop("popSim has no EA or HH level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  
  # set n if unset
  if(is.null(n)) {
    if(!is.null(hhPopDat)) {
      n = length(hhDatList)
    } else {
      n = length(eaDatList)
    }
  }
  
  if((length(clustpaList) > 1) && (length(clustpaList) != n)) {
    stop(paste0("length mismatch between n (", n, ") and length(clustpaList) (", length(clustpaList), ")"))
  }
  
  # if popDat has length 1 and n > 1, sample multiple surveys from the same population
  if(!is.null(eaPopDat)) {
    resamplePop = ifelse((length(eaPopDat) == 1) && (n > 1), TRUE, FALSE)
  } else if(!is.null(hhPopDat)) {
    resamplePop = ifelse((length(hhPopDat) == 1) && (n > 1), TRUE, FALSE)
  }
  
  
  print("Simulating surveys")
  surveys = list()
  for(i in 1:n) {
    print(paste0("Sampling survey ", i, "/", n))
    
    # get number of clusters to sample per area
    thisClustpa = clustpaList[[min(c(length(clustpaList), i))]]
    
    # get HH level population information
    if(!is.null(eaPopDat)) {
      # if no HH info available, first get the EA level info, then use to get HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        eaDat = eaPopDat[[1]]
      } else {
        eaDat = eaPopDat[[i]]
      }
      
      hhDat = getHHpop(list(eaDatList = list(eaDat)), fixPopPerHH = fixPopPerHH)[[1]]
      
    } else if(!is.null(hhPopDat)) {
      # we already have HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        hhDat = hhPopDat[[1]]
      } else {
        hhDat = hhPopDat[[i]]
      }
      eaDat = NULL
    }
    
    # obtain info about the EAs
    uniqueEAIs = sort(unique(hhDat$eaIs))
    eaUrbs = hhDat$urban[match(uniqueEAIs, hhDat$eaIs)]
    eaAreas = hhDat$area[match(uniqueEAIs, hhDat$eaIs)]
    
    if(eaSampleStrat == "pps") {
      # calculate number of HHs per EA
      aggHHs = aggregate(hhDat$nHH, by=list(hhDat$eaIs), FUN=sum)
      eaHHs = aggHHs$x
    }
    
    # sample EAs
    sampledEAIs = numeric(0)
    inclusionProbs = numeric(0)
    for(j in 1:nrow(thisClustpa)) {
      
      thisArea = thisClustpa$area[j]
      nUrbEA = thisClustpa$EAUrb[j]
      nRurEA = thisClustpa$EARur[j]
      thisEAIs = uniqueEAIs[eaAreas == thisArea]
      thisEAUrbs = eaUrbs[eaAreas == thisArea]
      thisEAIsUrb = thisEAIs[thisEAUrbs]
      thisEAIsRur = thisEAIs[!thisEAUrbs]
      
      if(eaSampleStrat == "pps") {
        thisEAhhs = eaHHs[eaAreas == thisArea]
        thisEAhhsUrb = thisEAhhs[thisEAUrbs]
        thisEAhhsRur = thisEAhhs[!thisEAUrbs]
      }
      
      # sample EAs and HHs for this area
      sampUrbEAIs = numeric(0)
      sampRurEAIs = numeric(0)
      inclusionProbsUrb = numeric(0)
      inclusionProbsRur = numeric(0)
      if(eaSampleStrat == "srs") {
        if(nUrbEA != 0) {
          sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F)
          inclusionProbsUrb = rep(nUrbEA/length(thisEAIsUrb), nUrbEA)
        }
        if(nRurEA != 0) {
          sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F)
          inclusionProbsRur = rep(nRurEA/length(thisEAIsRur), nRurEA)
        }
      } else if(eaSampleStrat == "pps") {
        require(sampling)
        if(nUrbEA != 0) {
          # sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F, prob=thisEAhhsUrb/sum(thisEAhhsUrb))
          inclusionProbsUrb = nUrbEA * thisEAhhsUrb/sum(thisEAhhsUrb)
          sampUrbEAIs = thisEAIsUrb[as.logical(UPmidzuno(inclusionProbsUrb))]
        }
        if(nRurEA != 0) {
          # sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F, prob=thisEAhhsRur/sum(thisEAhhsRur))
          inclusionProbsRur = nRurEA * thisEAhhsRur/sum(thisEAhhsRur)
          sampRurEAIs = thisEAIsRur[as.logical(UPmidzuno(inclusionProbsRur))]
        }
      } else {
        stop(paste0("eaSampleStrat '", eaSampleStrat, "' not supported"))
      }
      
      # concatenate urban and rural EAs sampled from this area
      thisSampEAIs = c(sampUrbEAIs, sampRurEAIs)
      thisInclusionProbs = c(inclusionProbsUrb, inclusionProbsRur)
      
      # concatenate to vector of all EAs sampled from all areas
      sampledEAIs = c(sampledEAIs, thisSampEAIs)
      inclusionProbs = c(inclusionProbs, thisInclusionProbs)
    }
    
    # subset hhDat to only EAs sampled
    hhSubdat = hhDat[hhDat$eaIs %in% sampledEAIs,]
    
    # sample HHs within chosen EAs
    hhIsTab = aggregate(1:nrow(hhSubdat), by=list(eaIs=hhSubdat$eaIs), function(x) {
      sample(x, HHperClust, replace=FALSE)
    })
    hhIs = sort(c(as.matrix(hhIsTab[,-1])))
    
    hhDatSample = hhSubdat[hhIs,]
    
    # aggregate HH level data for only the EAs sampled
    aggTab = lapply(1:ncol(hhSubdat), function(j) {
      varName = names(hhSubdat)[j]
      aggregate(hhSubdat[,j], by=list(hhSubdat$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhSubdat)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    eaDatSampled = aggTab
    
    # do the same for the actual sampled HHs
    aggTab = lapply(1:ncol(hhDatSample), function(j) {
      varName = names(hhDatSample)[j]
      aggregate(hhDatSample[,j], by=list(hhDatSample$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhDatSample)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    surveyDat = aggTab
    
    # calculate final sampling weight and number HHs in the full EA
    surveyDat$nHHsEA = eaDatSampled$nHH[match(surveyDat$eaIs, eaDatSampled$eaIs)]
    surveyDat$includeProbHH = surveyDat$nHH / surveyDat$nHHsEA
    surveyDat$samplingWeight = surveyDat$N / (surveyDat$includeProbHH * surveyDat$includeProbEA)
    
    
    # concatenate results
    surveys = c(surveys, list(surveyDat))
  }
  
  surveys
}

# stratOrder: if provided, will sort resulting strata to be in the same order as in 
#             stratOrder Otherwise, sorted alphabetically
getClustpaFromSurvey = function(survDat=ed, stratOrder=easpaNGA$area, stratName="area", HHperEA=25) {
  
  if("ns" %in% names(survDat)) {
    survDat$n = survDat$ns
  }
  
  nAreas = length(unique(survDat[[stratName]]))
  
  # first the number of sampled clusters
  nEAtabTmp = aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=length, drop=FALSE)
  nEAtabTmp[is.na(nEAtabTmp[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  nEAtab = data.frame(nEAtabTmp[1:nAreas, 1], EAUrb=nEAtabTmp[(nAreas+1):(2*nAreas), 3], EARur=nEAtabTmp[1:nAreas, 3])
  names(nEAtab)[1] = stratName
  nEAtab$EATotal = nEAtab$EAUrb + nEAtab$EARur
  
  # initialize clustpa
  clustpa = nEAtab
  
  # second the number of households
  clustpa$HHUrb = clustpa$EAUrb * HHperEA
  clustpa$HHRur = clustpa$EARur * HHperEA
  clustpa$HHTotal = clustpa$EATotal * HHperEA
  
  # third the number of people
  popTabTmp = aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=sum, drop=FALSE)
  popTabTmp[is.na(popTabTmp[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  popTab = data.frame(popTabTmp[1:nAreas, 1], popUrb=popTabTmp[(nAreas+1):(2*nAreas), 3], popRur=popTabTmp[1:nAreas, 3])
  names(popTab)[1] = stratName
  popTab$popTotal = popTab$popUrb + popTab$popRur
  
  # concatenate cluster level denominator info to clustpa info
  clustpa$popUrb = popTab$popUrb
  clustpa$popRur = popTab$popRur
  clustpa$popTotal = popTab$popTotal
  
  # sort if need be
  if(!is.null(stratOrder)) {
    ordering = order(stratOrder)
    clustpa = clustpa[ordering,]
  }
  
  clustpa
}

makeIntegrationPointsSim1 = function() {
  
  KMICS=100
  KDHSurb = 11 # 3 rings of 5 each
  JInnerUrban = 3
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  JInnerRural = 3
  JOuterRural = 1
  
  out = load("savedOutput/simStudy1/simSurveys.RData")
  
  for(i in 1:length(surveysDHS)) {
    thisEdDHS = surveysDHS[[i]]
    
    intPtsDHS = makeAllIntegrationPointsDHS(cbind(thisEdDHS$east, thisEdDHS$north), thisEdDHS$urban, 
                                            areaNames=thisEdDHS$subarea, popPrior=TRUE, 
                                            numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                            JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                            JOuterRural=JOuterRural, adminMap=adm2Full)
    
    save(intPtsDHS, file=paste0("intPtsDHS_simStudy1_", i, ".RData"))
  }
  
  invisible(NULL)
}

makeInputsSim1 = function() {
  
  out = load("savedOutput/simStudy1/simSurveys.RData")
  
  for(i in 1:length(surveysDHS)) {
    thisEdDHS = surveysDHS[[i]]
    thisEdMICS = surveysMICS[[i]]
    
    out = load(paste0("intPtsDHS_simStudy1_", i, ".RData"))
    
    
    
  }
  
  
}






