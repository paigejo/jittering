# script for collecting inputs for TMB for the M_DM model

makeInputsMDM = function(datDHS, datMCIS, intPtsMICS=NULL, intPtsDHS=NULL, 
                         KMICS=100,
                         KDHSurb = 11, # 3 rings of 5 each
                         JInnerUrban = 3,
                         KDHSrur = 16, # 3 inner + 1 outer rings of 5 each
                         JInnerRural = 3,
                         JOuterRural = 1, admMICS=admFinal, adm2DHS=adm2Full) {
  
  # make sure edMICS is sorted by Stratum column according to admMICS
  edMICS = sortByCol(edMICS, "Stratum", admMICS$NAME_FINAL)
  
  # make integration points if necessary
  if(is.null(intPtsMICS)) {
    intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
                                              numPtsRur=KMICS, numPtsUrb=KMICS, saveOutput=FALSE)
  }
  
  if(is.null(intPtsDHS)) {
    # make integration points if necessary
    intPtsDHS = makeAllIntegrationPointsDHS(cbind(datDHS$east, datDHS$north), datDHS$urban, 
                                            areaNames=datDHS$subarea, popPrior=TRUE, 
                                            numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                            JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                            JOuterRural=JOuterRural, adminMap=adm2DHS, saveOutput=FALSE)
  }
  
  intPtsMICS = straightenMICS(intPtsMICS)
  
  # AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admMICS$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  # ARurDHS = makeApointToArea(intPtsDHS$areasRural, admMICS$NAME_FINAL) # 41 x 810
  
  AUrbDHS = makeApointToArea(rep(adm2ToStratumMICS(datDHS$subarea[datDHS$urban]), times=KDHSurb), admMICS$NAME_FINAL) # 775 x 6259 nArea x nObsUrb
  ARurDHS = makeApointToArea(rep(adm2ToStratumMICS(datDHS$subarea[!datDHS$urban]), times=KDHSrur), admMICS$NAME_FINAL) # 775 x 12960
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(datMICS$Stratum, by=list(strat=datMICS$Stratum, urb=datMICS$urban), FUN=length, drop=FALSE)
  allNumPerStrat = straightenNumPerStrat(allNumPerStrat, admMICS$NAME_FINAL)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  # AUrbMICS = makeApointToArea(datMICS$Stratum[datMICS$urban], admMICS$NAME_FINAL)
  # TODO: EXTEND AMICS TO BE LARGER, INCLUDE DIFFERENT ROW FOR EACH INTEGRATION POINT AND OBSERVATION
  
  # numPerStratUrb = table(datMICS$Stratum[datMICS$urban])
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # startInds = seq(1, KMICS*length(admMICS@data$NAME_FINAL), by=KMICS)
  # getInds = function(intPtI = 1, numPerStrat) {
  #   unlist(mapply(rep, startInds+intPtI-1, each=numPerStrat))
  # }
  # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrb))
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  startStratInds = which(XUrb$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XUrb)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrb[x])})) # length nUrb, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
  nUrb = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIUrb = allAreaIs + (allIntIs-1)*nAreas
  XUrb = XUrb[transformIUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  AUrbMICS = makeApointToArea(adm2ToStratumMICS(XUrb$subarea), admMICS$NAME_FINAL)
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  # ARurMICS = makeApointToArea(datMICS$Stratum[!datMICS$urban], admMICS$NAME_FINAL)
  # numPerStratRur = table(datMICS$Stratum[!datMICS$urban])
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRur))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  startStratInds = which(XRur$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XRur)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRur[x])})) # length nRur, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
  nRur = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIRur = allAreaIs + (allIntIs-1)*nAreas
  XRur = XRur[transformIRur,]
  
  ARurMICS = makeApointToArea(adm2ToStratumMICS(XRur$subarea), admMICS$NAME_FINAL)
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(datMICS$Stratum, admMICS$NAME_FINAL)
  # datMICS = datMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = datMICS[datMICS$urban,]$ys
  nsUrbMICS = datMICS[datMICS$urban,]$ns
  ysRurMICS = datMICS[!datMICS$urban,]$ys
  nsRurMICS = datMICS[!datMICS$urban,]$ns
  
  ysUrbDHS = datDHS$y[datDHS$urban]
  ysRurDHS = datDHS$y[!datDHS$urban]
  nsUrbDHS = datDHS$n[datDHS$urban]
  nsRurDHS = datDHS$n[!datDHS$urban]
  
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
  # AUrbMICS = as(AUrbMICS, "sparseMatrix")
  # ARurMICS = as(ARurMICS, "sparseMatrix")
  # AUrbDHS = as(AUrbDHS, "sparseMatrix")
  # ARurDHS = as(ARurDHS, "sparseMatrix")
  AUrbMICS = as.matrix(AUrbMICS)
  ARurMICS = as.matrix(ARurMICS)
  AUrbDHS = as.matrix(AUrbDHS)
  ARurDHS = as.matrix(ARurDHS)
  
  areaidxlocUrbanMICS = apply(AUrbMICS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
  areaidxlocRuralMICS = apply(ARurMICS, 1, function(x) {match(1, x)}) - 1
  areaidxlocUrbanMICS = as.integer(areaidxlocUrbanMICS)
  areaidxlocRuralMICS = as.integer(areaidxlocRuralMICS)
  areaidxlocUrbanDHS = apply(AUrbDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
  areaidxlocRuralDHS = apply(ARurDHS, 1, function(x) {match(1, x)}) - 1
  areaidxlocUrbanDHS = as.integer(areaidxlocUrbanDHS)
  areaidxlocRuralDHS = as.integer(areaidxlocRuralDHS)
  
  list(AUrbMICS, ARurMICS, AUrbDHS, ARurDHS, intPtsDHS, intPtsMICS, 
       areaidxlocUrbanMICS, areaidxlocRuralMICS, 
       areaidxlocUrbanDHS, areaidxlocRuralDHS, 
       ysUrbMICS, nsUrbMICS, ysRurMICS, nsRurMICS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS)
}

makeInputsMdm = function(inputsMDM) {
  
  # set weights to be put entirely on a random integration point.
  inputsMDM$intPtsDHS$wUrban[,-1] = 0
  inputsMDM$intPtsDHS$wUrban[,1] = 1
  inputsMDM$intPtsDHS$wRural[,-1] = 0
  inputsMDM$intPtsDHS$wRural[,1] = 1
  
  inputsMDM$intPtsMICS$wUrban[,-1] = 0
  inputsMDM$intPtsMICS$wUrban[,1] = 1
  inputsMDM$intPtsMICS$wRural[,-1] = 0
  inputsMDM$intPtsMICS$wRural[,1] = 1
  
  inputsMDM
}






