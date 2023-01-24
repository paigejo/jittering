# functions for constructing jittering integration points for both DHS and MICS 
# data

##### Jittering for DHS ----

# constructs the points over which to integrate the likelihood 
# with respect to the jittering distribution, and their weights, 
# where the weights are related to the jittering distribution
# Arguments: 
# urban: whether or not to generate integration points from the 
#        urban or rural jittering distributions.
# numPoints: the number of integration points. 1 goes in the first 
#            ring, and the rest are evenly split among the other 
#            rings
# scalingFactor: adjust typical DHS jitter distance by a scaling factor
# JInner: the number of `inner' rings within 5km, including the first central ring
# JOuter: the number of `outer' rings beyond 5km. In urban case, only 
#         M=JInner+JOuter is used
# integrationPointType: either the integration point is set as the 
#                       center of mass within the integration area 
#                       or the midpoint of the radius and angle
# verbose: whether to run the function in verbose mode
getIntegrationPointsDHS = function(urban=TRUE, numPoints=ifelse(urban, 11, 16), 
                                   scalingFactor=1, 
                                   JInner=3, JOuter=ifelse(urban, 0, 1), 
                                   integrationPointType=c("mean", "midpoint"), 
                                   verbose=TRUE) {
  integrationPointType = match.arg(integrationPointType)
  
  M = JInner + JOuter
  ms = c(1, rep(floor((numPoints - 1) / (M - 1)), M - 1))
  msInner = ms[1:JInner]
  msOuter = ms[(JInner + 1):M]
  if(sum(ms) != numPoints) {
    stop("(numPoints - 1)/M must be a whole number")
  }
  
  if(urban) {
    # 2 is largest possible displacement distance
    rsInner = 2 * cumsum(ms) / sum(ms)
    rsOuter = NULL
    # r2 = (m2 + 1) * r3 / sum(ms)
    # r1 = r3 / (1 + m2 + m3)
  } else {
    rsInner = 5 * cumsum(msInner) / sum(msInner)
    rsOuter = 5 + 5 * cumsum(msOuter) / sum(msOuter)
    # r3 = 10 # largest possible displacement distance for 1% of data
    # r2 = 5 # largest possible displacement distance for 99% of data
    # r1 = r2 / (1 + m2)
  }
  # r1 = r2 / sqrt(m2 + 1) # for uniform distribution
  rs = c(rsInner, rsOuter)
  
  # areas for each ring segment
  As = pi * rs^2 # get areas of M discs of radius r
  As = As - c(0, As[1:(M-1)]) # subtract area of the smaller desk to get area of annulus
  As = As / ms # divide annulus area by number of integration areas per annulus
  # A1 = pi * r1^2
  # A2 = (pi * r2^2 - A1) / m2
  # A3 = (pi * r3^2 - pi * r2^2) / m3
  
  # these probabilities and weights were for a uniform distribution: 
  # p1 = p2 = 99/100 * 1/(pi*r2^2) + 1/100 * 1/(pi*r3^2)
  # p3 = 1/100 * 1/(pi*r3^2)
  # 
  # w1 = p1 * A1
  # w2 = p2 * A2
  # w3 = p3 * A3
  
  # helper function for calculting values of pi(s_i | s_ijk^*) for 
  # integration points in ring 1, 2, or 3
  densityFunTemp = function(d, Di=rs[M]) {
    out = 1 / (Di * 2 * pi * d)
    out[d >= Di] = 0
    out[d < 0] = NaN
    out
  }
  densityFunFinal = function(d) {
    if(urban) {
      densityFunTemp(d, rs[M])
    } else {
      densityFunTemp(d, rs[JInner]) * 99 / 100 + densityFunTemp(d, rs[M]) * 1 / 100
    }
  }
  
  # p1 = Inf
  # p2 = densityFunFinal(mean(c(r1, r2)))
  # p3 = densityFunFinal(mean(c(r2, r3)))
  rsIntegrationPointsMidpoint = c(0, rs[-M] + diff(rs) / 2) # midpoint solution
  if(integrationPointType == "mean") {
    aDiff = c(0, 2 * pi / ms[-1])
    shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
    shrinkFactor[1] = 1
    rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
    tooShrunk = rsIntegrationPoints < c(0, rs[-M])
    if(any(tooShrunk)) {
      warning(paste0("Center of mass is outside integration area for rings ", 
                     paste(tooShrunk, collapse=", "), 
                     ". Setting integration point to closest within the integration area"))
      rsIntegrationPoints[toShrunk] = c(0, rs[-M])[tooShrunk]
    }
  } else if(integrationPointType == "midpoint") {
    rsIntegrationPoints = rsIntegrationPointsMidpoint
  }
  ps = densityFunFinal(rsIntegrationPoints)
  
  # these weights are for a distribution uniform in radial distance and in angle, marginally:
  if(urban) {
    # ws[1] = rs[1] / rs[M]
    # w2 = (r2 - r1) / r3 * 1 / m2
    # w3 = (r3 - r2) / r3 * 1 / m3
    ws = (rs - c(0, rs[-M])) / rs[M] * 1 / ms
    
    if(verbose) {
      print("This should be 1 (the sum of the weights), and should all be equal:")
    }
  } else {
    # w1 = r1 / r2 * (99 / 100) + r1 / r3 / 100
    # w2 = (r2 - r1) / (r2 * m2) * (99 / 100) + (r2 - r1) / (r3 * m2) / 100
    # w3 = (r3 - r2) / (r3 * m3) / 100
    annulusWidths = diff(c(0, rs))
    annulusWidthsInner = annulusWidths[1:JInner]
    annulusWidthsOuter = annulusWidths[-(1:JInner)]
    wsInner = annulusWidthsInner * ( 
      (99 / 100) / (rsInner[JInner] * msInner) + 
        (1 / 100) / (rs[M] * msInner) )
    wsOuter = annulusWidthsOuter * (1 / 100) / (rs[M] * msOuter)
    ws = c(wsInner, wsOuter)
    
    if(verbose) {
      print(paste0("This should be 1 (the sum of the weights), and the first ", JInner, " (and last ", JOuter, ") should be equal:"))
    }
  }
  if(verbose) {
    print(sum(ms*ws))
    print(matrix(ws, nrow=1))
  }
  
  # pts1 = matrix(c(0, 0), nrow=1)
  # as2 = seq(0, 2*pi, l=m2+1)[-(m2+1)]
  # pts2 = (r1 + r2) / 2 * cbind(cos(as2), sin(as2))
  # as3 = seq(0, 2*pi, l=m3+1)[-(m3+1)] + pi/m3
  # pts3 = (r2 + r3) / 2 * cbind(cos(as3), sin(as3))
  pts = list()
  as = list()
  for(j in 1:M) {
    if(j == 1) {
      as = list(as = 0)
      pts = list(pts = matrix(c(0, 0), nrow=1))
    } else {
      thisas = seq(0, 2*pi, l=ms[j]+1)[-(ms[j]+1)]
      if(j > 1 && j %% 2 == 1) {
        thisas = thisas + pi / ms[j]
      }
      as = c(as, as=list(thisas))
      thisPts = rsIntegrationPoints[j] * cbind(cos(thisas), sin(thisas))
      pts = c(pts, pts=list(thisPts))
    }
  }
  names(as) = paste0("as", 1:M)
  names(pts) = paste0("pts", 1:M)
  
  if(scalingFactor != 1) {
    rs = rs * scalingFactor
    As = As * scalingFactor^2
    ps = ps / scalingFactor^2
    pts = lapply(pts, function(x){x * scalingFactor})
    densityFunFinalScaled = function(d) {
      (1 / scalingFactor^2) * densityFunFinal(d / scalingFactor)
    }
  } else {
    densityFunFinalScaled = densityFunFinal
  }
  
  # list(r1=r1, r2=r2, r3=r3, m1=m1, m2=m2, m3=m3, 
  #      A1=A1, A2=A2, A3=A3, 
  #      w1=w1, w2=w2, w3=w3, 
  #      p1=p1, p2=p2, p3=p3, 
  #      as2=as2, as3=as3, 
  #      pts1=pts1, pts2=pts2, pts3=pts3, 
  #      densityFun=densityFunFinal)
  list(rs=rs, ms=ms, As=As, ws=ws, ps=ps, 
       pts=pts, as=as, ptRs=rsIntegrationPoints, 
       densityFun=densityFunFinalScaled)
}

# construct integration points as well as weights. 
# Output: 3 pairs of matrices of dimension nCoordsInStratum x nIntegrationPointsInStratum, 
#         each pair contains one urban matrix and one equivalent rural matrix
#   x: x/easting coordinates
#   y: y/northing coordinates
#   w: integration weights
# Input: 
#   coords: 2 column matrix of observation easting/northing coordinates
#   urbanVals: vector of observation urbanicity classifications
#   numPointsUrban: number of urban numerical integration points
#   numPointsRural: number of rural numerical integration points
#   scalingFactor: factor by which to scale the jittering distribution. 
#                  1 corresponds to standard DHS jittering, larger than 1 
#                  corresponds to more jittering than DHS
#   JInnerUrban: number of inner integration rings for urban points
#   JOuterUrban: number of outer integration rings for urban points
#   JInnerRural: number of inner integration rings for rural points
#   JOuterRural: number of outer integration rings for rural points
#   integrationPointType: 'mean' is center of mass, 'midpoint' is the 
#                         median angle and median radius within the 
#                         integration area
makeAllIntegrationPointsDHS = function(coords, urbanVals, 
                                       numPointsUrban=11, numPointsRural=16, 
                                       scalingFactor=1, 
                                       JInnerUrban=3, JOuterUrban=0, 
                                       JInnerRural=3, JOuterRural=1, 
                                       integrationPointType=c("mean", "midpoint")) {
  
  # calculate integration points and weights relative to individual points
  outUrban = getIntegrationPoints(urban=TRUE, numPointsUrban, 
                                  scalingFactor, 
                                  JInnerUrban, JOuterUrban, 
                                  integrationPointType, 
                                  verbose=FALSE)
  
  outRural = getIntegrationPoints(urban=FALSE, numPointsRural, 
                                  scalingFactor, 
                                  JInnerRural, JOuterRural, 
                                  integrationPointType, 
                                  verbose=FALSE)
  
  # concatenate integration points and weights into a single vector
  xsUrbanVec = unlist(sapply(outUrban$pts, function(x) {x[,1]}))
  ysUrbanVec = unlist(sapply(outUrban$pts, function(x) {x[,2]}))
  wsUrbanVec = rep(outUrban$ws, outUrban$ms)
  xsRuralVec = unlist(sapply(outRural$pts, function(x) {x[,1]}))
  ysRuralVec = unlist(sapply(outRural$pts, function(x) {x[,2]}))
  wsRuralVec = rep(outRural$ws, outRural$ms)
  
  # separate coordinates in urban and rural
  coordsUrban = coords[urbanVals,]
  coordsRural = coords[!urbanVals,]
  nUrban = nrow(coordsUrban)
  nRural = nrow(coordsRural)
  
  # calculate the six matrices
  xUrban = outer(coordsUrban[,1], xsUrbanVec, "+")
  yUrban = outer(coordsUrban[,2], ysUrbanVec, "+")
  wUrban = matrix(wsUrbanVec, ncol=length(wsUrbanVec), nrow=nUrban, byrow=TRUE)
  xRural = outer(coordsRural[,1], xsRuralVec, "+")
  yRural = outer(coordsRural[,2], ysRuralVec, "+")
  wRural = matrix(wsRuralVec, ncol=length(wsRuralVec), nrow=nRural, byrow=TRUE)
  
  # return list of all matrices
  list(xUrban=xUrban, yUrban=yUrban, wUrban=wUrban, 
       xRural=xRural, yRural=yRural, wRural=wRural)
}

# Inputs:
#  integrationPointInfo: outfrom from makeAllIntegrationPoints
#  ys: observation vector
#  urbanicity: vector of TRUE/FALSE depending on urbanicity of the observation
#  ns: observation denominator vector
#  spdeMesh: spde triangular basis function mesh object
# Outputs: 
#  dat: data.frame containing y, n, east, north, w, urban for each integration point
#  AUrban: (nObs x nUrbanIntegrationPts) x nMesh sparse spatial projection matrix for 
#          urban observations. Every nObsUrban rows is the same observation but a new 
#          integration point
#  ARural: (nObs x nRuralIntegrationPts) x nMesh sparse spatial projection matrix for 
#          rural observations. Every nObsRural rows is the same observation but a new 
#          integration point
makeJitterDataForTMB = function(integrationPointInfo, ys, urbanicity, ns, spdeMesh) {
  # first extract the integration point information
  xUrban = integrationPointInfo$xUrban
  yUrban = integrationPointInfo$yUrban
  wUrban = integrationPointInfo$wUrban
  xRural = integrationPointInfo$xRural
  yRural = integrationPointInfo$yRural
  wRural = integrationPointInfo$wRural
  
  # get the long set of coordinates
  coordsUrban = cbind(c(xUrban), c(yUrban))
  coordsRural = cbind(c(xRural), c(yRural))
  
  # separate observations by urbanicity
  ysUrban = ys[urbanicity]
  ysRural = ys[!urbanicity]
  nsUrban = ns[urbanicity]
  nsRural = ns[!urbanicity]
  
  # # gather data into data frame
  # dat = data.frame(y = c(rep(ysUrban, ncol(xUrban)), rep(ysRural, ncol(xRural))), 
  #                  n = c(rep(nsUrban, ncol(xUrban)), rep(nsRural, ncol(xRural))), 
  #                  east = c(c(xUrban), c(xRural)), 
  #                  north = c(c(yUrban), c(yRural)), 
  #                  w = c(c(wUrban), c(wRural)),
  #                  urban = c(rep(TRUE, length(xUrban)), rep(FALSE, length(xRural))))
  
  # construct `A' matrices
  AUrban = inla.spde.make.A(mesh = spdeMesh,
                            loc = coordsUrban)
  ARural = inla.spde.make.A(mesh = spdeMesh,
                            loc = coordsRural)
  
  list(ysUrban=ysUrban, ysRural=ysRural, 
       nsUrban=nsUrban, nsRural=nsRural, 
       AUrban=AUrban, ARural=ARural)
}

##### Jittering for MICS ----

# constructs the points over which to integrate the likelihood 
# with respect to the jittering distribution, and their weights, 
# where the weights are related to the jittering distribution
# Arguments: 
# strat: MICS stratum name
# kmresFineStart: km resolution for the starting fine grid of 
#   integration points
# numPts: number of integration points for the area (produces fewer if there are 
#   fewer than that many fine scale grid points)
# propUrb: Proportion of area's population that is urban 
# proj: projection to use for the points
# projArea: projection to use for the areas
# spatialAsCovariate: whether to use easting/northing as a covariate
# lambda: spatial scaling coefficient. Default is 1 / priorSD for 
#   priorSD = (domainDiameter / 5) / 2
# domainDiameter: used for calculating default lambda. Default is 1463.733 km, 
#   Nigeria's diameter
getIntegrationPointsMICS = function(strat, kmresFineStart=2.5, numPtsUrb=25, numPtsRur=25, 
                                    stratumMICSMapDat=admFinalFull, stratumMICSNameVar="NAME_FINAL", 
                                    subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                    poppsub=poppsubNGAThresh, 
                                    normalized=TRUE, useThreshPopMat=TRUE, 
                                    proj=projNigeria, projArea=projNigeriaArea, 
                                    spatialAsCovariate=FALSE, 
                                    lambda=NULL, domainDiameter=NULL, 
                                    returnFineGrid=FALSE) {
  
  
  fineIntPtsTab = getFineIntPointsInfoMICS(stratumName=strat, kmresStart=kmresFineStart, 
                                           minPointsUrb=numPtsUrb, minPointsRur=numPtsRur, 
                                           stratumMICSMapDat=stratumMICSMapDat, stratumMICSNameVar=stratumMICSNameVar, 
                                           subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                           poppsub=poppsub, 
                                           normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                           proj=proj, projArea=projArea)
  
  X = fineIntPtsTab[,c(11:14, 16)] # don't use urbanicity to generate clustering
  pop = fineIntPtsTab$pop
  urb = fineIntPtsTab$urban
  ENCoords = cbind(fineIntPtsTab$east, fineIntPtsTab$north)
  
  # add spatial coordinates as covariates if user requests
  if(spatialAsCovariate) {
    if(is.null(domainDiameter)) {
      # adm0Boundary = gBoundary(adm0)
      # adm0BoundaryCoords = do.call("rbind", lapply(adm0Boundary@lines[[1]]@Lines, function(x) {x@coords}))
      # domainDiameter = max(rdist.earth(adm0BoundaryCoords, miles=FALSE))
      domainDiameter = 1463.733 # in km
    }
    
    # normalize spatial coordinates based on prior median effective range
    if(is.null(lambda)) {
      priorSD = (domainDiameter / 5) / 2
      lambda = 1 / priorSD
    }
    
    ENCoordsNorm = sweep(ENCoords, 2, colMeans(ENCoords), FUN="-")
    ENCoordsNorm = ENCoordsNorm * lambda
    X = cbind(X, ENCoordsNorm)
  }
  
  nPtsUrb = sum(urb)
  nPtsRur = sum(!urb)
  
  # if((nPtsUrb != numPtsUrb) || (nPtsRur != numPtsRur)) {
  #   browser()
  # }
  
  hasUrbPop = TRUE
  hasRurPop = TRUE
  if(nPtsUrb > 0) {
    # do the weighted K-medoids (technically PAM): 
    # Maechler, M., P. Rousseeuw, A. Struyf, M. Hubert and K. Hornik (2011).
    # cluster: Cluster Analysis Basics and Extensions. R package version 1.14.1
    XmatUrb = as.matrix(X)[urb,]
    distMatUrb = rdist(XmatUrb, XmatUrb)
    medoidsUrb = wcKMedoids(distMatUrb^2, numPtsUrb, weights=pop[urb], method="PAM")
    medoidIUrb = medoidsUrb$clustering
    
    # calculate weights of the medoids
    totalPopUrb = aggregate(pop[urb], by=list(medoid=medoidIUrb), FUN=sum)
    weightsUrb = totalPopUrb[,2] / sum(totalPopUrb[,2])
  } else {
    hasUrbPop = FALSE
    XmatUrb = NULL
    distMatUrb = NULL
    medoidsUrb = NULL
    medoidIUrb = NULL
    
    totalPopUrb = NULL
    weightsUrb = NULL
  }
  
  if(nPtsRur > 0) {
    XmatRur = as.matrix(X)[!urb,]
    distMatRur = rdist(XmatRur, XmatRur)
    medoidsRur = wcKMedoids(distMatRur^2, numPtsRur, weights=pop[!urb], method="PAM")
    medoidIRur = medoidsRur$clustering
    
    totalPopRur = aggregate(pop[!urb], by=list(medoid=medoidIRur), FUN=sum)
    weightsRur = totalPopRur[,2] / sum(totalPopRur[,2])
  } else {
    hasRurPop = FALSE
    XmatRur = NULL
    distMatRur = NULL
    medoidsRur = NULL
    medoidIRur = NULL
    
    totalPopRur = NULL
    weightsRur = NULL
  }
  
  # # assign fine grid points to each medoid
  # distToMedoidsUrb = distMatUrb[,medoidIUrb]
  # medoidAssignedUrb = apply(distToMedoidsUrb, 1, which.min)
  # pointIAssignedUrb = medoidIUrb[medoidAssignedUrb]
  # 
  # distToMedoidsRur = distMatRur[,medoidIRur]
  # medoidAssignedRur = apply(distToMedoidsRur, 1, which.min)
  # pointIAssignedRur = medoidIRur[medoidAssignedRur]
  
  # return results
  sortIUrb = which(fineIntPtsTab$urban)[sort(unique(medoidIUrb))]
  sortIRur = which(!fineIntPtsTab$urban)[sort(unique(medoidIRur))]
  out = list(strat=strat, hasUrbPop=hasUrbPop, hasRurPop=hasRurPop, 
       pts=fineIntPtsTab[c(sortIUrb, sortIRur),], weights=c(weightsUrb, weightsRur), 
       ptsUrb=fineIntPtsTab[sortIUrb,], weightsUrb=weightsUrb, 
       ptsRur=fineIntPtsTab[sortIRur,], weightsRur=weightsRur)
  
  if(returnFineGrid) {
    list(fineGrid = fineIntPtsTab, 
         intPoints = out, 
         medoidIRur=medoidIRur, 
         medoidIUrb=medoidIUrb)
  } else {
    out
  }
}

# Contructs all integration points for MICS data, over all areas. Calls 
# getIntegrationPointsMICS for each stratum
# area: shapefile of the area to transform
# kmresFine: km resolution for the fine grid used to produce the final set of 
#   integration points
# numPts: number of integration points for the area (produces fewer if there are 
#   fewer than that many fine scale grid points)
# propUrb: Proportion of area's population that is urban 
# proj: projection to use for the points
# projArea: projection to use for the areas
# spatialAsCovariate: whether to use easting/northing as a covariate
# lambda: spatial scaling coefficient. Default is 1 / priorSD for 
#   priorSD = (domainDiameter / 5) / 2
# domainDiameter: used for calculating default lambda. Default is 1463.733 km, 
#   Nigeria's diameter
makeAllIntegrationPointsMICS = function(datStrata=NULL, datUrb=NULL, kmresFineStart=2.5, 
                                        numPtsUrb=25, numPtsRur=25, 
                                        stratumMICSMapDat=admFinalFull, 
                                        stratumMICSNameVar="NAME_FINAL", 
                                        subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                        poppsub=poppsubNGAThresh, 
                                        normalized=TRUE, useThreshPopMat=TRUE, 
                                        proj=projNigeria, projArea=projNigeriaArea, 
                                        spatialAsCovariate=FALSE, 
                                        lambda=NULL, domainDiameter=NULL, 
                                        fileNameRoot="MICSintPts_", loadSavedIntPoints=TRUE) {
  
  if(is.null(domainDiameter)) {
    domainDiameter = 1463.733 # in km
  }
  
  # normalize spatial coordinates based on prior median effective range
  if(is.null(lambda)) {
    priorSD = (domainDiameter / 5) / 2
    lambda = 1 / priorSD
  }
  
  # set file name root for saving the results
  fileNameRoot = paste0(fileNameRoot, "_km", kmresFineStart, 
                        "_nPtU", numPtsUrb, "_nPtR", numPtsRur, 
                        "_norm", as.numeric(normalized), 
                        "_thresh", as.numeric(useThreshPopMat), 
                        "_spatCov", as.numeric(spatialAsCovariate), 
                        "_lam", round(lambda, 4))
  
  # For each stratum, generate the integration points
  allIntPts = list()
  strataMICS = stratumMICSMapDat@data[[stratumMICSNameVar]]
  for(i in 1:length(strataMICS)) {
    thisStrat = strataMICS[i]
    if(thisStrat == "Lake Chad") {
      # no people in Lake Chad
      next
    }
    
    print(paste0("Constructing integration points for MICS stratum ", thisStrat, 
                 " (", i, "/", length(strataMICS), ")"))
    
    # output (thisIntPoints) is a list with pts, weights, and area
    if(loadSavedIntPoints && file.exists(paste0(fileNameRoot, "_i", i, ".RData"))) {
      out = load(paste0(fileNameRoot, "_i", i, ".RData"))
    } else {
      thisIntPoints= getIntegrationPointsMICS(strat=thisStrat, kmresFineStart=kmresFineStart, 
                                              numPtsUrb=numPtsUrb, numPtsRur=numPtsRur, 
                                              stratumMICSMapDat=stratumMICSMapDat, 
                                              stratumMICSNameVar=stratumMICSNameVar, 
                                              subareaMapDat=subareaMapDat, 
                                              subareaNameVar=subareaNameVar, 
                                              poppsub=poppsub, 
                                              normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                              proj=proj, projArea=projArea, 
                                              spatialAsCovariate=spatialAsCovariate, 
                                              lambda=lambda, domainDiameter=domainDiameter)
      
      save(thisIntPoints, file=paste0(fileNameRoot, "_i", i, ".RData"))
    }
    
    allIntPts = c(allIntPts, list(thisIntPoints))
  }
  
  # concatenate integration points and weights into matrices
  getCovIUrb = function(i) {
    c(do.call("rbind", lapply(allIntPts, function(x) {x$ptsUrb[,i]})))
  }
  getCovIRur = function(i) {
    c(do.call("rbind", lapply(allIntPts, function(x) {x$ptsRur[,i]})))
  }
  allCovsUrb = data.frame(lapply(1:ncol(allIntPts[[1]]$ptsUrb), getCovIUrb))
  names(allCovsUrb) = names(allIntPts[[1]]$ptsUrb)
  allCovsRur = data.frame(lapply(1:ncol(allIntPts[[1]]$ptsRur), getCovIRur))
  names(allCovsRur) = names(allIntPts[[1]]$ptsRur)
  
  wsUrban = do.call("rbind", (lapply(allIntPts, function(x) {x$weightsUrb})))
  wsRural = do.call("rbind", (lapply(allIntPts, function(x) {x$weightsRur})))
  
  browser()
  if(is.null(datStrata)) {
    if(!is.null(datUrb)) {
      stop("datUrb provided but not datStrata")
    }
    
    # return list of all matrices at the stratum rather than observation level
    list(strataMICS=strataMICS, 
         XUrb=allCovsUrb, XRur=allCovsRur, 
         wUrban=wsUrban, wRural=wsRural)
  } else {
    if(is.null(datUrb)) {
      stop("non-stratified integration point construction not currently supported")
    }
    
    # we must expand the matrices to be at the observation rather than stratum level
    strataUrb = datStrata[datUrb]
    strataRur = datStrata[!datUrb]
    strataUrbU = sort(unique(strataUrb))
    strataRurU = sort(unique(strataRur))
    
    browser()
    allCovsUrb = cbind(stratI=rep(1:length(strataUrbU), numPtsUrb), 
                       intI=rep(1:numPtsUrb, each=length(strataUrbU)), 
                       allCovsUrb)
    allCovsUrbByObs = as.data.frame(obsI=rep(1:length(strataUrb), numPtsUrb), 
                                    stratI=rep(strataUrb, numPtsUrb), 
                                    intI=rep(1:numPtsUrb, each=length(strataUrb)))
    XUrbFull = merge(allCovsUrbByObs, allCovsUrb)
    
    allCovsRur = cbind(stratI=rep(1:length(strataRurU), numPtsRur), 
                       intI=rep(1:numPtsRur, each=length(strataRurU)), 
                       allCovsRur)
    allCovsRurByObs = as.data.frame(obsI=rep(1:length(strataRur), numPtsRur), 
                                    stratI=rep(strataRur, numPtsRur), 
                                    intI=rep(1:numPtsRur, each=length(strataRur)))
    XRurFull = merge(allCovsRurByObs, allCovsRur)
    
    # return list of all matrices at the observation level
    list(strataMICS=strataMICS, 
         XUrb=XUrbFull, XRur=XRurFull, 
         wUrban=wsUrban, wRural=wsRural)
  }
}

# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
getFineIntPointsInfoMICS = function(stratumName, kmresStart=2.5, minPointsUrb=20, minPointsRur=20, 
                                    stratumMICSMapDat=admFinalFull, stratumMICSNameVar="NAME_FINAL", 
                                    subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                    poppsub=poppsubNGAThresh, 
                                    normalized=TRUE, useThreshPopMat=TRUE, 
                                    proj=projNigeria, projArea=projNigeriaArea) {
  
  # get saved urbanicity and population integration matrix
  if(useThreshPopMat) {
    out = load("savedOutput/global/popMatNGAThresh.RData")
    popMat = popMatNGAThresh
  } else {
    out = load("savedOutput/global/popMatNGA.RData")
    popMat = popMatNGA
  }
  stratumPopMat = popMat[popMat$strat==stratumName,]
  
  # obtain the SpatialPolygon associated with this MICS stratum
  stratPolygonI = which(stratumMICSMapDat@data[[stratumMICSNameVar]] == stratumName)
  stratumMICSPolygons = SpatialPolygons(stratumMICSMapDat@polygons)
  stratumMICSSpatialPolygon = stratumMICSPolygons[stratPolygonI]
  stratumMICSPolygons@proj4string = stratumMICSMapDat@proj4string
  stratumMICSSpatialPolygon@proj4string = stratumMICSMapDat@proj4string
  
  # get all subareas associated with this MICS stratum
  allStrataMICS = getMICSstratumNigeria(poppsub$subarea, poppsub$area)
  thisI = allStrataMICS == stratumName
  thisPoppsub = poppsub[thisI,]
  uniqueSubareas = thisPoppsub$subarea
  
  totalUrbPop = sum(thisPoppsub$popUrb)
  totalRurPop = sum(thisPoppsub$popRur)
  
  # subset subareaMapDat to subareas in relevant MICS stratum only
  subareaPolygons = SpatialPolygons(subareaMapDat@polygons)
  subareaPolygonI = which(subareaMapDat@data[[subareaNameVar]] %in% uniqueSubareas)
  thisSubareaPolygons = subareaPolygons[subareaPolygonI]
  thisData = subareaMapDat@data[subareaPolygonI,]
  thisSubareaMapDat = SpatialPolygonsDataFrame(thisSubareaPolygons, thisData)
  
  # subset grid so it's in the domain
  # inDomain = in.poly(lonLatGrid, domainPoly)
  # determine version of PROJ
  ver = terra::gdal(lib="proj")
  PROJ6 <- as.numeric(substr(ver, 1, 1)) >= 6
  
  # determine whether subareas have custom grid points
  customGridPoints = logical(length(uniqueSubareas))
  for(i in 1:length(uniqueSubareas)) {
    thisSubarea = uniqueSubareas[i]
    thisPopMat = stratumPopMat[stratumPopMat$subarea == thisSubarea,]
    
    # get number of urban and rural points in the subarea
    thisNUrb = sum(thisPopMat$urban)
    thisNRur = sum(!thisPopMat$urban)
    
    customGridPoints[i] = FALSE
    if((thisNUrb == 1) && (thisNRur == 1)) {
      urbPt = thisPopMat[thisPopMat$urban,]
      rurPt = thisPopMat[!thisPopMat$urban,]
      
      if((urbPt$east == rurPt$east) && (urbPt$north == rurPt$north)) {
        customGridPoints[i] = TRUE
      }
    }
  }
  
  # from lon/lat coords to easting/northing
  bbox = projNigeriaBBox(stratumMICSSpatialPolygon@bbox)
  eastRange = sort(bbox[1,])
  northRange = sort(bbox[2,])
  kmres = kmresStart * 2
  npUrb = 0
  npRur = 0
  while((npUrb < minPointsUrb) && (npRur < minPointsRur)) {
    kmres = kmres/2
    print(paste0("Constructing fine grid with resolution ", kmres))
    eastPoints = seq(eastRange[1], eastRange[2], by=kmres)
    northPoints = seq(northRange[1], northRange[2], by=kmres)
    allPointsEN = make.surface.grid(list(east=eastPoints, north=northPoints))
    allPointsLL = proj(allPointsEN, inverse=TRUE)
    
    # keep only points in the area of interest
    if(!PROJ6) {
      allPointsLLsp = sp::SpatialPoints(allPointsLL, sp::CRS("+proj=longlat"))
    } else {
      allPointsLLsp = sp::SpatialPoints(allPointsLL, sp::CRS(SRS_string="EPSG:4326"))
    }
    inArea = sp::over(allPointsLLsp, stratumMICSSpatialPolygon)
    inArea = !is.na(inArea)
    allPointsEN = allPointsEN[inArea,]
    allPointsLL = allPointsLL[inArea,]
    
    # get subareas associated with the points
    theseSubareas = getRegion2(allPointsLL, mapDat=thisSubareaMapDat, nameVar=subareaNameVar)
    theseSubareas = theseSubareas$regionNames
    naSubs = is.na(theseSubareas)
    
    allPointsEN = allPointsEN[!naSubs,]
    allPointsLL = allPointsLL[!naSubs,]
    theseSubareas = theseSubareas[!naSubs]
    
    uniqueSubareas = sort(unique(theseSubareas))
    
    # popMatIs = sapply(1:nrow(allPointsEN), function(i) {
    #   pts = allPointsEN[i,]
    #   thisSubarea = theseSubareas[i]
    #   thisPopMatI = which(stratumPopMat$subarea == thisSubarea)
    #   thisPopMat = stratumPopMat[thisPopMatI,]
    #   thisSubareaI = match(thisSubarea, uniqueSubareas)
    #   
    #   dists = rdist(matrix(pts, ncol=2), cbind(thisPopMat$east, thisPopMat$north))
    #   closeI = which.min(dists)
    #   thisPopMatI[which(closeI)]
    # })
    getPopMatIs = function(i) {
      pts = allPointsEN[i,]
      thisSubarea = theseSubareas[i]
      thisPopMatI = which(stratumPopMat$subarea == thisSubarea)
      thisPopMat = stratumPopMat[thisPopMatI,]
      thisSubareaI = match(thisSubarea, uniqueSubareas)
      
      # dists = rdist(matrix(pts, ncol=2, cbind(thisPopMat$east, thisPopMat$north))
      
      # fullPopMat is 5km resolution, so must be within 2.5 km in easting and northing directions
      closeE = (pts[1] > thisPopMat$east - 2.5) & (pts[1] <= thisPopMat$east + 2.5)
      closeN = (pts[2] > thisPopMat$north - 2.5) & (pts[2] <= thisPopMat$north + 2.5)
      closeI = closeE & closeN
      
      # there should be exactly 1 point we're closest to within this subarea
      if(sum(closeI > 1)) {
        stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      } else if(sum(closeI) == 0) {
        # this case should only happen at the edges, but just take closest point then
        dists = rdist(rbind(pts), cbind(thisPopMat$east, thisPopMat$north))
        # warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
        return(thisPopMatI[which.min(dists)])
      } else {
        return(thisPopMatI[which(closeI)])
      }
    }
    
    popMatIs = unlist(sapply(1:nrow(allPointsEN), getPopMatIs))
    urbVals = stratumPopMat$urban[popMatIs]
    
    npUrb = sum(urbVals)
    npRur = sum(!urbVals)
    
    print(paste0("Fine grid has ", npUrb, " and ", npRur,  " urban and rural points respectively"))
    
    if(totalUrbPop == 0) {
      npUrb = minPointsUrb
    }
    if(totalRurPop == 0) {
      npRur = minPointsRur
    }
  }
  
  adm2Vals = stratumPopMat$subarea[popMatIs]
  
  fineGridCoordsEN = allPointsEN
  fineGridCoordsLL = allPointsLL
  if(!PROJ6) {
    fineGridCoordsLLsp = sp::SpatialPoints(fineGridCoordsLL, sp::CRS("+proj=longlat"))
  } else {
    fineGridCoordsLLsp = sp::SpatialPoints(fineGridCoordsLL, sp::CRS(SRS_string="EPSG:4326"))
  }
  
  popVals = terra::extract(pop, fineGridCoordsLLsp, method="simple")
  popVals[is.na(popVals)] = 0
  
  # recalibrate populations in urban/rural parts of the area
  urbText = sapply(urbVals, function(x) {ifelse(x, "U", "R")})
  adm2TimesUR = paste(adm2Vals, urbText, sep=",")
  regions = sort(unique(adm2Vals))
  thisPoppsub = poppsub[poppsub$subarea %in% regions,]
  thisPoppsub = thisPoppsub[order(thisPoppsub$subarea),]
  regionTotals = c(thisPoppsub$popRur, thisPoppsub$popUrb)
  regionsUR = c(paste(regions, "R", sep=","), paste(regions, "U", sep=","))
  zeroPops = regionTotals == 0
  regionTotals = regionTotals[!zeroPops]
  regionsUR = regionsUR[!zeroPops]
  
  finalPopVals = SUMMER::calibrateByRegion(pointTotals=popVals, pointRegions=adm2TimesUR, 
                                           regions=regionsUR, regionTotals=regionTotals)
  finalPopMat = data.frame(east=fineGridCoordsEN[,1], north=fineGridCoordsEN[,2], 
                           lon=fineGridCoordsLL[,1], lat=fineGridCoordsLL[,2], 
                           pop=finalPopVals, urban=urbVals, area=stratumPopMat$area[popMatIs], 
                           subarea=adm2Vals, strat=stratumPopMat$strat[popMatIs], 
                           popMatIs=which(popMat$strat==stratumName)[popMatIs])
  
  if(normalized) {
    out = load("savedOutput/global/covariatesNorm.RData")
    # out = load("savedOutput/global/covariates.RData")
    inf = sessionInfo()
    if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
      urb@file@name = "~/git/jittering/savedOutput/global/urb.tif"
      accessNorm@file@name = "~/git/jittering/savedOutput/global/accessNorm.tif"
      elevNorm@file@name = "~/git/jittering/savedOutput/global/elevNorm.tif"
      minDistRiverLakesNorm@file@name = "~/git/jittering/savedOutput/global/minDistRiverLakesNorm.tif"
    }
    
    urbanicityVals = terra::extract(urb, fineGridCoordsLL, method="bilinear") # don't normalize urbanicity
    accessVals = terra::extract(accessNorm, fineGridCoordsLL, method="bilinear")
    elevVals = terra::extract(elevNorm, fineGridCoordsLL, method="bilinear")
    distVals = terra::extract(minDistRiverLakesNorm, fineGridCoordsLL, method="bilinear")
  } else {
    out = load("savedOutput/global/covariates.RData")
    
    inf = sessionInfo()
    if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
      urb@file@name = "~/git/jittering/savedOutput/global/urb.tif"
      access@file@name = "~/git/jittering/savedOutput/global/access.tif"
      elev@file@name = "~/git/jittering/savedOutput/global/elev.tif"
      dist@file@name = "~/git/jittering/savedOutput/global/dist.tif"
    }
    
    urbanicityVals = extract(urb, fineGridCoordsLL, method="bilinear") # don't normalize urbanicity
    accessVals = extract(access, fineGridCoordsLL, method="bilinear")
    elevVals = extract(elev, fineGridCoordsLL, method="bilinear")
    distVals = extract(minDistRiverLakes, fineGridCoordsLL, method="bilinear")
  }
  
  # remove NA covariate points and points with zero population, and, if there 
  # aren't enough points, restart this function with finer resolution
  naCovRowIs = is.na(urbanicityVals) | is.na(accessVals) | is.na(elevVals) | 
    is.na(distVals) | (finalPopVals == 0)
  
  if(((sum(finalPopMat$urb[!naCovRowIs]) < minPointsUrb) && totalUrbPop > 0) || ((sum(!finalPopMat$urb[!naCovRowIs]) < minPointsRur) && totalRurPop > 0)) {
    warning("NA covariates and zero pop points reduced number of urban and rural points to below minimum. Increasing resolution...")
    
    out = getFineIntPointsInfoMICS(stratumName=stratumName, kmresStart=kmres/2, minPointsUrb=minPointsUrb, minPointsRur=minPointsRur, 
                                   stratumMICSMapDat=stratumMICSMapDat, stratumMICSNameVar=stratumMICSNameVar, 
                                   subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                   poppsub=poppsub, 
                                   normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                   proj=proj, projArea=projArea)
  } else {
    out = cbind(finalPopMat, int=1, access=accessVals, elev=elevVals, distRiversLakes=distVals, urbanicity=urbanicityVals)[!naCovRowIs,]
    
    # recalibrate final population densities and include log population density as a covariate
    trueFinalUrbVals = finalPopMat$urban[!naCovRowIs]
    urbText = sapply(trueFinalUrbVals, function(x) {ifelse(x, "U", "R")})
    adm2TimesUR = paste(finalPopMat$subarea[!naCovRowIs], urbText, sep=",")
    trueFinalPopVals = SUMMER::calibrateByRegion(pointTotals=finalPopMat$pop[!naCovRowIs], pointRegions=adm2TimesUR, 
                                                 regions=regionsUR, regionTotals=regionTotals)
    
    out = cbind(out, log1pPop=log1p(trueFinalPopVals))
  }
  
  out
}

