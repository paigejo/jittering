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

# constructs the sub-integration points for each integration point. Used 
# to adjust the weights associated with each integration point.
# Arguments: 
# integrationPoints: output of getIntegrationPoints
# adminPoly: polygon of Admin area
# nSubAPerPoint: number of unique angles of sub-integration 
#                points within any given integration point area
# nSubRPerPoint: number of unique radii of sub-integration 
#                points within any given integration point area
getSubIntegrationPointsDHS = function(integrationPoints, centerCoords=cbind(0, 0), 
                                   nSubAPerPoint=10, nSubRPerPoint=10) {
  rs = integrationPoints$rs
  ms = integrationPoints$ms
  As = integrationPoints$As
  ws = integrationPoints$ws
  pts = integrationPoints$pts
  as = integrationPoints$as
  ptRs = integrationPoints$ptRs
  densityFun = integrationPoints$densityFun
  
  centerCoords = pts[[1]]
  
  # generates sub-integration points for any point
  getSubIntegrationPointsForOnePoint = function(minR, maxR, minA, maxA, theseCenterCoords) {
    widthR = (maxR - minR)/(nSubRPerPoint)
    widthA = (maxA - minA)/(nSubAPerPoint)
    
    # calculate radial coordinates of the sub-integration points
    
    # get angular coordinates of the centers of mass
    if(minA <= maxA) {
      theseAs = seq(minA + widthA/2, maxA - widthA/2, by=widthA)
    } else {
      theseAs = seq(minA + widthA/2 - 2*pi, maxA - widthA/2, by=widthA)
      theseAs[theseAs < 0] = theseAs[theseAs < 0] + 2*pi
    }
    
    # now get radial centers of mass
    theseRs = seq(minR + widthR/2, maxR - widthR/2, by=widthR)
    rsIntegrationPointsMidpoint = theseRs # midpoint solution
    
    aDiff = widthA
    shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
    shrinkFactor[1] = 1
    rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
    tooShrunk = rsIntegrationPoints < (rsIntegrationPointsMidpoint - widthR/2)
    if(any(tooShrunk)) {
      warning(paste0("Center of mass is outside integration area for rings ", 
                     paste(tooShrunk, collapse=", "), 
                     ". Setting integration point to closest within the integration area"))
      rsIntegrationPoints[toShrunk] = rsIntegrationPointsMidpoint[toShrunk] - widthR/2
    }
    theseRs = rsIntegrationPoints
    
    thesePointsRadial = make.surface.grid(list(rs=theseRs, as=theseAs))
    
    # convert to Euclidean coordinates
    thesePointsEuclidean = cbind(thesePointsRadial[,1]*cos(thesePointsRadial[,2]), 
                                 thesePointsRadial[,1]*sin(thesePointsRadial[,2]))
    
    # translate coordinates based on the jittered observation coordinates
    sweep(thesePointsEuclidean, 2, c(centerCoords), "+")
  }
  
  # for every ring:
  #   for every point:
  #     get sub-integration points
  subWs = list()
  subPts = list()
  for(i in 1:length(rs)) {
    thisIntW = ws[i]
    theseSubWs = rep(thisIntW/(nSubRPerPoint*nSubAPerPoint), 
                     each=nSubRPerPoint*nSubAPerPoint*nrow(pts[[i]]))
    
    theseas = as[[i]]
    theseMinR = ifelse(i==1, 0, rs[i-1])
    theseMaxR = rs[i]
    if(length(theseas) != 1) {
      aWidth = theseas[2] - theseas[1]
    } else {
      aWidth = 2*pi
    }
    theseSubPts = c()
    
    for(j in 1:length(theseas)) {
      # determine boundaries of this integration area
      thisMinA = theseas[j] - aWidth/2
      thisMaxA = theseas[j] + aWidth/2
      thisMinR = theseMinR
      thisMaxR = theseMaxR
      
      # obtain sub-integration points
      thisSubPts = getSubIntegrationPointsForOnePoint(minR=thisMinR, maxR=thisMaxR, 
                                                      minA=thisMinA, maxA=thisMaxA)
      theseSubPts = rbind(theseSubPts, thisSubPts)
    }
    
    subPts = c(subPts, list(theseSubPts))
    subWs = c(subWs, list(theseSubWs))
  }
  
  list(subPts=subPts, subWs=subWs)
}

updateWeightsByAdminArea = function(coords, 
                                    urbanVals, adminMap, 
                                    integrationPointsUrban, 
                                    integrationPointsRural, areas=NULL, 
                                    nSubAPerPoint=10, nSubRPerPoint=10, 
                                    testMode=FALSE, areaNameVar="NAME_FINAL", proj=projNigeria) {
  
  adminMapPoly = as.SpatialPolygons.PolygonsList(adminMap@polygons, adminMap@proj4string)
  
  # calculate set of typical sub-integration points for urban and rural clusters
  subIntegrationPointsUrban = getSubIntegrationPointsDHS(integrationPoints=integrationPointsUrban, 
                                                      nSubAPerPoint=nSubAPerPoint, 
                                                      nSubRPerPoint=nSubRPerPoint)
  
  subIntegrationPointsRural = getSubIntegrationPointsDHS(integrationPoints=integrationPointsRural, 
                                                      nSubAPerPoint=nSubAPerPoint, 
                                                      nSubRPerPoint=nSubRPerPoint)
  
  # get admin areas associated with coordinates
  coordsLonLat = proj(coords, inverse=TRUE)
  spCoordsLonLat = SpatialPoints(coordsLonLat, proj4string=adminMap@proj4string, bbox = NULL)
  out = over(spCoordsLonLat, adminMap, returnList = FALSE)
  adminNames = out[[areaNameVar]]
  # adminIDs = out$OBJECTID
  
  # get the ID of the admin each point is associated with, making sure to take the closest if a point isn't in any
  temp = over(spCoordsLonLat, adminMap, returnList=FALSE)
  nas = is.na(temp[[areaNameVar]])
  
  # calculate distances to admin boundaries for unknown points
  adminMapPolygons = as.SpatialPolygons.PolygonsList(adminMap@polygons, adminMap@proj4string)
  require(geosphere)
  if(is.null(areas)) {
    naClosestIDs = sapply(which(nas), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons)[4]})
    adminID = rep(1, nrow(temp))
    adminID[nas] = naClosestIDs
    adminID[!nas] = match(temp[[areaNameVar]][!nas], adminMap[[areaNameVar]])
  } else {
    adminID = match(areas, adminMap[[areaNameVar]])
  }
  
  
  # for each jittered coordinate:
  #   for each integration point:
  #     get associated sub-integration points
  #     get proportion in correct admin area
  #     update integration point weight
  wsUrban = matrix(nrow=sum(urbanVals), ncol=sum(sapply(integrationPointsUrban$pts, function(x) {nrow(x)})))
  wsRural = matrix(nrow=sum(!urbanVals), ncol=sum(sapply(integrationPointsRural$pts, function(x) {nrow(x)})))
  iUrban = 1
  iRural = 1
  for(i in 1:nrow(coords)) {
    # time1 = proc.time()[3]
    theseCoords = matrix(coords[i,], nrow=1)
    thisArea = adminNames[i]
    thisAreaID = adminID[i]
    thisPoly = adminMapPoly[thisAreaID]
    
    # get sub-integration points
    isUrban = urbanVals[i]
    if(isUrban) {
      thisSubOut = subIntegrationPointsUrban
    } else {
      thisSubOut = subIntegrationPointsRural
    }
    thisSubWs = thisSubOut$subWs
    thisSubPtsEN = lapply(thisSubOut$subPts, function(x) {sweep(x, 2, theseCoords, "+")})
    thisSubPtsLL = lapply(thisSubPtsEN, proj, inverse=TRUE)
    
    # project subPts to correct projection
    thisSubPtsSPLonLat = lapply(thisSubPtsLL, function(x) {SpatialPoints(x, proj4string=adminMap@proj4string)})
    
    # determine if each sub-integration point is in correct admin area
    goodAreas <- lapply(thisSubPtsSPLonLat, function(x) {!is.na(over(x, thisPoly, returnList=FALSE))})
    
    
    # update weights for sub-integration points
    updatedSubWs = thisSubWs
    updatedSubWs = lapply(1:length(updatedSubWs), function(x) {
      temp = updatedSubWs[[x]]
      thisNotGoodAreas = !goodAreas[[x]]
      if(!all(thisNotGoodAreas)) {
        temp[thisNotGoodAreas] = 0
      }
      else {
        warning(paste0("point ", i, " (", coords[i,1], ", ", coords[i,2], ") outside of assigned area with ",
                       "no integration points in assigned area. No integration weights set to zero."))
      }
      temp
    })
    
    # sum sub-integration weights to get new (unnormalized) integration weights
    nSubPts = nSubRPerPoint * nSubAPerPoint
    tempWs = lapply(updatedSubWs, function(x) {
      nIntPts = length(x) / nSubPts
      aggIDs = rep(1:nIntPts, each=nSubPts)
      aggregate(x, by=list(aggIDs), FUN=sum)$x
    })
    tempWs = unlist(tempWs)
    
    # normalize new integration weights to sum to 1
    finalWs = tempWs/sum(tempWs)
    
    # update weights matrix
    if(isUrban) {
      wsUrban[iUrban,] = finalWs
      iUrban = iUrban + 1
    } else {
      wsRural[iRural,] = finalWs
      iRural = iRural + 1
    }
    
    # time2 = proc.time()[3]
    # print(paste0("Iteration ", i, "/", nrow(coords), " took ", round(time2-time1, 2), " seconds"))
  }
  
  if(!testMode) {
    list(wUrban=wsUrban, wRural=wsRural)
  } else {
    list(wUrban=wsUrban, wRural=wsRural, 
         subPts=thisSubPts, goodPts=goodAreas, updatedSubWs=updatedSubWs)
  }
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
#   popPrior: use population density as prior for updating weights if need be.
#   setMissingToAvg: sets NA covariates to 0 if they are pop or urban or normalized
makeAllIntegrationPointsDHSold = function(coords, urbanVals, 
                                    numPointsUrban=11, numPointsRural=16, 
                                    scalingFactor=1, 
                                    JInnerUrban=3, JOuterUrban=0, 
                                    JInnerRural=3, JOuterRural=1, 
                                    integrationPointType=c("mean", "midpoint"), 
                                    adminMap=admFinalFull, areaNameVar="NAME_FINAL", nSubAPerPoint=10, nSubRPerPoint=10, 
                                    popPrior=TRUE, testMode=FALSE, proj=projNigeria, 
                                    outFile="savedOutput/global/intPtsDHS.RData", 
                                    getCovariates=TRUE, normalized=TRUE, useThreshPopMat=TRUE, 
                                    extractMethod="bilinear", setMissingToAvg=TRUE) {
  
  # calculate integration points and weights relative to individual points
  outUrban = getIntegrationPointsDHS(urban=TRUE, numPointsUrban, 
                                  scalingFactor, 
                                  JInnerUrban, JOuterUrban, 
                                  integrationPointType, 
                                  verbose=FALSE)
  
  outRural = getIntegrationPointsDHS(urban=FALSE, numPointsRural, 
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
  
  if(!is.null(adminMap)) {
    # if adminMap is input, integration weights will be adjusted using a much 
    # finer "sub" integration grid
    
    # first subset jittered points that are close enough to the border to have 
    # a possibility of being adjusted
    
    # determine if points are within max distance of admin area:
    # first calculate max distance from admin area
    maxUrbanDistance = 2 * scalingFactor
    maxRuralDistance = 10 * scalingFactor
    maxDist = rep(maxRuralDistance, nrow(coords))
    maxDist[urbanVals] = maxUrbanDistance
    
    # determine what admin area each point is in
    coordsLonLat = proj(coords, inverse=TRUE)
    spCoordsLonLat = SpatialPoints(coordsLonLat, adminMap@proj4string)
    temp = over(spCoordsLonLat, adminMap, returnList=FALSE)
    nas = is.na(temp[[areaNameVar]])
    
    # calculate distances to admin boundaries for unknown points
    adminMapPolygons = as.SpatialPolygons.PolygonsList(adminMap@polygons, adminMap@proj4string)
    require(geosphere)
    naClosestIDs = sapply(which(nas), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons)[4]})
    adminID = rep(1, nrow(temp))
    adminID[nas] = naClosestIDs
    adminID[!nas] = match(temp[[areaNameVar]][!nas], adminMap[[areaNameVar]])
    areas = adminMap[[areaNameVar]][adminID]
    
    # calculate distances to admin boundaries
    dists = sapply(1:nrow(coords), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons[adminID[ind]])[1]}) * (1/1000)
    
    # set whether or not to update weights based on distance to admin boundaries
    updateI = dists < maxDist
    urbanUpdateI = updateI[urbanVals]
    ruralUpdateI = updateI[!urbanVals]
    
    # calculate updated weights for the integration points near the borders
    tempCoords = coords[updateI,]
    tempUrbanVals = urbanVals[updateI]
    require(fields)
    
    if(testMode) {
      # in this case, we take only one set of coords that is very close to border, 
      # and adjust its weight, saving relevant results for plotting
      tempDists = dists[updateI]
      closeCoords = which.min(tempDists)
      smallTempCoords = matrix(tempCoords[closeCoords,], nrow=1)
      smallTempUrbanVals = tempUrbanVals[closeCoords]
      
      tempNewWsForPlotting = updateWeightsByAdminArea(coords=smallTempCoords, urbanVals=smallTempUrbanVals, 
                                                      adminMap=adminMap, 
                                                      integrationPointsUrban=outUrban, 
                                                      integrationPointsRural=outRural, 
                                                      nSubAPerPoint=nSubAPerPoint, 
                                                      nSubRPerPoint=nSubRPerPoint, 
                                                      testMode=testMode)
      
      subPts = tempNewWsForPlotting$subPts
      goodSubPts = tempNewWsForPlotting$goodPts
      subWs = tempNewWsForPlotting$updatedSubWs
      
      if(smallTempUrbanVals) {
        thisIntPts = outUrban
        thisIntWs = tempNewWsForPlotting$wUrban
      } else {
        thisIntPts = outRural
        thisIntWs = tempNewWsForPlotting$wRural
      }
      
      return(list(centerCoords=smallTempCoords, isUrban=smallTempUrbanVals, 
                  intPts=thisIntPts, intWs=thisIntWs, 
                  subPts=subPts, goodSubPts=goodSubPts, subWs=subWs))
    }
    
    tempNewWs = updateWeightsByAdminArea(coords=tempCoords, urbanVals=tempUrbanVals, 
                                         adminMap=adminMap, 
                                         integrationPointsUrban=outUrban, 
                                         integrationPointsRural=outRural, 
                                         nSubAPerPoint=nSubAPerPoint, 
                                         nSubRPerPoint=nSubRPerPoint)
    
    # update the weights with the new values
    wUrban[urbanUpdateI,] = tempNewWs$wUrban
    wRural[ruralUpdateI,] = tempNewWs$wRural
  }
  else {
    area = NULL
  }
  
  if(popPrior || getCovariates) {
    out = load("savedOutput/global/covariates.RData")
    
    # multiply weights by their associated populations and renormalize after
    # xUrban and yUrban so they are [K x nObs] x 1 after c()
    urbanPts = proj(cbind(c(xUrban), c(yUrban)), inverse=TRUE)
    ruralPts = proj(cbind(c(xRural), c(yRural)), inverse=TRUE)
    
    spCoordsUrban = SpatialPoints(urbanPts, adminMap@proj4string)
    spCoordsRural = SpatialPoints(ruralPts, adminMap@proj4string)
    popValsUrb = extract(pop, spCoordsUrban, method=extractMethod)
    popValsRur = extract(pop, spCoordsRural, method=extractMethod)
    popValsUrb[is.na(popValsUrb)] = 0
    popValsRur[is.na(popValsRur)] = 0
  }
  
  # return list of all matrices
  intPtsDHS = list(xUrban=xUrban, yUrban=yUrban, wUrban=wUrban, areasUrban=areas[urbanVals], 
             xRural=xRural, yRural=yRural, wRural=wRural, areasRural=areas[!urbanVals])
  
  popValsUrbOrig = popValsUrb
  popValsRurOrig = popValsRur
  if(getCovariates) {
    # coordsLonLatUrb = proj(urbanPts, inverse=TRUE)
    # coordsLonLatRur = proj(ruralPts, inverse=TRUE)
    # spCoordsLonLatUrb = SpatialPoints(coordsLonLatUrb, adminMap@proj4string)
    # spCoordsLonLatRur = SpatialPoints(coordsLonLatRur, adminMap@proj4string)
    
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
      
      urbanicityValsUrb = rep(TRUE, nrow(coordsUrban))
      accessValsUrb = terra::extract(accessNorm, spCoordsUrban, method=extractMethod)
      elevValsUrb = terra::extract(elevNorm, spCoordsUrban, method=extractMethod)
      distValsUrb = terra::extract(minDistRiverLakesNorm, spCoordsUrban, method=extractMethod)
      
      urbanicityValsRur = rep(FALSE, nrow(coordsRural))
      accessValsRur = terra::extract(accessNorm, spCoordsRural, method=extractMethod)
      elevValsRur = terra::extract(elevNorm, spCoordsRural, method=extractMethod)
      distValsRur = terra::extract(minDistRiverLakesNorm, spCoordsRural, method=extractMethod)
      
      load("savedOutput/global/popMeanSDCal.RData")
      popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
      popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
      popValsUrb = (log1p(popValsUrb) - popMean) * (1/popSD)
      popValsRur = (log1p(popValsRur) - popMean) * (1/popSD)
    } else {
      out = load("savedOutput/global/covariates.RData")
      warning("normalized set to FALSE, but normalization is good practice here")
      
      inf = sessionInfo()
      if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
        urb@file@name = "~/git/jittering/savedOutput/global/urb.tif"
        access@file@name = "~/git/jittering/savedOutput/global/access.tif"
        elev@file@name = "~/git/jittering/savedOutput/global/elev.tif"
        dist@file@name = "~/git/jittering/savedOutput/global/dist.tif"
      }
      
      urbanicityValsUrb = rep(TRUE, nrow(coordsUrban))
      accessValsUrb = terra::extract(access, spCoordsUrban, method=extractMethod)
      elevValsUrb = terra::extract(elev, spCoordsUrban, method=extractMethod)
      distValsUrb = terra::extract(minDistRiverLakes, spCoordsUrban, method=extractMethod)
      
      urbanicityValsRur = rep(FALSE, nrow(coordsRural))
      accessValsRur = terra::extract(access, spCoordsRural, method=extractMethod)
      elevValsRur = terra::extract(elev, spCoordsRural, method=extractMethod)
      distValsRur = terra::extract(minDistRiverLakes, spCoordsRural, method=extractMethod)
    }
    
    if(setMissingToAvg) {
      # urban covariates first
      if(normalized) {
        accessValsUrb[is.na(accessValsUrb)] = 0
        elevValsUrb[is.na(elevValsUrb)] = 0
        distValsUrb[is.na(distValsUrb)] = 0
        
        accessValsRur[is.na(accessValsRur)] = 0
        elevValsRur[is.na(elevValsRur)] = 0
        distValsRur[is.na(distValsRur)] = 0
      }
    }
    
    # set NA covariate points and points with zero population to have 0 weight 
    # unless all integration points for a given cluster would get 0 weight
    naCovRowsUrb = is.na(urbanicityValsUrb) | is.na(accessValsUrb) | is.na(elevValsUrb) | 
      is.na(distValsUrb) | (popValsUrbOrig == 0)
    naCovRowsRur = is.na(urbanicityValsRur) | is.na(accessValsRur) | is.na(elevValsRur) | 
      is.na(distValsRur) | (popValsRurOrig == 0)
    
    # combine everything together into final [K x nObs] x nPar matrix of covariates
    # the following is [nObs x k] x nPar by construction of spCoordsRural and spCoordsUrban
    covsUrb = cbind(int=1, urban=1, access=accessValsUrb, elev=elevValsUrb, distRiversLakes=distValsUrb, pop=popValsUrb)
    covsRur = cbind(int=1, urban=0, access=accessValsRur, elev=elevValsRur, distRiversLakes=distValsRur, pop=popValsRur)
    
    # make sure there are no NAs. They will be given zero weight anyway
    covsUrb[is.na(covsUrb)] = 0
    covsRur[is.na(covsRur)] = 0
  }
  else {
    covsUrb = covsRur = NULL
  }
  
  # adjust weights so that points with any NA covariates are given 0 weight
  if(getCovariates) {
    # make sure NA covariate values get 0 weight
    notnaCovMatUrb = matrix(!naCovRowsUrb, nrow=nrow(xUrban), ncol=ncol(xUrban))
    notnaCovMatRur = matrix(!naCovRowsRur, nrow=nrow(xRural), ncol=ncol(xRural))
    mode(notnaCovMatUrb) = "numeric"
    mode(notnaCovMatRur) = "numeric"
    wUrban = wUrban * notnaCovMatUrb
    wRural = wRural * notnaCovMatRur
  }
  
  # multiply weights by population density prior and renormalize
  if(popPrior) {
    popMatUrban = matrix(popValsUrbOrig, nrow=nrow(xUrban), ncol=ncol(xUrban))
    popMatRural = matrix(popValsRurOrig, nrow=nrow(xRural), ncol=ncol(xRural))
    wUrbanTemp = wUrban * popMatUrban
    wRuralTemp = wRural * popMatRural
    rowSumsUrban = rowSums(wUrbanTemp)
    rowSumsRural = rowSums(wRuralTemp)
    
    # just in case, only multiply weights by population density if it is nonzero for at least 1 int point
    nonzeroRowsUrban = rowSumsUrban != 0
    nonzeroRowsRural = rowSumsRural != 0
    wUrban[nonzeroRowsUrban,] = sweep(wUrbanTemp[nonzeroRowsUrban,], 1, 1/rowSumsUrban[nonzeroRowsUrban], "*")
    wRural[nonzeroRowsRural,] = sweep(wRuralTemp[nonzeroRowsRural,], 1, 1/rowSumsRural[nonzeroRowsRural], "*")
    
    intPtsDHS$wUrban = wUrban
    intPtsDHS$wRural = wRural
  }
  
  intPtsDHS$covsUrb = covsUrb
  intPtsDHS$covsRur = covsRur
  intPtsDHS$wUrban = wUrban
  intPtsDHS$wRural = wRural
  
  save(intPtsDHS, file=outFile)
  
  intPtsDHS
}

# construct integration points as well as weights. Same as makeAllIntegrationPointsDHSold, 
# except makes covariates (population density in particular) in the same way as MICS
# Output: 3 pairs of matrices of dimension nCoordsInStratum x nIntegrationPointsInStratum, 
#         each pair contains one urban matrix and one equivalent rural matrix
#   x: x/easting coordinates
#   y: y/northing coordinates
#   w: integration weights
# Input: 
#   coords: 2 column matrix of observation easting/northing coordinates
#   urbanVals: vector of observation urbanicity classifications
#   areaNames: vector of areas associated with points (can only jitter within the same areas)
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
#   popPrior: use population density as prior for updating weights if need be.
#   setMissingToAvg: sets NA covariates to 0 if they are pop or urban or normalized
makeAllIntegrationPointsDHS = function(coords, urbanVals, areaNames=NULL, 
                                       numPointsUrban=11, numPointsRural=16, 
                                       scalingFactor=1, 
                                       JInnerUrban=3, JOuterUrban=0, 
                                       JInnerRural=3, JOuterRural=1, 
                                       integrationPointType=c("mean", "midpoint"), 
                                       adminMap=adm2Full, areaNameVar="NAME_2", nSubAPerPoint=10, nSubRPerPoint=10, 
                                       popPrior=TRUE, testMode=FALSE, proj=projNigeria, 
                                       outFile="savedOutput/global/intPtsDHS.RData", 
                                       getCovariates=TRUE, normalized=TRUE, useThreshPopMat=TRUE, 
                                       extractMethod="bilinear", setMissingToAvg=TRUE) {
  
  # calculate integration points and weights relative to individual points
  outUrban = getIntegrationPointsDHS(urban=TRUE, numPointsUrban, 
                                     scalingFactor, 
                                     JInnerUrban, JOuterUrban, 
                                     integrationPointType, 
                                     verbose=FALSE)
  
  outRural = getIntegrationPointsDHS(urban=FALSE, numPointsRural, 
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
  
  if(!is.null(adminMap)) {
    # if adminMap is input, integration weights will be adjusted using a much 
    # finer "sub" integration grid
    
    # first subset jittered points that are close enough to the border to have 
    # a possibility of being adjusted
    
    # determine if points are within max distance of admin area:
    # first calculate max distance from admin area
    maxUrbanDistance = 2 * scalingFactor
    maxRuralDistance = 10 * scalingFactor
    maxDist = rep(maxRuralDistance, nrow(coords))
    maxDist[urbanVals] = maxUrbanDistance
    
    # get all points and map data in sp format
    coordsLonLat = proj(coords, inverse=TRUE)
    spCoordsLonLat = SpatialPoints(coordsLonLat, adminMap@proj4string)
    adminMapPolygons = as.SpatialPolygons.PolygonsList(adminMap@polygons, adminMap@proj4string)
    
    # determine what admin area each point is in
    temp = over(spCoordsLonLat, adminMap, returnList=FALSE)
    nas = is.na(temp[[areaNameVar]])
    
    # calculate distances to admin boundaries for unknown points
    require(geosphere)
    naClosestIDs = sapply(which(nas), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons)[4]})
    adminID = rep(1, nrow(temp))
    adminID[nas] = naClosestIDs
    adminID[!nas] = match(temp[[areaNameVar]][!nas], adminMap[[areaNameVar]])
    
    # set names of what area each point is in if need be
    if(is.null(areaNames)) {
      areas = adminMap[[areaNameVar]][adminID]
    } else {
      areas = areaNames
    }
    
    # calculate distances to admin boundaries
    dists = sapply(1:nrow(coords), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons[adminID[ind]])[1]}) * (1/1000)
    
    # set whether or not to update weights based on distance to admin boundaries
    updateI = dists < maxDist
    urbanUpdateI = updateI[urbanVals]
    ruralUpdateI = updateI[!urbanVals]
    
    # calculate updated weights for the integration points near the borders
    tempCoords = coords[updateI,]
    tempUrbanVals = urbanVals[updateI]
    require(fields)
    
    if(testMode) {
      # in this case, we take only one set of coords that is very close to border, 
      # and adjust its weight, saving relevant results for plotting
      tempDists = dists[updateI]
      closeCoords = which.min(tempDists)
      smallTempCoords = matrix(tempCoords[closeCoords,], nrow=1)
      smallTempUrbanVals = tempUrbanVals[closeCoords]
      
      tempNewWsForPlotting = updateWeightsByAdminArea(coords=smallTempCoords, urbanVals=smallTempUrbanVals, 
                                                      adminMap=adminMap, areas=areas[which(updateI)[closeCoords]], 
                                                      areaNameVar=areaNameVar, 
                                                      integrationPointsUrban=outUrban, 
                                                      integrationPointsRural=outRural, 
                                                      nSubAPerPoint=nSubAPerPoint, 
                                                      nSubRPerPoint=nSubRPerPoint, 
                                                      testMode=testMode)
      
      subPts = tempNewWsForPlotting$subPts
      goodSubPts = tempNewWsForPlotting$goodPts
      subWs = tempNewWsForPlotting$updatedSubWs
      
      if(smallTempUrbanVals) {
        thisIntPts = outUrban
        thisIntWs = tempNewWsForPlotting$wUrban
      } else {
        thisIntPts = outRural
        thisIntWs = tempNewWsForPlotting$wRural
      }
      
      return(list(centerCoords=smallTempCoords, isUrban=smallTempUrbanVals, 
                  intPts=thisIntPts, intWs=thisIntWs, 
                  subPts=subPts, goodSubPts=goodSubPts, subWs=subWs))
    }
    
    tempNewWs = updateWeightsByAdminArea(coords=tempCoords, urbanVals=tempUrbanVals, 
                                         adminMap=adminMap, areas=areas[updateI], 
                                         areaNameVar=areaNameVar, 
                                         integrationPointsUrban=outUrban, 
                                         integrationPointsRural=outRural, 
                                         nSubAPerPoint=nSubAPerPoint, 
                                         nSubRPerPoint=nSubRPerPoint)
    
    # update the weights with the new values
    wUrban[urbanUpdateI,] = tempNewWs$wUrban
    wRural[ruralUpdateI,] = tempNewWs$wRural
  }
  else {
    area = NULL
  }
  
  if(popPrior || getCovariates) {
    out = load("savedOutput/global/covariates.RData")
    
    # multiply weights by their associated populations and renormalize after
    # xUrban and yUrban so they are [K x nObs] x 1 after c()
    urbanPts = proj(cbind(c(xUrban), c(yUrban)), inverse=TRUE)
    ruralPts = proj(cbind(c(xRural), c(yRural)), inverse=TRUE)
    
    spCoordsUrban = SpatialPoints(urbanPts, adminMap@proj4string)
    spCoordsRural = SpatialPoints(ruralPts, adminMap@proj4string)
    popValsUrb = terra::extract(pop, spCoordsUrban, method=extractMethod)
    popValsRur = terra::extract(pop, spCoordsRural, method=extractMethod)
    popValsUrb[is.na(popValsUrb)] = 0
    popValsRur[is.na(popValsRur)] = 0
  }
  
  # return list of all matrices
  intPtsDHS = list(xUrban=xUrban, yUrban=yUrban, wUrban=wUrban, areasUrban=areas[urbanVals], 
                   xRural=xRural, yRural=yRural, wRural=wRural, areasRural=areas[!urbanVals])
  
  popValsUrbOrig = popValsUrb
  popValsRurOrig = popValsRur
  if(getCovariates) {
    covsUrb = getDesignMatPopNorm(spCoordsUrban@coords, useThreshPopMat=useThreshPopMat, 
                                   proj=proj, testMode=testMode, setMissingToAvg=TRUE)
    covsRur = getDesignMatPopNorm(spCoordsRural@coords, useThreshPopMat=useThreshPopMat, 
                                   proj=proj, testMode=testMode, setMissingToAvg=TRUE)
    
    covsUrb = covsUrb[,c(1, 3, 4, 5, 6, 2)]
    names(covsUrb)[2] = "urban"
    
    covsRur = covsRur[,c(1, 3, 4, 5, 6, 2)]
    names(covsRur)[2] = "urban"
    
    # figure out what rows have NAs
    naCovRowsUrb = apply(covsUrb, 1, anyNA)
    naCovRowsRur = apply(covsRur, 1, anyNA)
    
    # make sure integration points correspond to the same urbanicity as the corresponding data
    covsUrb[,2] = 1
    covsRur[,2] = 0
    
    # make sure there are no NAs. They will be given zero weight anyway
    covsUrb[is.na(covsUrb)] = 0
    covsRur[is.na(covsRur)] = 0
  }
  else {
    covsUrb = covsRur = NULL
  }
  
  # adjust weights so that points with any NA covariates are given 0 weight
  if(getCovariates) {
    # make sure NA covariate values get 0 weight
    notnaCovMatUrb = matrix(!naCovRowsUrb, nrow=nrow(xUrban), ncol=ncol(xUrban))
    notnaCovMatRur = matrix(!naCovRowsRur, nrow=nrow(xRural), ncol=ncol(xRural))
    mode(notnaCovMatUrb) = "numeric"
    mode(notnaCovMatRur) = "numeric"
    wUrban = wUrban * notnaCovMatUrb
    wRural = wRural * notnaCovMatRur
  }
  
  # multiply weights by population density prior and renormalize
  if(popPrior) {
    popMatUrban = matrix(popValsUrbOrig, nrow=nrow(xUrban), ncol=ncol(xUrban))
    popMatRural = matrix(popValsRurOrig, nrow=nrow(xRural), ncol=ncol(xRural))
    wUrbanTemp = wUrban * popMatUrban
    wRuralTemp = wRural * popMatRural
    rowSumsUrban = rowSums(wUrbanTemp)
    rowSumsRural = rowSums(wRuralTemp)
    
    # just in case, only multiply weights by population density if it is nonzero for at least 1 int point
    nonzeroRowsUrban = rowSumsUrban != 0
    nonzeroRowsRural = rowSumsRural != 0
    wUrban[nonzeroRowsUrban,] = sweep(wUrbanTemp[nonzeroRowsUrban,], 1, 1/rowSumsUrban[nonzeroRowsUrban], "*")
    wRural[nonzeroRowsRural,] = sweep(wRuralTemp[nonzeroRowsRural,], 1, 1/rowSumsRural[nonzeroRowsRural], "*")
    
    intPtsDHS$wUrban = wUrban
    intPtsDHS$wRural = wRural
  }
  
  intPtsDHS$covsUrb = covsUrb
  intPtsDHS$covsRur = covsRur
  intPtsDHS$wUrban = wUrban
  intPtsDHS$wRural = wRural
  
  save(intPtsDHS, file=outFile)
  
  intPtsDHS
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
# Details: FOR THE SPDE MODEL!!!!
makeJitterDataForTMB_SPDE = function(integrationPointInfo, ys, urbanicity, ns, spdeMesh) {
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
# Details: FOR THE SPDE MODEL!!!!
makeJitterDataForTMB_BYM2 = function(intPtInfoDHS, ysDHS, urbanicityDHS, nsDHS, 
                                     intPtInfoMICS, ysMICS, urbanicityMICS, nsMICS, 
                                     mapDat=admFinalFull, proj=projNigeria) {
  
  # first extract the integration point information
  xUrbanDHS = intPtInfoDHS$xUrban
  yUrbanDHS = intPtInfoDHS$yUrban
  wUrbanDHS = intPtInfoDHS$wUrban
  xRuralDHS = intPtInfoDHS$xRural
  yRuralDHS = intPtInfoDHS$yRural
  wRuralDHS = intPtInfoDHS$wRural
  
  # get the long set of coordinates
  coordsUrbanDHS = cbind(c(xUrbanDHS), c(yUrbanDHS))
  coordsRuralDHS = cbind(c(xRuralDHS), c(yRuralDHS))
  
  # separate observations by urbanicity
  ysUrbanDHS = ysDHS[urbanicity]
  ysRuralDHS = ysDHS[!urbanicity]
  nsUrbanDHS = nsDHS[urbanicity]
  nsRuralDHS = nsDHS[!urbanicity]
  
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
# lambda: spatial scaling coefficient. Default is 1 / priorRange for 
#   priorRange = (domainDiameter / 5) / 2
# domainDiameter: used for calculating default lambda. Default is 1463.733 km, 
#   Nigeria's diameter
getIntegrationPointsMICS = function(strat, kmresFineStart=2.5, numPtsUrb=25, numPtsRur=25, 
                                    stratumMICSMapDat=admFinalFull, stratumMICSNameVar="NAME_FINAL", 
                                    subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                    poppsub=poppsubNGAThresh, 
                                    normalized=TRUE, useThreshPopMat=TRUE, 
                                    proj=projNigeria, projArea=projNigeriaArea, 
                                    spatialAsCovariate=FALSE, adm2AsCovariate=FALSE, 
                                    lambda=NULL, domainDiameter=NULL, 
                                    returnFineGrid=FALSE, testMode=FALSE, extractMethod="bilinear") {
  
  
  fineIntPtsTab = getFineIntPointsInfoMICS(stratumName=strat, kmresStart=kmresFineStart, 
                                           minPointsUrb=numPtsUrb, minPointsRur=numPtsRur, 
                                           stratumMICSMapDat=stratumMICSMapDat, stratumMICSNameVar=stratumMICSNameVar, 
                                           subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                           poppsub=poppsub, 
                                           normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                           proj=proj, projArea=projArea, testMode=testMode, extractMethod=extractMethod)
  
  X = fineIntPtsTab[,c(11:14, 16)] # don't use urbanicity to generate clustering
  pop = fineIntPtsTab$pop
  urb = fineIntPtsTab$urban
  ENCoords = cbind(fineIntPtsTab$east, fineIntPtsTab$north)
  
  if(spatialAsCovariate) {
    # add spatial coordinates as covariates if user requests
    if(is.null(domainDiameter)) {
      # adm0Boundary = gBoundary(adm0)
      # adm0BoundaryCoords = do.call("rbind", lapply(adm0Boundary@lines[[1]]@Lines, function(x) {x@coords}))
      # domainDiameter = max(rdist.earth(adm0BoundaryCoords, miles=FALSE))
      domainDiameter = 1463.733 # in km
    }
    
    # normalize spatial coordinates based on prior median effective range
    if(is.null(lambda)) {
      priorRange = (domainDiameter / 5) / 2
      lambda = 1 / priorRange
    }
    
    ENCoordsNorm = sweep(ENCoords, 2, colMeans(ENCoords), FUN="-")
    ENCoordsNorm = ENCoordsNorm * lambda
    
    X = cbind(X, ENCoordsNorm)
  }
  
  subareaIDs = as.numeric(factor(fineIntPtsTab$subarea))
  if(adm2AsCovariate) {
    
    # get new covariate: centers of subareas
    meanXs = aggregate(ENCoords[,1], by=list(subareas=subareaIDs), mean)
    meanYs = aggregate(ENCoords[,2], by=list(subareas=subareaIDs), mean)
    thisXs = meanXs$x[match(subareaIDs, meanXs$subareas)]
    thisYs = meanYs$x[match(subareaIDs, meanYs$subareas)]
    
    # test if stratum spans more distance vertically or horizontally. Order 
    # subareas along longest dimension by their mean coordinate (DONT DO THIS)
    # newSubareaIDs = order(meanCoords$x)
    
    # scale and center subareaIDs
    # subareaCov = (newSubareaIDs - mean(newSubareaIDs))/sd(newSubareaIDs)
    
    if(is.null(lambda)) {
      lambda = 1
    }
    
    # scale and center new coordinates
    nSubareas = length(unique(subareaIDs))
    thisXs = lambda * nSubareas * (thisXs - mean(thisXs))/sd(thisXs)
    thisYs = lambda * nSubareas * (thisYs - mean(thisYs))/sd(thisYs)
    subareaCov = cbind(thisXs, thisYs)
    X = cbind(X, subareaCov)
  }
  
  nPtsUrb = sum(urb, na.rm=TRUE)
  nPtsRur = sum(!urb, na.rm=TRUE)
  
  # if((nPtsUrb != numPtsUrb) || (nPtsRur != numPtsRur)) {
  #   browser()
  # }
  
  hasUrbPop = TRUE
  hasRurPop = TRUE
  if(nPtsUrb > 0) {
    # do the weighted K-medoids (technically PAM): 
    # Maechler, M., P. Rousseeuw, A. Struyf, M. Hubert and K. Hornik (2011).
    # cluster: Cluster Analysis Basics and Extensions. R package version 1.14.1
    XmatUrb = matrix(unlist(X[urb,]), ncol=ncol(X))
    # if(all(rep(t(XmatUrb)[,1], times=nrow(XmatUrb)) == c(t(XmatUrb)))) {
    #   XmatUrb = cbind(XmatUrb, ENCoordsNorm[urb,])
    # }
      
    distMatUrb = rdist(XmatUrb, XmatUrb)
    medoidsUrb = wcKMedoids(distMatUrb^2, numPtsUrb, weights=pop[urb], method="PAM")
    medoidIUrb = medoidsUrb$clustering
    
    # calculate weights of the medoids
    totalPopUrb = aggregate(pop[urb], by=list(medoid=medoidIUrb), FUN=sum)
    weightsUrb = totalPopUrb[,2] / sum(totalPopUrb[,2])
    
    adm2Wurb = matrix(as.numeric(outer(subareaIDs[urb], sort(unique(subareaIDs)), FUN="==")), ncol=length(unique(subareaIDs)))
    XmatUrb = cbind(XmatUrb, adm2Wurb)
    
    # calculate averages of covariates, both fine scale and at the integration points
    finePopWeights = pop[urb]
    finePopWeights[is.na(finePopWeights)] = 0
    finePopWeights = finePopWeights/sum(finePopWeights)
    fineAvgsUrb = colSums(sweep(XmatUrb, 1, finePopWeights, "*"))
    
    intAvgsUrb = colSums(sweep(XmatUrb[sort(unique(medoidIUrb)),], 1, weightsUrb, "*"))
    
    # calculate ranges of covariates, both fine scale and at the integration points 
    fineSDUrb = apply(XmatUrb, 2, wtdSD, weights=finePopWeights)
    intSDUrb = apply(XmatUrb[sort(unique(medoidIUrb)),], 2, wtdSD, weights=weightsUrb)
    
    # get weights of admin2 areas
    fineAdm2Wurb = colSums(sweep(adm2Wurb, 1, finePopWeights, "*"))
    intAdm2Wurb = colSums(sweep(adm2Wurb[sort(unique(medoidIUrb)),], 1, weightsUrb, "*"))
  } else {
    hasUrbPop = FALSE
    XmatUrb = NULL
    distMatUrb = NULL
    medoidsUrb = NULL
    medoidIUrb = NULL
    
    totalPopUrb = NULL
    weightsUrb = NULL
    
    fineAvgsUrb = NULL
    intAvgsUrb = NULL
    fineSDUrb = NULL
    intSDUrb = NULL
    
    fineAdm2Wrur = NULL
    intAdm2Wrur = NULL
  }
  
  if(nPtsRur > 0) {
    XmatRur = matrix(unlist(X[!urb,]), ncol=ncol(X))
    # if(all(rep(t(XmatRur)[,1], times=nrow(XmatRur)) == c(t(XmatRur)))) {
    #   XmatRur = cbind(XmatRur, ENCoordsNorm[!urb,])
    # }
    
    distMatRur = rdist(XmatRur, XmatRur)
    medoidsRur = wcKMedoids(distMatRur^2, numPtsRur, weights=pop[!urb], method="PAM")
    medoidIRur = medoidsRur$clustering
    
    totalPopRur = aggregate(pop[!urb], by=list(medoid=medoidIRur), FUN=sum)
    weightsRur = totalPopRur[,2] / sum(totalPopRur[,2])
    
    adm2Wrur = matrix(as.numeric(outer(subareaIDs[!urb], sort(unique(subareaIDs)), FUN="==")), ncol=length(unique(subareaIDs)))
    XmatRur = cbind(XmatRur, adm2Wrur)
    
    # calculate averages of covariates, both fine scale and at the integration points
    finePopWeights = pop[!urb]
    finePopWeights[is.na(finePopWeights)] = 0
    finePopWeights = finePopWeights/sum(finePopWeights)
    fineAvgsRur = colSums(sweep(XmatRur, 1, finePopWeights, "*"))
    
    intAvgsRur = colSums(sweep(XmatRur[sort(unique(medoidIRur)),], 1, weightsRur, "*"))
    
    # calculate ranges of covariates, both fine scale and at the integration points 
    fineSDRur = apply(XmatRur, 2, wtdSD, weights=finePopWeights)
    intSDRur = apply(XmatRur[sort(unique(medoidIRur)),], 2, wtdSD, weights=weightsRur)
    
    # get weights of admin2 areas
    fineAdm2Wrur = colSums(sweep(adm2Wrur, 1, finePopWeights, "*"))
    intAdm2Wrur = colSums(sweep(adm2Wrur[sort(unique(medoidIRur)),], 1, weightsRur, "*"))
  } else {
    hasRurPop = FALSE
    XmatRur = NULL
    distMatRur = NULL
    medoidsRur = NULL
    medoidIRur = NULL
    
    totalPopRur = NULL
    weightsRur = NULL
    
    fineAvgsRur = NULL
    intAvgsRur = NULL
    fineSDRur = NULL
    intSDRur = NULL
    
    fineAdm2Wrur = NULL
    intAdm2Wrur = NULL
  }
  
  # # assign fine grid points to each medoid
  # distToMedoidsUrb = distMatUrb[,medoidIUrb]
  # medoidAssignedUrb = apply(distToMedoidsUrb, 1, which.min)
  # pointIAssignedUrb = medoidIUrb[medoidAssignedUrb]
  # 
  # distToMedoidsRur = distMatRur[,medoidIRur]
  # medoidAssignedRur = apply(distToMedoidsRur, 1, which.min)
  # pointIAssignedRur = medoidIRur[medoidAssignedRur]
  
  # return results, and calculate pop averages of covariates compared to with 
  # the fine integration points
  sortIUrb = which(fineIntPtsTab$urban)[sort(unique(medoidIUrb))]
  sortIRur = which(!fineIntPtsTab$urban)[sort(unique(medoidIRur))]
  
  out = list(strat=strat, hasUrbPop=hasUrbPop, hasRurPop=hasRurPop, 
       pts=fineIntPtsTab[c(sortIUrb, sortIRur),], weights=c(weightsUrb, weightsRur), 
       ptsUrb=fineIntPtsTab[sortIUrb,], weightsUrb=weightsUrb, 
       ptsRur=fineIntPtsTab[sortIRur,], weightsRur=weightsRur, 
       fineAvgsUrb=fineAvgsUrb, intAvgsUrb=intAvgsUrb, 
       fineAvgsRur=fineAvgsRur, intAvgsRur=intAvgsRur, 
       fineSDUrb=fineSDUrb, intSDUrb=intSDUrb, 
       fineSDRur=fineSDRur, intSDRur=intSDRur, 
       fineAdm2Wurb = fineAdm2Wurb, 
       intAdm2Wurb = intAdm2Wurb, 
       fineAdm2Wrur = fineAdm2Wrur, 
       intAdm2Wrur = intAdm2Wrur)
  
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
# lambda: spatial scaling coefficient. Default is 1 / priorRange for 
#   priorRange = (domainDiameter / 5) / 2
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
                                        spatialAsCovariate=FALSE, adm2AsCovariate=FALSE, 
                                        lambda=NULL, domainDiameter=NULL, 
                                        fileNameRoot="MICSintPts_", loadSavedIntPoints=TRUE, 
                                        extractMethod="bilinear", outFile=NULL) {
  
  # normalize spatial coordinates based on prior median effective range
  if(is.null(lambda)) {
    if(!adm2AsCovariate) {
      priorRange = (domainDiameter / 5) / 2
      lambda = 1 / priorRange
    } else {
      lambda = 1
    }
  }
  
  if(is.null(outFile)) {
    lambdaText = ifelse(lambda == 1, "", paste0("_lam", round(lambda, 2)))
    adm2CovText = ""
    if(adm2AsCovariate) {
      adm2CovText = paste0("_adm2Cov", lambdaText)
    }
    if((numPtsUrb == numPtsRur) && (numPtsUrb == 25)) {
      outFile = paste0("savedOutput/global/intPtsMICS", adm2CovText, ".RData")
    } else if(numPtsUrb == numPtsRur) {
      outFile = paste0("savedOutput/global/intPtsMICS_", numPtsUrb, adm2CovText, ".RData")
    } else {
      outFile = paste0("savedOutput/global/intPtsMICS_u", numPtsUrb, "_r", numPtsRur, adm2CovText, ".RData")
    }
  }
  
  if(is.null(domainDiameter)) {
    domainDiameter = 1463.733 # in km
  }
  
  # set file name root for saving the results
  fileNameRoot = paste0(fileNameRoot, "_km", kmresFineStart, 
                        "_nPtU", numPtsUrb, "_nPtR", numPtsRur, 
                        "_norm", as.numeric(normalized), 
                        "_thresh", as.numeric(useThreshPopMat), 
                        "_spatCov", as.numeric(spatialAsCovariate), 
                        "_lam", round(lambda, 4), 
                        "_adm2Cov", as.numeric(adm2AsCovariate))
  
  # For each stratum, generate the integration points
  allIntPts = list()
  strataMICS = sort(stratumMICSMapDat@data[[stratumMICSNameVar]])
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
                                              poppsub=poppsub, adm2AsCovariate=adm2AsCovariate, 
                                              normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                              proj=proj, projArea=projArea, 
                                              spatialAsCovariate=spatialAsCovariate, 
                                              lambda=lambda, domainDiameter=domainDiameter, 
                                              extractMethod=extractMethod)
      
      save(thisIntPoints, file=paste0(fileNameRoot, "_i", i, ".RData"))
    }
    
    print(paste0("fineIntPts cov averages (urban): ", paste(thisIntPoints$fineAvgsUrb, collapse= " ")))
    print(paste0("intPts cov averages (urban): ", paste(thisIntPoints$intAvgsUrb, collapse= " ")))
    print(paste0("fineIntPts cov averages (rural): ", paste(thisIntPoints$fineAvgsRur, collapse= " ")))
    print(paste0("intPts cov averages (rural): ", paste(thisIntPoints$intAvgsRur, collapse= " ")))
    print(paste0("fineIntPts cov SDs (urban): ", paste(thisIntPoints$fineSDUrb, collapse= " ")))
    print(paste0("intPts cov SDs (urban): ", paste(thisIntPoints$intSDUrb, collapse= " ")))
    print(paste0("fineIntPts cov SDs (rural): ", paste(thisIntPoints$fineSDRur, collapse= " ")))
    print(paste0("intPts cov SDs (rural): ", paste(thisIntPoints$intSDRur, collapse= " ")))
    
    if(length(thisIntPoints$weightsUrb) == 0) {
      browser()
    }
    
    allIntPts = c(allIntPts, list(c(thisIntPoints, list(strat=thisStrat))))
  }
  
  # concatenate integration points and weights into matrices
  allStrat = sapply(allIntPts, function(x) {x$strat})
  getCovIUrb = function(i) {
    thisPtsUrbI = lapply(allIntPts, function(x) {x$ptsUrb[,i]})
    lengths = sapply(thisPtsUrbI, length)
    thisLength = max(lengths)
    thisMissing = which(lengths == 0)
    for(j in thisMissing) {
      thisPtsUrbI[[thisMissing]] = rep(99, thisLength)
    }
    c(do.call("rbind", thisPtsUrbI))
  }
  getCovIRur = function(i) {
    thisPtsRurI = lapply(allIntPts, function(x) {x$ptsRur[,i]})
    lengths = sapply(thisPtsRurI, length)
    thisLength = max(lengths)
    thisMissing = which(lengths == 0)
    for(j in thisMissing) {
      thisPtsRurI[[thisMissing]] = rep(99, thisLength)
    }
    c(do.call("rbind", thisPtsRurI))
  }
  allCovsUrb = data.frame(lapply(1:ncol(allIntPts[[1]]$ptsUrb), getCovIUrb))
  names(allCovsUrb) = names(allIntPts[[1]]$ptsUrb)
  allCovsRur = data.frame(lapply(1:ncol(allIntPts[[1]]$ptsRur), getCovIRur))
  names(allCovsRur) = names(allIntPts[[1]]$ptsRur)
  
  wsUrban = do.call("rbind", (lapply(allIntPts, function(x) {x$weightsUrb})))
  wsRural = do.call("rbind", (lapply(allIntPts, function(x) {x$weightsRur})))
  
  # get average value of covariates
  fineIntPtAvgsUrb = do.call("rbind", lapply(allIntPts, function(x) {x$fineAvgsUrb}))
  fineIntPtAvgsRur = do.call("rbind", lapply(allIntPts, function(x) {x$fineAvgsRur}))
  intPtAvgsUrb = do.call("rbind", lapply(allIntPts, function(x) {x$intAvgsUrb}))
  intPtAvgsRur = do.call("rbind", lapply(allIntPts, function(x) {x$intAvgsRur}))
  errorUrb = intPtAvgsUrb - fineIntPtAvgsUrb
  errorRur = intPtAvgsRur - fineIntPtAvgsRur
  error = rbind(errorUrb, errorRur)
  absPctErrorUrb = abs(100*(intPtAvgsUrb - fineIntPtAvgsUrb)/fineIntPtAvgsUrb)
  absPctErrorRur = abs(100*(intPtAvgsRur - fineIntPtAvgsRur)/fineIntPtAvgsRur)
  absPctError = rbind(absPctErrorUrb, 
                      absPctErrorRur)
  absMeanPctErrorUrb = colMeans(absPctErrorUrb)
  absMeanPctErrorRur = colMeans(absPctErrorRur)
  absMeanPctError = colMeans(absPctError)
  absMaxPctErrorUrb = apply(absPctErrorUrb, 2, max)
  absMaxPctErrorRur = apply(absPctErrorRur, 2, max)
  absMaxPctError = apply(absPctError, 2, max)
  
  fineIntPtSDUrb = do.call("rbind", lapply(allIntPts, function(x) {x$fineSDUrb}))
  fineIntPtSDRur = do.call("rbind", lapply(allIntPts, function(x) {x$fineSDRur}))
  intPtSDUrb = do.call("rbind", lapply(allIntPts, function(x) {x$intSDUrb}))
  intPtSDRur = do.call("rbind", lapply(allIntPts, function(x) {x$intSDRur}))
  absPctErrorSDUrb = abs(100*(intPtSDUrb - fineIntPtSDUrb)/fineIntPtSDUrb)
  absPctErrorSDRur = abs(100*(intPtSDRur - fineIntPtSDRur)/fineIntPtSDRur)
  absPctErrorSD = rbind(absPctErrorSDUrb, 
                      absPctErrorSDRur)
  absMeanPctErrorSDUrb = colMeans(absPctErrorSDUrb)
  absMeanPctErrorSDRur = colMeans(absPctErrorSDRur)
  absMeanPctErrorSD = colMeans(absPctErrorSD)
  absMaxPctErrorSDUrb = apply(absPctErrorSDUrb, 2, max)
  absMaxPctErrorSDRur = apply(absPctErrorSDRur, 2, max)
  absMaxPctErrorSD = apply(absPctErrorSD, 2, max)
  
  # make sure to fill in gaps where there is no urban or rural population
  hasUrbPop = sapply(allIntPts, function(x) {x$hasUrbPop})
  hasRurPop = sapply(allIntPts, function(x) {x$hasRurPop})
  hasUrbPopI = which(hasUrbPop)
  hasRurPopI = which(hasRurPop)
  
  # first expand weight matrices to include the missing rows
  if(any(!hasUrbPop)) {
    # first fix wsUrb (extra rows get zero weight)
    wsUrbanTemp = matrix(0, nrow=length(allIntPts), ncol=ncol(wsUrban))
    wsUrbanTemp[hasUrbPop,] = wsUrban
    wsUrban = wsUrbanTemp
    
    # # then fix allCovsUrb. Replace missing rows with the urban rows. 
    # # They will get zero weight anyway
    # numMissingRows = numPtsUrb*length(allIntPts) - nrow(allCovsUrb)
    # tempRows = allCovsUrb[1:numMissingRows,]
    # allCovsUrbTemp = rbind(allCovsUrb, 
    #                        tempRows)
    # hasUrbPopIndsAll = c(unlist(sapply(1:length(hasUrbPopI), function(i) {
    #   startInd = (hasUrbPopI[i] - 1) * numPtsUrb + 1
    #   endInd = startInd + numPtsUrb - 1
    #   startInd:endInd
    # })))
    # 
    # allCovsUrbTemp[hasUrbPopIndsAll,] = allCovsUrb
  } else {
    # allCovsUrbTemp = allCovsUrb
  }
  
  if(any(!hasRurPop)) {
    # first fix wsRural
    wsRuralTemp = matrix(0, nrow=length(allIntPts), ncol=ncol(wsRural))
    wsRuralTemp[hasRurPop,] = wsRural
    wsRural = wsRuralTemp
    
    # # then fix allCovsRur. Replace missing rows with the urban rows. 
    # # They will get zero weight anyway
    # numMissingRows = numPtsRur*length(allIntPts) - nrow(allCovsRur)
    # tempRows = allCovsRur[1:numMissingRows,]
    # 
    # tempRows$strat = allStrat[!hasRurPop]
    # allCovsRurTemp = rbind(allCovsRur, 
    #                        tempRows)
    # hasRurPopIndsAll = c(unlist(sapply(1:length(hasRurPopI), function(i) {
    #   startInd = (hasRurPopI[i] - 1) * numPtsRur + 1
    #   endInd = startInd + numPtsRur - 1
    #   startInd:endInd
    # })))
    # noRurPopIndsAll = setdiff(1:(numPtsRur*length(allIntPts)), hasRurPopIndsAll)
    # 
    # allCovsRurTemp[hasRurPopIndsAll,] = allCovsRur
    # allCovsRurTemp[noRurPopIndsAll,] = 4
  } else {
    # allCovsRurTemp = allCovsRur
  }
  
  # now fill in the missing strata and areas with elements from the opposite stratum
  if(any(!hasUrbPop) || any(!hasRurPop)) {
    if(numPtsRur != numPtsUrb) {
      stop("numPtsRur != numPtsUrb not yet implemented")
    }
    # fill in everything but the covariates of interest, urbanicity
    allCovsUrb[allCovsUrb$strat == 99, c(1:5, 7:9)] = allCovsRur[allCovsUrb$strat == 99, c(1:5, 7:9)]
    allCovsRur[allCovsRur$strat == 99, c(1:5, 7:9)] = allCovsUrb[allCovsRur$strat == 99, c(1:5, 7:9)]
    allCovsUrb$urban = 1
    allCovsRur$urban = 0
  }
  
  if(is.null(datStrata)) {
    if(!is.null(datUrb)) {
      stop("datUrb provided but not datStrata")
    }
    
    # return list of all matrices at the stratum rather than observation level
    intPtsMICS = list(strataMICS=strataMICS, 
         XUrb=allCovsUrb, XRur=allCovsRur, 
         wUrban=wsUrban, wRural=wsRural, 
         absMeanPctErrorUrb=absMeanPctErrorUrb, 
         absMeanPctErrorRur=absMeanPctErrorRur, 
         absMeanPctError=absMeanPctError, 
         absMaxPctErrorUrb=absMaxPctErrorUrb, 
         absMaxPctErrorRur=absMaxPctErrorRur, 
         absMaxPctError=absMaxPctError, 
         absMeanPctErrorSDUrb=absMeanPctErrorSDUrb, 
         absMeanPctErrorSDRur=absMeanPctErrorSDRur, 
         absMeanPctErrorSD=absMeanPctErrorSD, 
         absMaxPctErrorSDUrb=absMaxPctErrorSDUrb, 
         absMaxPctErrorSDRur=absMaxPctErrorSDRur, 
         absMaxPctErrorSD=absMaxPctErrorSD, 
         errorUrb=errorUrb, 
         errorRur=errorRur, 
         error=error, 
         fineIntPtAvgsUrb=fineIntPtAvgsUrb, 
         fineIntPtAvgsRur=fineIntPtAvgsRur, 
         intPtAvgsUrb=intPtAvgsUrb,
         intPtAvgsRur=intPtAvgsRur)
  } else {
    if(is.null(datUrb)) {
      stop("non-stratified integration point construction not currently supported")
    }
    
    # we must expand the matrices to be at the observation rather than stratum level
    strataUrb = datStrata[datUrb]
    strataRur = datStrata[!datUrb]
    strataUrbU = sort(unique(strataUrb))
    strataRurU = sort(unique(strataRur))
    
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
    
    # return list of all matrices and diagnostics at the observation level
    intPtsMICS = list(strataMICS=strataMICS, 
         XUrb=XUrbFull, XRur=XRurFull, 
         wUrban=wsUrban, wRural=wsRural, 
         absMeanPctErrorUrb=absMeanPctErrorUrb, 
         absMeanPctErrorRur=absMeanPctErrorRur, 
         absMeanPctError=absMeanPctError, 
         absMaxPctErrorUrb=absMaxPctErrorUrb, 
         absMaxPctErrorRur=absMaxPctErrorRur, 
         absMaxPctError=absMaxPctError, 
         absMeanPctErrorSDUrb=absMeanPctErrorSDUrb, 
         absMeanPctErrorSDRur=absMeanPctErrorSDRur, 
         absMeanPctErrorSD=absMeanPctErrorSD, 
         absMaxPctErrorSDUrb=absMaxPctErrorSDUrb, 
         absMaxPctErrorSDRur=absMaxPctErrorSDRur, 
         absMaxPctErrorSD=absMaxPctErrorSD, 
         errorUrb=errorUrb, 
         errorRur=errorRur, 
         error=error, 
         fineIntPtAvgsUrb=fineIntPtAvgsUrb, 
         fineIntPtAvgsRur=fineIntPtAvgsRur, 
         intPtAvgsUrb=intPtAvgsUrb,
         intPtAvgsRur=intPtAvgsRur)
  }
  
  save(intPtsMICS, file=outFile)
  
  intPtsMICS
}



# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
# Inputs:
# setMissingToAvg: if TRUE and normalized is TRUE, sets NA covariates to 0
getFineIntPointsInfoMICSold = function(stratumName, kmresStart=2.5, minPointsUrb=20, minPointsRur=20, 
                                    stratumMICSMapDat=admFinalFull, stratumMICSNameVar="NAME_FINAL", 
                                    subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                    poppsub=poppsubNGAThresh, 
                                    normalized=TRUE, useThreshPopMat=TRUE, 
                                    proj=projNigeria, projArea=projNigeriaArea, 
                                    testMode=FALSE, extractMethod="bilinear", 
                                    setMissingToAvg=TRUE) {
  
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
    
    # get subareas associated with the points (considering only subareas within 
    # relevant MICS stratum)
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
      # thisSubarea = theseSubareas[i]
      # thisPopMatI = which(stratumPopMat$subarea == thisSubarea)
      thisPopMat = stratumPopMat #[thisPopMatI,]
      # thisSubareaI = match(thisSubarea, uniqueSubareas)
      
      # dists = rdist(matrix(pts, ncol=2, cbind(thisPopMat$east, thisPopMat$north))
      
      # fullPopMat is 5km resolution, so must be within 2.5 km in easting and northing directions
      closeE = (pts[1] > thisPopMat$east - 2.5) & (pts[1] <= thisPopMat$east + 2.5)
      closeN = (pts[2] > thisPopMat$north - 2.5) & (pts[2] <= thisPopMat$north + 2.5)
      closeI = closeE & closeN
      
      # there should be exactly 1 point we're closest to within this subarea
      if(sum(closeI > 1)) {
        stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      } else if(sum(closeI) == 0) {
        # warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
        
        if(testMode) {
          return(NA)
        }
        # this case should only happen at the edges, but just take closest point then
        dists = rdist(rbind(pts), cbind(thisPopMat$east, thisPopMat$north))
        # return(thisPopMatI[which.min(dists)])
        which.min(dists)
      } else {
        # return(thisPopMatI[which(closeI)])
        which(closeI)
      }
    }
    
    popMatIs = unlist(sapply(1:nrow(allPointsEN), getPopMatIs))
    urbVals = stratumPopMat$urban[popMatIs]
    
    npUrb = sum(urbVals, na.rm=TRUE)
    npRur = sum(!urbVals, na.rm=TRUE)
    
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
  
  popVals = terra::extract(pop, fineGridCoordsLLsp, method=extractMethod)
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
    
    urbanicityVals = terra::extract(urb, fineGridCoordsLL, method=extractMethod) # don't normalize urbanicity
    accessVals = terra::extract(accessNorm, fineGridCoordsLL, method=extractMethod)
    elevVals = terra::extract(elevNorm, fineGridCoordsLL, method=extractMethod)
    distVals = terra::extract(minDistRiverLakesNorm, fineGridCoordsLL, method=extractMethod)
    
    if(setMissingToAvg) {
      accessVals[is.na(accessVals)] = 0
      elevVals[is.na(elevVals)] = 0
      distVals[is.na(distVals)] = 0
    }
  } else {
    out = load("savedOutput/global/covariates.RData")
    
    inf = sessionInfo()
    if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
      urb@file@name = "~/git/jittering/savedOutput/global/urb.tif"
      access@file@name = "~/git/jittering/savedOutput/global/access.tif"
      elev@file@name = "~/git/jittering/savedOutput/global/elev.tif"
      dist@file@name = "~/git/jittering/savedOutput/global/dist.tif"
    }
    
    urbanicityVals = extract(urb, fineGridCoordsLL, method=extractMethod) # don't normalize urbanicity
    accessVals = extract(access, fineGridCoordsLL, method=extractMethod)
    elevVals = extract(elev, fineGridCoordsLL, method=extractMethod)
    distVals = extract(minDistRiverLakes, fineGridCoordsLL, method=extractMethod)
  }
  
  # remove NA covariate points and points with zero population, and, if there 
  # aren't enough points, restart this function with finer resolution
  naCovRowIs = is.na(urbanicityVals) | is.na(accessVals) | is.na(elevVals) | 
    is.na(distVals) | (finalPopVals == 0)
  if(testMode) {
    naCovRowIs[naCovRowIs] = FALSE
    # just return the entire set of results, don't continue to make the grid finer
  }
  
  if(((sum(finalPopMat$urb[!naCovRowIs], na.rm=TRUE) < minPointsUrb) && totalUrbPop > 0) || ((sum(!finalPopMat$urb[!naCovRowIs], na.rm=TRUE) < minPointsRur) && totalRurPop > 0)) {
    warning("NA covariates and zero pop points reduced number of urban and rural points to below minimum. Increasing resolution...")
    
    out = getFineIntPointsInfoMICS(stratumName=stratumName, kmresStart=kmres/2, minPointsUrb=minPointsUrb, minPointsRur=minPointsRur, 
                                   stratumMICSMapDat=stratumMICSMapDat, stratumMICSNameVar=stratumMICSNameVar, 
                                   subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                   poppsub=poppsub, 
                                   normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                   proj=proj, projArea=projArea, extractMethod=extractMethod)
  } else {
    out = cbind(finalPopMat, int=1, access=accessVals, elev=elevVals, distRiversLakes=distVals, urbanicity=urbanicityVals)[!naCovRowIs,]
    
    # recalibrate final population densities and include log population density as a covariate
    trueFinalUrbVals = finalPopMat$urban[!naCovRowIs]
    urbText = sapply(trueFinalUrbVals, function(x) {ifelse(x, "U", "R")})
    adm2TimesUR = paste(finalPopMat$subarea[!naCovRowIs], urbText, sep=",")
    trueFinalPopVals = SUMMER::calibrateByRegion(pointTotals=finalPopMat$pop[!naCovRowIs], pointRegions=adm2TimesUR, 
                                                 regions=regionsUR, regionTotals=regionTotals)
    
    load("savedOutput/global/popMeanSDCal.RData")
    popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
    popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
    out = cbind(out, normPop=(log1p(trueFinalPopVals)-popMeanCal)/popSDCal)
  }
  
  out
}

# second version of getFineIntPointsInfoMICSold using same covariates scheme as DHS integration points
# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
# Inputs:
# setMissingToAvg: if TRUE and normalized is TRUE, sets NA covariates to 0
getFineIntPointsInfoMICS = function(stratumName, kmresStart=2.5, minPointsUrb=20, minPointsRur=20, 
                                    stratumMICSMapDat=admFinalFull, stratumMICSNameVar="NAME_FINAL", 
                                    subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                    poppsub=poppsubNGAThresh, 
                                    normalized=TRUE, useThreshPopMat=TRUE, 
                                    proj=projNigeria, projArea=projNigeriaArea, 
                                    testMode=FALSE, extractMethod="bilinear", 
                                    setMissingToAvg=TRUE) {
  
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
    
    # get subareas associated with the points (considering only subareas within 
    # relevant MICS stratum)
    theseSubareas = getRegion2(allPointsLL, mapDat=thisSubareaMapDat, nameVar=subareaNameVar)
    theseSubareas = theseSubareas$regionNames
    naSubs = is.na(theseSubareas)
    
    allPointsEN = allPointsEN[!naSubs,]
    allPointsLL = allPointsLL[!naSubs,]
    theseSubareas = theseSubareas[!naSubs]
    
    uniqueSubareas = sort(unique(theseSubareas))
    
    getPopMatIs = function(i) {
      pts = allPointsEN[i,]
      # thisSubarea = theseSubareas[i]
      # thisPopMatI = which(stratumPopMat$subarea == thisSubarea)
      thisPopMat = stratumPopMat #[thisPopMatI,]
      # thisSubareaI = match(thisSubarea, uniqueSubareas)
      
      # dists = rdist(matrix(pts, ncol=2, cbind(thisPopMat$east, thisPopMat$north))
      
      # fullPopMat is 5km resolution, so must be within 2.5 km in easting and northing directions
      closeE = (pts[1] > thisPopMat$east - 2.5) & (pts[1] <= thisPopMat$east + 2.5)
      closeN = (pts[2] > thisPopMat$north - 2.5) & (pts[2] <= thisPopMat$north + 2.5)
      closeI = closeE & closeN
      
      # there should be exactly 1 point we're closest to within this subarea
      if(sum(closeI > 1)) {
        stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      } else if(sum(closeI) == 0) {
        # warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
        
        if(testMode) {
          return(NA)
        }
        # this case should only happen at the edges, but just take closest point then
        dists = rdist(rbind(pts), cbind(thisPopMat$east, thisPopMat$north))
        # return(thisPopMatI[which.min(dists)])
        which.min(dists)
      } else {
        # return(thisPopMatI[which(closeI)])
        which(closeI)
      }
    }
    
    popMatIs = unlist(sapply(1:nrow(allPointsEN), getPopMatIs))
    urbVals = stratumPopMat$urban[popMatIs]
    
    npUrb = sum(urbVals, na.rm=TRUE)
    npRur = sum(!urbVals, na.rm=TRUE)
    
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
  
  # NOTE: getDesignMat gets urbanicity from nearest popMat pixel, whereas urbanicity 
  # in this case should be set based on nearest pixel IN THE STRATUM. So we must correct it
  fineIntPtInfo = getDesignMat(fineGridCoordsLL, useThreshPopMat=useThreshPopMat, normalized=normalized, 
                               proj=proj, testMode=testMode, setMissingToAvg=setMissingToAvg)[,-1]
  fineIntPtInfo[,2] = urbVals
  
  # calibrate population to sum to correct totals within subareas x urban/rural 
  # (but only the population used for calculating aggregation weights)
  pointUrbText = sapply((fineIntPtInfo[,2]==1), function(x) {ifelse(x, "U", "R")})
  pointSubareaUR = paste(adm2Vals, pointUrbText, sep="")
  thisPoppsub = poppsub[adm2ToStratumMICS(poppsub$subarea) == stratumName,]
  regionUR = c(paste(thisPoppsub$subarea, "U", sep=""), paste(thisPoppsub$subarea, "R", sep=""))
  regionURpop = c(thisPoppsub$popUrb, thisPoppsub$popRur)
  calPop = SUMMER::calibrateByRegion(fineIntPtInfo[,1], pointSubareaUR, regionUR, regionURpop)
  
  fineIntPtInfo[,1] = calPop
  
  finalPopMat = data.frame(east=fineGridCoordsEN[,1], north=fineGridCoordsEN[,2], 
                           lon=fineGridCoordsLL[,1], lat=fineGridCoordsLL[,2], 
                           pop=fineIntPtInfo[,1], urban=(fineIntPtInfo[,2]==1), 
                           area=stratumPopMat$area[popMatIs], 
                           subarea=adm2Vals, strat=stratumPopMat$strat[popMatIs], 
                           popMatIs=which(popMat$strat==stratumName)[popMatIs])
  
  # remove NA covariate points and points with zero population, and, if there 
  # aren't enough points, restart this function with finer resolution
  naCovRowIs = apply(fineIntPtInfo, 1, anyNA) | (fineIntPtInfo[,1] == 0)
  if(testMode) {
    naCovRowIs[naCovRowIs] = FALSE
    # just return the entire set of results, don't continue to make the grid finer
  }
  
  # browser()
  if(((sum(finalPopMat$urban[!naCovRowIs], na.rm=TRUE) < minPointsUrb) && totalUrbPop > 0) || ((sum(!finalPopMat$urban[!naCovRowIs], na.rm=TRUE) < minPointsRur) && totalRurPop > 0)) {
    warning("NA covariates and zero pop points reduced number of urban and rural points to below minimum. Increasing resolution...")
    
    out = getFineIntPointsInfoMICS(stratumName=stratumName, kmresStart=kmres/2, minPointsUrb=minPointsUrb, minPointsRur=minPointsRur, 
                                   stratumMICSMapDat=stratumMICSMapDat, stratumMICSNameVar=stratumMICSNameVar, 
                                   subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                   poppsub=poppsub, 
                                   normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                   proj=proj, projArea=projArea, extractMethod=extractMethod)
  } else {
    out = cbind(finalPopMat, int=1, access=fineIntPtInfo[,3], elev=fineIntPtInfo[,4], distRiversLakes=fineIntPtInfo[,5], urbanicity=fineIntPtInfo[,6])[!naCovRowIs,]
    
    # get normalized population density 
    load("savedOutput/global/covariatesNorm.RData")
    popVals = extract(pop, fineGridCoordsLLsp, method="bilinear")
    
    load("savedOutput/global/popMeanSDCal.RData")
    popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
    popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
    out = cbind(out, normPop=(log1p(popVals[!naCovRowIs])-popMeanCal)/popSDCal)
  }
  
  out
}

# function for simulating fake MICS points using pop density
simMICSlocs = function(nsim=20, popMat=popMatNGAThresh, targetPopMat=popMatNGAThresh, 
                       poppsub=poppsubNGAThresh, saveFile="savedOutput/validation/simEdMICS.RData", 
                       seed=123) {
  set.seed(seed)
  
  # load in MICS data and integration points
  out = load("savedOutput/global/edMICS.RData")
  out=load("~/git/jittering/savedOutput/validation/edMICSval.RData")
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/global/intPtsMICS.RData")
  
  # we need the following format for the faux easpa to input to SUMMER and 
  # simulate the population. Only EAUrb, EARur, and EATotal matter for what we 
  # want, we just need to make sure the other values are reasonable so the 
  # simulation function doesn't break.
  #              area EAUrb EARur EATotal HHUrb  HHRur HHTotal popUrb  popRur popTotal
  # 1         Baringo   216  1754    1970 16322  94327  110649  71588  547153   618741
  clustpaMICS = aggregate(edMICSval$Stratum, by=list(edMICSval$Stratum, edMICSval$urban), FUN=length, drop=FALSE)
  clustpaMICS[is.na(clustpaMICS)] = 0
  rurInds = 1:41
  urbInds = 42:82
  clustpaMICS = cbind(clustpaMICS[urbInds, c(1, 3)], clustpaMICS[rurInds, 3])
  names(clustpaMICS) = c("area", "EAUrb", "EARur")
  clustpaMICS$EATotal = clustpaMICS$EAUrb + clustpaMICS$EARur
  clustpaMICS$HHUrb = clustpaMICS$EAUrb * 25 + 5
  clustpaMICS$HHRur = clustpaMICS$EARur * 25 + 5
  clustpaMICS$HHTotal = clustpaMICS$HHUrb + clustpaMICS$HHRur
  clustpaMICS$popUrb = clustpaMICS$HHUrb
  clustpaMICS$popRur = clustpaMICS$HHRur
  clustpaMICS$popTotal = clustpaMICS$popUrb + clustpaMICS$popRur
  
  # use MICS stratum as the "area" level in poppsub
  poppsub$origAdm1 = poppsub$area
  poppsub$area = adm2ToStratumMICS(poppsub$subarea)
  
  popMat$area = popMat$stratumMICS
  targetPopMat$area = targetPopMat$stratumMICS
  out = simPopCustom(logitRiskDraws=matrix(rep(0, nsim*nrow(targetPopMat)), ncol=nsim), 
                     sigmaEpsilonDraws=rep(0, nsim), easpa=clustpaMICS, 
                     popMat=popMat, targetPopMat=targetPopMat, 
                     stratifyByUrban=TRUE, validationPixelI=NULL, validationClusterI=NULL, 
                     clustersPerPixel=NULL, 
                     doFineScaleRisk=FALSE, doSmoothRisk=FALSE, 
                     doSmoothRiskLogisticApprox=TRUE, 
                     poppsub=poppsub, subareaLevel=TRUE, 
                     min1PerSubarea=FALSE, 
                     returnEAinfo=TRUE, epsc=NULL)
  
  eaDatList = out$eaDatList
  
  # extract the simulated locations from the simulated datasets and add in covariate 
  # information
  MICSlocs = lapply(eaDatList, function(x) {
    # extract coords from sim data
    temp = x[c("lon", "lat", "east", "north", "area", "subarea", "urban")]
    names(temp)[5] = "Stratum"
    
    # get covariate info
    lonLatCoords = temp[,1:2]
    xNorm = getDesignMatPopNorm(lonLatCoords)
    
    # concatenate all info together (except for urbanicity, which we're not including)
    temp = cbind(temp, xNorm[,-ncol(xNorm)])
    
    temp
    })
  browser()
  if(FALSE) {
    # test points to make sure they are correct
    dat = MICSlocs[[1]]
    
    tempPoints = dat[dat$Stratum == "Federal Capital Territory",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum %in% c("Kano North", "Kano Central", "Kano South"),]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum == "Kano Central",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum %in% c("Lagos East", "Lagos Central", "Lagos West"),]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum == "Lagos West",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
  }
  
  # give the MICS data the simulated locations. In order to merge we need to sort the 
  # tables to be in the same order, then combine, then return the ordering to how it was
  edMICSurb = sapply(edMICSval$urban, function(x) {ifelse(x, "U", "R")})
  edMICSstratUrb = sapply(1:length(edMICSval$Stratum), function(i) {paste(edMICSval$Stratum[i], edMICSurb[i], sep=",")})
  orderEdMICS = order(edMICSstratUrb)
  revOrderEdMICS = match(1:length(orderEdMICS), orderEdMICS)
  
  edMICSordered = edMICSval[orderEdMICS,]
  
  # combine simulated MICS locations with edMICS data. Also add on covariates at 
  # the locations
  edMICSlist = lapply(MICSlocs, function(x) {
    thisUrb = sapply(x$urban, function(x) {ifelse(x, "U", "R")})
    thisStratUrb = sapply(1:length(x$Stratum), function(i) {paste(x$Stratum[i], thisUrb[i], sep=",")})
    thisOrder = order(thisStratUrb)
    
    xOrdered = x[thisOrder,]
    
    thisEdMICSordered = cbind(edMICSordered, xOrdered[,-which(colnames(xOrdered) %in% c("Stratum", "urban"))])
    thisEdMICSordered[revOrderEdMICS,]
  })
  
  if(FALSE) {
    # TODO: Figure out why the wrong number of EAs are simulated in each stratum. 
    #       clustpaMICS appears to be correct, so problem must be in the 
    #       simulation code
    cbind(table(MICSlocs[[1]]$Stratum), table(edMICSval$Stratum), clustpaMICS$area, clustpaMICS$EATotal)
    # [,1] [,2] [,3]                        [,4]
    # ...
    # Ekiti                     "55" "51" "Ekiti"                     "51"
    
    # test points to make sure they are correct
    dat = data.frame(lon=edMICSlist[[1]]$lon, lat=edMICSlist[[1]]$lat, Stratum=edMICSlist[[1]]$Stratum)
    
    tempPoints = dat[dat$Stratum == "Federal Capital Territory",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum %in% c("Kano North", "Kano Central", "Kano South"),]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum == "Kano Central",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum %in% c("Lagos East", "Lagos Central", "Lagos West"),]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
    
    tempPoints = dat[dat$Stratum == "Lagos West",]
    plotMapDat(adm1)
    points(tempPoints[,1:2], pch=19, col="blue", cex=.2)
  }
  
  # get covariates 
  
  simEdMICS = edMICSlist
  save(simEdMICS, file="savedOutput/validation/simEdMICS.RData")
  
  invisible(simEdMICS)
}

# by default, finds resolution within 1% of fine scale average
# tol: tolerance in percent
getBestResMICS = function(startRes=25, tolMean=.01, tolMax=.05) {
  
  lastRes = startRes
  thisRes = startRes
  lastAbsMeanPctError = NULL
  finished = FALSE
  allRes = startRes
  allAbsMeanPctErrorUrb = c()
  allAbsMeanPctErrorRur = c()
  allAbsMeanPctError = c()
  allAbsMaxPctErrorUrb = c()
  allAbsMaxPctErrorRur = c()
  allAbsMaxPctError = c()
  allAbsMeanPctErrorSDUrb = c()
  allAbsMeanPctErrorSDRur = c()
  allAbsMeanPctErrorSD = c()
  allAbsMaxPctErrorSDUrb = c()
  allAbsMaxPctErrorSDRur = c()
  allAbsMaxPctErrorSD = c()
  while(!finished) {
    thisOutfile = paste0("savedOutput/global/intPtsMICS_", thisRes, ".RData")
    if(!file.exists(thisOutfile)) {
      micsPts = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
                                             numPtsUrb = thisRes, numPtsRur=thisRes, 
                                             outFile=thisOutfile)
    } else {
      out = load(thisOutfile)
      micsPts = intPtsMICS
    }
    
    # get diagnostics
    absMeanPctErrorUrb = micsPts$absMeanPctErrorUrb
    absMeanPctErrorRur = micsPts$absMeanPctErrorRur
    absMeanPctError = micsPts$absMeanPctError
    absMaxPctErrorUrb = micsPts$absMaxPctErrorUrb
    absMaxPctErrorRur = micsPts$absMaxPctErrorRur
    absMaxPctError = micsPts$absMaxPctError
    absMeanPctErrorSDUrb = micsPts$absMeanPctErrorSDUrb
    absMeanPctErrorSDRur = micsPts$absMeanPctErrorSDRur
    absMeanPctErrorSD = micsPts$absMeanPctErrorSD
    absMaxPctErrorSDUrb = micsPts$absMaxPctErrorSDUrb
    absMaxPctErrorSDRur = micsPts$absMaxPctErrorSDRur
    absMaxPctErrorSD = micsPts$absMaxPctErrorSD
    
    allAbsMeanPctErrorUrb = c(allAbsMeanPctErrorUrb, absMeanPctErrorUrb)
    allAbsMeanPctErrorRur = c(allAbsMeanPctErrorRur, absMeanPctErrorRur)
    allAbsMeanPctError = c(allAbsMeanPctError, absMeanPctError)
    allAbsMaxPctErrorUrb = c(allAbsMaxPctErrorUrb, absMaxPctErrorUrb)
    allAbsMaxPctErrorRur = c(allAbsMaxPctErrorRur, absMaxPctErrorRur)
    allAbsMaxPctError = c(allAbsMaxPctError, absMaxPctError)
    allAbsMeanPctErrorSDUrb = c(allAbsMeanPctErrorSDUrb, absMeanPctErrorSDUrb)
    allAbsMeanPctErrorSDRur = c(allAbsMeanPctErrorSDRur, absMeanPctErrorSDRur)
    allAbsMeanPctErrorSD = c(allAbsMeanPctErrorSD, absMeanPctErrorSD)
    allAbsMaxPctErrorSDUrb = c(allAbsMaxPctErrorSDUrb, absMaxPctErrorSDUrb)
    allAbsMaxPctErrorSDRur = c(allAbsMaxPctErrorSDRur, absMaxPctErrorSDRur)
    allAbsMaxPctErrorSD = c(allAbsMaxPctErrorSD, absMaxPctErrorSD)
    
    fineIntPtAvgsUrb = micsPts$fineIntPtAvgsUrb
    fineIntPtAvgsRur = micsPts$fineIntPtAvgsRur
    intPtAvgsUrb = micsPts$intPtAvgsUrb
    intPtAvgsRur = micsPts$intPtAvgsRur
    error = micsPts$error
    
    meanErrors = colMeans(abs(error))
    maxErrors = apply(abs(error), 2, max)
    
    browser()
    if((max(meanErrors) < tolMean) && (max(maxErrors) < tolMax)) {
      finished = TRUE
    } else {
      allRes = c(allRes, thisRes)
      
      # Reimann integrals converge at rate 1/n
      Cmean = max(meanErrors)*thisRes
      Cmax = max(maxErrors)*thisRes
      
      # C/n = tolMean
      meanRes = ceiling(Cmean/tolMean)
      maxRes = ceiling(Cmax/tolMax)
      
      lasRes = thisRes
      thisRes = max(c(meanRes, maxRes))
      
      if(thisRes == lastRes) {
        stop("why???")
      }
      
      print(paste0("meanError: ", max(meanErrors), ", maxError: ", max(maxErrors)))
      print(paste0("newRes: ", thisRes, ", meanRes: ", meanRes, ", maxRes: ", maxRes))
    }
  }
  
  browser()
  
  list(finalMICSpts=micsPts, allRes=allRes, 
       allAbsMeanPctErrorUrb=allAbsMeanPctErrorUrb, 
       allAbsMeanPctErrorRur=allAbsMeanPctErrorRur, 
       allAbsMeanPctError=allAbsMeanPctError, 
       allAbsMaxPctErrorUrb=allAbsMaxPctErrorUrb, 
       allAbsMaxPctErrorRur=allAbsMaxPctErrorRur, 
       allAbsMaxPctError=allAbsMaxPctError, 
       allAbsMeanPctErrorSDUrb=allAbsMeanPctErrorSDUrb, 
       allAbsMeanPctErrorSDRur=allAbsMeanPctErrorSDRur, 
       allAbsMeanPctErrorSD=allAbsMeanPctErrorSD, 
       allAbsMaxPctErrorSDUrb=allAbsMaxPctErrorSDUrb, 
       allAbsMaxPctErrorSDRur=allAbsMaxPctErrorSDRur, 
       allAbsMaxPctErrorSD=allAbsMaxPctErrorSD)
}




