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
getIntegrationPointsMICS = function(urban=TRUE, integrationResolution=25, 
                                   integrationPointType=c("mean", "midpoint"), 
                                   verbose=TRUE) {
}



