# make figures for illustrating how the integration points 
# are constructed and weighted
# 
# notation:
# group 1: the central point (the jittered cluster location)
# group 2: the first ring of points, over the smaller radius
# group 3: the second ring of points, over the larger radius
# ri: radius containing group i
# mi: number of integration points in group i
# Ai: area of section associated with an individual integration point in group i
makeIntegrationIllustration = function(numPointsUrban=11, MInnerUrban=3, MOuterUrban=0, 
                                       numPointsRural=16, MInnerRural=3, MOuterRural=1) {
  require(latex2exp)
  require(xtable)
  
  outUrban = getIntegrationPoints(numPoints=numPointsUrban, 
                                  MInner=MInnerUrban, MOuter=MOuterUrban)
  outRural = getIntegrationPoints(urban=FALSE, numPoints=numPointsRural, 
                                  MInner=MInnerRural, MOuter=MOuterRural)
  
  # quantities that will be the same for both urban and rural
  as = outUrban$as
  ms = outUrban$ms
  rs = outUrban$rs
  pts = outUrban$pts
  
  # ranges of variables that will be different for urban and rural
  densityFunUrban = outUrban$densityFun
  densityFunRural = outRural$densityFun
  pRange = range(c(densityFunUrban(.1), outUrban$ps[-1], 
                   densityFunRural(.1), outRural$ps[-1]))
  wRange = range(c(outUrban$ws, outRural$ws))
  
  # urban specific values
  ms = outUrban$ms
  ws = outUrban$ws
  rs = outUrban$rs
  pts = outUrban$pts
  as = outUrban$as
  
  colScale = makeGreenSequentialColors(64)
  
  pdf(file=paste0("figures/integrationWeightsUrban_numPoints", numPointsUrban, ".pdf"), 
      width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  vals = ws
  cols = getColorsFromScale(vals, valRange=c(0, max(wRange)), colScale)
  plot(c(-10, 10), c(-10, 10), type="n", asp=1, 
       xlab="East Displacement (km)", ylab="North Displacement (km)", 
       main=paste0("Integration Weights (Urban Points)"))
  for(i in 1:length(rs)) {
    draw.circle(0, 0, rs[i])
  }
  
  for(i in 2:length(rs)) {
    thisas = as[[i]] + pi/ms[i]
    for(j in 1:ms[i]) {
      thisA = thisas[j]
      thisx1 = rs[i-1]*cos(thisA)
      thisy1 = rs[i-1]*sin(thisA)
      thisx2 = rs[i]*cos(thisA)
      thisy2 = rs[i]*sin(thisA)
      lines(c(thisx1, thisx2), c(thisy1, thisy2))
    }
  }
  
  for(i in 1:length(rs)) {
    points(pts[[i]], pch=21, bg=cols[i], cex=.9, col="black")
  }
  
  # image.plot(zlim=c(0, max(pRange)), nlevel=length(colScale), legend.only=TRUE, horizontal=FALSE,
  #            col=colScale, add = TRUE)
  
  dev.off()
  
  pdf(file=paste0("figures/integrationWeightsUrban_numPoints", numPointsUrban, "DiffScale.pdf"), 
      width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  vals = ws
  cols = getColorsFromScale(vals, valRange=c(0, max(wRange)), colScale)
  plot(c(-2, 2), c(-2, 2), type="n", asp=1, 
       xlab="East Displacement (km)", ylab="North Displacement (km)", 
       main=paste0("Integration Weights (Urban Points)"))
  
  for(i in 1:length(rs)) {
    draw.circle(0, 0, rs[i])
  }
  
  for(i in 2:length(rs)) {
    thisas = as[[i]] + pi/ms[i]
    for(j in 1:ms[i]) {
      thisA = thisas[j]
      thisx1 = rs[i-1]*cos(thisA)
      thisy1 = rs[i-1]*sin(thisA)
      thisx2 = rs[i]*cos(thisA)
      thisy2 = rs[i]*sin(thisA)
      lines(c(thisx1, thisx2), c(thisy1, thisy2))
    }
  }
  
  for(i in 1:length(rs)) {
    points(pts[[i]], pch=21, bg=cols[i], cex=.9, col="black")
  }
  
  # image.plot(zlim=c(0, max(pRange)), nlevel=length(colScale), legend.only=TRUE, horizontal=FALSE,
  #            col=colScale, add = TRUE)
  
  dev.off()
  
  # rural specific values
  ms = outRural$ms
  ws = outRural$ws
  rs = outRural$rs
  as = outRural$as
  pts = outRural$pts
  densityFunRural = outRural$densityFun
  
  pdf(file=paste0("figures/integrationWeightsRural_numPoints", numPointsRural, ".pdf"), 
      width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  vals = ws
  cols = getColorsFromScale(vals, valRange=c(0, max(wRange)), colScale)
  plot(c(-10, 10), c(-10, 10), type="n", asp=1, 
       xlab="East Displacement (km)", ylab="North Displacement (km)", 
       main=paste0("Integration Weights (Rural Points)"))
  for(i in 1:length(rs)) {
    draw.circle(0, 0, rs[i])
  }
  
  for(i in 2:length(rs)) {
    thisas = as[[i]] + pi/ms[i]
    for(j in 1:ms[i]) {
      thisA = thisas[j]
      thisx1 = rs[i-1]*cos(thisA)
      thisy1 = rs[i-1]*sin(thisA)
      thisx2 = rs[i]*cos(thisA)
      thisy2 = rs[i]*sin(thisA)
      lines(c(thisx1, thisx2), c(thisy1, thisy2))
    }
  }
  
  for(i in 1:length(rs)) {
    points(pts[[i]], pch=21, bg=cols[i], cex=2, col="black")
  }
  
  image.plot(zlim=c(0, max(wRange)), nlevel=length(colScale), legend.only=TRUE, horizontal=FALSE,
             col=colScale, add = TRUE)
  
  dev.off()
  
  ##### Make simplified illustration versus radius
  pdf(file=paste0("figures/jitteringPriorRadial_numPointsUrban", numPointsUrban, "_numPointsRural", numPointsRural, ".pdf"), 
      width=5, height=5)
  
  yRange = range(c(pRange, wRange))
  
  plot(c(0, 10), yRange, type="n", 
       xlab=TeX("Displacement (km)"), ylab="", 
       main="Jittering Density and Integration Weights", 
       log="y")
  
  # urban specific values
  ms = outUrban$ms
  ps = outUrban$ps
  ws = outUrban$ws
  rs = outUrban$rs
  ptRs = outUrban$ptRs
  
  # lines(c(0, r2), c(p2, p2), col="blue")
  # lines(c(r2, r3), c(p3, p3), col="blue")
  rseq = seq(0, 2, l=100)
  ps = densityFunUrban(rseq)
  lines(rseq, ps, col="blue")
  
  points(ptRs, ws, col="blue", pch="x")
  
  # rural specific values
  ms = outRural$ms
  ps = outRural$ps
  ws = outRural$ws
  rs = outRural$rs
  ptRs = outRural$ptRs
  
  # lines(c(0, r2), c(p2, p2), col="green", lty=2)
  # lines(c(r2, r3), c(p3, p3), col="green", lty=2)
  rseq = seq(0, 5, l=100)[-100]
  ps = densityFunRural(rseq)
  lines(rseq, ps, col="green", lty=2)
  rseq = seq(5, 10, l=100)
  ps = densityFunRural(rseq)
  lines(rseq, ps, col="green", lty=2)
  
  points(ptRs, ws, col="green", pch="+")
  
  legend("topright", 
         c("Urban density", "Rural density", 
           "Urban weights", "Rural weights"), 
         lty=c(1, 2, NA, NA), 
         pch=c(NA, NA, "x", "+"), 
         col=rep(c("blue", "green"), 2))
  
  dev.off()
  
  # make prior and integration weight table
  tab = cbind(c(outUrban$ptRs, outRural$ptRs), 
              c(outUrban$ms, outRural$ms), 
              c(outUrban$ws, outRural$ws))
  # tab = rbind(c(0, outUrban$m1, outUrban$w1), 
  #             c(mean(c(outUrban$r1, outUrban$r2)), outUrban$m2, outUrban$w2), 
  #             c(mean(c(outUrban$r2, outUrban$r3)), outUrban$m3, outUrban$w3), 
  #             c(0, outRural$m1, outRural$w1), 
  #             c(mean(c(outRural$r1, outRural$r2)), outRural$m2, outRural$w2), 
  #             c(mean(c(outRural$r2, outRural$r3)), outRural$m3, outRural$w3))
  
  row.names(tab) = c("Urban  Ring 1", 
                     "Urban  Ring 2", 
                     "Urban  Ring 3", 
                     "Rural  Ring 1", 
                     "Rural  Ring 2", 
                     "Rural  Ring 3", 
                     "Rural  Ring 4")
  colnames(tab) = c("Displacement (km)", "Number of Points", 
                    "Integration Weights")
  # format(tab, digits=3)
  print(xtable(tab, digits=c(0, 2, 0, 3), 
               display=c("s", "f", "d", "f")))
  
  # print(xtable(t(tab), digits=3))
}

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
getIntegrationPoints = function(urban=TRUE, numPoints=ifelse(urban, 11, 16), 
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
makeAllIntegrationPoints = function(coords, urbanVals, 
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

# Long story short, use this function when using a log, logit, or some other fancy 
# transformation of a color scale.
# given data to be assigned colors from scale (vals), the range of the values, the center 
# of the value range (only for for diverging scale), the color scale as a vector 
# of colors (cols), the scale of the color scale (e.g. identity, log, logit), 
# and whether or not to force the colors into valRange, this function returns the 
# colors associated with each value in vals.
# Example:
# # construct data:
# testX = exp(rnorm(100))
# # construct centered color scale on log scale
# test = centerColorScale(64, testX, center=1, colScale=makeRedBlueDivergingColors, scaleFun=log)
# # get the colors associated with each value
# testCols = getColorsFromScale(testX, center=1, cols=test, scaleFun=log)
# # plot the result
# plot(testX, col=testCols, pch=19)
getColorsFromScale = function(vals, valRange=NULL, center=NULL, cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  
  if(!is.null(center)) {
    center = scaleFun(center)
    n = length(cols)
    
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = cols
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
    
    cols = tempColors
  }
  
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}

makeGreenSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev, n_interp=n)
}