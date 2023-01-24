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
makeIntegrationIllustrationDHS = function(numPointsUrban=11, MInnerUrban=3, MOuterUrban=0, 
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

makeIntegrationIllustrationMICS = function(stratum="Federal Capital Territory", 
                                           kmresFineStart=1, 
                                           numPtsUrb=3, numPtsRur=3, 
                                           stratumMICSMapDatFull=admFinalFull, 
                                           stratumMICSMapDat=admFinal, 
                                           stratumMICSNameVar="NAME_FINAL", 
                                           subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                           poppsub=poppsubNGAThresh, 
                                           normalized=TRUE, useThreshPopMat=TRUE, 
                                           proj=projNigeria, projArea=projNigeriaArea, 
                                           spatialAsCovariate=FALSE, 
                                           lambda=NULL, domainDiameter=NULL, 
                                           fileNameRoot="MICSintPts_", loadSavedIntPoints=FALSE) {
  
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
  i = which(strataMICS == stratum)
  
  thisStrat = strataMICS[i]
  if(thisStrat == "Lake Chad") {
    # no people in Lake Chad
    stop("No integration points generated for Lake Chad")
  }
  
  print(paste0("Constructing integration points for MICS stratum ", thisStrat, 
               " (", i, "/", length(strataMICS), ")"))
  
  # output (thisIntPoints) is a list with pts, weights, and area
  if(loadSavedIntPoints && file.exists(paste0(fileNameRoot, "_i", i, ".RData"))) {
    out = load(paste0(fileNameRoot, "_i", i, ".RData"))
  } else {
    thisIntPoints= getIntegrationPointsMICS(strat=thisStrat, kmresFineStart=kmresFineStart, 
                                            numPtsUrb=numPtsUrb, numPtsRur=numPtsRur, 
                                            stratumMICSMapDat=stratumMICSMapDatFull, 
                                            stratumMICSNameVar=stratumMICSNameVar, 
                                            subareaMapDat=subareaMapDat, 
                                            subareaNameVar=subareaNameVar, 
                                            poppsub=poppsub, 
                                            normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                            proj=proj, projArea=projArea, 
                                            spatialAsCovariate=spatialAsCovariate, 
                                            lambda=lambda, domainDiameter=domainDiameter, 
                                            returnFineGrid=TRUE)
  }
  
  # plot the fine grid of integration points, circling the selected points
  intPoints = thisIntPoints$intPoints
  fineGrid = thisIntPoints$fineGrid
  medoidIUrb = thisIntPoints$medoidIUrb
  medoidIRur = thisIntPoints$medoidIRur
  sortIUrb = sort(unique(medoidIUrb))
  sortIRur = sort(unique(medoidIRur))
  
  # obtain the SpatialPolygon associated with this MICS stratum
  stratPolygonI = which(stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum)
  stratumMICSPolygons = SpatialPolygons(stratumMICSMapDat@polygons)
  stratumMICSSpatialPolygon = stratumMICSPolygons[stratPolygonI]
  stratumMICSPolygons@proj4string = stratumMICSMapDat@proj4string
  stratumMICSSpatialPolygon@proj4string = stratumMICSMapDat@proj4string
  
  bbox = stratumMICSSpatialPolygon@bbox
  
  # plot colored by cluster
  pdf("figures/integration/integrationMICSbyCluster.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  nClustersUrb = nrow(intPoints$ptsUrb)
  nClustersRur = nrow(intPoints$ptsRur)
  cols = rainbow(nClustersUrb + nClustersRur)
  cols[1:nClustersUrb] = sample(cols[1:nClustersUrb], nClustersUrb)
  cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)] = sample(cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], nClustersRur)
  urbCols = cols[match(medoidIUrb, sort(unique(medoidIUrb)))]
  rurCols = cols[nClustersUrb + match(medoidIRur, sort(unique(medoidIRur)))]
  
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Integration points and associated locations"), leaveRoomForLegend = TRUE, 
             legend.mar=3, xlab="")
  plotWithColor(fineGrid$lon[fineGrid$urban], fineGrid$lat[fineGrid$urban], 
                match(medoidIUrb, sort(unique(medoidIUrb))), pch=15, cex=.3, 
                colScale=cols[1:nClustersUrb], add=TRUE, addColorBar=TRUE, 
                legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                                legend.cex=1, smallplot=c(.9,.93,.55,.9), 
                                legend.lab="urban", legend.line=0.5))
  plotWithColor(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], 
                nClustersUrb + match(medoidIRur, sort(unique(medoidIRur))), 
                pch=15, cex=.3, colScale=cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], 
                add=TRUE, addColorBar=TRUE, 
                legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                                legend.cex=1, smallplot=c(.9,.93,.15,.5), 
                                legend.lab="rural", legend.line=0.5))
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col=cols[1:nClustersUrb], pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col=cols[(nClustersUrb+1):length(cols)], pch=1)
  points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, pch=0)
  points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, pch=1)
  dev.off()
  
  browser()
  
  pdf("figures/integration/integrationMICSbyClusterQuilt.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  nClustersUrb = nrow(intPoints$ptsUrb)
  nClustersRur = nrow(intPoints$ptsRur)
  cols = rainbow(nClustersUrb + nClustersRur)
  cols[1:nClustersUrb] = sample(cols[1:nClustersUrb], nClustersUrb)
  cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)] = sample(cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], nClustersRur)
  urbCols = cols[match(medoidIUrb, sort(unique(medoidIUrb)))]
  rurCols = cols[nClustersUrb + match(medoidIRur, sort(unique(medoidIRur)))]
  
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Integration points and associated locations"), leaveRoomForLegend = TRUE, 
             legend.mar=3, xlab="")
  nx = 95
  myQuiltPlot(x=fineGrid$lon[fineGrid$urban], y=fineGrid$lat[fineGrid$urban], 
              z=match(medoidIUrb, sort(unique(medoidIUrb))), nx=nx, asp=1, 
              xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
              colScale=cols[1:nClustersUrb], add=TRUE, addColorBar=TRUE, 
              legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                              legend.cex=1, smallplot=c(.9,.93,.55,.9), 
                              legend.lab="urban", legend.line=0.5))
  myQuiltPlot(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], 
              nClustersUrb + match(medoidIRur, sort(unique(medoidIRur))), 
              colScale=cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], 
              xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
              add=TRUE, addColorBar=TRUE, nx=nx, asp=1, 
              legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                              legend.cex=1, smallplot=c(.9,.93,.15,.5), 
                              legend.lab="rural", legend.line=0.5))
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col=cols[1:nClustersUrb], pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col=cols[(nClustersUrb+1):length(cols)], pch=1)
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, pch=1)
  outUrb = myQuiltPlot(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, 
                       1:length(intPoints$ptsUrb$lon), 
                       colScale=rep("black", length(intPoints$ptsUrb$lon)), 
                       xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
                       add=TRUE, addColorBar=FALSE, nx=nx, asp=1, plot=FALSE)
  lonUrb = rev(outUrb$x[outUrb$ind[,1]])[-(1:4)]
  latUrb = rev(outUrb$y[outUrb$ind[,2]])[-(1:4)]
  outRur = myQuiltPlot(intPoints$ptsRur$lon, intPoints$ptsRur$lat, 
                       1:length(intPoints$ptsRur$lon), 
                       colScale=rep("black", length(intPoints$ptsRur$lon)), 
                       xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
                       add=TRUE, addColorBar=FALSE, nx=nx, asp=1, plot=FALSE)
  lonRur = rev(outRur$x[outRur$ind[,1]])[-(1:4)]
  latRur = rev(outRur$y[outRur$ind[,2]])[-(1:4)]
  points(lonUrb, latUrb, cex=.5, pch=0)
  points(lonRur, latRur, cex=.5, pch=1)
  
  dev.off()
  
  # plot covariates
  popCols = makeBlueSequentialColors(64)
  covCols = makePurpleYellowSequentialColors(64)
  
  # urban
  pdf("figures/integration/integrationMICSbyUrban.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Urban/rural classification"), leaveRoomForLegend=TRUE, 
             legend.mar=3, xlab="")
  points(fineGrid$lon[fineGrid$urban], fineGrid$lat[fineGrid$urban], pch=".", col="blue")
  points(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], pch=".", col="green")
  points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col="darkblue", pch=0)
  points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col="darkgreen", pch=1)
  image.plot(x=c(1, 2, 1, 2), y=c(1, 1, 2, 2), z=1:4, breaks=c(0, 2), 
             legend.only=TRUE, add=TRUE, col="blue", 
             axis.args=list(labels=FALSE, tick=FALSE), 
             legend.cex=1, smallplot=c(.9,.93,.55,.9), 
             legend.lab="urban", legend.line=0.5)
  image.plot(x=c(1, 2, 1, 2), y=c(1, 1, 2, 2), z=1:4, breaks=c(0, 2), 
             legend.only=TRUE, add=TRUE, col="green", 
             axis.args=list(labels=FALSE, tick=FALSE), 
             legend.cex=1, smallplot=c(.9,.93,.15,.5), 
             legend.lab="rural", legend.line=0.5)
  dev.off()
  
  # log1pPop
  pdf("figures/integration/integrationMICSbyLog1pPop.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Log1p population"), legend.mar=3, xlab="", ylab="")
  zlim = range(fineGrid$log1pPop)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$log1pPop, colScale=popCols, 
                add=TRUE, pch=16, cex=.3, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  plotWithColor(intPoints$pts$lon, intPoints$pts$lat, intPoints$pts$log1pPop, colScale=popCols, 
                add=TRUE, pch=1, cex=.5, zlim=zlim, addColorBar=FALSE)
  dev.off()
  
  # access
  pdf("figures/integration/integrationMICSbyAccess.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Access to healthcare"), legend.mar=3)
  zlim = range(fineGrid$access)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$access, colScale=covCols, 
                add=TRUE, pch=16, cex=.3, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  plotWithColor(intPoints$pts$lon, intPoints$pts$lat, intPoints$pts$access, colScale=covCols, 
                add=TRUE, pch=1, cex=.5, zlim=zlim, addColorBar=FALSE)
  dev.off()
  
  # elev
  pdf("figures/integration/integrationMICSbyElev.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Elevation"), legend.mar=3, ylab="")
  zlim = range(fineGrid$elev)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$elev, colScale=covCols, 
                add=TRUE, pch=16, cex=.3, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  plotWithColor(intPoints$pts$lon, intPoints$pts$lat, intPoints$pts$elev, colScale=covCols, 
                add=TRUE, pch=1, cex=.5, zlim=zlim, addColorBar=FALSE)
  dev.off()
  
  # distRiversLakes
  pdf("figures/integration/integrationMICSbyDist.pdf", width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Distance to rivers and lakes"), legend.mar=3, xlab="", ylab="")
  zlim = range(fineGrid$distRiversLakes)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$distRiversLakes, colScale=covCols, 
                add=TRUE, pch=16, cex=.3, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  plotWithColor(intPoints$pts$lon, intPoints$pts$lat, intPoints$pts$distRiversLakes, colScale=covCols, 
                add=TRUE, pch=1, cex=.5, zlim=zlim, addColorBar=FALSE)
  dev.off()
  
  browser()
}

# plot rasters and Nigeria DHS and MICS data
plotDatasets = function(kmres=5) {
  browser()
  
  # plot rasters ----
  
  # project area to easting/northing in km, and make a fine easting/northing
  # grid of points
  projArea = projNigeriaArea(adm0)
  xCoords = seq(projArea@bbox[1,1], projArea@bbox[1,2], by=kmres)
  yCoords = seq(projArea@bbox[2,1], projArea@bbox[2,2], by=kmres)
  ENCoords = make.surface.grid(list(x=xCoords, y=yCoords))
  
  # convert the grid of points back to longitude/latitude and get covariate
  # values
  LLCoords = projNigeria(ENCoords, inverse=TRUE)
  X = getDesignMat(LLCoords, normalized=FALSE)
  
  
  quilt.plot(1:2, 1:2, rep(NA, 2), nx=2, ny=2, add.legend = FALSE, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
             ylab="Latitude", main="Population density")
  plot(pop, add=TRUE)
  plotMapDat(mapDat=adm0, border=rgb(1,1,1,0))
  plotMapDat(mapDat=adm1)
  plotMapDat(mapDat=adm2, lwd=.2, border=rgb(.5, .5, .5, .5), add=TRUE)
  
  
  
  # plot DHS data ----
  
  # plot MICS data ----
  
}







