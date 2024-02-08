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

makeAllIntegrationIllustrationsMICS = function(stratum="Federal Capital Territory", 
                                               kmresFineStart=1, 
                                               Kvals=c(3, 6, 12, 24), 
                                               stratumMICSMapDatFull=admFinalFull, 
                                               stratumMICSMapDat=admFinal, 
                                               stratumMICSNameVar="NAME_FINAL", 
                                               subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                               poppsub=poppsubNGAThresh, 
                                               normalized=TRUE, useThreshPopMat=TRUE, 
                                               proj=projNigeria, projArea=projNigeriaArea, 
                                               spatialAsCovariate=FALSE, 
                                               lambda=NULL, domainDiameter=NULL, 
                                               fileNameRoot="MICSintPts_", loadSavedIntPoints=FALSE, 
                                               plotLongitude=c(12, 24), plotLatitude=c(3, 12)) {
  
  for(k in Kvals) {
    print(paste0("Generating illustration for K=", k))
    subDir = paste0("K", k, "/")
    thisPlotLon = k %in% plotLongitude
    thisPlotLat = k %in% plotLatitude
    makeIntegrationIllustrationMICS(stratum=stratum, 
                                    kmresFineStart=kmresFineStart, 
                                    numPtsUrb=k, numPtsRur=k, 
                                    stratumMICSMapDatFull=stratumMICSMapDatFull, 
                                    stratumMICSMapDat=stratumMICSMapDat, 
                                    stratumMICSNameVar=stratumMICSNameVar, 
                                    subareaMapDat=subareaMapDat, subareaNameVar=subareaNameVar, 
                                    poppsub=poppsub, 
                                    normalized=normalized, useThreshPopMat=useThreshPopMat, 
                                    proj=proj, projArea=projArea, 
                                    spatialAsCovariate=spatialAsCovariate, 
                                    lambda=lambda, domainDiameter=domainDiameter, 
                                    fileNameRoot=fileNameRoot, 
                                    loadSavedIntPoints=loadSavedIntPoints, 
                                    subDir=subDir, plotLongitude=thisPlotLon, 
                                    plotLatitude=thisPlotLat)
  }
}

makeIntegrationIllustrationMICS = function(stratum="Federal Capital Territory", 
                                           kmresFineStart=1, 
                                           numPtsUrb=3, numPtsRur=numPtsUrb, 
                                           stratumMICSMapDatFull=admFinalFull, 
                                           stratumMICSMapDat=admFinal, 
                                           stratumMICSNameVar="NAME_FINAL", 
                                           subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                                           poppsub=poppsubNGAThresh, 
                                           normalized=TRUE, useThreshPopMat=TRUE, 
                                           proj=projNigeria, projArea=projNigeriaArea, 
                                           spatialAsCovariate=FALSE, 
                                           lambda=NULL, domainDiameter=NULL, 
                                           fileNameRoot="MICSintPts_", loadSavedIntPoints=FALSE, 
                                           subDir="", plotLongitude=TRUE, plotLatitude=TRUE) {
  
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
  
  if(FALSE) {
    thisIntPointsTest = getIntegrationPointsMICS(strat=thisStrat, kmresFineStart=kmresFineStart, 
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
                                            returnFineGrid=TRUE, testMode=TRUE)
    
    intPointsTest = thisIntPointsTest$intPoints
    fineGridTest = thisIntPointsTest$fineGrid
    medoidIUrbTest = thisIntPointsTest$medoidIUrb
    medoidIRurTest = thisIntPointsTest$medoidIRur
    sortIUrbTest = sort(unique(medoidIUrbTest))
    sortIRurTest = sort(unique(medoidIRurTest))
  }
  
  # obtain the SpatialPolygon associated with this MICS stratum
  stratPolygonI = which(stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum)
  stratumMICSPolygons = SpatialPolygons(stratumMICSMapDat@polygons)
  stratumMICSSpatialPolygon = stratumMICSPolygons[stratPolygonI]
  stratumMICSPolygons@proj4string = stratumMICSMapDat@proj4string
  stratumMICSSpatialPolygon@proj4string = stratumMICSMapDat@proj4string
  
  bbox = stratumMICSSpatialPolygon@bbox
  
  # plot colored by cluster
  require(RColorBrewer)
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyCluster.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  nClustersUrb = nrow(intPoints$ptsUrb)
  nClustersRur = nrow(intPoints$ptsRur)
  cols = rainbow(nClustersUrb + nClustersRur)
  # cols = brewer.pal(n = nClustersUrb + nClustersRur, name = "Dark2")
  # cols[1:nClustersUrb] = sample(cols[1:nClustersUrb], nClustersUrb)
  # cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)] = sample(cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], nClustersRur)
  tempCols = cols
  cols[1:nClustersUrb] = tempCols[seq(1, length(cols), by=2)]
  cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)] = tempCols[seq(2, length(cols), by=2)]
  urbCols = cols[match(medoidIUrb, sort(unique(medoidIUrb)))]
  rurCols = cols[nClustersUrb + match(medoidIRur, sort(unique(medoidIRur)))]
  
  xlab = ifelse(plotLongitude, "Longitude", "")
  ylab = ifelse(plotLongitude, "Latitude", "")
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=bquote(Integration~points~and~associated~locations~(K^{MICS} == .(numPtsUrb))), 
             leaveRoomForLegend = TRUE, legend.mar=3, xlab=xlab, ylab=ylab)
  plotWithColor(fineGrid$lon[fineGrid$urban], fineGrid$lat[fineGrid$urban], 
                match(medoidIUrb, sort(unique(medoidIUrb))), pch=15, cex=.35, 
                colScale=cols[1:nClustersUrb], add=TRUE, addColorBar=TRUE, 
                legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                                legend.cex=1, smallplot=c(.9,.93,.55,.9), 
                                legend.lab="urban", legend.line=0.5))
  plotWithColor(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], 
                nClustersUrb + match(medoidIRur, sort(unique(medoidIRur))), 
                pch=15, cex=.35, colScale=cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], 
                add=TRUE, addColorBar=TRUE, 
                legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                                legend.cex=1, smallplot=c(.9,.93,.15,.5), 
                                legend.lab="rural", legend.line=0.5))
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col=cols[1:nClustersUrb], pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col=cols[(nClustersUrb+1):length(cols)], pch=1)
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, pch=1)
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  dev.off()
  
  # pdf("figures/integration/integrationMICSbyClusterQuilt.pdf", width=5, height=5)
  # par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # nClustersUrb = nrow(intPoints$ptsUrb)
  # nClustersRur = nrow(intPoints$ptsRur)
  # cols = rainbow(nClustersUrb + nClustersRur)
  # cols[1:nClustersUrb] = sample(cols[1:nClustersUrb], nClustersUrb)
  # cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)] = sample(cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], nClustersRur)
  # urbCols = cols[match(medoidIUrb, sort(unique(medoidIUrb)))]
  # rurCols = cols[nClustersUrb + match(medoidIRur, sort(unique(medoidIRur)))]
  # 
  # # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  # plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
  #            main=paste0("Integration points and associated locations"), leaveRoomForLegend = TRUE, 
  #            legend.mar=3, xlab="")
  # nx = 95
  # myQuiltPlot(x=fineGrid$lon[fineGrid$urban], y=fineGrid$lat[fineGrid$urban], 
  #             z=match(medoidIUrb, sort(unique(medoidIUrb))), nx=nx, asp=1, 
  #             xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
  #             colScale=cols[1:nClustersUrb], add=TRUE, addColorBar=TRUE, 
  #             legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
  #                             legend.cex=1, smallplot=c(.9,.93,.55,.9), 
  #                             legend.lab="urban", legend.line=0.5))
  # myQuiltPlot(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], 
  #             nClustersUrb + match(medoidIRur, sort(unique(medoidIRur))), 
  #             colScale=cols[(nClustersUrb+1):(nClustersUrb+nClustersRur)], 
  #             xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
  #             add=TRUE, addColorBar=TRUE, nx=nx, asp=1, 
  #             legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
  #                             legend.cex=1, smallplot=c(.9,.93,.15,.5), 
  #                             legend.lab="rural", legend.line=0.5))
  # # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col=cols[1:nClustersUrb], pch=0)
  # # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col=cols[(nClustersUrb+1):length(cols)], pch=1)
  # # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, pch=0)
  # # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, pch=1)
  # outUrb = myQuiltPlot(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, 
  #                      1:length(intPoints$ptsUrb$lon), 
  #                      colScale=rep("black", length(intPoints$ptsUrb$lon)), 
  #                      xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
  #                      add=TRUE, addColorBar=FALSE, nx=nx, asp=1, plot=FALSE)
  # lonUrb = rev(outUrb$x[outUrb$ind[,1]])[-(1:4)]
  # latUrb = rev(outUrb$y[outUrb$ind[,2]])[-(1:4)]
  # outRur = myQuiltPlot(intPoints$ptsRur$lon, intPoints$ptsRur$lat, 
  #                      1:length(intPoints$ptsRur$lon), 
  #                      colScale=rep("black", length(intPoints$ptsRur$lon)), 
  #                      xlim=range(fineGrid$lon), ylim=range(fineGrid$lat), 
  #                      add=TRUE, addColorBar=FALSE, nx=nx, asp=1, plot=FALSE)
  # lonRur = rev(outRur$x[outRur$ind[,1]])[-(1:4)]
  # latRur = rev(outRur$y[outRur$ind[,2]])[-(1:4)]
  # # points(lonUrb, latUrb, cex=.5, pch=0)
  # # points(lonRur, latRur, cex=.5, pch=1)
  # points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  # 
  # dev.off()
  
  # plot covariates
  popCols = makeBlueSequentialColors(64)
  covCols = makePurpleYellowSequentialColors(64)
  
  if(FALSE) {
    pdf(paste0("figures/integration/", subDir, "integrationMICSbyUrbanTest.pdf"), width=5, height=5)
    par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
    # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
    plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
               main=paste0("Urban/rural classification"), leaveRoomForLegend=TRUE, 
               legend.mar=3, xlab="")
    points(fineGridTest$lon[fineGridTest$urban], fineGridTest$lat[fineGridTest$urban], pch=15, cex=.35, col="blue")
    points(fineGridTest$lon[!fineGridTest$urban], fineGridTest$lat[!fineGridTest$urban], pch=15, cex=.35, col="green")
    points(fineGridTest$lon[is.na(fineGridTest$urban)], fineGridTest$lat[is.na(fineGridTest$urban)], pch=15, cex=.35, col="red")
    # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col="darkblue", pch=0)
    # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col="darkgreen", pch=1)
    points(popMatNGAThresh$lon, popMatNGAThresh$lat, cex=.5, pch=0)
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
  }
  
  # urban
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyUrban.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Urban/rural classification"), leaveRoomForLegend=TRUE, 
             legend.mar=3, xlab="")
  points(fineGrid$lon[fineGrid$urban], fineGrid$lat[fineGrid$urban], pch=15, cex=.35, col="blue")
  points(fineGrid$lon[!fineGrid$urban], fineGrid$lat[!fineGrid$urban], pch=15, cex=.35, col="green")
  # points(intPoints$ptsUrb$lon, intPoints$ptsUrb$lat, cex=.5, col="darkblue", pch=0)
  # points(intPoints$ptsRur$lon, intPoints$ptsRur$lat, cex=.5, col="darkgreen", pch=1)
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
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
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyLog1pPop.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Log1p population density"), legend.mar=3, xlab="", ylab="")
  zlim = range(fineGrid$normPop)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$normPop, colScale=popCols, 
                add=TRUE, pch=15, cex=.35, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  # plotWithColor(intPoints$pts$lon, intPoints$pts$lat, intPoints$pts$log1pPop, colScale=popCols, 
  #               add=TRUE, pch=1, cex=.5, zlim=zlim, addColorBar=FALSE)
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  dev.off()
  
  # access
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyAccess.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Healthcare Inaccessibility"), legend.mar=3)
  zlim = range(fineGrid$access)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$access, colScale=covCols, 
                add=TRUE, pch=15, cex=.35, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  dev.off()
  
  # elev
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyElev.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Elevation"), legend.mar=3, ylab="")
  zlim = range(fineGrid$elev)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$elev, colScale=covCols, 
                add=TRUE, pch=15, cex=.35, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  dev.off()
  
  # distRiversLakes
  pdf(paste0("figures/integration/", subDir, "integrationMICSbyDist.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 2, 0), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  # plot(stratumMICSMapDat[stratumMICSMapDat@data[[stratumMICSNameVar]] == stratum,])
  plotMapDat(mapDat=stratumMICSMapDatFull, xlim=bbox[1,], ylim=bbox[2,], 
             main=paste0("Distance to rivers and lakes"), legend.mar=3, xlab="", ylab="")
  zlim = range(fineGrid$distRiversLakes)
  plotWithColor(fineGrid$lon, fineGrid$lat, fineGrid$distRiversLakes, colScale=covCols, 
                add=TRUE, pch=15, cex=.35, zlim=zlim, addColorBar=TRUE, legend.mar=0, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.15,.9)))
  points(intPoints$pts$lon, intPoints$pts$lat, cex=.5, pch=0)
  dev.off()
}

# plot rasters and Nigeria DHS and MICS data
plotDatasets = function(kmres=5) {
  
  # plot rasters ----
  
  # project area to easting/northing in km, and make a fine easting/northing
  # grid of points
  projArea = projNigeriaArea(adm0Full)
  xCoords = seq(projArea@bbox[1,1], projArea@bbox[1,2], by=kmres)
  yCoords = seq(projArea@bbox[2,1], projArea@bbox[2,2], by=kmres)
  ENCoords = make.surface.grid(list(x=xCoords, y=yCoords))
  
  # convert the grid of points back to longitude/latitude and get covariate
  # values
  LLCoords = projNigeria(ENCoords, inverse=TRUE)
  inArea = over(SpatialPoints(LLCoords, CRS(SRS_string="EPSG:4326")), adm0Full)
  inArea = !is.na(inArea$GID_0)
  
  ENCoords = ENCoords[inArea,]
  LLCoords = LLCoords[inArea,]
  X = getDesignMat(LLCoords, normalized=FALSE)
  Xtest = getDesignMat(LLCoords, normalized=FALSE, testMode=TRUE)
  
  out = load("savedOutput/global/covariates.RData")
  
  
  # pdf("figures/data/pop.pdf", width=5, height=5)
  # par(mar=c(3.5, 4, 2.5, 5))
  # # quilt.plot(1:2, 1:2, rep(NA, 2), nx=2, ny=2, add.legend = FALSE, 
  # #            xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
  # #            ylab="Latitude", main="Population density", legend.mar=7)
  # myQuiltPlot(1:2, 1:2, 1:2, nx=2, ny=2, colScale=rep(rgb(0,0,0,0), 2), 
  #            xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
  #            ylab="Latitude", main="Population density", 
  #            addColorBar=FALSE)
  # plotMapDat(mapDat=adm0, border=rgb(1,1,1,1), new=FALSE)
  # plotMapDat(mapDat=adm1, new=FALSE)
  # plotMapDat(mapDat=adm2, lwd=.2, border=rgb(.5, .5, .5, .5), new=FALSE)
  # plot(pop, add=TRUE)
  # dev.off()
  
  popColScale = makeBlueSequentialColors(64)
  
  if(FALSE) {
    pdf("figures/data/popQuiltTest.pdf", width=5, height=3.8)
    par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
    # myQuiltPlot(1:2, 1:2, 1:2, nx=2, ny=2, colScale=rep(rgb(0,0,0,0), 2), 
    #             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
    #             ylab="Latitude", main="Population density", 
    #             addColorBar=FALSE)
    nas = is.na(Xtest[,7])
    myQuiltPlot(LLCoords[!nas,1], LLCoords[!nas,2], Xtest[!nas,2], nx=256, colScale=popColScale,
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude",
                ylab="Latitude", main="Population density", legend.mar=4.7,
                addColorBar=TRUE, scaleFun=log1p, 
                ticks=c(0, 1, 10, 100, 1000, 10000, 100000))
    myQuiltPlot(LLCoords[nas,1], LLCoords[nas,2], 1:sum(nas), nx=256, colScale=c("green", "green"),
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, legend.mar=4.7,
                addColorBar=FALSE, add=TRUE)
    # plotMapDat(mapDat=adm2, lwd=0.5, new=FALSE, 
    #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
    # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
    plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
               border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
    plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
    # plot(pop, add=TRUE)
    dev.off()
  }
  
  pdf("figures/data/popQuilt.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  # myQuiltPlot(1:2, 1:2, 1:2, nx=2, ny=2, colScale=rep(rgb(0,0,0,0), 2), 
  #             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
  #             ylab="Latitude", main="Population density", 
  #             addColorBar=FALSE)
  myQuiltPlot(LLCoords[,1], LLCoords[,2], X[,2], nx=256, colScale=popColScale,
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude",
              ylab="Latitude", main="Population density", legend.mar=4.7,
              addColorBar=TRUE, scaleFun=log1p, 
              ticks=c(0, 1, 10, 100, 1000, 10000, 100000))
  # plotMapDat(mapDat=adm2, lwd=0.5, new=FALSE, 
  #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  # plot(pop, add=TRUE)
  dev.off()
  
  pdf("figures/data/urbQuilt.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  myQuiltPlot(LLCoords[X[,3]==1,1], LLCoords[X[,3]==1,2], X[X[,3]==1,3] + runif(sum(X[,3]==1)), nx=256, colScale=c("blue", "blue"), 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Urban classification", legend.mar=4.7,
              addColorBar=FALSE, leaveRoomForLegend=TRUE)
  myQuiltPlot(LLCoords[,1], LLCoords[,2], X[,3], nx=256, colScale=c("blue", "blue"), 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Urban classification", legend.mar=4.7,
              addColorBar=TRUE, leaveRoomForLegend=TRUE, legendOnly=TRUE, 
              legendArgs=list(axis.args=list(labels=FALSE, tick=FALSE), 
                              legend.cex=1, smallplot=c(.825,.86,.19,.85), 
                              legend.lab="urban", legend.line=0.5))
  # plotMapDat(mapDat=adm2, lwd=0.5, new=FALSE, 
  #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  # plotMapDat(mapDat=adm0, border=rgb(0,0,0,1), new=FALSE, lwd=2)
  # plot(pop, add=TRUE)
  dev.off()
  
  covColScale = makeBlueGreenYellowSequentialColors(64)
  
  pdf("figures/data/accessQuilt.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  myQuiltPlot(LLCoords[,1], LLCoords[,2], X[,4], nx=256, colScale=covColScale, 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Healthcare inaccessibility", legend.mar=4.7,
              addColorBar=TRUE, scaleFun=log1p, ticks=c(0, 5, 50, 500))
  # plotMapDat(mapDat=adm2, lwd=.5, new=FALSE, 
  #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  # plot(pop, add=TRUE)
  dev.off()
  
  pdf("figures/data/elevQuilt.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  myQuiltPlot(LLCoords[,1], LLCoords[,2], X[,5], nx=256, colScale=covColScale, 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Elevation", legend.mar=4.7, zlim=c(0, 2040), 
              addColorBar=TRUE, scaleFun=function(x){sign(x)*sqrt(abs(x))}, 
              ticks=c(0, 100, 500, 1000, 1500, 2000))
  # plotMapDat(mapDat=adm2, lwd=.5, new=FALSE, 
  #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  # plot(pop, add=TRUE)
  dev.off()
  
  pdf("figures/data/distQuilt.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  myQuiltPlot(LLCoords[,1], LLCoords[,2], X[,6], nx=256, colScale=covColScale, 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Distance to rivers and lakes", legend.mar=4.7,
              addColorBar=TRUE)
  # plotMapDat(mapDat=adm2, lwd=.5, new=FALSE, 
  #            border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  # plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  # plot(pop, add=TRUE)
  dev.off()
  
  browser()
  
  # plot DHS data ----
  # ed
  
  pdf("figures/data/edDHS.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  plotWithColor(ed$lon, ed$lat, ed$y/ed$n, colScale=covColScale, 
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
                ylab="Latitude", main="Women's secondary education prevalence (DHS)", legend.mar=4.7,
                addColorBar=TRUE, pch=19, cex=.3, ordering="increasing")
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  dev.off()
  
  # plot MICS data ----
  # edMICS
  
  covColScale = makeBlueGreenYellowSequentialColors(64)
  
  # first average at the MICS stratum level
  ys = aggregate(edMICS$ys, by=list(edMICS$Stratum), FUN=sum)
  ns = aggregate(edMICS$ns, by=list(edMICS$Stratum), FUN=sum)
  prevs = ys[,2]/ns[,2]
  
  pdf("figures/data/edMICS.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  plotMapDat(mapDat=admFinal, cols=covColScale, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
             ylab="Latitude", main="Women's secondary education prevalence (MICS)", legend.mar=4.7,
             addColorBar=TRUE, plotVar=prevs, varAreas=ys[,1], regionNames=admFinal$NAME_FINAL)
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  dev.off()
}

# Application section ----

# Plots:
#  1)-3) For Md, M_D, Mdm, and M_DM models, plot posterior mean prevalence and 95% 
#        CI width at pixel, Admin2, and Admin1 levels
#  4)-6) For M_D, M_M, and M_DM models, plot posterior mean prevalence and 95% CI 
#        width at pixel, Admin2, and Admin1 levels
#  7)-9) Plot posterior mean prevalence and 95% CI width for M_DM model, and 
#        percent diff in mean and CI widths to Md and Mdm models at pixel, 
#        Admin2, and Admin1 levels
# 10-12) Same as 7)-9) but plot percent diffs on vert axis versus M_DM values on 
#        horiz axis
plotAppPreds = function(adm2Model=FALSE) {
  sep = TRUE
  sepRepar = TRUE
  adm2Text = ifelse(adm2Model, "2", "1")
  
  # load predictions
  if(!adm2Model) {
    # grid preds
    out = load("savedOutput/ed/gridPredsM_DMSepRepar.RData")
    gridPredsM_DM = gridPreds
    # out = load("savedOutput/ed/gridPredsMdmSepRepar.RData")
    # gridPredsMdm = gridPreds
    out = load("savedOutput/ed/gridPredsM_DSepRepar.RData")
    gridPredsM_D = gridPreds
    # out = load("savedOutput/ed/gridPredsM_MSepRepar.RData")
    # gridPredsM_M = gridPreds
    out = load("savedOutput/ed/gridPredsMdSepRepar.RData")
    gridPredsMd = gridPreds
    
    # admin2 preds
    out = load("savedOutput/ed/admin2PredsM_DMSepRepar.RData")
    admin2PredsM_DM = admin2Preds
    # out = load("savedOutput/ed/admin2PredsMdmSepRepar.RData")
    # admin2PredsMdm = admin2Preds
    out = load("savedOutput/ed/admin2PredsM_DSepRepar.RData")
    admin2PredsM_D = admin2Preds
    # out = load("savedOutput/ed/admin2PredsM_MSepRepar.RData")
    # admin2PredsM_M = admin2Preds
    out = load("savedOutput/ed/admin2PredsMdSepRepar.RData")
    admin2PredsMd = admin2Preds
    
    # admin1 preds
    out = load("savedOutput/ed/admin1PredsM_DMSepRepar.RData")
    admin1PredsM_DM = admin1Preds
    # out = load("savedOutput/ed/admin1PredsMdmSepRepar.RData")
    # admin1PredsMdm = admin1Preds
    out = load("savedOutput/ed/admin1PredsM_DSepRepar.RData")
    admin1PredsM_D = admin1Preds
    # out = load("savedOutput/ed/admin1PredsM_MSepRepar.RData")
    # admin1PredsM_M = admin1Preds
    out = load("savedOutput/ed/admin1PredsMdSepRepar.RData")
    admin1PredsMd = admin1Preds
  } else {
    # grid preds
    out = load("savedOutput/ed/gridPredsM_DM2SepRepar.RData")
    gridPredsM_DM = gridPreds
    # out = load("savedOutput/ed/gridPredsMdm2SepRepar.RData")
    # gridPredsMdm = gridPreds
    out = load("savedOutput/ed/gridPredsM_D2SepRepar.RData")
    gridPredsM_D = gridPreds
    # out = load("savedOutput/ed/gridPredsM_M2SepRepar.RData")
    # gridPredsM_M = gridPreds
    out = load("savedOutput/ed/gridPredsMd2SepRepar.RData")
    gridPredsMd = gridPreds
    
    # admin2 preds
    out = load("savedOutput/ed/admin2PredsM_DM2SepRepar.RData")
    admin2PredsM_DM = admin2Preds
    # out = load("savedOutput/ed/admin2PredsMdm2SepRepar.RData")
    # admin2PredsMdm = admin2Preds
    out = load("savedOutput/ed/admin2PredsM_D2SepRepar.RData")
    admin2PredsM_D = admin2Preds
    # out = load("savedOutput/ed/admin2PredsM_M2SepRepar.RData")
    # admin2PredsM_M = admin2Preds
    out = load("savedOutput/ed/admin2PredsMd2SepRepar.RData")
    admin2PredsMd = admin2Preds
    
    # admin1 preds
    out = load("savedOutput/ed/admin1PredsM_DM2SepRepar.RData")
    admin1PredsM_DM = admin1Preds
    # out = load("savedOutput/ed/admin1PredsMdm2SepRepar.RData")
    # admin1PredsMdm = admin1Preds
    out = load("savedOutput/ed/admin1PredsM_D2SepRepar.RData")
    admin1PredsM_D = admin1Preds
    # out = load("savedOutput/ed/admin1PredsM_M2SepRepar.RData")
    # admin1PredsM_M = admin1Preds
    out = load("savedOutput/ed/admin1PredsMd2SepRepar.RData")
    admin1PredsMd = admin1Preds
  }
  
  LLCoords = cbind(popMatNGAThresh$lon, popMatNGAThresh$lat)
  
  CIwidth = function(x, signif=.95) {
    alpha = 1-signif
    diff(quantile(prob=c(alpha/2, 1-alpha/2), x, na.rm=TRUE))
  }
  
  predCols = makeBlueGreenYellowSequentialColors(64)
  widthCols = makePurpleYellowSequentialColors(64)
  
  if(adm2Model) {
    models = c(expression(M[d]^2), expression(M[D]^2), expression(M[DM]^2))
  } else {
    models = c(expression(M[d]^1), expression(M[D]^1), expression(M[DM]^1))
  }
  
  # lwd = ifelse(doAdm2, .6, 1)
  # border = ifelse(doAdm2, rgb(.6, .6, .6), "black")
  
  lwd = .6
  border = rgb(.6, .6, .6)
  
  # 1) Md, M_D, M_DM models preds and CI widths at pixel level
  predsM_DM = rowMeans(gridPredsM_DM$gridDraws)
  predsM_D = rowMeans(gridPredsM_D$gridDraws)
  predsMd = rowMeans(gridPredsMd$gridDraws)
  widthM_DM = apply(gridPredsM_DM$gridDraws, 1, CIwidth)
  widthM_D = apply(gridPredsM_D$gridDraws, 1, CIwidth)
  widthMd = apply(gridPredsMd$gridDraws, 1, CIwidth)
  predsLim = range(c(predsM_DM, predsM_D, predsMd))
  widthLim = range(c(widthM_DM, widthM_D, widthMd))
  
  pdf(paste0("figures/ed/main", adm2Text, "PredsPixel.pdf"), width=8, height=4.5)
  par(mar=c(1, 1, 1, 2), mgp=c(1.7, .5, 0), mfrow=c(2,3), oma=c(2.1,5,2,3))
  
  # preds
  myQuiltPlot(LLCoords[,1], LLCoords[,2], predsMd, nx=256, colScale=predCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=predsLim, asp=1, xlab="", 
              ylab="", main="", legend.mar=0,
              addColorBar=FALSE, leaveRoomForLegend=FALSE)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("Estimates", side=2, line=4)
  mtext(models[1], side=3, line=.7)
  
  myQuiltPlot(LLCoords[,1], LLCoords[,2], predsM_D, nx=256, colScale=predCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=predsLim, asp=1, xlab="", 
              ylab="", main="", legend.mar=0,
              addColorBar=FALSE, leaveRoomForLegend=FALSE)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  
  mtext(models[2], side=3, line=1)
  
  myQuiltPlot(LLCoords[,1], LLCoords[,2], predsM_DM, nx=256, colScale=predCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=predsLim, asp=1, xlab="", 
              ylab="", main="", addColorBar=TRUE, leaveRoomForLegend=TRUE, 
              legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                              legend.cex=1, smallplot= c(.96,1,.1,.9)), 
              legend.width=3, legend.mar=0)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  
  mtext(models[3], side=3, line=1)
  
  # widths
  myQuiltPlot(LLCoords[,1], LLCoords[,2], widthMd, nx=256, colScale=widthCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=widthLim, asp=1, xlab="", 
              ylab="", main="", legend.mar=0,
              addColorBar=FALSE, leaveRoomForLegend=FALSE)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("95% CI Width", side=2, line=4)
  mtext("Longitude", side=1, line=2, cex=.8)
  
  
  myQuiltPlot(LLCoords[,1], LLCoords[,2], widthM_D, nx=256, colScale=widthCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=widthLim, asp=1, xlab="", 
              ylab="", main="", legend.mar=0,
              addColorBar=FALSE, leaveRoomForLegend=FALSE)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  myQuiltPlot(LLCoords[,1], LLCoords[,2], widthM_DM, nx=256, colScale=widthCols, 
              xlim=lonLimNGA, ylim=latLimNGA, zlim=widthLim, asp=1, xlab="", 
              ylab="", main="", addColorBar=TRUE, leaveRoomForLegend=TRUE, 
              legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                              legend.cex=1, smallplot= c(.96,1,.1,.9)), 
              legend.width=3, legend.mar=0)
  plotMapDat(mapDat=adm2, new=FALSE, lwd=lwd, border=border)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  dev.off()
  
  
  
  
  
  # 2) Md, M_D, M_DM models preds and CI widths at Admin2 level
  predsM_DM = rowMeans(admin2PredsM_DM$aggregationResults$p, na.rm=TRUE)
  predsM_D = rowMeans(admin2PredsM_D$aggregationResults$p, na.rm=TRUE)
  predsMd = rowMeans(admin2PredsMd$aggregationResults$p, na.rm=TRUE)
  widthM_DM = apply(admin2PredsM_DM$aggregationResults$p, 1, CIwidth)
  widthM_D = apply(admin2PredsM_D$aggregationResults$p, 1, CIwidth)
  widthMd = apply(admin2PredsMd$aggregationResults$p, 1, CIwidth)
  predsLim = range(c(predsM_DM, predsM_D, predsMd), na.rm=TRUE)
  widthLim = range(c(widthM_DM, widthM_D, widthMd), na.rm=TRUE)
  
  pdf(paste0("figures/ed/main", adm2Text, "PredsAdmin2.pdf"), width=8, height=4.5)
  par(mar=c(1, 1, 1, 2), mgp=c(1.7, .5, 0), mfrow=c(2,3), oma=c(2.1,5,2,3))
  
  # preds
  plotMapDat(adm2, predsMd, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=predCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=predsLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("Estimates", side=2, line=4)
  mtext(models[1], side=3, line=.7)
  
  plotMapDat(adm2, predsM_D, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=predCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=predsLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  plotMapDat(admFinal, new=FALSE)
  
  mtext(models[2], side=3, line=1)
  
  plotMapDat(adm2, predsM_DM, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=predCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, zlim=predsLim, asp=1, xlab="", 
             addColorBar=TRUE, leaveRoomForLegend=TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                             legend.cex=1, smallplot= c(.96,1,.1,.9)), 
             legend.width=3, legend.mar=0)
  plotMapDat(admFinal, new=FALSE)
  
  mtext(models[3], side=3, line=1)
  
  # widths
  plotMapDat(adm2, widthMd, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=widthCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=widthLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  plotMapDat(admFinal, new=FALSE)
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("95% CI Width", side=2, line=4)
  mtext("Longitude", side=1, line=2, cex=.8)
  
  
  plotMapDat(adm2, widthM_D, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=widthCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=widthLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  plotMapDat(adm2, widthM_DM, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             cols=widthCols, crosshatchNADensity=30, lwd=lwd, border=border, 
             xlim=lonLimNGA, ylim=latLimNGA, zlim=widthLim, asp=1, xlab="", 
             addColorBar=TRUE, leaveRoomForLegend=TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                             legend.cex=1, smallplot= c(.96,1,.1,.9)), 
             legend.width=3, legend.mar=0)
  plotMapDat(admFinal, new=FALSE)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  dev.off()
  
  
  
  
  
  # 3) Md, M_D, M_DM models preds and CI widths at Admin1 level
  predsM_DM = rowMeans(admin1PredsM_DM$aggregationResults$p, na.rm=TRUE)
  predsM_D = rowMeans(admin1PredsM_D$aggregationResults$p, na.rm=TRUE)
  predsMd = rowMeans(admin1PredsMd$aggregationResults$p, na.rm=TRUE)
  widthM_DM = apply(admin1PredsM_DM$aggregationResults$p, 1, CIwidth)
  widthM_D = apply(admin1PredsM_D$aggregationResults$p, 1, CIwidth)
  widthMd = apply(admin1PredsMd$aggregationResults$p, 1, CIwidth)
  predsLim = range(c(predsM_DM, predsM_D, predsMd), na.rm=TRUE)
  widthLim = range(c(widthM_DM, widthM_D, widthMd), na.rm=TRUE)
  
  pdf(paste0("figures/ed/main", adm2Text, "PredsAdmin1.pdf"), width=8, height=4.5)
  par(mar=c(1, 1, 1, 2), mgp=c(1.7, .5, 0), mfrow=c(2,3), oma=c(2.1,5,2,3))
  
  # preds
  plotMapDat(adm1, predsMd, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=predCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=predsLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("Estimates", side=2, line=4)
  mtext(models[1], side=3, line=.7)
  
  plotMapDat(adm1, predsM_D, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=predCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=predsLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  
  mtext(models[2], side=3, line=1)
  
  plotMapDat(adm1, predsM_DM, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=predCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, zlim=predsLim, asp=1, xlab="", 
             addColorBar=TRUE, leaveRoomForLegend=TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                             legend.cex=1, smallplot= c(.96,1,.1,.9)), 
             legend.width=3, legend.mar=0)
  
  mtext(models[3], side=3, line=1)
  
  # widths
  plotMapDat(adm1, widthMd, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=widthCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=widthLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  
  mtext("Latitude", side=2, line=2, cex=.8)
  mtext("95% CI Width", side=2, line=4)
  mtext("Longitude", side=1, line=2, cex=.8)
  
  
  plotMapDat(adm1, widthM_D, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=widthCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="", 
             ylab="", main="", legend.mar=0, zlim=widthLim, addColorBar=FALSE, 
             leaveRoomForLegend=FALSE)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  plotMapDat(adm1, widthM_DM, varAreas=adm1@data$NAME_1, regionNames=adm1@data$NAME_1, 
             cols=widthCols, crosshatchNADensity=30, 
             xlim=lonLimNGA, ylim=latLimNGA, zlim=widthLim, asp=1, xlab="", 
             addColorBar=TRUE, leaveRoomForLegend=TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=0), 
                             legend.cex=1, smallplot= c(.96,1,.1,.9)), 
             legend.width=3, legend.mar=0)
  
  mtext("Longitude", side=1, line=2, cex=.8)
  
  dev.off()
}











