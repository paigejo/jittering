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
  
  # first average at the MICS stratum level
  ys = aggregate(edMICS$secondaryEd, by=list(edMICS$Stratum), FUN=sum)
  ns = aggregate(edMICS$secondaryEd, by=list(edMICS$Stratum), FUN=length)
  prevs = ys[,2]/ns[,2]
  
  pdf("figures/data/edMICS.pdf", width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  plotMapDat(mapDat=admFinal, cols=covColScale, 
             xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
             ylab="Latitude", main="Women's secondary education prevalence (MICS)", legend.mar=4.7,
             addColorBar=TRUE, plotVar=prevs, varAreas=ys[,1])
  plotMapDat(mapDat=sen, lwd=.5, new=FALSE,
             border=do.call("rgb", c(as.list(rep(0, 3)), list(1))))
  plotMapDat(mapDat=adm1, new=FALSE, lwd=1.5)
  dev.off()
}







