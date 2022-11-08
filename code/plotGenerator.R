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

# plot rasters and Nigeria DHS and MICS data
plotDatasets = function(kmres=5) {
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







