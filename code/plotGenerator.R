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

makeRedBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
  else
    scale_colour_continuous_sequential(h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3, n_interp=n)
}

makeGreenBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  else
    scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
}

makeGreenBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE, p1=1) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  # else
  #   scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev)
    else
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeGreenBlueDivergingColors(totalColors, rev=rev, p1=p1)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n, mid=center)
    }
  }
}

makePurpleYellowSequentialColors = function(n, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev)
  else
    scale_colour_continuous_sequential(h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev, n_interp=n)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    if(!ggplot)
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev, mid=center)
    }
  }
}

makeBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE, n_interp=n)
}

makeGreenSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE, n_interp=n)
}

makeBlueYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n)
}

makeYellowRedSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2)
  else
    scale_colour_continuous_sequential(h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2, n_interp=n)
}

makeRedGreenDivergingColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  else
    scale_colour_continuous_sequential(h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5, n_interp=n)
}

# given continuous color scale and range, chooses colors based on a set of values
getColorsFromScale = function(vals, valRange=range(vals), cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}







