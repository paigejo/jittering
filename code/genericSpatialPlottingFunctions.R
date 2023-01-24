# helpful spatial plotting functions using base R and fields:

# continuous plotting ----
# x, y: horizontal and vertical spatial coordinates
# z: response
# zlim: range of the response
# cols: color vector representing the color scale
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# legend.mar, n.ticks, legend.width, legendArgs (as legend.args): see ?image.plot
# min.n: approximate number of ticks in color scale. See ?pretty
# orderI: a specific ordering to plot the points in
# ordering: in what order to plot the points, where the ordering is based on the response z
plotWithColor = function(x, y, z, zlim=NULL, colScale=tim.colors(), 
                         legend.mar=7, add=FALSE, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                         n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                         legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                         ordering=c("none", "increasing", "decreasing"), colorName = c("col", "bg"), ...) {
  ordering = match.arg(ordering)
  colorName = match.arg(colorName)
  
  # remove NA points
  nas = is.na(x) | is.na(y) | is.na(z)
  if(any(nas)) {
    warning("Removing NAs")
    x = x[!nas]
    y = y[!nas]
    z = z[!nas]
  }
  
  # do setup for plotting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(z))
    zlim = range(z[!nas])
  }
  
  # order the plotting of the points
  if(is.null(orderI)) {
    if(ordering == "increasing") {
      orderI = sort(z, index.return=TRUE)$ix
    } else if(ordering == "decreasing") {
      orderI = sort(z, decreasing=TRUE, index.return=TRUE)$ix
    } else {
      orderI = 1:length(z)
    }
  }
  x = x[orderI]
  y = y[orderI]
  z = z[orderI]
  
  # if(forceColorsInRange) {
  #   z[z > zlim[2]] = zlim[2]
  #   z[z < zlim[1]] = zlim[1]
  # }
  
  # get colors of points
  cols = getColorsFromScale(z, zlim, cols=colScale, scaleFun=scaleFun, 
                            forceValuesInRange=forceColorsInRange)
  
  # generate new plot if necessary
  # browser()
  if(!add) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    if(colorName == "col") {
      do.call("plot", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("plot", c(list(x=x, y=y, bg=cols), list(...)))
    }
  } else {
    if(colorName == "col") {
      do.call("points", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("points", c(list(x=x, y=y, bg=cols), list(...)))
    }
  }
  
  if(addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(tickLabels))
      setTickLabels = TRUE
    
    if(is.null(ticks)) {
      if(setTickLabels)
        tickLabels = pretty(zlim, n=n.ticks, min.n=min.n)
      ticks = scaleFun(tickLabels)
    }
    else {
      if(setTickLabels)
        tickLabels = ticks
      ticks = scaleFun(ticks)
    }
    if(setTickLabels)
      tickLabels = tickLabels[is.finite(ticks)]
    ticks = ticks[is.finite(ticks)]
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    
    legendArgs$zlim=scaleFun(zlim)
    legendArgs$nlevel=length(colScale)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=colScale
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      if(is.null(legendArgs$axis.args$at)) {
        legendArgs$axis.args$at=ticks
      }
      if(is.null(legendArgs$axis.args$labels)) {
        legendArgs$axis.args$labels=tickLabels
      }
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# x, y: horizontal and vertical spatial coordinates
# z: response
# zlim: range of the response
# cols: color vector representing the color scale
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# legend.mar, n.ticks, legend.width, legendArgs (as legend.args): see ?image.plot
# min.n: approximate number of ticks in color scale. See ?pretty
# orderI: a specific ordering to plot the points in
# ordering: in what order to plot the points, where the ordering is based on the response z
# setGridResMethod: if nx and ny are both NULL, this is the method used to set them automatically.
#     Rice: From Rice 1944. Set number of bins within the domain of the data to be 
#           2*n^(1/3)
#     Sturges: From Sturges 1926. Set number of bins within the domain of the data to be 
#              ceiling(log2(n)) + 1. Decide nx, ny proportionally to xlim, ylim
#     minDist: Best for if the data is already on a grid. Calculates the minimum distance 
#              between x values and between y values and sets nx and ny accordingly. 
#              Recommended for data on an exact grid.
# xyResRatio: x resolution divided by y resolution (e.g. if 2, twice as many rows per 
#             unit length as columns). Used if Defaults to 1, and used only if 
#             setGridResMethod != "minDist" or one of nx or ny has been specified
# minDistTol: for setGridResMethod == "minDist", removes distances lower than the tolerance 
#             times the x and y limit widths when checking for the lowest distances.
myQuiltPlot = function(x, y, z, zlim=NULL, colScale=tim.colors(), nx=64, ny=NULL, 
                       xlim=NULL, ylim=NULL, legend.mar=7, add=FALSE, plot=TRUE, 
                       scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                       n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                       legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                       FUN=function(x){mean(x, na.rm=TRUE)}, na.rm=FALSE, setGridResMethod=c("Rice", "Sturges", "minDist"), 
                       xyResRatio=1, minDistTol=1e-3, ...) {
  setGridResMethod = match.arg(setGridResMethod)
  
  ## code from fields:::quilt.plot
  x <- as.matrix(x)
  if (ncol(x) == 2) {
    z <- y
  }
  if (ncol(x) == 1) {
    x <- cbind(x, y)
  }
  if (ncol(x) == 3) {
    z <- x[, 3]
    x <- x[, 1:2]
  }
  
  ## modifications to fields code:
  
  if(is.null(xlim)) {
    xlim = range(x[,1], na.rm=TRUE)
  }
  if(is.null(ylim)) {
    ylim = range(x[,2], na.rm=TRUE)
  }
  
  if(is.null(nx) || is.null(ny)) {
    ensureNxyProp = setGridResMethod != "minDist"
  }
  else {
    ensureNxyProp = FALSE
  }
  
  # if nx or ny is set and ensureNxyProp is TRUE, set the unset of nx or ny accordingly
  if(ensureNxyProp) {
    if(!is.null(nx) && is.null(ny)) {
      LyLxRatio = diff(ylim)/diff(xlim)
      ny = round(nx * LyLxRatio/xyResRatio)
    }
    else if(is.null(nx) && !is.null(ny)) {
      LyLxRatio = diff(ylim)/diff(xlim)
      nx = round(xyResRatio * ny/LyLxRatio)
    }
    
    print(paste0("xyResRatio auto set nx and ny to ", nx, " and ", ny))
  }
  
  if(setGridResMethod == "minDist") {
    #     minDist: Best for if the data is already on a grid. Calculates the minimum distance 
    #              between x values and between y values and sets nx and ny accordingly. 
    #              Recommended for data on an exact grid.
    
    if(is.null(nx)) {
      minDistTolX = diff(xlim) * minDistTol
      dists = diff(sort(unique(x[,1])))
      minDistX = min(dists[dists > minDistTolX])
      nx = round(diff(xlim) / minDistX)
    }
    if(is.null(ny)) {
      minDistTolY = diff(ylim) * minDistTol
      dists = diff(sort(unique(x[,2])))
      minDistY = min(dists[dists > minDistTolY])
      ny = round(diff(ylim) / minDistY)
    }
    
    print(paste0("method ", setGridResMethod, " auto set nx and ny to ", nx, " and ", ny))
  }
  
  ## first determine nx, ny automatically if they are null
  if(is.null(nx) && is.null(ny)) {
    # setGridResMethod=c("Rice", "Sturges", "minDist")
    #     Rice: From Rice 1944. Set number of bins within the domain of the data to be 
    #           
    #     Sturges: From Sturges 1926. Set number of bins within the domain of the data to be 
    #              ceiling(log2(n)) + 1. Decide nx, ny proportionally to xlim, ylim
    
    n = nrow(x)
    if(setGridResMethod == "Rice") {
      nBins = 2*n^(1/3)
    }
    else if(setGridResMethod == "Sturges") {
      nBins = ceiling(log2(n)) + 1
    }
    LxShort = diff(range(x[,1]))
    LyShort = diff(range(x[,2]))
    nxShort = sqrt(LxShort/LyShort * xyResRatio * nBins)
    nyShort = nxShort * (LyShort / LxShort) / xyResRatio
    nx = round(nxShort * diff(xlim)/LxShort)
    ny = round(nyShort * diff(ylim)/LyShort)
    
    print(paste0("method ", setGridResMethod, " auto set nx and ny to ", nx, " and ", ny))
  }
  
  # artificially add NA points at boundaries so grids match up when adding multiple grids
  z = c(z, rep(NA, 4))
  x = rbind(x, 
            c(xlim[1], ylim[1]), 
            c(xlim[1], ylim[2]), 
            c(xlim[2], ylim[1]), 
            c(xlim[2], ylim[2]))
  out.p <- as.image(z, x = x, nx = nx, ny = ny, grid = NULL, 
                   FUN = FUN, na.rm = na.rm, boundary.grid = FALSE)
  
  # do setup for plotting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(out.p$z))
    zlim = range(out.p$z[!nas])
  }
  
  if(forceColorsInRange) {
    out.p$z[out.p$z < zlim[1]] = zlim[1]
    out.p$z[out.p$z > zlim[2]] = zlim[2]
  }
  
  
    # if (add.legend) {
    #   image.plot(out.p, nlevel = nlevel, col = col, add = add, 
    #              ...)
    # }
    # else {
    #   image(out.p, col = col, add = add, ...)
    # }
  
  out.p$z = scaleFun(out.p$z)
  
  if(plot) {
    image.args = c(list(out.p, col=colScale, add=add, xlim=xlim, 
                        ylim=ylim, zlim=scaleFun(zlim)), list(...))
    do.call("image", image.args)
    
    if(addColorBar) {
      # add legend
      # par( oma=c(0,0,0,2))
      if(is.null(tickLabels))
        setTickLabels = TRUE
      
      if(is.null(ticks)) {
        if(setTickLabels)
          tickLabels = pretty(zlim, n=n.ticks, min.n=min.n)
        ticks = scaleFun(tickLabels)
      }
      else {
        if(setTickLabels)
          tickLabels = ticks
        ticks = scaleFun(ticks)
      }
      if(setTickLabels)
        tickLabels = tickLabels[is.finite(ticks)]
      ticks = ticks[is.finite(ticks)]
      
      # par( oma=c( 0,0,0,3))
      
      # set list of arguments to image.plot
      
      legendArgs$zlim=scaleFun(zlim)
      legendArgs$nlevel=length(colScale)
      legendArgs$legend.only=TRUE
      legendArgs$horizontal=FALSE
      legendArgs$col=colScale
      legendArgs$add = TRUE
      if(is.null(legendArgs$axis.args))
        legendArgs$axis.args=list(at=ticks, labels=tickLabels)
      else {
        if(is.null(legendArgs$axis.args$at)) {
          legendArgs$axis.args$at=ticks
        }
        if(is.null(legendArgs$axis.args$labels)) {
          legendArgs$axis.args$labels=tickLabels
        }
      }
      legendArgs$legend.mar=legend.mar
      legendArgs$legend.width=legend.width
      
      do.call("image.plot", legendArgs)
      
      # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
      #            col=cols, add = TRUE)
    }
  }
  
  
  invisible(out.p)
}

# TODO: make nice version of quilt.plot including scales for color scale
# base it on this example. Make sure color scale range is centered correctly for
# diverging color scales
# quilt.plot(popMat$lon, popMat$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
#            zlim=log(meanRangePixel), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
#            add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
#            xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
# plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
# points(mort$lon, mort$lat, pch=19, cex=.1)
# 
# meanTicksPixel = getLogScaleTicks(meanRangePixel)
# meanTickLabelsPixel = as.character(meanTicksPixel)
# meanTicks = pretty(meanRange, n=5)
# meanTickLabels = as.character(meanTicks)
# image.plot(zlim=range(log(meanRangePixel)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
#            col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
#            legend.mar = 0, legend.cex=2, legend.width=3, smallplot=c(.88,.91,.1,.9))

# discrete plotting ----

# for plotting areal values assuming plotVar is in alphabetical order of the area names. 
# mapDat: a SpatialPolygonsDataFrame object, where each polygon may have a value you wish to plot
# plotVar: if null, plot the areal boundaries only. Else plot plotVar values for each area
# varAreas: area names associated with plotVar
# zlim: range of the response
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using myProjection function.  This can be used when plotting the projected `easting' 
#          and `northing' variables for instance.
# cols: color vector representing the color scale
# legend.mar, legend.args, n.ticks: see ?image.plot
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# crosshatchNADensity: Adds crosshatching for areas with NA values. See ?polygon density argument.
# # min.n: approximate number of ticks in color scale. See ?pretty
# myProjection: a map projection function taking a 2 column matrix of coordinates 
#   and projects them.
# ...: arguments to polygon function
plotMapDat = function(mapDat, plotVar=NULL, varAreas, regionNames=sort(unique(varAreas)), zlim=NULL, project=FALSE, cols=tim.colors(), 
                      legend.mar=7, new=TRUE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                      ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, addColorBar=TRUE, 
                      legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, 
                      crosshatchNADensity=10, myProjection=NULL, ...) {
  require(fields)
  
  # do setup for plotting data by area if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    if(forceColorsInRange) {
      plotVar[plotVar > zlim[2]] = zlim[2]
      plotVar[plotVar < zlim[1]] = zlim[1]
    }
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(project) {
      if(is.null(xlab))
        xlab = "East (km)"
      if(is.null(xlim))
        xlim = eastLim
      if(is.null(ylab))
        ylab = "North (km)"
      if(is.null(ylim))
        ylim = northLim
    }
    else {
      if(is.null(xlab))
        xlab = "Longitude"
      if(is.null(ylab))
        ylab = "Latitude"
      if(is.null(xlim)) {
        xlim = mapDat@bbox[1,]
      }
      if(is.null(ylim)) {
        ylim = mapDat@bbox[2,]
      }
    }
    if(is.null(main))
      main = ""
    
    if(is.null(plotArgs)) {
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=asp)
    } else {
      plotArgs$main = main
      plotArgs$xlab = xlab
      plotArgs$ylab = ylab
      plotArgs$xlim = xlim
      plotArgs$ylim = ylim
      plotArgs$asp = asp
    }
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotArea = function(i) {
    areaPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      if(!project)
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(areaPolys[[x]]@coords), list(...)))})
      else
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(myProjection(areaPolys[[x]]@coords)), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this area
      thisI = which(varAreas == regionNames[i])
      
      # if there is no matching region name, do nothing
      if(length(thisI) == 0) {
        return(NULL)
      }
      
      # get color to plot
      vals = c(zlim, scaleFun(plotVar[thisI]))
      vals = vals-zlim[1]
      vals = vals/(zlim[2] - zlim[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      if(is.na(vals[3])) {
        thisDensity = crosshatchNADensity
      } else {
        thisDensity = NULL
      }
      
      if(!project)
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(areaPolys[[x]]@coords, col=col, density=thisDensity), list(...)))})
      else
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(myProjection(areaPolys[[x]]@coords), col=col, density=thisDensity), list(...)))})
    }
    
  }
  
  sapply(1:length(polys), plotArea)
  
  if(!is.null(plotVar) && addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(scaleFunInverse(zlim), n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=zlim
    legendArgs$nlevel=length(cols)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=cols
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# for plotting administration area names (for use with plotMapDat)
# varAreas: the names of the areas we wish to plot
# mapDat: spatial polygon file containing information on the areas
# ...: additional arguments to the text function
addMapLabels = function(varAreas, mapDat, offsets=NULL, areaVarName, ...) {
  
  # determine the names of the mapDat areas
  regionNames = as.character(mapDat@data[[areaVarName]])
  
  # plot map labels
  xs = coordinates(mapDat)
  includeI = match(varAreas, regionNames)
  xs = xs[includeI,]
  
  if(!is.null(offsets)) {
    xs = xs + offsets
  }
  
  text(xs, as.character(varAreas), ...)
  
  invisible(NULL)
}

# color scales ----
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
    else
      scale_colour_continuous_diverging(n_interp=n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
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
      
      if(propUp >= propDown && totalMissingColors > 0)
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

makeGreenSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev, n_interp=n)
}

makeBlueGreenYellowSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, rev=rev)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n, rev=rev)
}

makeRedYellowBlueColors = function(n, ggplot=FALSE) {
  if(!ggplot)
    divergingx_hcl(n, palette="RdYlBu")
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

# color scale utility functions ----
# combine two color scale functions (that return vector of colors given number of colors), 
# given the number of colors in the scale desired
combineTwoScales = function(n, scale1, scale2, args1, args2) {
  if(n %% 2 == 0) {
    n1 = n2 = n/2
  } else {
    n1 = ceiling(n/2)
    n2 = floor(n/2)
  }
  
  c(do.call(scale1, c(args1, list(n=n1))), 
    do.call(scale2, c(args2, list(n=n2))))
}

# convert a single color sequential scale into a diverging scale
makeDivergingScale = function(n, scale, ...) {
  do.call("combineTwoScales", list(n=n, scale1=scale, scale2=scale, args1=list(...), args2=list(...)))
}

# centers a color scale at its midpoint. Returns vector of the centered color scale. 
# Useful when using diverging scales centered at 0 for data with asymmetric range
# colScale a function taking 'n' as input and return a color scale centered in the middle
# n: number of colors
# vals: values to plot (unscaled)
# valRange: the range of vals you want to plot
# center: the center of the range of values (this is the center of the diverging color scale)
# colScale: function return a vector of colors taking inputs 'n' and '...', e.g. makeRedBlueDivergingColors
# scaleFun: scales the color scale to be on, e.g., log scale
# compressed: rather than cutting off a color scale to make the center at the desired point, 
#             compress the scale on the shorter side so it changes color faster.
# ...: other arguments to colScale, such as 'rev'
# Example:
# # construct data:
# testX = exp(rnorm(100))
# # construct centered color scale on log scale
# test = centerColorScale(64, testX, center=1, colScale=makeRedBlueDivergingColors, scaleFun=log)
# # get the colors associated with each value
# testCols = getColorsFromScale(testX, center=1, cols=test, scaleFun=log)
# # plot the result
# plot(testX, col=testCols, pch=19)
centerColorScale = function(n, vals=NULL, valRange=NULL, center, colScale, scaleFun=function(x) {x}, 
                            compressed=FALSE, ...) {
  require("colorspace")
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  valRange = scaleFun(valRange)
  center = scaleFun(center)
  
  propUp = (valRange[2] - center) / diff(valRange)
  propDown = 1 - propUp
  totalColors = ceiling(2 * max(propUp, propDown) * n)
  tempColors = do.call(colScale, c(list(totalColors), list(...)))
  totalMissingColors = totalColors - n
  
  if(!compressed) {
    if(propUp >= propDown && totalMissingColors > 0)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
  } else {
    if(propUp >= propDown && totalMissingColors > 0) {
      lowIndices = round(seq(1, round(totalColors/2), l=round(totalColors/2)-totalMissingColors))
      highIndices = round(seq(round(totalColors/2)+1, totalColors, l=round(totalColors/2)))
      tempColors[c(lowIndices, highIndices)]
    }
    else {
      lowIndices = round(seq(1, round(totalColors/2), l=round(totalColors/2)))
      highIndices = round(seq(round(totalColors/2)+1, totalColors, l=round(totalColors/2)-totalMissingColors))
      tempColors[c(lowIndices, highIndices)]
    }
  }
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

# axis scale ticks ---- 
# x: numbers on probability scale
getLogitScaleTicks = function(x, nint=3, add.5=FALSE) {
  minX = min(x)
  maxX = max(x)
  
  # check if range contains .5. If so, make sure .5 is in range of data
  if(minX <= .5 && maxX >= .5) {
    x = c(x, .5)
  }
  
  # first generate ticks below .5, then flip about .5
  rx = x
  rx[rx > .5] = 1 - rx[rx > .5]
  
  # now add log scale ticks
  lowerTicks = axisTicks(range(log10(rx)), log=TRUE, nint=nint)
  upperTicks = rev(1 - lowerTicks)
  
  if(add.5) {
    c(lowerTicks, .5, upperTicks)
  } else {
    c(lowerTicks, upperTicks)
  }
}

# x: numbers on positive real scale
getLogScaleTicks = function(x, nint=5) {
  axisTicks(range(log10(x)), log=TRUE, nint=nint)
}

## Functions from jamba package ----
getCustomScaleTicks = function(usr, scaleFun=sqrt, nint=5, log=FALSE) {
  # axisTicks(range(scaleFun(x)), log=log, nint=nint)
  axp <- unlist(.axisPars(range(scaleFun(usr)), log = log, nintLog = nint), 
                use.names = FALSE)
  .Call(C_R_CreateAtVector, axp, usr, nint, log)
}


#' Determine square root axis tick mark positions
#'
#' Determine square root axis tick mark positions
#'
#' This function calculates positions for tick marks for data
#' that has been transformed with `sqrt()`, specifically a directional
#' transformation like `sqrt(abs(x)) * sign(x)`.
#'
#' The main goal of this function is to provide reasonably placed
#' tick marks using integer values.
#'
#' @return
#' Invisibly returns a numeric vector of axis tick positions,
#' named by the display label.
#' The axis values are in square root space while the labels represent
#' the normal space values.
#'
#' @family jam plot functions
#'
#' @param side integer value indicating the axis position, as used
#'    by `axis()`, 1=bottom, 2=left, 3=top, 4=right.
#' @param x optional numeric vector representing the numeric range
#'    to be labeled.
#' @param pretty.n numeric value indicating the number of desired
#'    tick marks, passed to `pretty()`.
#' @param u5.bias numeric value passed to `pretty()` to influence the
#'    frequency of intermediate tick marks.
#' @param big.mark character value passed to `format()` which helps
#'    visually distinguish numbers larger than 1000.
#' @param plot logical indicating whether to plot the axis tick
#'    marks and labels.
#' @param las,cex.axis numeric values passed to `axis()` when drawing
#'    the axis, by default `las=2` plots labels rotated
#'    perpendicular to the axis.
#' @param ... additional parameters are passed to `pretty()`.
#'
#' @export
sqrtAxis <- function
(side=1,
 x=NULL,
 pretty.n=10,
 u5.bias=0,
 big.mark=",",
 plot=TRUE,
 las=2,
 cex.axis=0.6,
 ...)
{
  ## Purpose is to generate a set of tick marks for sqrt
  ## transformed data axes.  It assumes data is already sqrt-transformed,
  ## and that negative values have been treated like:
  ## sqrt(abs(x))*sign(x)
  if (length(side) > 2) {
    x <- side;
    side <- 0;
  }
  if (length(side) == 0) {
    side <- 0;
  }
  if (1 %in% side) {
    xRange <- par("usr")[1:2];
  } else if (2 %in% side) {
    xRange <- par("usr")[3:4];
  } else if (length(x) > 0) {
    xRange <- range(x, na.rm=TRUE);
  }
  
  subdivideSqrt <- function(atPretty1, n=pretty.n, ...) {
    ## Purpose is to take x in form of 0,x1,
    ## and subdivide using pretty()
    atPretty1a <- unique(sort(abs(atPretty1)));
    atPretty1b <- tail(atPretty1a, -2);
    atPretty2a <- pretty(head(atPretty1a,2), n=n, ...);
    return(unique(sort(c(atPretty2a, atPretty1b))));
  }
  
  ## Determine tick positions
  nSubFactor <- 2.44;
  
  atPretty1 <- pretty(xRange^2*sign(xRange),
                      u5.bias=u5.bias,
                      n=(pretty.n)^(1/nSubFactor));
  
  atPretty1old <- atPretty1;
  while (length(atPretty1) <= pretty.n) {
    atPretty1new <- subdivideSqrt(atPretty1,
                                  n=noiseFloor(minimum=2, (pretty.n)^(1/nSubFactor)));
    atPretty1 <- atPretty1new[atPretty1new <= max(abs(xRange^2))];
    atPretty1old <- atPretty1;
  }
  atPretty3 <- unique(sort(
    rep(atPretty1,
        each=length(unique(sign(xRange)))) * sign(xRange)));
  atPretty <- atPretty3[
    (atPretty3 >= head(xRange,1)^2*sign(head(xRange,1)) &
       atPretty3 <= tail(xRange, 1)^2*sign(tail(xRange, 1)))];
  
  xLabel <- sapply(atPretty, function(i){
    format(i,
           trim=TRUE,
           digits=2,
           big.mark=big.mark);
  });
  ## Transform to square root space
  atSqrt <- sqrt(abs(atPretty))*sign(atPretty);
  if (plot) {
    axis(side=side,
         at=atSqrt,
         labels=xLabel,
         las=las,
         cex.axis=cex.axis,
         ...);
  }
  
  invisible(nameVector(atSqrt, xLabel));
}

#' Apply noise floor and ceiling to numeric vector
#'
#' Apply noise floor and ceiling to numeric vector
#'
#' A noise floor is useful when detected numeric values are sometimes below
#' a clear noise threshold, and where some downstream ratio may be calculated
#' using these values. Applying a noise floor ensures the ratios and not
#' artificially higher, especially in cases where the values involved are
#' least reliable. This procedure is expected to produce more conservative
#' and appropriate ratios in that scenario.
#'
#' A ceiling is similar, values above the ceiling are set to the ceiling,
#' which is practical when values above a certain threshold are conceptually
#' similar to those at the threshold. One clear example is plotting
#' `-log10(Pvalue)` when the range of P-values might approach 1e-1000.
#' In this case, setting a ceiling of 50 conceptually equates P-values
#' below 1e-50, while also restricting the axis range of a plot.
#'
#' The ability to set values at the floor to a different value, using
#' `newValue` different from `minimum`, is intended to allow separation
#' of numeric values from the floor for illustrative purposes.
#'
#' @return
#' A numeric vector or matrix, matching the input type `x` where numeric
#' values are fixed to the `minimum` and `ceiling` values as defined
#' by `newValue` and `newCeiling`, respectively.
#'
#' @family jam numeric functions
#'
#' @param x numeric vector or matrix
#' @param minimum numeric floor value
#' @param newValue numeric, by default the same as the floor value. Sometimes
#'    it can be useful to define a different value, one example is to define
#'    values as NA, or another distinct number away from the floor.
#' @param adjustNA logical whether to change NA values to the newValue.
#' @param ceiling numeric value, optionally a ceiling. If defined, then values
#'    above the ceiling value are set to newCeiling.
#' @param newCeiling numeric value when ceiling is defined, values above the
#'    ceiling are set to this numeric value.
#' @param ... additional parameters are ignored.
#'
#' @examples
#' # start with some random data
#' n <- 2000;
#' x1 <- rnorm(n);
#' y1 <- rnorm(n);
#'
#' # apply noise floor and ceiling
#' x2 <- noiseFloor(x1, minimum=-2, ceiling=2);
#' y2 <- noiseFloor(y1, minimum=-2, ceiling=2);
#'
#' # apply noise floor and ceiling with custom replacement values
#' xm <- cbind(x=x1, y=y1);
#' xm3 <- noiseFloor(xm,
#'    minimum=-2, newValue=-3,
#'    ceiling=2, newCeiling=3);
#'
#' parMfrow <- par("mfrow");
#' par("mfrow"=c(2,2));
#' plotSmoothScatter(x1, y1);
#' plotSmoothScatter(x2, y2);
#' plotSmoothScatter(xm3);
#' par("mfrow"=parMfrow);
#'
#' @export
noiseFloor <- function
(x,
 minimum=0,
 newValue=minimum,
 adjustNA=FALSE,
 ceiling=NULL,
 newCeiling=ceiling,
 ...)
{
  ## Purpose is to apply a noise floor, that is, to set all values
  ## to be at least 'minimum' amount.
  ## This function performs no scaling or normalization.
  if (length(x) == 0) {
    return(x);
  }
  if (length(minimum) > 0) {
    if (adjustNA) {
      x[is.na(x) | (!is.na(x) & x < minimum)] <- newValue;
    } else {
      x[!is.na(x) & x < minimum] <- newValue;
    }
  }
  if (length(ceiling) > 0) {
    x[!is.na(x) & x > ceiling] <- newCeiling;
  }
  return(x);
}

#' assign unique names for a vector
#'
#' assign unique names for a vector
#'
#' This function assigns unique names to a vector, if necessary it runs
#' \code{\link{makeNames}} to create unique names. It differs from
#' \code{\link[stats]{setNames}} in that it ensures names are unique,
#' and when no names are supplied, it uses the vector itself to define
#' names. It is helpful to run this function inside an \code{\link[base]{lapply}}
#' function call, which by default maintains names, but does not assign
#' names if the input data did not already have them.
#'
#' When used with a data.frame, it is particularly convenient to pull out
#' a named vector of values. For example, log2 fold changes by gene, where
#' the gene symbols are the name of the vector.
#'
#' \code{nameVector(genedata[,c("Gene","log2FC")])}
#'
#' @return vector with names defined
#'
#' @family jam string functions
#'
#' @param x vector input, or data.frame, matrix, or tibble with two columns,
#'    the second column is used to name values in the first column.
#' @param y NULL or character vector of names. If NULL then x is used.
#'    Note that y is recycled to the length of x, prior to being sent
#'    to the makeNamesFunc.
#'    In fringe cases, y can be a matrix, data.frame, or tibble, in which
#'    case \code{\link{pasteByRow}} will be used to create a character string
#'    to be used for vector names. Note this case is activated only when x
#'    is not a two column matrix, data.frame, or tibble.
#' @param makeNamesFunc function to make names unique, by default
#'    \code{\link{makeNames}} which ensures names are unique.
#' @param ... passed to \code{\link{makeNamesFunc}}, or to
#'    \code{\link{pasteByRow}} if y is a two column data.frame, matrix, or
#'    tibble. Thus, \code{sep} can be defined here as a delimiter between
#'    column values.
#'
#' @examples
#' # it generally just creates names from the vector values
#' nameVector(LETTERS[1:5]);
#'
#' # if values are replicated, the makeNames() function makes them unique
#' V <- rep(LETTERS[1:5], each=3);
#' nameVector(V);
#'
#' # for a two-column data.frame, it creates a named vector using
#' # the values in the first column, and names in the second column.
#' df <- data.frame(seq_along(V), V);
#' df;
#' nameVector(df);
#'
#' # Lastly, admittedly a fringe case, it can take a multi-column data.frame
#' # to generate labels:
#' nameVector(V, df);
#'
#' @export
nameVector <- function
(x,
 y=NULL,
 makeNamesFunc=makeNames,
 ...)
{
  ## Purpose is to name a vector with its own values,
  ## useful for lapply which only names output if the input
  ## vector has names.
  ##
  ## A neat trick is to use the _v# naming scheme in makeNames to
  ## create unique names based upon a single label, e.g.
  ## set1colors <- nameVector(brewer.pal(15, "Set1"), "Set1");
  ##   Set1_v1   Set1_v2   Set1_v3   Set1_v4   Set1_v5   Set1_v6   Set1_v7   Set1_v8   Set1_v9
  ## "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
  ##
  ## Added bonus, if given a 2-column table, it'll use them as x and y
  if (igrepHas("dataframe", class(x))) {
    x <- as.data.frame(x);
  }
  if (igrepHas("data.frame", class(x)) && ncol(x) == 2) {
    y <- x[[2]];
    x <- x[[1]];
  } else if (igrepHas("matrix", class(x)) && ncol(x) == 2) {
    y <- x[,2];
    x <- x[,1];
  }
  if (length(y) > 0) {
    if (igrepHas("data.frame|matrix", class(y))) {
      ## If given a data.frame use pasteByRow() to create a string
      y <- pasteByRow(y, ...);
    }
    names(x) <- makeNamesFunc(rep(y, length.out=length(x)), ...);
  } else {
    names(x) <- makeNamesFunc(x, ...);
  }
  return(x);
}

#' vector contains any case-insensitive grep match
#'
#' vector contains any case-insensitive grep match
#'
#' This function checks the input vector for any elements matching the
#' grep pattern. The grep is performed case-insensitive (igrep). This function
#' is particularly useful when checking function arguments or object class,
#' where the class(a) might return multiple values, or where the name of
#' the class might be slightly different than expected, e.g. data.frame,
#' data_frame, DataFrame.
#'
#' @param pattern the grep pattern to use with `base::grep()`
#' @param x vector to use in the grep
#' @param ignore.case logical default TRUE, meaning the grep will be performed
#'    in case-insensitive mode.
#' @param minCount integer minimum number of matches required to return TRUE.
#' @param naToBlank logical whether to convert NA to blank, instead of
#'    allowing grep to handle NA values as-is.
#'
#' @return logical indicating whether the grep match criteria were met,
#'    TRUE indicates the grep pattern was present in minCount or more
#'    number of entries.
#'
#' @seealso `base::grep()`
#'
#' @examples
#' a <- c("data.frame","data_frame","tibble","tbl");
#' igrepHas("Data.*Frame", a);
#' igrepHas("matrix", a);
#'
#' @family jam grep functions
#'
#' @export
igrepHas <- function
(pattern, x=NULL, ignore.case=TRUE,
 minCount=1, naToBlank=FALSE,
 ...)
{
  ## Purpose is a quick check for greppable substring, for if() statements
  ##
  ## naToBlank=TRUE will convert NA values to "" prior to running grep
  ##
  ## The special case where minCount is negative (minCount == -1) or larger
  ## than length(x), it will be set to length(x) and therefore
  ## requires all elements of x to meet the grep criteria
  if (minCount < 0 || minCount > length(x)) {
    minCount <- length(x);
  }
  if (length(x) == 0) {
    return(FALSE);
  } else {
    if (naToBlank && any(is.na(x))) {
      x[is.na(x)] <- "";
    }
    length(grep(pattern=pattern,
                x=x,
                ignore.case=ignore.case,
                ...)) >= as.integer(minCount);
  }
}

#' make unique vector names
#'
#' make unique vector names
#'
#' This function extends the basic goal from \code{\link[base]{make.names}}
#' which is intended to make syntactically valid names from a character vector.
#' This makeNames function makes names unique, and offers configurable methods
#' to handle duplicate names. By default, any duplicated entries receive a
#' suffix _v# where # is s running count of entries observed, starting at 1.
#' The \code{\link[base]{make.names}} function, by contrast, renames the
#' second observed entry starting at .1, leaving the original entry
#' unchanged. Optionally, makeNames can rename all entries with a numeric
#' suffix, for consistency.
#'
#' For example:
#' \code{A, A, A, B, B, C}
#' becomes:
#' \code{A_v1, A_v2, A_v3, B_v1, B_v2, C}
#'
#' Also, makeNames always allows "_".
#'
#' This makeNames function is similar to \code{\link[base]{make.unique}}
#' which also converts a vector into a unique vector by adding suffix values,
#' however the \code{\link[base]{make.unique}} function intends to allow
#' repeated operations which recognize duplicated entries and continually
#' increment the suffix number. This makeNames function currently does not
#' handle repeat operations. The recommended approach to workaround having
#' pre-existing versioned names would be to remove suffix values prior to
#' running this function. One small distinction from
#' \code{\link[base]{make.unique}} is that makeNames does version the first
#' entry in a set.
#'
#' @return character vector of unique names
#'
#' @family jam string functions
#'
#' @param x character vector to be used when defining names. All other
#'    vector types will be coerced to character prior to use.
#' @param unique argument which is ignored, included only for
#'    compatibility with `base::make.names`. All results from
#'    `makeNames()` are unique.
#' @param suffix character separator between the original entry and the
#'    version, if necessary.
#' @param renameOnes logical whether to rename single, unduplicated, entries.
#' @param doPadInteger logical whether to pad integer values to a consistent
#'    number of digits, based upon all suffix values needed. This output
#'    allows for more consistent sorting of names. To define a fixed number
#'    of digits, use the useNchar parameter.
#' @param useNchar integer or NULL, number of digits to use when padding
#'    integer values with leading zero, only relevant when usePadInteger=TRUE.
#' @param startN integer number used when numberStyle is "number", this integer
#'    is used for the first entry to be renamed. You can use this value to
#'    make zero-based suffix values, for example.
#' @param numberStyle character style for version numbering
#'    \describe{
#'       \item{"number"}{Use integer numbers to represent each duplicated
#'          entry.}
#'       \item{"letters"}{Use lowercase letters to represent each duplicated
#'          entry. The 27th entry uses the pattern "aa" to represent two
#'          26-base digits. When doPadInteger=TRUE, a zero is still used
#'          to pad the resulting version numbers, again to allow easy sorting
#'          of text values, but also because there is no letter equivalent
#'          for the number zero.
#'          It is usually best to change the suffix to "_" or "" when using
#'          "letters".}
#'       \item{"LETTERS"}{Use uppercase letters to represent each duplicated
#'          entry, with the same rules as applied to "letters".}
#'    }
#' @param renameFirst logical whether to rename the first entry in a set of
#'    duplicated entries. If FALSE then the first entry in a set will not
#'    be versioned, even when renameOnes=TRUE.
#' @param keepNA logical whether to retain NA values using the string "NA".
#'    If keepNA is FALSE, then NA values will remain NA, thus causing some
#'    names to become `<NA>`, which can cause problems with some downstream
#'    functions which assume all names are either NULL or non-NA.
#'
#' @examples
#' V <- rep(LETTERS[1:3], c(2,3,1));
#' makeNames(V);
#' makeNames(V, renameOnes=TRUE);
#' makeNames(V, renameFirst=FALSE);
#' exons <- makeNames(rep("exon", 3), suffix="");
#' makeNames(rep(exons, c(2,3,1)), numberStyle="letters", suffix="");
#'
#' @export
makeNames <- function
(x,
 unique=TRUE,
 suffix="_v",
 renameOnes=FALSE,
 doPadInteger=FALSE,
 startN=1,
 numberStyle=c("number","letters","LETTERS"),
 useNchar=NULL,
 renameFirst=TRUE,
 keepNA=TRUE,
 ...)
{
  ## Purpose is to make unique names without the R mangling that comes
  ## with make.names().
  ## By default, unique entries are not renamed, and entries with two or
  ## more replicates are renamed to NAME_v1, NAME_v2, NAME_v3, etc.
  ##
  ## if renameOnes=TRUE, it will rename singlets to NAME_v1 even if there
  ## is only one entry.
  ##
  ## renameFirst=TRUE will rename each duplicated entry NAME_v1, NAME_v2,
  ## NAME_v3, etc.
  ## renameFirst=FALSE will not rename the first in a set of duplicated
  ## entries, e.g. NAME, NAME_v1, NAME_v2, etc.
  ##
  ## The distinction between renameOnes and renameFirst:
  ## renameOnes=TRUE will rename all singlets and duplicated entries,
  ## starting with the first entry.
  ## renameOnes=FALSE will not rename singlet entries.
  ## renameFirst=TRUE will only rename duplicated entries, starting with
  ## the first entry.
  ## renameFirst=FALSE will not rename the first entry in a set of
  ## duplicated entries.
  ##
  ## the suffix can be changed, e.g. "_r" will name names NAME_r1,
  ## NAME_r2, NAME_r3, etc.
  ##
  ## numberStyle="number" uses integers as the suffix
  ## numberStyle="letters" uses lowercase letters as digits, similar to Excel column names
  ## numberStyle="LETTERS" uses uppercase letters as digits, similar to Excel column names
  ## Be aware that letters can only go to roughly 18,000 entries, given the current implementation
  ## of colNum2excelName
  ##
  ## When useNchar is numeric, it sets doPadInteger=TRUE, and will use at least
  ## that many digits in padding the integer.
  ##
  ##
  ## TODO:
  ## Update logic to be analogous to using make.unique(), which intends
  ## to maintain previous versioning of names without appending deeper
  ## suffices as appropriate.
  ## E.g.     c("",    "",    "",    "_v1", "_v2", "_v3")
  ## becomes  c("_v4", "_v5", "_v6", "_v1", "_v2", "_v3")
  ## instead of
  ##          c("_v1", "_v2", "_v3", "_v1", "_v2", "_v3")
  ## or
  ##          c("_v1_v1", "_v2_v1", "_v3_v1", "_v1_v2", "_v2_v2", "_v3_v2")
  ##
  if (length(x) == 0) {
    return(x);
  }
  numberStyle <- match.arg(numberStyle);
  if (!is.null(useNchar)) {
    useNchar <- as.integer(useNchar);
    doPadInteger=TRUE;
  }
  if (any(c("factor", "ordered") %in% class(x))) {
    x <- as.character(x);
  }
  if (keepNA && any(is.na(x))) {
    x <- rmNA(x,
              naValue="NA");
  }
  
  ## First check for duplicates using anyDuplicated()
  ## version 0.0.35.900, this change speeds assignment
  ## in large vectors when most entries are not duplicated.
  dupes <- duplicated(x);
  
  ## Convert entries to a named count of occurences of each entry
  if (any(dupes)) {
    xSubDupes <- table(x[dupes]) + 1;
    maxCt <- max(c(1,xSubDupes));
    xSubOnes <- setNames(rep(1, sum(!dupes)), x[!dupes]);
    xSub <- c(xSubOnes, xSubDupes);
  } else {
    xSub <- setNames(rep(1, sum(!dupes)), x[!dupes]);
    maxCt <- 1;
  }
  ## version 0.0.34.900 and previous used the method below
  #xSub <- table(as.character(x));
  
  ## Vector of counts to be used
  versionsV <- as.integer(renameFirst):maxCt + startN - 1;
  
  ## If using letters, define the set of letter upfront to save processing
  if (igrepHas("letters", numberStyle)) {
    if (numberStyle %in% "letters") {
      useLetters <- letters[1:26];
      zeroVal <- "A";
    } else {
      useLetters <- LETTERS[1:26];
      zeroVal <- "a";
    }
    num2letters <- colNum2excelName(versionsV,
                                    useLetters=useLetters,
                                    zeroVal=zeroVal,
                                    ...);
    versionsV <- num2letters;
  }
  if (doPadInteger) {
    versionsV <- padInteger(versionsV,
                            useNchar=useNchar,
                            ...);
  }
  
  ## If no duplicated entries
  if (max(xSub) %in% c(-Inf,1)) {
    ## If not renaming the singlet entries, send the same list back
    if (!renameOnes) {
      return(x);
    } else {
      ## If renaming singlets, simply paste the suffix and first entry
      return(paste0(x, suffix, head(versionsV, 1)));
    }
  }
  
  if (renameOnes) {
    xUse <- 1:length(x);
  } else {
    xUse <- (x %in% names(xSub)[xSub > 1]);
  }
  xSub1 <- x[xUse];
  
  ## Preserve the original order
  names(xSub1) <- padInteger(seq_along(xSub1));
  ## Split the vector into a list of vectors, by name
  xSub2 <- split(xSub1, xSub1);
  names(xSub2) <- NULL;
  
  ## Optionally pad the integer to facilitate sorting
  ## Note: This padding only pads integers within each name,
  ## not across all names.
  xSub3 <- lapply(xSub2, function(i){
    versionsV[seq_along(i)];
  });
  
  ## Now simply paste the value, the suffix, and the new version
  xSub1v <- paste0(unlist(xSub2), suffix, unlist(xSub3));
  
  ## Re-order the vector using the original order
  names(xSub1v) <- names(unlist(xSub2));
  xSub1v2 <- xSub1v[names(xSub1)];
  
  ## Assign only the entries we versioned
  x[xUse] <- xSub1v2;
  
  ## Last check for renameFirst=FALSE, in which case we remove the first
  ## versioned entry
  if (!renameFirst) {
    escapeRegex <- function(string){
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string);
    }
    firstVer <- paste0(escapeRegex(paste0(suffix, head(versionsV, 1))), "$");
    x[xUse] <- gsub(firstVer, "", x[xUse]);
  }
  
  return(x);
}