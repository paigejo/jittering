# this script contains some miscellaneous, but useful functions

logit <- function(x) {
  log(x/(1-x))
}

expit2 <- function(x) {
  res = exp(x)/(1+exp(x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}

expit <- function(x) {
  res = 1/(1+exp(-x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}

# inla.seed: seed input to inla.qsample. 0L sets seed intelligently, > 0 sets a 
#            specific seed, < 0 keeps existing RNG
simSPDE = function(coords, nsim=1, mesh=NULL, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, kenya=FALSE, inla.seed=0L) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    if(kenya)
      mesh = getSPDEMeshKenya(coords, doPlot = FALSE)
    else
      mesh = getSPDEMesh(doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = inla.qsample(nsim, Q, seed = inla.seed)
  simDat = as.matrix(A %*% simField)
  
  simDat
}

meanSegmentLength = function(mesh, filterLargerThan=NULL) {
  t.sub = 1:nrow(mesh$graph$tv)
  idx = cbind(mesh$graph$tv[t.sub, c(1:3, 1), drop = FALSE], NA)
  x = mesh$loc[t(idx), 1]
  y = mesh$loc[t(idx), 2]
  indices = 1:4 + rep(seq(from=0, to=length(x)-5, by=5), each=4)
  segx = x[indices]
  segy = y[indices]
  coords = cbind(segx, segy)
  dists = rdist.vec(coords[1:(length(segx) - 1),], coords[2:length(segx),])
  dists = dists[-seq(from=4, to=length(dists), by=4)]
  
  if(!is.null(filterLargerThan))
    dists = dists[dists <= filterLargerThan]
  
  mean(dists)
}
# meanSegmentLength(getSPDEMeshKenya(), 50)
# meanSegmentLength(getSPDEMesh(), 0.01+1e-6)

rpcvar = function(n, alpha=.01, u=1) {
  1 / inla.pc.rprec(n, alpha=alpha, u=u)
}

dpcvar = function(x, alpha=.01, u=1) {
  inla.pc.dprec(1 / x, alpha=alpha, u=u) / x^2
}

qpcvar = function(p, alpha=.01, u=1) {
  1 / inla.pc.qprec(1-p, alpha=alpha, u=u)
}

ppcvar = function(q, alpha=.01, u=1, tol = 1e-10) {
  fun = function(x) {dpcvar(x, alpha=alpha, u=u)}
  integrate(fun, lower = tol, upper=q)$value
}

# generate simulations on rectangular domain (e.g. unit square) from a 
# Matern Gaussian Process with zero mean
genSimsMatern = function(xRange=c(0,1), yRange=c(0,1), n=100, nsim=100, 
                         beta=(xRange[2]-xRange[1])/10, nu=1.5, sigmaSq=1, tauSq=sigmaSq/10) {
  # generate locations of data
  xs = matrix(runif(n*nsim)*(xRange[2] - xRange[1]) + xRange[1], ncol=nsim)
  ys = matrix(runif(n*nsim)*(yRange[2] - yRange[1]) + yRange[1], ncol=nsim)
  
  # simulate from standard normal
  zSims = matrix(rnorm(n*nsim), ncol=nsim)
  errs = matrix(rnorm(n*nsim, sd=sqrt(tauSq)), ncol=nsim)
  
  # generate a simulation
  genSim = function(i) {
    # print progress
    if(i %% 10 == 0)
      print(paste0("iteration ", i, "/", nsim))
    
    # use Cholesky decomposition of covariance matrix to simulate
    Sigma = stationary.cov(cbind(xs[,i], ys[,i]), Covariance="MaternLR", theta=beta, phi=sigmaSq, nu=nu)
    L = t(chol(Sigma))
    L %*% zSims[,i]
  }
  trueMat = sapply(1:nsim, genSim)
  obsMat = trueMat + errs
  
  list(xs=xs, ys=ys, trueMat=trueMat, obsMat=obsMat)
}

# modified version of fields packages Matern function to code LR2011
# parameterization of the Matern covariance
# range = beta (this is the effective range for 99% correlation, roughly)
# all other parameters are the same as in ?Matern from fields package
MaternLR = function (d, range = 1, beta=range, alpha = 1/beta, smoothness = 0.5, nu = smoothness, 
                     phi = 1) {
  Matern(sqrt(8*nu)*d*alpha, range, alpha, smoothness, nu, phi)
}

# construct numerical integration matrix, constructing mx*my aggregation regions by 
# diving xRange and yRange into mx x my grid of regions
# predPts: Ideally, predPoints should be a regular grid of locations, since they are 
# all weighed equally when being aggregated in each region
makeNumericalIntegralMat = function(predPts, xRange=c(-1, 1), yRange=c(-1, 1), mx=3, my=3) {
  # construct aggregation matrix for predictions by testing which prediction locations 
  # are in which aggregation regions
  xRegionGrid = seq(xRange[1], xRange[2], l=mx + 1)[-1]
  yRegionGrid = seq(yRange[1], yRange[2], l=my + 1)[-1]
  xRegion = function(x) {
    match(TRUE, x <= xRegionGrid)
  }
  yRegion = function(y) {
    match(TRUE, y <= yRegionGrid)
  }
  xRegionI = sapply(predPts[,1], xRegion)
  yRegionI = sapply(predPts[,2], yRegion)
  regionI = (yRegionI-1)*mx + xRegionI
  getARow = function(ai) {regionI == ai}
  
  A = t(sapply(1:(mx*my), getARow))
  A = sweep(A, 1, rowSums(A), "/")
  
  A
}

# for computing what administrative 1 regions the given points are in
# project: project to longitude/latitude coordinates
getRegion = function(pts, project=FALSE) {
  getCounty(pts, project)
}

# for computing what general administrative regions the given points are in
# project: project to longitude/latitude coordinates
getRegion2 = function(pts, project=FALSE, mapDat=adm1compressed, nameVar="NAME_1") {
  
  # project pts to lon/lat coordinate system if user specifies
  if(project)
    pts = projNigeria(pts, inverse=TRUE)
  
  regionNames = mapDat@data[[nameVar]]
  
  # make sure county names are consistent for mapDat == adm1
  # regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
  # regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  
  # get region map polygons and set helper function for testing if pts are in the regions
  polys = mapDat@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {in.poly(pts, countyPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  regionID = apply(out, 1, function(vals) {match(1, vals != 0)})
  regionNameVec = regionNames[regionID]
  list(regionID=regionID, regionNames=regionNameVec, multipleRegs=multipleRegs)
}

getRegionRobust = function(pts, mapDat=adm1, regionNameVar="NAME_1") {
  spCoordsLonLat = SpatialPoints(pts, mapDat@proj4string)
  temp = over(spCoordsLonLat, mapDat, returnList=FALSE)
  nas = is.na(temp[[regionNameVar]])
  
  # calculate which distances to admin boundaries for unknown points are closest
  mapDatPolygons = as.SpatialPolygons.PolygonsList(mapDat@polygons, mapDat@proj4string)
  require(geosphere)
  naClosestIDs = sapply(which(nas), function(ind) {dist2Line(spCoordsLonLat[ind], mapDatPolygons)[4]})
  
  # set associated region to the closest one
  regionID = rep(1, nrow(temp))
  regionID[nas] = naClosestIDs
  regionID[!nas] = match(temp[[regionNameVar]][!nas], mapDat[[regionNameVar]])
  
  # get the region names associated with the region IDs
  regionNames = mapDat@data[[regionNameVar]][regionID]
  
  list(regionID=regionID, regionNames=regionNames)
}

# for computing subareas given what areas the given points are in. If a point 
# isn't in any subarea, finds the closest subarea in its area
# project: project to longitude/latitude coordinates
# warningDist: if distance (in degrees) is above this value, produce a warning
getSubareaRobust = function(pts, areas, subareaMapDat=adm2, 
                            subareaNameVar="NAME_2", 
                            areaNameVar="NAME_1", warningDist=.2) {
  
  subareaNames = subareaMapDat@data[[subareaNameVar]]
  areaNames = subareaMapDat@data[[areaNameVar]]
  
  # make sure county names are consistent for mapDat == adm1
  # regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
  # regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  
  # get region map polygons and set helper function for testing if pts are in the regions
  polys = subareaMapDat@polygons
  inSubarea = function(i) {
    subareaPolys = polys[[i]]@Polygons
    inside = sapply(1:length(subareaPolys), function(x) {in.poly(pts, subareaPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inSubarea)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  subareaID = apply(out, 1, function(vals) {match(1, vals != 0)})
  subareaNameVec = subareaNames[subareaID]
  
  # check if any points are in the wrong area or not in any subareas. If there 
  # are, match them to the closest possible subarea within the correct area.
  outAreas = areaNames[match(subareaNameVec, subareaNames)]
  badPts = is.na(subareaNameVec) | (outAreas != areas)
  badPtDists = rep(0, sum(badPts))
  if(any(badPts)) {
    badIs = which(badPts)
    
    for(i in 1:length(badIs)) {
      ind = badIs[i]
      
      # get the point and convert to SpatialPoint
      thisPtSP = SpatialPoints(matrix(pts[ind,], nrow=1), proj4string=subareaMapDat@proj4string)
      
      # get the area associate with this point and the subareas it contains. 
      # Subset the polygons to only be the ones associated with the correct area
      thisArea = areas[ind]
      subareaIs = which(areaNames == thisArea)
      thisSubareas = subareaNames[subareaIs]
      thisPolys = subareaMapDat@polygons[subareaIs]
      
      # determine which polygon is closest
      dists = gDistance(thisPtSP, subareaMapDat, byid=TRUE)[subareaIs]
      closestI = which.min(dists)
      closestSubarea = thisSubareas[closestI]
      badPtDists[i] = dists[closestI]
      
      # update outputs
      subareaNameVec[ind] = closestSubarea
      subareaID[ind] = which(subareaNames == closestSubarea)
      multipleRegs[ind] = FALSE
    }
  }
  
  # produce warnings
  warnIDs = which(badPtDists > warningDist)
  if(length(warnIDs) > 0) {
    warning(paste0("points ", paste(warnIDs, collapse=", "), " had dists ", paste(badPtDists[warnIDs], collapse=", "), " from nearest subareas"))
  }
  
  list(regionID=subareaID, regionNames=subareaNameVec, multipleRegs=multipleRegs)
}

# for computing what administrative regions the given points are in
# project: project to longitude/latitude coordinates
getProvince = function(pts, project=FALSE) {
  out = getCounty(pts, project)
  countyToRegion(out)
}

# for computing what administrative regions the given points are in
# project: project to longitude/latitude coordinates
getCounty = function(pts, project=FALSE) {
  out = getConstituency(pts, project)
  constituencies = out$constituencyNames
  constituencyToCounty(constituencies)
}

getUrbanRural = function(utmGrid) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("../U5MR/Kenya2014Pop/pop.RData")
  load("../U5MR/adminMapData.RData")
  load("../U5MR/Kenya2014Pop/pop.RData")
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid, FALSE)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  threshes = setThresholds()
  popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
  badPoints = which(sapply(popThreshes, length) > 1)
  if(length(badPoints) >= 1) {
    for(i in 1:length(badPoints)) {
      warnings(paste0("Point ", badPoints[i], " not in Kenya"))
      popThreshes[[badPoints[i]]] = NA
    }
  }
  popThreshes = unlist(popThreshes)
  urban = newPop$popOrig > popThreshes
  
  urban
}

# subareas: Admin2 names
# areas: Admin1 names
getMICSstratumNigeria = function(subareas, areas) {
  # load table giving what Admin2 areas are in what senatorial district and make 
  # sure the Admin2 names match the GADM names in Kano and Lagos
  adm2ToSen = read.csv2("data/admin2ToSen.csv")
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Kano Municipal"] = "Kano"
  adm2ToSen$admin2RefName[(adm2ToSen$admin2RefName == "Nasarawa") & (adm2ToSen$admin1Name_en == "Kano")] = "Nasarawa,Kano"
  
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Lagos Mainland"] = "Mainland"
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Ajeromi-Ifelodun"] = "Ajeromi/Ifelodun"
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Amuwo-Odofin"] = "Amuwo Odofin"
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Ibeju/Lekki"] = "Ibeju-Lekki"
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Ifako-Ijaye"] = "Ifako/Ijaye"
  adm2ToSen$admin2RefName[adm2ToSen$admin2RefName == "Oshodi-Isolo"] = "Oshodi/Isolo"
  adm2ToSen$admin2RefName[(adm2ToSen$admin2RefName == "Surulere") & (adm2ToSen$admin1Name_en == "Lagos")] = "Surulere,Lagos"
  
  kanoSens = c("Kano South", "Kano Central", "Kano North")
  lagosSens = c("Lagos West", "Lagos Central", "Lagos East")
  allSens = c(kanoSens, lagosSens)
  
  # match Admin2 area names to row of adm2ToSen table to get associated 
  # senatorial district. Stratum is Admin1 area unless Admin2 area is in Kano
  rowIs = match(subareas, adm2ToSen$admin2RefName)
  senDistrict = adm2ToSen$SenDist_en[rowIs]
  goodSen = senDistrict %in% allSens
  stratumMICS = areas
  stratumMICS[goodSen] = senDistrict[goodSen]
  
  stratumMICS
}

# set thresholds within each county based on percent population urban
setThresholds = function() {
  require(raster)
  
  load("../U5MR/kenyaPopProj.RData")
  load("../U5MR/adminMapData.RData")
  load("../U5MR/poppc.RData")
  
  getCountyThresh = function(countyName) {
    # if Nairobi or Mombasa, always urban
    if((countyName == "Nairobi") || (countyName == "Mombasa"))
      return(-Inf)
    
    # do the setup
    thisCounty = as.character(kenyaPop$admin1) == countyName
    thisPop = kenyaPop$popOrig[thisCounty]
    thisTot = sum(thisPop)
    pctUrb = poppc$pctUrb[poppc$County == countyName]/100
    pctRural = 1 - pctUrb
    
    # objective function to minimize
    # objFun = function(thresh) {
    #   curPctUrb = sum(thisPop[thisPop > thresh])/thisTot
    #   (curPctUrb - pctUrb)^2
    # }
    
    # do optimization
    # out = optim(10, objFun)
    # thresh = out$par
    # out = optimize(objFun, c(.01, 50000))
    # thresh = out$par
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    thresh = sortedPop[threshI]
    
    # print(paste0("pctUrb: ", pctUrb, "; resPctUrb: ", sum(thisPop[thisPop > thresh])/thisTot, "; thresh: ", thresh, "; obj: ", out$objective))
    thresh
  }
  
  # compute threshold for each county
  counties = poppc$County
  threshes = sapply(counties, getCountyThresh)
  
  list(counties=counties, threshes=threshes)
}

# generate the population density surface along with urbanicity estimates. 
# If any constituencies have no pixel centroids in them, they are added as 
# separate areal units in the pixellated grid.
# delta, mean.neighbor: argument passed to fields.rdist.near via getConstituency
# poppcon: if not NULL, renormalizes population grid values to have these 
#          populations per constituency.
# mapDat: SpatialPolygons object
# polygonSubsetI: index in mapDat 
makeInterpPopGrid = function(kmRes=5, adjustPopSurface=FALSE, targetPop=c("children", "women"), 
                             mean.neighbor=50, delta=.1, conMap=adm2, poppcon=NULL, 
                             mapDat=NULL, polygonSubsetI=NULL) {
  thresholdUrbanBy = ifelse(is.null(poppcon), "county", "constituency")
  
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("../U5MR/Kenya2014Pop/pop.RData")
  load("../U5MR/lims.RData")
  load(paste0(globalDirectory, "adminMapData.RData"))
  kenyaMap = adm0
  countyMap = adm1
  constituencyMap = conMap
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  
  if(!is.null(polygonSubsetI)) {
    # get range of the grid that we actually need
    temp = mapDat@polygons[[polygonSubsetI]]
    allPolygons = temp@Polygons
    eastNorthCoords = do.call("rbind", lapply(1:length(allPolygons), function(i) {projKenya(allPolygons[[i]]@coords)}))
    eastSubRange = range(eastNorthCoords[,1])
    northSubRange = range(eastNorthCoords[,2])
    
    # subset grid to the range we need
    eastGrid = eastGrid[eastGrid >= eastSubRange[1]]
    eastGrid = eastGrid[eastGrid <= eastSubRange[2]]
    northGrid = northGrid[northGrid >= northSubRange[1]]
    northGrid = northGrid[northGrid <= northSubRange[2]]
  }
  
  utmGrid = matrix(make.surface.grid(list(east=eastGrid, north=northGrid)), ncol=2)
  
  # project coordinates into lat/lon
  if(length(utmGrid) > 0) {
    lonLatGrid = matrix(projKenya(utmGrid, inverse=TRUE), ncol=2)
  } else {
    warning(paste0("no grid cell centroid are in the areas of interest. ", 
                   "Integration grid will be composed entirely of custom ", 
                   "integration points at the centroids of the areas of interest"))
    lonLatGrid = utmGrid
  }
  
  # subset grid so it's in Kenya
  polys = kenyaMap@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(lonLatGrid, kenyaPoly)
  utmGrid = matrix(utmGrid[inKenya,], ncol=2)
  lonLatGrid = matrix(lonLatGrid[inKenya,], ncol=2)
  
  # compute counties associated with locations
  # counties = getRegion(lonLatGrid)$regionNames
  # constituencies = getConstituency(lonLatGrid, countyNames=counties)$constituencyNames
  # counties = getCounty(lonLatGrid)
  # provinces = getProvince(lonLatGrid)
  if(length(lonLatGrid) > 0) {
    constituencies = getConstituency(lonLatGrid, mean.neighbor=mean.neighbor, delta=delta)$constituencyNames
    counties = constituencyToCounty(constituencies)
    provinces = countyToRegion(counties)
  } else {
    constituencies = character(0)
    counties = character(0)
    provinces = character(0)
  }
  
  # check to make sure every constituency has at least 2 pixels
  constituenciesFactor = factor(constituencies, levels=sort(unique(conMap@data$CONSTITUEN)))
  if(length(lonLatGrid) > 0) {
    out = aggregate(constituencies, by=list(constituency=constituenciesFactor), FUN=length, drop=FALSE)
  } else {
    out = data.frame(constituency=sort(unique(conMap@data$CONSTITUEN)), 
                     x=NA)
  }
  noPixels = is.na(out$x)
  onePixel = out$x == 1
  onePixel[is.na(onePixel)] = FALSE
  onePixelNames = out$constituency[onePixel]
  badConstituencies = noPixels | onePixel
  badConstituencyNames = as.character(out$constituency[badConstituencies])
  
  if(any(badConstituencies)) {
    if(is.null(poppcon)) {
      # make sure we know the population of each constituency
      stop("If any constituency has < 2 pixels, must specify poppcon")
    } else {
      # get centroids of the constituencies (or it's single pixel coordinates)
      thisSpatialPolyList = as.SpatialPolygons.PolygonsList(conMap@polygons)
      centroidsLonLat = matrix(ncol=2, nrow=nrow(conMap@data))
      
      for(i in 1:nrow(conMap@data)) {
        thisCon = conMap@data$CONSTITUEN[i]
        if(thisCon %in% onePixelNames) {
          thisCentroid = lonLatGrid[constituencies == thisCon,]
        } else {
          thisCentroid = coordinates(thisSpatialPolyList[i])
        }
        
        centroidsLonLat[i,] = thisCentroid
      }
      
      # sort to match results of aggregate (alphabetical order)
      sortI = sort(conMap@data$CONSTITUEN, index.return=TRUE)$ix
      centroidsLonLat = centroidsLonLat[sortI,]
      
      # remove the one pixel for constituencies with only one pixel 
      # (we will add it in again later, twice: both urban and rural)
      onePixel = which(constituencies %in% onePixelNames)
      if(length(onePixel) > 0) {
        lonLatGrid = lonLatGrid[-onePixel,]
        utmGrid = utmGrid[-onePixel,]
        constituencies = constituencies[-onePixel]
        counties = counties[-onePixel]
        provinces = provinces[-onePixel]
      }
      
      # add centroids of only the bad constituencies
      centroidsLonLat = centroidsLonLat[badConstituencies,]
      
      # convert to east/north
      centroidsEastNorth = projKenya(centroidsLonLat[,1], centroidsLonLat[,2])
      
      # only add centroid if bad constituencies have any population in the stratum
      hasUrbanPop = (poppcon$popUrb > 0)[badConstituencies]
      hasRuralPop = (poppcon$popRur > 0)[badConstituencies]
      # centroidsLonLatUrban = centroidsLonLat[hasUrbanPop,]
      # centroidsLonLatRural = centroidsLonLat[hasRuralPop,]
      
      # add centroids to the matrices of pixellated grid coordinates. 
      # Add them twice: once for urban, once for rural
      lonLatGrid = rbind(lonLatGrid, centroidsLonLat[hasUrbanPop,], centroidsLonLat[hasRuralPop,])
      utmGrid = rbind(utmGrid, centroidsEastNorth[hasUrbanPop,], centroidsEastNorth[hasRuralPop,])
      
      # add associated consituencies, counties, provinces to respective vectors
      newConstituencies = c(badConstituencyNames[hasUrbanPop], badConstituencyNames[hasRuralPop])
      newCounties = constituencyToCounty(newConstituencies)
      newProvinces = countyToRegion(newCounties)
      constituencies = c(constituencies, newConstituencies)
      counties = c(counties, newCounties)
      provinces = c(provinces, newProvinces)
    }
  }
  
  # get population density at those coordinates
  if(!PROJ6) {
    interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  } else {
    proj4string(pop) = CRS(SRS_string="EPSG:4326")
    pop@file@name = "/home/ahomei/j/johnpai/git/U5MR/Kenya2014Pop/worldpop_total_1y_2014_00_00.tif"
    interpPopVals = extract(pop, SpatialPoints(lonLatGrid, proj4string=CRS(SRS_string="EPSG:4326")), method="bilinear")
    
    # interpPopVals = extract(pop, SpatialPoints(lonLatGrid, proj4string=CRS("+init=epsg:4326")), method="bilinear")
  }
  
  if(any(badConstituencies)) {
    # make sure population densities in the bad constituencies 
    # are slightly different so one will be classified as urban and 
    # the other as rural
    nUnits = length(interpPopVals)
    nNewRural = sum(hasRuralPop)
    if(nNewRural >= 1) {
      interpPopVals[(nUnits-nNewRural + 1):nUnits] = interpPopVals[(nUnits-nNewRural + 1):nUnits] / 2
    }
  }
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, region=provinces, admin1=counties, admin2=constituencies))
  if(thresholdUrbanBy == "county") {
    threshes = setThresholds()
    popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
    urban = newPop$popOrig >= unlist(popThreshes)
    newPop$urban = urban
  } else {
    if(is.null(poppcon)) {
      stop("if thresholdUrbanBy == 'constituency' must provide poppcon")
    }
    tempCon = poppcon
    tempPop = newPop
    names(tempCon) = c("subarea", "area", "popUrb", "popRur", "popTotal")
    names(tempPop)[c(3, 5:6)] = c("pop", "area", "subareas")
    allSubareas = sort(unique(tempPop$subarea))
    tempCon = tempCon[tempCon$subarea %in% allSubareas,]
    threshes = SUMMER::setThresholdsSubarea(tempPop, poppsub = tempCon)
    popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$subareas == newPop$admin2[i]]})
    urban = newPop$popOrig >= unlist(popThreshes)
    newPop$urban = urban
  }
  
  newPop$east = utmGrid[,1]
  newPop$north = utmGrid[,2]
  
  # if necessary, renormalize population values within constituencies crossed with 
  # urban/rural to be the correct value
  if(!is.null(poppcon)) {
    for(i in 1:nrow(poppcon)) {
      thisCon = poppcon$Constituency[i]
      substratumUrban = (constituencies == thisCon) & urban
      substratumRural = (constituencies == thisCon) & !urban
      factorUrban = poppcon$popUrb[i] / sum(newPop$popOrig[substratumUrban])
      factorRural = poppcon$popRur[i] / sum(newPop$popOrig[substratumRural])
      if(is.finite(factorUrban)) {
        newPop$popOrig[substratumUrban] = newPop$popOrig[substratumUrban] * factorUrban
      } else if(!is.nan(factorUrban) && factorUrban == Inf) {
        newPop$popOrig[substratumUrban] = poppcon$popUrb[i] / sum(substratumUrban)
      }
      if(is.finite(factorRural)) {
        newPop$popOrig[substratumRural] = newPop$popOrig[substratumRural] * factorRural
      } else if(!is.nan(factorRural) && factorRural == Inf) {
        newPop$popOrig[substratumRural] = poppcon$popRur[i] / sum(substratumRural)
      }
    }
  }
  
  # if necessary, adjust the population surface so that it better represents the the child population density 
  # rather than the total population density
  if(adjustPopSurface) {
    targetPop = match.arg(targetPop)
    
    # sort easpc by county name alphabetically
    counties=sort(unique(poppc$County))
    
    sortI = sort(easpc$County, index.return=TRUE)$ix
    temp = easpc[sortI,]
    
    # calculate the number of children per stratum using true total eas and empirical children per ea from census data
    load("../U5MR/empiricalDistributions.RData")
    if(targetPop == "children") {
      # targetPopPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
      #   ecdfExpectation(empiricalDistributions$childrenUrban)
      # targetPopPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
      #   ecdfExpectation(empiricalDistributions$childrenRural)
      targetPopPerStratumUrban = temp$HHUrb * ecdfExpectation(empiricalDistributions$mothersUrban) * 
        ecdfExpectation(empiricalDistributions$childrenUrban)
      targetPopPerStratumRural = temp$HHRur * ecdfExpectation(empiricalDistributions$mothersRural) * 
        ecdfExpectation(empiricalDistributions$childrenRural)
    }
    else {
      # targetPopPerStratumUrban = temp$HHUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
      # targetPopPerStratumRural = temp$HHRur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
      targetPopPerStratumUrban = temp$HHUrb * ecdfExpectation(empiricalDistributions$womenUrban)
      targetPopPerStratumRural = temp$HHRur * ecdfExpectation(empiricalDistributions$womenRural)
    }
    
    # generate 2 47 x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
    getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
      counties = as.character(counties)
      
      mat = t(sapply(counties, function(countyName) {newPop$admin1 == countyName}))
      mat = sweep(mat, 2, newPop$popOrig, "*")
      sweep(mat, 2, newPop$urban == getUrban, "*")
    }
    urbanIntegrationMat = getCountyStratumIntegrationMatrix()
    ruralIntegrationMat = getCountyStratumIntegrationMatrix(FALSE)
    
    # calculate number of people per stratum by integrating the population density surface
    urbanPopulations = rowSums(urbanIntegrationMat)
    ruralPopulations = rowSums(ruralIntegrationMat)
    
    # adjust each row of the integration matrices to get the correct expected number of children per stratum
    urbanIntegrationMat = sweep(urbanIntegrationMat, 1, targetPopPerStratumUrban / urbanPopulations, "*")
    ruralIntegrationMat = sweep(ruralIntegrationMat, 1, targetPopPerStratumRural / ruralPopulations, "*")
    ruralIntegrationMat[ruralPopulations == 0,] = 0
    
    # the column sums of the matrices give the correct modified population densities
    newPop$popOrig = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  }
  
  newPop
}

# generate the population density surface along with urbanicity estimates
# poppcAdjusted: either at the county or constituency level depending on adjustLevel
adjustPopGrid = function(popMat, poppcAdjusted=NULL, adjustLevel=c("County", "Constituency")) {
  adjustLevel = match.arg(adjustLevel)
  areaNamePoppc = "area"
  areaNamePopMat = "area"
  if(adjustLevel == "County") {
    if("County" %in% names(poppcAdjusted)) {
      areaNamePoppc = "County"
    } else if("admin1" %in% names(poppcAdjusted)) {
      areaNamePoppc = "admin1"
    }
    
    if("admin1" %in% names(popMat)) {
      areaNamePopMat = "admin1"
    }
  } else if(adjustLevel == "Constituency") {
    if("Constituency" %in% names(poppcAdjusted)) {
      areaNamePoppc = "Constituency"
    } else if("admin2" %in% names(poppcAdjusted)) {
      areaNamePoppc = "admin2"
    }
    
    if("admin2" %in% names(popMat)) {
      areaNamePopMat = "admin2"
    }
  }
  # sort get population per stratum from poppcAdjusted
  # if("County" %in% names(poppcAdjusted)) {
  #   sortI = sort(unique(poppcAdjusted$County), index.return=TRUE)$ix
  #   areas=poppcAdjusted$County[sortI]
  # } else if("area" %in% names(poppcAdjusted)) {
  #   sortI = sort(unique(poppcAdjusted$area), index.return=TRUE)$ix
  #   areas=poppcAdjusted$area[sortI]
  # }
  sortI = sort(unique(as.character(poppcAdjusted[[areaNamePoppc]])), index.return=TRUE)$ix
  areas=poppcAdjusted[[areaNamePoppc]][sortI]
  
  targetPopPerStratumUrban = poppcAdjusted$popUrb[sortI]
  targetPopPerStratumRural = poppcAdjusted$popRur[sortI]
  targetPopPerStratumTotal = poppcAdjusted$popTotal[sortI]
  
  # generate 2 47 x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
  getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
    areas = as.character(areas)
    
    mat = t(sapply(areas, function(areaName) {
      # if("admin1" %in% names(popMat)) {
      #   popMat$admin1 == areaName
      # } else {
      #   popMat$area == areaName
      # }
      as.character(popMat[[areaNamePopMat]]) == as.character(areaName)
    }))
    if("popOrig" %in% names(popMat)) {
      mat = sweep(mat, 2, popMat$popOrig, "*")
    } else {
      mat = sweep(mat, 2, popMat$pop, "*")
    }
    
    sweep(mat, 2, popMat$urban == getUrban, "*")
  }
  urbanIntegrationMat = getCountyStratumIntegrationMatrix()
  ruralIntegrationMat = getCountyStratumIntegrationMatrix(FALSE)
  
  # calculate number of people per stratum by integrating the population density surface
  urbanPopulations = rowSums(urbanIntegrationMat)
  ruralPopulations = rowSums(ruralIntegrationMat)
  
  # adjust each row of the integration matrices to get the correct expected number of children per stratum
  urbanIntegrationMat = sweep(urbanIntegrationMat, 1, targetPopPerStratumUrban / urbanPopulations, "*")
  ruralIntegrationMat = sweep(ruralIntegrationMat, 1, targetPopPerStratumRural / ruralPopulations, "*")
  urbanIntegrationMat[urbanPopulations == 0,] = 0
  ruralIntegrationMat[ruralPopulations == 0,] = 0
  
  # the column sums of the matrices give the correct modified population densities
  if("popOrig" %in% names(popMat)) {
    popMat$popOrig = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  } else {
    popMat$pop = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  }
  
  popMat
}

# add in binomial variation to the probability sampling matrix
addBinomialVar = function(probMatrix, ns) {
  simulatedObservations = matrix(rbinom(n=length(probMatrix), size=rep(ns, ncol(probMatrix)), prob=c(as.matrix(probMatrix))), nrow=nrow(probMatrix))
  sweep(simulatedObservations, 1, 1/ns, "*")
}

getArea = function(thisMap = NULL, nameVar="NAME_1", sortAreas=FALSE, project=FALSE) {
  require(shapefiles)
  
  getOneArea = function(poly) {
    if("Polygons" %in% slotNames(poly)) {
      areas = sapply(poly@Polygons, function(x) {x@area})
      
      correctPolyI = which.max(areas) # sum or max?
      poly = poly@Polygons[[correctPolyI]]
      # thisCentroid = centroid(poly@coords)
      # allPoints = rbind(thisCentroid, 
      #                   poly@coords)
      allPoints = poly@coords
    } else {
      areas = sapply(poly, function(x) {x@area})
      
      correctPolyI = which.max(areas) # sum or max?
      poly = poly[[correctPolyI]]
      # thisCentroid = centroid(poly@coords)
      # allPoints = rbind(thisCentroid, 
      #                   poly@coords)
      allPoints = poly@coords
    }
    
    if(project) {
      allProjected = projNigeria(allPoints[,1], allPoints[,2], inverse=FALSE)
    } else {
      allProjected = allPoints
    }
    
    # downsample spatial polygons as necessary
    if(nrow(allProjected) >= 6000) {
      # simplify the polygon
      newProjected = allProjected
      tolerance = .05
      while(nrow(newProjected) >= 6000) {
        newProjected = do.call("cbind", dp(list(x=allProjected[,1], y=allProjected[,2]), tolerance=tolerance))
        tolerance = tolerance * 2
      }
      allProjected = newProjected
    }
    
    maps:::area.polygon(allProjected)
  }
  # calculate areas in km^2
  areas = sapply(thisMap@polygons, getOneArea)
  
  areaNames = as.character(thisMap@data[[nameVar]])
  if(sortAreas) {
    sortI = sort(areaNames, index.return=TRUE)$ix
  } else {
    sortI = 1:length(areaNames)
  }
  areas = areas[sortI]
  names(areas) = areaNames[sortI]
  areas
}

getAreaPerObservation = function(areaLevel=c("Region", "County"), dataType=c("ed", "mort"), doLog = TRUE) {
  areaLevel = match.arg(areaLevel)
  areas = getArea(areaLevel)
  
  # load the dataset
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    dat = mort
  }
  else {
    dat = ed
  }
  
  # get area names
  allCounties = as.character(dat$admin1)
  counties = sort(unique(as.character(dat$admin1)))
  allRegions = countyToRegion(allCounties)
  regions = sort(unique(allRegions))
  if(areaLevel == "Region") {
    allAreas = allRegions
    uniqueAreas = regions
  } else if(areaLevel == "County") {
    allAreas = allCounties
    uniqueAreas = counties
  }
  
  # count the number of observations per area
  counts = table(allAreas)
  sortI = sort(names(counts), index.return=TRUE)$ix
  counts = counts[sortI]
  
  if(!doLog) {
    areas / counts
  } else {
    log(areas) - log(counts)
  }
}

getAreaPerObservationTicks = function(areaLevel=c("Region", "County"), dataType=c("ed", "mort")) {
  areaLevel = match.arg(areaLevel)
  dataType = match.arg(dataType)
  
  if(dataType != "ed")
    stop("Only education dataset currently supported")
  
  if(areaLevel == "Region") {
    c(50, 100, 200, 400, 800, 1200)
  } else {
    c(10, 25, 50, 100, 200, 400, 800, 1600, 2400)
  }
}



getRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    # thisCentroid = centroid(poly@coords)
    # allPoints = rbind(thisCentroid, 
    #                   poly@coords)
    allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    # downsample spatial polygons as necessary
    if(nrow(allProjected) >= 6000) {
      # simplify the polygon
      newProjected = allProjected
      tolerance = .05
      while(nrow(newProjected) >= 6000) {
        newProjected = do.call("cbind", dp(list(x=allProjected[,1], y=allProjected[,2]), tolerance=tolerance))
        tolerance = tolerance * 2
      }
      allProjected = newProjected
    }
    max(rdist(allProjected)) / 2
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}

getMeanRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    thisCentroid = centroid(poly@coords)
    allPoints = rbind(thisCentroid,
                      poly@coords)
    # allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    mean(rdist.vec(matrix(allProjected[1,], nrow=nrow(allProjected)-1, ncol=2, byrow=TRUE), allProjected[-1,]))
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}

getPredictionDistance = function(doLog=TRUE, dataType=c("ed", "mort")) {
  # load the pixel/prediction grid
  out = load("../U5MR/popGrid.RData")
  predictionPoints = cbind(popGrid$east, popGrid$north)
  
  # load the dataset
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    dat = mort
  }
  else {
    dat = ed
  }
  observationPoints = cbind(dat$east, dat$north)
  
  # calculate distance of nearest observation to each prediction point
  distances = rdist(observationPoints, predictionPoints)
  distances = apply(distances, 2, min)
  
  if(doLog)
    log(distances)
  else
    distances
}

getPredictionDistanceTicks = function(dataType=c("ed", "mort")) {
  dataType = match.arg(dataType)
  
  if(dataType != "ed")
    stop("Only education dataset currently supported")
  
  c(.1, 1, 5, 10, 20, 40, 80)
}

# draw random numbers from an ecdf object
recdf = function(n, distribution) {
  probs = runif(n)
  # quantile(distribution, probs, type=1) # type=1 signifies inverse ecdf
  if(!("function" %in% class(distribution)) && !is.null(distribution$rfun))
    distribution$rfun(n)
  else
    ecdfQuantile(probs, distribution)
}

# draw random numbers from a composition of ecdf objects. I.e. we want children per EA and 
# we have households/EA, mothers/household, and children/mother
recdfComposed = function(n, distributions) {
  results = rep(1, n)
  for(i in 1:length(distributions)) {
    results = sapply(results, function(x) {sum(recdf(x, distributions[[i]]))})
  }
  
  results
}

# draw random numbers from an ecdf object (this could be implemented efficiently)
decdf = function(x, distribution) {
  distributionKnots = knots(distribution)
  masses = diff(c(0, evalq(y, environment(distribution))))
  i = match(x, distributionKnots)
  if(length(i) == 1 && is.na(i))
    return(0)
  else {
    out = masses[i]
    out[is.na(i)] = 0
    return(0)
  }
}

# get expectation of an ecdf object
ecdfExpectation = function(distribution) {
  distributionKnots = knots(distribution)
  distributionKnots = c(distributionKnots[1] - 1, distributionKnots)
  probs = distribution(distributionKnots[2:length(distributionKnots)]) - distribution(distributionKnots[1:(length(distributionKnots) - 1)])
  sum(distributionKnots[2:length(distributionKnots)] * probs)
}
# ecdfExpectation(empiricalDistributions$households) * ecdfExpectation(empiricalDistributions$mothers) * ecdfExpectation(empiricalDistributions$children)

# Becuase stats:::quantile.ecdf is terribly inefficient...
ecdfQuantile = function(p, distribution) {
  distributionKnots = knots(distribution)
  cumulativeMass = evalq(y, environment(distribution))
  indices= sapply(p, function(x) {match(TRUE, x <= cumulativeMass)})
  distributionKnots[indices]
}

ecdf2edfun = function(distribution) {
  samples = evalq(rep.int(x, diff(c(0, round(nobs * y)))), environment(distribution))
  edfun(samples)
}

# takes the poppc table, containing the proportion of population that is urban and rural in each stratum, and
# adjusts it to be representative of the children in urban and rural areas per stratum based on census data
adjustPopulationPerCountyTable = function(dataType=c("children", "women")) {
  dataType = match.arg(dataType)
  
  # calculate the number of childrenor women per stratum using true total eas and empirical children per ea from census data
  load(paste0(globalDirectory, "empiricalDistributions.RData"))
  if(dataType == "children") {
    # targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
    #   ecdfExpectation(empiricalDistributions$childrenUrban)
    # targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
    #   ecdfExpectation(empiricalDistributions$childrenRural)
    targetPopPerStratumUrban = easpc$HHUrb * ecdfExpectation(empiricalDistributions$mothersUrban) * 
      ecdfExpectation(empiricalDistributions$childrenUrban)
    targetPopPerStratumRural = easpc$HHRur * ecdfExpectation(empiricalDistributions$mothersRural) * 
      ecdfExpectation(empiricalDistributions$childrenRural)
  }
  else {
    # targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
    # targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
    targetPopPerStratumUrban = easpc$HHUrb * ecdfExpectation(empiricalDistributions$womenUrban)
    targetPopPerStratumRural = easpc$HHRur * ecdfExpectation(empiricalDistributions$womenRural)
  }
  
  
  # adjust poppc table to be representative of the number of children per stratum
  newPopTable = poppc
  targetPopPerCounty = targetPopPerStratumUrban + targetPopPerStratumRural
  newPopTable$popUrb = targetPopPerStratumUrban
  newPopTable$popRur = targetPopPerStratumRural
  newPopTable$popTotal = targetPopPerCounty
  newPopTable$pctUrb = newPopTable$popUrb / targetPopPerCounty  * 100
  newPopTable$pctTotal = newPopTable$popTotal/sum(newPopTable$popTotal) * 100
  
  # return results
  newPopTable
}

# adapted from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat has a mean and standard deviation 
# on the logit scale
logitNormMean = function(muSigmaMat, parClust=NULL, logisticApproximation=TRUE, splineApproximation=FALSE, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust) && !splineApproximation) {
      apply(muSigmaMat, 1, logitNormMean, logisticApproximation=logisticApproximation, ...)
    }
    else if(splineApproximation) {
      # not parallel and using spline approximation
      uniqueSigmas = sort(unique(muSigmaMat[,2]))
      muSigmaIndsList = lapply(uniqueSigmas, function(s) {
        inds = which(muSigmaMat[,2] == s)
        list(mu=muSigmaMat[inds,1], inds=inds, sigma=s)
      })
      outList = lapply(muSigmaIndsList, function(l) {
        mus = l$mu
        sigma = l$sigma
        logitNormMeanSplineApprox(mus, sigma, ...)
      })
      outVals = 1:nrow(muSigmaMat)
      for(i in 1:length(muSigmaIndsList)) {
        inds = muSigmaIndsList[[i]]$inds
        outVals[inds] = outList[[i]]$vals
      }
      outVals
    } else {
      parApply(parClust, muSigmaMat, 1, logitNormMean, logisticApproximation=logisticApproximation, ...)
    }
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      expit(mu)
    else {
      if(any(is.na(c(mu, sigma))))
        NA
      else if(!logisticApproximation) {
        # numerically calculate the mean
        fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
        integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        expit(mu / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
}

# Approximates logitNormMean at a single value of sigma and many mus using spline 
# on a logit scale. npts determines number of values of mu in the range of mu 
# over which the monotonic spline function is generated.
# Note: Uses a monotonic cubic spline. See ?splinefun for method="Hyman", and:
# Hyman, J. M. (1983). Accurate monotonicity preserving cubic interpolation. 
# SIAM Journal on Scientific and Statistical Computing, 4, 645–654. 
# doi:10.1137/0904045.
logitNormMeanSplineApprox = function(mus, sigma, npts=250, ...) {
  
  rangeExpitMu = expit(range(mus))
  seqMus = logit(seq(rangeExpitMu[1], rangeExpitMu[2], l=npts))
  
  muSigmaMat = cbind(seqMus, sigma)
  vals = logit(logitNormMean(muSigmaMat, logisticApproximation=FALSE, splineApproximation=FALSE))
  
  spFun = splinefun(seqMus, vals)
  
  outVals = expit(spFun(mus))
  
  list(vals = outVals, fun=spFun, range=logit(rangeExpitMu))
}

getPixelIndex = function(eastNorth, popMat=NULL, clusterAreas=NULL, enforceSameArea=TRUE) {
  
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  
  # construct distance matrix
  distMat = rdist(eastNorth, cbind(popMat$east, popMat$north))
  
  # For each observation location, get closest pixel (that is in the same area, if necessary)
  if(enforceSameArea) {
    if(is.null(clusterAreas))
      stop("clusterAreas must be included if enforceSameArea is set to TRUE")
    
    # set distance between clusters and pixels that are in different areas to be infinite
    if("admin1" %in% names(popMat)) {
      differentArea = outer(clusterAreas, popMat$admin1, FUN=function(x, y) {x != y})
    } else if("area" %in% names(popMat)) {
      differentArea = outer(clusterAreas, popMat$area, FUN=function(x, y) {x != y})
    } else {
      stop("popMat must have either 'admin1' or 'area' as a variable")
    }
    
    distMat[differentArea] = Inf
  }
  
  apply(distMat, 1, which.min)
}

# make sure borders of polygons in mapDat are within that of border
makeInBorder = function(mapDat, border=adm0) {
  out = intersect(mapDat, border)
  
  if(FALSE) {
    plot(border)
    for(i in 1:length(out)) {
      # if(i %% 10 == 0) {
      #   browser()
      # }
      plot(out[i,], add=TRUE)
    }
  }
  
  out
}

myin.poly = function (xd, xp, convex.hull = FALSE, inflation = 1e-07) 
{
  if (convex.hull) {
    xp <- xp[chull(xp), ]
  }
  nd <- as.integer(nrow(xd))
  np <- as.integer(nrow(xp))
  if (convex.hull) {
    xm <- matrix(c(mean(xp[, 1]), mean(xp[, 2])), nrow = np, 
                 ncol = 2, byrow = TRUE)
    xp <- (xp - xm) * (1 + inflation) + xm
  }
  ind <- .Fortran("inpoly", PACKAGE = "fields", nd = as.integer(nd), 
                  as.single(matrix(xd[, 1], ncol=1)), as.single(matrix(xd[, 2], ncol=1)), np = np, as.single(xp[, 
                                                                                                                1]), as.single(matrix(xp[, 2], ncol=1)), ind = as.integer(rep(-1, 
                                                                                                                                                                              nd)))$ind
  as.logical(ind)
}

getPopMatResolution = function(popMat) {
  out = aggregate(popMat$east, by=list(east=popMat$east), FUN=length)
  whichI = which.max(out[,2])
  temp = sort(popMat$north[popMat$east == out[whichI,1]])
  min(diff(temp))
}

# calculate true object size in bytes
trueObjectSize = function(x) {
  length(serialize(x,NULL))
}

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

getCustomScaleTicks = function(usr, scaleFun=sqrt, nint=5, log=FALSE) {
  # axisTicks(range(scaleFun(x)), log=log, nint=nint)
  axp <- unlist(.axisPars(range(scaleFun(usr)), log = log, nintLog = nint), 
                use.names = FALSE)
  .Call(C_R_CreateAtVector, axp, usr, nint, log)
}

# i: either the first index of the two input parameters, or the job index if rev==TRUE
# j: the second index of the two input parameters
# maxJ: maximum of j for each i if rev==FALSE (can be a vector of length maxI)
# rev: if TRUE, performs inverse operation of if rev==FALSE given the job index i
getJobIndices = function(i=1, j=1:100, maxJ=100, rev=FALSE) {
  
  if(rev) {
    if(length(maxJ) == 1) {
      outI = ((i-1) %/% maxJ) + 1
      outJ = ((i-1) %% maxJ) + 1
    } else {
      cumJs = c(0, cumsum(maxJ))
      firstCumJbiggerI = match(TRUE, cumJs >= i)
      outI = firstCumJbiggerI-1
      outJ = i - cumJs[firstCumJbiggerI-1]
    }
    return(cbind(i=outI, j=outJ))
  } else {
    if(length(maxJ) == 1) {
      c(jobIndex=(i-1)*maxJ + j)
    } else {
      if(any(j > maxJ[i])) {
        stop(paste0("maxJ for this i is ", maxJ[i], " but maximum j is ", max(j)))
      }
      startI = sum(c(0, maxJ)[1:i])
      c(jobIndex=startI + j)
    }
  }
}

# iNames: vector of names corresponding to i job index from getJobIndices
# jNames: vector of names corresponding to j job index from getJobIndices
# NOTE: iNames and jNames must be names that uniquely characterize files 
#       with the i and j index respectively
generateJobList = function(workDir="savedOutput/validation/folds/", 
                           filePrefix="scores", fileSuffix=".RData", 
                           iNames=c("Md2", "M_D2", "Mdm2", "M_DM2"), 
                           jNames=paste(paste("fold", 1:20, sep=""), ".", sep=""), 
                           extensiveCheck=FALSE, extensiveCheckNULL=FALSE, 
                           excludeFiles=c("d_", "D_", "dm_", "DM_")) {
  # save current directory to return to later. Set directory to job file locations
  thisDir = getwd()
  setwd(workDir)
  
  # get indices of relevant characters in the file names
  out = system(paste0("ls ", filePrefix, "*", fileSuffix), intern=TRUE)
  
  # remove excluded files if necessary
  if(length(excludeFiles) > 1) {
    badFiles = matrix(sapply(excludeFiles, function(x) {
      grepl(x, out)
    }), ncol=length(excludeFiles))
    badFiles = apply(badFiles, 1, any)
    out = out[!badFiles]
  }
  
  
  # get index of character in x representing the start of the iName or jName. 
  # Assumes x has at most 1 iName or jName
  getIJNameInds = function(x, type=c("i", "j")) {
    type=match.arg(type)
    thisNames = iNames
    if(type == "j") {
      thisNames = jNames
    }
    myGregexpr = function(y, z) {
      gregexpr(y, z, fixed=TRUE)[[1]][1]
    }
    out = which(sapply(thisNames, myGregexpr, z=x) != -1)
    
    if(length(out) > 1) {
      stop(paste("multiple matches for type", type, "and names", paste(thisNames[out], collapse=" "), "in file", x, sep=" "))
    }
    out
  }
  iIs = sapply(out, getIJNameInds, type="i")
  jIs = sapply(out, getIJNameInds, type="j")
  
  # determine whether the relevant file exists or not
  fileExists = matrix(FALSE, nrow=length(iNames), ncol=length(jNames))
  for(i in 1:length(out)) {
    thisFilename = out[i]
    iI = iIs[i]
    jI = jIs[i]
    
    if(!extensiveCheck) {
      fileExists[iI, jI] = TRUE
    } else if(extensiveCheck || extensiveCheckNULL) {
      # check to make sure we can actually load the relevant file
      
      canLoad <<- TRUE
      tmp = tryCatch(load(thisFilename), error= function(e) {canLoad <<- FALSE})
      
      if(canLoad && extensiveCheckNULL) {
        canLoad <<- !is.null(preds)
      }
      
      fileExists[iI, jI] = canLoad
    }
  }
  missingJobInds = (1:length(fileExists))[!c(t(fileExists))]
  
  if(length(missingJobInds) == 0) {
    setwd(thisDir)
    print("No missing jobs")
    return(invsible(NULL))
  }else if(length(missingJobInds) == 1) {
    setwd(thisDir)
    return(missingJobInds)
  } else {
    print(paste0("Total missing job results: ", length(missingJobInds)))
  }
  
  # shorten the string to make it readable:
  
  # add a fake missing job index at the end to make sure we get the final missing jobs
  missingJobInds = c(missingJobInds, length(fileExists) + 100)
  
  # look for groups of consecutive missing job indices
  lastInd = missingJobInds[1]
  firstJobInd = lastInd
  jobList = ""
  for(i in 2:length(missingJobInds)) {
    thisInd = missingJobInds[i]
    thisDiff = thisInd - lastInd
    commaStr = ifelse(firstJobInd == missingJobInds[1], "", ",")
    
    if(thisDiff != 1) {
      # this index and last index aren't in the same group, must add 
      # last index (or its associated group) to the job list string
      
      if(firstJobInd == lastInd) {
        # last index was a singleton
        jobList = paste0(jobList, commaStr, lastInd)
      } else {
        # last index was the last consecutive job in the group
        jobList = paste0(jobList, commaStr, "[", firstJobInd, "-", lastInd, "]")
      }
      
      firstJobInd = thisInd
    } else {
      # this index and last index are in the same group. Don't need to 
      # do anything here
    }
    lastInd = thisInd
  }
  
  setwd(thisDir)
  
  return(jobList)
}

shrinkpdf<-function(pdf,maxsize=5,suffix="_small",verbose=T){  
  require(multicore)  
  wd=getwd()  
  td=paste(tempdir(),"/pdf",sep="")  
  if(!file.exists(td)) dir.create(td)  
  if(verbose) print("Performing initial compression")  
  system(paste("ps2pdf ",pdf," ",td,"/test.pdf",sep=""))  
  setwd(td)  
  system(paste("pdftk ",td,"/test.pdf burst",sep=""))  
  files=list.files(pattern="pg_")  
  sizes=sapply(files,function(x) file.info(x)$size)*0.000001 #get sizes of individual pages  
  toobig=sizes>=maxsize  
  if(verbose)  print(paste("Resizing ",sum(toobig)," pages:  (",paste(files[toobig],collapse=","),")",sep=""))  
  
  mclapply(files[toobig],function(i){  
    system(paste("gs -dBATCH -dTextAlphaBits=4 -dNOPAUSE -r300 -q -sDEVICE=png16m -sOutputFile=",i,".png ",i,sep=""))  
    system(paste("convert -quality 100 -density 300 ",i,".png ",strsplit(i,".",fixed=T)[[1]][1],".pdf ",sep=""))  
    if(verbose) print(paste("Finished page ",i))  
    return()  
  })  
  if(verbose) print("Compiling the final pdf")  
  file.remove("test.pdf")  
  file.remove(list.files(pattern="png"))  
  setwd(wd)  
  browser()
  system(paste('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="',strsplit(pdf,".",fixed=T)[[1]][1],suffix,'.pdf" ',td,"/*.pdf",sep=""))  
  file.remove(list.files(td,full=T))  
  if(verbose) print("Finished!!")  
}

# finds the S3 method for an object.
# examples:
# findMethod(print, Sys.time())
## [1] "print.POSIXct"
findMethod <- function(generic, ...) {
  ch <- deparse(substitute(generic))
  f <- X <- function(x, ...) UseMethod("X")
  for(m in methods(ch)) assign(sub(ch, "X", m, fixed = TRUE), "body<-"(f, value = m))
  X(...)
}

# integrate out Gaussian noise from the square of a logit normal
logitNormSqMean = function(muSigmaMat, ...) 
{
  if (length(muSigmaMat) > 2) {
    apply(muSigmaMat, 1, logitNormSqMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if (sigma == 0) 
      SUMMER::expit(mu)^2
    else {
      if (any(is.na(c(mu, sigma)))) 
        NA
      else {
        fExp <- function(x) exp(stats::plogis(x, log.p = TRUE) * 2 + 
                                  stats::dnorm(x, mean = mu, sd = sigma, log = TRUE))
        stats::integrate(fExp, mu - 10 * sigma, mu + 
                           10 * sigma, abs.tol = 0, ...)$value
      }
    }
  }
}

projNigeria = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  # determine version of PROJ
  ver = rgdal::rgdal_extSoftVersion()
  theseNames = names(ver)
  thisI = which(grepl("PROJ", theseNames))
  PROJ6 <- as.numeric(substr(ver[thisI], 1, 1)) >= 6
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    if(!PROJ6) {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS("+proj=longlat"))
    } else {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS(SRS_string="EPSG:4326"))
    }
    coordsEN = sp::spTransform(lonLatCoords, sp::CRS("+init=epsg:26391 +units=m"))
    
    out = attr(coordsEN, "coords")
    colnames(out) = c("east", "north")
    
    # convert coordinates from m to km
    out = out/1000
  }
  else {
    # from easting/northing coords to lon/lat
    
    # first convert from km to m
    east = lon*1000
    north = lat*1000
    
    coordsEN = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS("+init=epsg:26391 +units=m"))
    if(!PROJ6) {
      lonLatCoords = sp::spTransform(coordsEN, sp::CRS("+proj=longlat"))
    } else {
      lonLatCoords = sp::spTransform(coordsEN, sp::CRS(SRS_string="EPSG:4326"))
    }
    
    out = attr(lonLatCoords, "coords")
    colnames(out) = c("lon", "lat")
  }
  
  out
}

# area: shapefile of the area to transform
projNigeriaArea = function(area, inverse=FALSE) {
  
  # determine version of PROJ
  ver = rgdal::rgdal_extSoftVersion()
  theseNames = names(ver)
  thisI = which(grepl("PROJ", theseNames))
  PROJ6 <- as.numeric(substr(ver[thisI], 1, 1)) >= 6
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    projArea = sp::spTransform(area, CRS("+init=epsg:26391 +units=m"))
    
    # convert coordinates from m to km
    projArea@bbox = projArea@bbox/1000
    projArea@polygons = lapply(projArea@polygons, function(x) {lapply(x@Polygons, function(x) {
      out = x; 
      # out@coords = projNigeria(x@coords)
      out@area = out@area / 1000^2
      out
      })})
  }
  else {
    # from easting/northing coords to lon/lat
    
    # first convert from km to m
    area@bbox = area@bbox*1000
    area@polygons = lapply(area@polygons, function(x) {lapply(x@Polygons, function(x) {
      out = x; 
      # out@coords = projNigeria(x@coords, inverse=TRUE)
      out@area = out@area * 1000^2
      out
    })})
    
    if(!PROJ6) {
      projArea = spTransform(area, CRS("+proj=longlat"))
    } else {
      projArea = spTransform(area, CRS(SRS_string="EPSG:4326"))
    }
  }
  
  # warning("projArea does not project polygon coordinates")
  projArea
}

projNigeriaBBox = function(bbox, inverse=FALSE, nlon=100, nlat=100) {
  lonRange = sort(bbox[1,])
  latRange = sort(bbox[2,])
  
  lonPts = seq(lonRange[1], lonRange[2], l=nlon)
  latPts = seq(latRange[1], latRange[2], l=nlat)
  
  allPointsLL = rbind(cbind(lonPts, latPts[1]), 
                      cbind(lonPts, latPts[length(latPts)]), 
                      cbind(lonPts[1], latPts), 
                      cbind(lonPts[length(lonPts)], latPts))
  
  projPts = projNigeria(allPointsLL, inverse=inverse)
  
  eastRange = range(projPts[,1])
  northRange = range(projPts[,2])
  
  rbind(eastRange, 
        northRange)
}

# this is specifically for Nigeria
adm2ToStratumMICS = function(adm2Names) {
  adm2ToSen = read.csv2("data/admin2ToSen.csv")
  
  # make sure adm2 names in Kano and Lagos match with GADM names
  klNamesTab = c(sort(adm2ToSen$admin2Name_en[adm2ToSen$admin1Name_en == "Kano"]), 
                 sort(adm2ToSen$admin2Name_en[adm2ToSen$admin1Name_en == "Lagos"]))
  klNamesGADM = c(sort(adm2@data$NAME_2[adm2@data$NAME_1 == "Kano"]), 
                  sort(adm2@data$NAME_2[adm2@data$NAME_1 == "Lagos"]))
  adm2ToSen$admin2Name_en[match(klNamesTab, adm2ToSen$admin2Name_en)] = klNamesGADM
  
  # kanoSens = c("Kano South", "Kano Central", "Kano North")
  # lagosSens = c("Lagos West", "Lagos Central", "Lagos East")
  
  # first get the GADM admin1 name from poppsubNGA
  poppsubIs = match(adm2Names, poppsubNGAThresh$subarea)
  strata = poppsubNGAThresh$area[poppsubIs]
  useSen = strata %in% c("Kano", "Lagos")
  
  # then select the senatorial district if need be
  adm2ToSenIs = match(adm2Names[useSen], adm2ToSen$admin2Name_en)
  strata[useSen] = adm2ToSen$SenDist_en[adm2ToSenIs]
  
  strata
}

