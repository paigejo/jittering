# Functions for handling covariate rasters

# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
getDesignMat = function(lonLatCoords, normalized=TRUE, stratifyPop=TRUE, 
                        useThreshPopMat=TRUE, proj=projNigeria) {
  
  # get saved urbanicity and population integration matrix
  if(useThreshPopMat) {
    out = load("savedOutput/global/popMatNGAThresh.RData")
    popMat = popMatNGAThresh
  } else {
    out = load("savedOutput/global/popMatNGA.RData")
    popMat = popMatNGA
  }
  
  # get associated points in popMat
  ENCoords = proj(lonLatCoords)
  popMatIs = apply(ENCoords, function(pts) {
    # popMat is 5km resolution, so must be within 2.5 km in easting and northing directions
    closeE = (pts[1] > popMat$east - 2.5) & (pts[1] <= popMat$east + 2.5)
    closeN = (pts[2] > popMat$north - 2.5) & (pts[2] <= popMat$north + 2.5)
    closeI = closeE & closeN
    
    # there should be exactly 1 point we're closest to
    if(sum(closeI > 1)) {
      stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
    } else if(sum(closeI) == 0) {
      # this case shouldn't happen, but just take closest point then
      dists = rdist(rbind(pts), cbind(popMat$east, popMat$north))
      warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      return(which.min(dists))
    } else {
      return(which(closeI))
    }
  })
  
  # get urban classification of each point (and don't normalize urban classification)
  urbVals = popMat$urban[popMatIs]
  
  if(normalized) {
    out = load("savedOutput/global/covariatesNorm.RData")
    out = load("savedOutput/global/covariates.RData")
    
    popVals = extract(pop, lonLatCoords, method="bilinear")
    urbanicityVals = extract(urb, lonLatCoords, method="bilinear") # don't normalize urbanicity or population
    accessVals = extract(accessNorm, lonLatCoords, method="bilinear")
    elevVals = extract(elevNorm, lonLatCoords, method="bilinear")
    distVals = extract(minDistRiverLakesNorm, lonLatCoords, method="bilinear")
  } else {
    out = load("savedOutput/global/covariates.RData")
    
    popVals = extract(pop, lonLatCoords, method="bilinear")
    urbanicityVals = extract(urb, lonLatCoords, method="bilinear")
    accessVals = extract(access, lonLatCoords, method="bilinear")
    elevVals = extract(elev, lonLatCoords, method="bilinear")
    distVals = extract(minDistRiverLakes, lonLatCoords, method="bilinear")
  }
  
  cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, 
        distRiversLakes=distVals, urbanicity=urbanicityVals)
}