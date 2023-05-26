# Functions for handling covariate rasters

# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
#   urbanicity
# Inputs:
# setMissingToAvg: if TRUE, sets missing covariates average value in NGA (or to 0 if normalized==TRUE)
getDesignMat = function(lonLatCoords, normalized=TRUE, 
                        useThreshPopMat=TRUE, proj=projNigeria, testMode=FALSE, 
                        setMissingToAvg=TRUE) {
  
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
  popMatIs = apply(ENCoords, 1, function(pts) {
    # popMat is 5km resolution, so must be within 2.5 km in easting and northing directions
    closeE = (pts[1] > popMat$east - 2.5) & (pts[1] <= popMat$east + 2.5)
    closeN = (pts[2] > popMat$north - 2.5) & (pts[2] <= popMat$north + 2.5)
    closeI = closeE & closeN
    
    # there should be exactly 1 point we're closest to
    if(sum(closeI > 1)) {
      stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
    } else if(sum(closeI) == 0) {
      if(testMode) {
        return(NA)
      }
      # this case shouldn't happen, but just take closest point then
      dists = rdist(rbind(pts), cbind(popMat$east, popMat$north))
      # warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      minI = which.min(dists)
      if(dists[minI] > 10) {
        warning(paste("nearest grid pt to: (", paste(pts, collapse=", "), ") is ", dists[minI], " km away.", collapse="", sep=""))
      }
      return(minI)
    } else {
      return(which(closeI))
    }
  })
  
  # get urban classification of each point (and don't normalize urban classification)
  urbVals = popMat$urban[popMatIs]
  
  if(normalized) {
    out = load("savedOutput/global/covariatesNorm.RData")
    # out = load("savedOutput/global/covariates.RData")
    
    popVals = terra::extract(pop, lonLatCoords, method="bilinear")
    urbanicityVals = terra::extract(urb, lonLatCoords, method="bilinear") # don't normalize urbanicity or population
    accessVals = terra::extract(accessNorm, lonLatCoords, method="bilinear")
    elevVals = terra::extract(elevNorm, lonLatCoords, method="bilinear")
    distVals = terra::extract(minDistRiverLakesNorm, lonLatCoords, method="bilinear")
    
    if(setMissingToAvg) {
      popNAs = is.na(popVals)
      urbNAs = is.na(urbanicityVals)
      accessNAs = is.na(accessVals)
      elevNAs = is.na(elevVals)
      distNAs = is.na(distVals)
      
      if(any(popNAs)) {
        popVals[popNAs] = 0
      }
      if(any(urbNAs)) {
        urbVals[urbNAs] = 0
      }
      if(any(accessNAs)) {
        accessVals[accessNAs] = 0
      }
      if(any(elevNAs)) {
        elevVals[elevNAs] = 0
      }
      if(any(distNAs)) {
        distVals[distNAs] = 0
      }
    }
  } else {
    out = load("savedOutput/global/covariates.RData")
    
    popVals = terra::extract(pop, lonLatCoords, method="bilinear")
    urbanicityVals = terra::extract(urb, lonLatCoords, method="bilinear")
    accessVals = terra::extract(access, lonLatCoords, method="bilinear")
    elevVals = terra::extract(elev, lonLatCoords, method="bilinear")
    distVals = terra::extract(minDistRiverLakes, lonLatCoords, method="bilinear")
    
    if(setMissingToAvg) {
      popNAs = is.na(popVals)
      urbNAs = is.na(urbanicityVals)
      # accessNAs = is.na(accessVals)
      # elevNAs = is.na(elevVals)
      # distNAs = is.na(distVals)
      
      if(any(popNAs)) {
        popVals[popNAs] = 0
      }
      if(any(urbNAs)) {
        urbVals[urbNAs] = 0
      }
    }
  }
  
  cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, 
        distRiversLakes=distVals, urbanicity=urbanicityVals)
}

# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
#   urbanicity
# Inputs:
# setMissingToAvg: if TRUE, sets missing covariates average value in NGA (or to 0 if normalized==TRUE)
getDesignMatPopNorm = function(lonLatCoords, 
                        useThreshPopMat=TRUE, proj=projNigeria, testMode=FALSE, 
                        setMissingToAvg=TRUE) {
  
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
  popMatIs = apply(ENCoords, 1, function(pts) {
    # popMat is 5km resolution, so must be within 2.5 km in easting and northing directions
    closeE = (pts[1] > popMat$east - 2.5) & (pts[1] <= popMat$east + 2.5)
    closeN = (pts[2] > popMat$north - 2.5) & (pts[2] <= popMat$north + 2.5)
    closeI = closeE & closeN
    
    # there should be exactly 1 point we're closest to
    if(sum(closeI > 1)) {
      stop(paste("close to multiple grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
    } else if(sum(closeI) == 0) {
      if(testMode) {
        return(NA)
      }
      # this case shouldn't happen, but just take closest point then
      dists = rdist(rbind(pts), cbind(popMat$east, popMat$north))
      # warning(paste("no close grid pts: (", paste(pts, collapse=", "), ")", collapse="", sep=""))
      minI = which.min(dists)
      if(dists[minI] > 10) {
        warning(paste("nearest grid pt to: (", paste(pts, collapse=", "), ") is ", dists[minI], " km away.", collapse="", sep=""))
      }
      return(minI)
    } else {
      return(which(closeI))
    }
  })
  
  # get urban classification of each point (and don't normalize urban classification)
  urbVals = popMat$urban[popMatIs]
  
  out = load("savedOutput/global/covariatesNorm.RData")
  # out = load("savedOutput/global/covariates.RData")
  
  popVals = terra::extract(pop, lonLatCoords, method="bilinear")
  urbanicityVals = terra::extract(urb, lonLatCoords, method="bilinear") # don't normalize urbanicity or population
  accessVals = terra::extract(accessNorm, lonLatCoords, method="bilinear")
  elevVals = terra::extract(elevNorm, lonLatCoords, method="bilinear")
  distVals = terra::extract(minDistRiverLakesNorm, lonLatCoords, method="bilinear")
  
  if(setMissingToAvg) {
    popNAs = is.na(popVals)
    urbNAs = is.na(urbanicityVals)
    accessNAs = is.na(accessVals)
    elevNAs = is.na(elevVals)
    distNAs = is.na(distVals)
    
    if(any(popNAs)) {
      popVals[popNAs] = 0
    }
    if(any(urbNAs)) {
      urbVals[urbNAs] = 0
    }
    if(any(accessNAs)) {
      accessVals[accessNAs] = 0
    }
    if(any(elevNAs)) {
      elevVals[elevNAs] = 0
    }
    if(any(distNAs)) {
      distVals[distNAs] = 0
    }
    
    # normalize the population densities
    out = load("savedOutput/global/popMeanSDCal.RData")
    popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
    popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
    popVals = (log1p(popVals) - popMean)*(1/popSD)
  }
  
  cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, 
        distRiversLakes=distVals, urbanicity=urbanicityVals)
}