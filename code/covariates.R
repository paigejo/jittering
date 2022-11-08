# Functions for handling covariate rasters

# covariates include:
#   intercept
#   pop
#   urb
#   access
#   elev
#   minDistRiverLakes
getDesignMat = function(lonLatCoords, normalized=TRUE, stratifyPop=TRUE) {
  
  if(normalized) {
    out = load("savedOutput/global/covariatesNorm.RData")
    out = load("savedOutput/global/covariates.RData")
    
    popVals = extract(pop, lonLatCoords, method="bilinear")
    urbVals = extract(urb, lonLatCoords, method="bilinear") # don't normalize urbanicity or population
    accessVals = extract(accessNorm, lonLatCoords, method="bilinear")
    elevVals = extract(elevNorm, lonLatCoords, method="bilinear")
    distVals = extract(minDistRiverLakesNorm, lonLatCoords, method="bilinear")
  } else {
    out = load("savedOutput/global/covariates.RData")
    
    popVals = extract(pop, lonLatCoords, method="bilinear")
    urbVals = extract(urb, lonLatCoords, method="bilinear")
    accessVals = extract(access, lonLatCoords, method="bilinear")
    elevVals = extract(elev, lonLatCoords, method="bilinear")
    distVals = extract(minDistRiverLakes, lonLatCoords, method="bilinear")
  }
  
  cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, distRiversLakes=distVals)
}