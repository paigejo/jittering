#' Aggregates population from the 
#' pixel level to the level of the area of interest.
#' 
#' @param popNumerators nIntegrationPoint x nsim matrix of simulated response (population numerators) for each pixel and sample
#' @param popDenominators nIntegrationPoint x nsim matrix of simulated counts (population denominators) for each pixel and sample
#' @param areas Either a nIntegrationPoint length character vector of areas (or subareas) or a nIntegrationPoint row x 
#' nDraws column matrix whose columns are vectors of areas associated with each value being aggregated.
#' @param urban nIntegrationPoint length vector of indicators specifying whether or not pixels are urban or rural. If areas is a 
#' matrix, this should be a matrix with equal dimension
#' @param targetPopMat same as in \code{\link{simPopCustom}}
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param normalize if TRUE, pixel level aggregation weights within specified area are normalized to sum to 1. This produces an 
#' average of the values in popNumerators rather than a sum. In general, should only be set to TRUE for smooth integrals of risk, say over 
#' target population density (i.e. if popDenominators is set to the target population density and popNumerators is set to the risk).
aggPreds = function(popNumerators, popDenominators, areas, urban=targetPopMat$urban, targetPopMat=NULL, 
                    stratifyByUrban=TRUE, normalize=FALSE, orderedAreas=sort(unique(areas)), estZeroPopAreas=FALSE) {
  
  # if(is.matrix(areas)) {
  #   if(normalize) {
  #     stop("normalize must be set to FALSE if areas is a matrix")
  #   }
  #   if(stratifyByUrban && !is.matrix(urban)) {
  #     stop("Either both area and urban must be matrices, or both must be vectors")
  #   } else if(!stratifyByUrban) {
  #     urban = NULL
  #   }
  #   return(aggPredsVariablePerArea(popNumerators=popNumerators, popDenominators=popDenominators, 
  #                                  areaMat=areas, urbanMat=urban))
  # }
  predsUrban = urban
  predsArea = areas
  
  # function to aggregate predictions to the given areal level. Use the 
  # following function to get numerical integration matrix for a given 
  # level of areal aggregation. Returned matrices have dimension 
  # length(unique(areaNames)) x length(areaNames)
  # areaNames: length(popNumerators) length vector of area names associated with each value 
  #            being aggregated
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    
    equalDensities = rep(1, nrow(popNumerators))
    densities = equalDensities
    
    uniqueNames = orderedAreas
    getMatrixHelper = function(i, thisUrban=NULL, thisNormalize=normalize) {
      areaI = areaNames == uniqueNames[i]
      
      theseDensities = equalDensities
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban != thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && thisNormalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!stratifyByUrban) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      rownames(integrationMatrix) = uniqueNames
      rownames(integrationMatrixUrban) = uniqueNames
      rownames(integrationMatrixRural) = uniqueNames
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!stratifyByUrban) {
      
      if(estZeroPopAreas) {
        # if there are any zero pop areas, take an equal weighted average
        zeroPopAreas = apply(A, 1, function(x) {all(popDenominators[x != 0] == 0, na.rm=TRUE)})
        if(any(zeroPopAreas)) {
          popDenominators[A[zeroPopAreas,] > 0] = .001
          warning(paste0(paste(uniqueNames[zeroPopAreas], collapse=", "), " areas have zero population. Taking equal weighted average"))
        }
        
        # set NAs and pixels without any sample size to 0
        popDenominators[is.na(popDenominators)] = 0
        popNumerators[popDenominators == 0] = 0
      }
      
      ZAggregated = A %*% popNumerators
      NAggregated = A %*% popDenominators
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      aggregationResults = list(p=pAggregated, Z=ZAggregated, N=NAggregated)
      aggregationMatrices = list(A=A, AUrban=NULL, ARural=NULL)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      if(estZeroPopAreas) {
        # if there are any zero pop areas, take an equal weighted average
        zeroPopAreasUrban = apply(AUrban, 1, function(x) {all(popDenominators[x != 0] == 0, na.rm=TRUE)})
        zeroPopAreasRural = apply(ARural, 1, function(x) {all(popDenominators[x != 0] == 0, na.rm=TRUE)})
        zeroPopAreas = zeroPopAreasUrban & zeroPopAreasRural
        if(any(zeroPopAreas)) {
          popDenominators[A[zeroPopAreas,] > 0] = .001
          warning(paste0(paste(orderedAreas[zeroPopAreas], collapse=", "), " areas have zero population. Taking equal weighted average"))
        }
        
        # set NAs and pixels without any sample size to 0
        popDenominators[is.na(popDenominators)] = 0
        popNumerators[popDenominators == 0] = 0
      }
      
      # first aggregate the numerator. The denominator will depend on the aggregation method
      ZAggregated = A %*% popNumerators
      ZAggregatedUrban = AUrban %*% popNumerators
      ZAggregatedRural = ARural %*% popNumerators
      
      # we must also aggregate the denominator to calculate 
      # the aggregated empirical proportions
      NAggregated = A %*% popDenominators
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      NAggregatedUrban = AUrban %*% popDenominators
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = NA
      
      NAggregatedRural = ARural %*% popDenominators
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = NA
      
      aggregationResults = list(region=orderedAreas, p=pAggregated, Z=ZAggregated, N=NAggregated, 
                                pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                                pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
      aggregationMatrices = list(A=A, AUrban=AUrban, ARural=ARural)
    }
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  areaResults = getIntegratedPredictions(predsArea)
  
  # return results
  areaResults
}