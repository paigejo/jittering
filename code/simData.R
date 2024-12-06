# simulate data from an SPDE model for simulation study 1

simData1 = function(nsim=100, margVar=1, effRange=NULL, gamma=NULL, 
                    mesh=getSPDEMesh(), easpaDat=NULL, 
                    popMat=popMatNGAThresh, targetPopMat=popMatNGAThresh, 
                    poppsub=poppsubNGAThresh, nHHSampled=16, seed=123) {
  set.seed(seed)
  
  if(is.null(easpaDat)) {
    # TODO
  }
  if(is.null(effRange)) {
    # by default, 348.7043 km
    effRange = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))/5
  }
  
  simData1 = 
    SUMMER::simPopSPDE(nsim=nsim, easpa=easpa, popMat=popMat, targetPopMat=targetPopMat, 
                       poppsub=poppsub, spdeMesh=spdeMesh, 
                       margVar=margVar, sigmaEpsilon=sigmaEpsilon, effRange=effRange, 
                       gamma=gamma, beta0=beta0, seed=NULL, nHHSampled=nHHSampled, 
                       stratifyByUrban=TRUE, subareaLevel=TRUE, 
                       doFineScaleRisk=FALSE, doSmoothRisk=FALSE, min1PerSubarea=TRUE
    )
  browser()
  
  save(simPops1, file="savedOutput/simStudy1/simPops1.RData")
  
  invisible(simPops1)
}