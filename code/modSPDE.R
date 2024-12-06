
# get default SPDE mesh
getSPDEMesh = function(locs=cbind(ed$east, ed$north), n=3000, max.n=5000, doPlot=TRUE, max.edge=c(25, 250), 
                       offset=-.08, cutoff=15) {
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    plotMapDat(mapDat=adm1, project=TRUE, border="blue", new=FALSE, 
               xlim=eastLimNGA, ylim=northLimNGA, myProjection=projNigeria)
  }
  
  mesh
}

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# or use inla.spde2.pcmatern (possibly allow (1/4,4) variance here rather than (1/2,2))
getSPDEPrior = function(mesh, sigma0=1, strictPrior=FALSE) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # 1277.237
  # range0 <- size/5
  # kappa0 <- sqrt(8)/range0
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
  #                           B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
  #                           theta.prior.prec = c(0.1, 1))
  
  range0 <- size/5
  if(!strictPrior)
    spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(1, 0.01))
  else
    spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(.15, 0.01))
  spde
}

# simulate from an SPDE model
simSPDE = function(coords, nsim=1, mesh=NULL, effRange=NULL, margVar=1) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(coords, doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  if(is.null(effRange)) {
    effRange = meshSize/5
  }
  
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
  simField = inla.qsample(nsim, Q)
  simDat = as.matrix(A %*% simField)
  
  simDat
}