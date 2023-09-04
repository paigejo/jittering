# BYM2 model

# Construct nPts x nAreas observation projection matrix of 0's and 1's
makeApointToArea = function(ptAreas, uniqueAreas) {
  t(sapply(uniqueAreas, function(a) {
    ptAreas == a
  }))
}

# Construct generalized precision matrix for Besag model
# See https://journals.sagepub.com/doi/pdf/10.1177/0962280216660421 and 
# ?inla.scale.model examples for details.
# constr: constraint for inla.scale.model
makeQBesag = function(graphObj, constr=TRUE, scale.model=TRUE, matrixType=c("spam", "TsparseMatrix")) {
  matrixType = match.arg(matrixType)
  
  Q = -inla.graph2matrix(graphObj)
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  
  if(scale.model) {
    n = dim(Q)[1]
    if(constr) {
      Q.scaled = inla.scale.model(Q, constr=list(A = matrix(1, 1, n), e=0))
    } else {
      Q.scaled = inla.scale.model(Q, constr=NULL)
    }
    
    Q = Q.scaled
  }
  
  if(matrixType == "spam") {
    Q = as.spam.dgCMatrix(as(Q, "dgCMatrix"))
  }
  else {
    Q = as(Q, "TsparseMatrix")
  }
  
  Q
}

# Construct generalized precision matrix for w = (w_1^T w_2^T)^T vector, where 
# w_1 is BYM2 effect and w_2 is BYM2 structured (Besag) effect. 
# See https://journals.sagepub.com/doi/pdf/10.1177/0962280216660421 and 
# ?inla.scale.model examples for details.
# constr: constraint for inla.scale.model
makeQBYM2repar = function(graphObj, constr=TRUE, scale.model=TRUE) {
  Q = -inla.graph2matrix(graphObj)
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  
  if(scale.model) {
    Q.scaled = inla.scale.model(Q, constr=constr)
    Q = Q.scaled
  }
  
  Q = as.spam.dgCMatrix(as(Q, "dgCMatrix"))
  Q
}

prepareBYM2argumentsForTMB = function(graphObj, u=0.5, alpha=2/3, constr=TRUE, scale.model=TRUE, tol=1e-8, matrixType=c("spam", "TsparseMatrix")) {
  matrixType = match.arg(matrixType)
  
  # make Q
  Q = makeQBesag(graphObj, constr, scale.model, matrixType=matrixType)
  
  # take sparse eigendecomposition
  if(matrixType == "TsparseMatrix") {
    tempQ = as.spam.dgCMatrix(as(Q, "dgCMatrix"))
  } else {
    tempQ = Q
  }
  out = eigen.spam(tempQ, symmetric=TRUE)
  gammas = out$values
  # gammaTildes = Schur(Q, vectors=FALSE)$EValues
  gammaTildes = 1/gammas
  gammaTildes[abs(gammas) < tol] = 0
  gammaTildesm1 = gammaTildes - 1
  
  # calculate trace(Q^-1)
  tr = sum(gammaTildes)
  
  lambda = getLambdaPCphi(u, alpha, Q, gammaTildesm1, tr, tol=tol)
  
  list(lambda=lambda, Q=Q, tr=tr, gammaTildesm1=gammaTildesm1, V=out$vectors)
}

dBYM2phiPC = function(phi, lambda, logitPhi=FALSE, Q=NULL, gammaTildesm1=NULL, tr=NULL, doLog=FALSE, tol=1e-8) {
  if(!logitPhi && (length(phi) == 1) && ((phi == 0) || (phi == 1))) {
    if(doLog) {
      return(-Inf)
    } else {
      return(0)
    }
  } else if((length(phi) == 1) && ((phi == -Inf) || (phi == Inf))) {
    if(doLog) {
      return(-Inf)
    } else {
      return(0)
    }
  }

  if(is.null(gammaTildesm1)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
    gammaTildesm1 = 1/gammas - 1
    gammaTildesm1[abs(gammas) < tol] = -1
  }
  
  if(is.null(tr)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    
    # calculate tr(Q^-1)
    gammaTildes = gammaTildesm1 + 1
    gammas = 1/gammaTildes
    gammaTildes[abs(gammas) < tol] = 0
    tr = sum(gammaTildes)
  }
  
  if(length(phi) > 1) {
    ldensity = sapply(phi, dBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr, doLog=TRUE, tol=tol)
  } else {
    if(logitPhi) {
      lPhi = phi
      phi = expit(phi)
    }
    
    # calculate log determinant, KLD(phi), and d(phi)
    logDet = sum(log1p(phi*gammaTildesm1))
    
    n = nrow(Q)
    KLD = 0.5 * (phi * tr - phi * n - logDet)
    d = sqrt(2*KLD)
    
    lexpDensity = log(lambda) - lambda * d
    ljacobian = - log(d) - log(2) + log(abs(tr - n - sum(gammaTildesm1/(1 + gammaTildesm1 * phi))))
    
    if(logitPhi) {
      ljacobian = ljacobian - lPhi - 2*log1p(exp(-lPhi))
    }
    
    ldensity = lexpDensity + ljacobian
  }
  
  if(!doLog) {
    exp(ldensity)
  } else {
    ldensity
  }
}

# ...: arguments to pexp
pBYM2phiPC = function(phi, lambda, Q=NULL, gammaTildesm1=NULL, tr=NULL, tol=1e-8, ...) {
  if(is.null(gammaTildesm1)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
    gammaTildesm1 = 1/gammas - 1
    gammaTildesm1[abs(gammas) < tol] = -1
  }
  
  if(is.null(tr)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    
    # calculate tr(Q^-1)
    gammaTildes = gammaTildesm1 + 1
    tr = sum(gammaTildes)
  }
  
  if(length(phi) > 1) {
    res = sapply(phi, pBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr, doLog=TRUE, tol=tol)
  } else {
    # calculate log determinant, KLD(phi), and d(phi)
    logDet = sum(log1p(phi*gammaTildesm1))
    
    n = nrow(Q)
    KLD = 0.5 * (phi * tr - phi * n - logDet)
    d = sqrt(2*KLD)
    
    # use exponential cdf wrt rate lambda
    pexp(d, rate=lambda, ...)
  }
}

getLambdaPCphi = function(u=0.5, alpha=2/3, Q=NULL, gammaTildesm1=NULL, tr=NULL, tol=1e-8) {
  if(is.null(gammaTildesm1)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
    gammaTildesm1 = 1/gammas - 1
    gammaTildesm1[abs(gammas) < tol] = -1
  }
  
  if(is.null(tr)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    
    # calculate tr(Q^-1)
    gammaTildes = gammaTildesm1 + 1
    tr = sum(gammaTildes)
  }
  
  # calculate log determinant, KLD(u), and d(u)
  logDet = sum(log1p(u*gammaTildesm1))
  
  n = nrow(Q)
  KLD = 0.5 * (u * tr - u * n - logDet)
  d = sqrt(2*KLD)
  
  -log(1 - alpha) / d
}

dPCprec = function(tau, lambda=NULL, u=NULL, alpha=NULL, doLog=FALSE) {
  if(is.null(lambda)) {
    lambda = getLambdaPCprec(u=u, alpha=alpha)
  }
  
  ldensity = log(lambda) - log(2) -(3/2)*log(tau) - lambda/sqrt(tau)
  
  if(doLog) {
    return(ldensity)
  } else {
    exp(ldensity)
  }
}

getLambdaPCprec = function(u=0.5, alpha=2/3) {
  -log(alpha) / u
}

testpBYM2phiPC = function(u=0.5, alpha=2/3, graphFile = "~/git/U5MR/Kenyaadm1.graph") {
  g = inla.read.graph(graphFile)
  
  pcArgs = prepareBYM2argumentsForTMB(g, u, alpha)
  lambda = pcArgs$lambda
  Q = pcArgs$Q
  tr = pcArgs$tr
  gammaTildesm1 = pcArgs$gammaTildesm1
  
  Qtest = as.dgCMatrix.spam(Q)
  
  browser()
  print(paste0("P(phi < ", u, ") = ", pBYM2phiPC(u, lambda=lambda, Q=Q, tr=tr, gammaTildesm1=gammaTildesm1)))
  print(paste0("P(phi < ", u, ") = ", pBYM2phiPCnumerical(u, lambda=lambda, Q=Q, tr=tr, gammaTildesm1=gammaTildesm1)))
  print(paste0("P(phi < ", u, ") = ", pBYM2phiPCnumericalINLA(logit(u), lambda=lambda, Q=Qtest, tr=tr, gammaTildesm1=gammaTildesm1)))
  
  phis = expit(seq(-7, 7, l=1000))
  cdfsVals = sapply(phis, pBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr)
  densityVals = sapply(phis, dBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr)
  logDensityVals = sapply(phis, dBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr, doLog=TRUE)
  
  
  gammas = 1/(gammaTildesm1+1)
  gammas[gammaTildesm1 == -1] = 0
  out = INLA:::inla.pc.bym.phi(Q=Qtest, lambda=lambda, rankdef=1)
  # out = INLA:::inla.pc.bym.phi(Q=Qtest, lambda=lambda, eigenvalues=gammas, rankdef=1, marginal.variances=gammas)
  
  pdf(file="~/git/jittering/figures/test/pcPhi.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  plot(phis, cdfsVals, type="l", col="blue", 
       main=expression(paste("CDF PC prior ", phi)), 
       xlab=expression(phi), ylab="CDF")
  plot(phis, densityVals, type="l", col="blue", 
       main=expression(paste("PC prior ", phi)), 
       xlab=expression(phi), ylab="Density")
  plot(phis, logDensityVals, type="l", col="blue", 
       main=expression(paste("Log PC prior ", phi)), 
       xlab=expression(phi), ylab="Log density")
  plot(phis, out(phis), type="l", col="blue", 
       main=expression(paste("INLA log PC prior ", phi)), 
       xlab=expression(phi), ylab="Log density")
  dev.off()
  
  dBYM2phiPC(0.5, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr)
  # 0.593995
  # log values:
  # expDensity
  # [1] -2.016251
  # jacobian
  # [1] 1.495367
  dBYM2phiPC(0, logitPhi=TRUE, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr)
  # [1] 0.1484987
  exp(out(0.5))
  # 0.33875
  debugonce(out)
  debugonce(INLA:::inla.pc.bym.phi)
  # log(abs(ff.d(log(1), deriv = 1)))
  # [1] 0.108295
  # log(lambda) - lambda * f.d(0.5)
  # [1] -2.016251
  
  browser()
  print(cbind(analytical=analytical, numerical=numerical))
}

testdPCprec = function() {
  browser()
  
  print(dPCprec(1, u=1, alpha=.95))
  print(inla.pc.dprec(1, u=1, alpha=.95))
}




pBYM2phiPCnumerical = function(phi, lambda, Q=NULL, gammaTildesm1=NULL, tr=NULL, low=expit(-15), tol=1e-8) {
  if(is.null(gammaTildesm1)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
    gammaTildesm1 = 1/gammas - 1
    gammaTildesm1[abs(gammas) < tol] = 0
  }
  
  if(is.null(tr)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    
    # calculate tr(Q^-1)
    gammaTildes = gammaTildesm1 + 1
    tr = sum(gammaTildes)
  }
  
  integrate(dBYM2phiPC, low, phi, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr, tol=tol)$value
}

pBYM2phiPCnumericalINLA = function(logitPhi, lambda, Q=NULL, gammaTildesm1=NULL, tr=NULL, low=-15, tol=1e-8) {
  if(is.null(gammaTildesm1)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
    gammaTildesm1 = 1/gammas - 1
    gammaTildesm1[abs(gammas) < tol] = 0
  }
  
  if(is.null(tr)) {
    if(is.null(Q)) {
      stop("must either provide Q or gammaTildesm1 and tr")
    }
    
    # calculate tr(Q^-1)
    gammaTildes = gammaTildesm1 + 1
    tr = sum(gammaTildes)
  }
  
  INLAdensity = INLA:::inla.pc.bym.phi(Q=Q, lambda=lambda, rankdef=1)
  
  integrate(function(x) {exp(INLAdensity(expit(x)))}, low, logitPhi)$value
}


# run the jittering model
runModBYM2jitter = function(modelName = c("")) {
  
  ys = dat$ys
  ns = dat$ns
  
  
}

plotPreds = function(SD0, tmbObj=NULL, popMat=popMatNGAThresh, gridPreds=NULL, 
                     arealPreds=NULL, normalized=TRUE, extractMethod="bilinear", 
                     nsim=1000, quantiles=c(0.025, 0.1, 0.9, 0.975), 
                     plotNameRoot="edFusion", plotNameRootAreal="Strat") {
  
  # get grid level predictions if need be
  if(is.null(gridPreds)) {
    gridPreds = predGrid(SD0, popMat, normalized, 
                         extractMethod, nsim, quantiles)
  }
  
  # extract variables from gridPreds
  popMat = gridPreds$popMat
  gridDraws = gridPreds$gridDraws
  epsDraws = gridPreds$epsDraws
  alphaDraws = gridPreds$alphaDraws
  betaDraws = gridPreds$betaDraws
  sigmaSqDraws = gridPreds$sigmaSqDraws
  phiDraws = gridPreds$phiDraws
  logitGridDrawsNoNug = gridPreds$logitGridDrawsNoNug
  sigmaEpsSqDraws = gridPreds$sigmaEpsSqDraws
  quants = gridPreds$quants
  pdHess = gridPreds$pdHess
  
  # print parameter summary table
  if(pdHess) {
    fixedMat = rbind(alphaDraws, 
                     betaDraws, 
                     sigmaSqDraws, 
                     phiDraws)
    betaNames = rep("beta", nrow(betaDraws))
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi")
  } else {
    fixedMat = matrix(c(alphaDraws, 
                     betaDraws, 
                     sigmaSqDraws, 
                     phiDraws), ncol=1)
    betaDraws = matrix(betaDraws, ncol=1)
    betaNames = rep("beta", nrow(betaDraws))
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi")
  }
  
  
  hasNugget = !is.null(sigmaEpsSqDraws)
  if(hasNugget) {
    fixedMat = rbind(fixedMat, 
                     sigmaEpsSqDraws)
    row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
  }
  
  # Make parameter summary tables
  parMeans = rowMeans(fixedMat)
  parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
  parSummary = cbind(parMeans, parQuants)
  colnames(parSummary)[1] = "Est"
  colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
  print(xtable(parSummary, digits=2))
  
  if(pdHess) {
    preds = rowMeans(gridDraws, na.rm=TRUE)
  } else {
    preds = gridDraws
  }
  
  
  predCols = makeBlueGreenYellowSequentialColors(64)
  quantCols = makePurpleYellowSequentialColors(64)
  
  pdf(paste0("figures/ed/", plotNameRoot, "PredGrid.pdf"), width=5, height=3.8)
  par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  myQuiltPlot(popMat$lon, popMat$lat, preds, colScale=predCols, nx=259, 
              xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
              ylab="Latitude", main="Predictions", legend.mar=4.7, zlim=c(0,1))
  plotMapDat(admFinal, new=FALSE)
  dev.off()
  
  print(paste0("mean predicted urban prob: ", mean(preds[popMat$urban], na.rm=TRUE)))
  print(paste0("mean predicted rural prob: ", mean(preds[!popMat$urban], na.rm=TRUE)))
  
  if(!is.null(tmbObj)) {
    temp = tryCatch({
      repOut = tmbObj$env$report(c(SD0$par.fixed, SD0$par.random))
      print(paste0("mean data urban prob: ", mean(expit(c(repOut$latentFieldUrbMICS, repOut$latentFieldUrbDHS)), na.rm=TRUE)))
      print(paste0("mean data rural prob: ", mean(expit(c(repOut$latentFieldRurMICS, repOut$latentFieldRurDHS)), na.rm=TRUE)))
    }, error=function(e) {print(e)})
  }
  
  if(!is.null(arealPreds)) {
    # plot areal predictions
    adm = arealPreds$adm
    orderedAreas = arealPreds$aggregationResults$region
    arealDraws = arealPreds$aggregationResults$p
    arealMean = rowMeans(arealDraws)
    arealQuants = apply(arealDraws, 1, quantile, prob=quantiles, na.rm=TRUE)
    arealQuants[!is.finite(arealQuants)] = NA
    
    pdf(paste0("figures/ed/", plotNameRoot, "Pred", plotNameRootAreal, ".pdf"), width=5, height=3.8)
    par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
    plotMapDat(adm, arealMean, varAreas=orderedAreas, regionNames=orderedAreas, 
               cols=predCols, crosshatchNADensity=30, 
               xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
               ylab="Latitude", main="Predictions", legend.mar=4.7)
    if(length(arealMean) > 41) {
      plotMapDat(admFinal, new=FALSE)
    }
    dev.off()
  }
  
  if(!SD0$pdHess) {
    warning("SD0 hessian not PD. Only predictions plotted")
    return(invisible(NULL))
  }
  
  # quantiles
  for(i in 1:length(quantiles)) {
    thisQuant = quantiles[i]
    pdf(paste0("figures/ed/", plotNameRoot, "Quant", thisQuant, "Grid.pdf"), width=5, height=3.8)
    par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
    myQuiltPlot(popMat$lon, popMat$lat, quants[i,], colScale=quantCols, nx=259, 
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
                ylab="Latitude", main=paste0(thisQuant, "th quantile"), 
                legend.mar=4.7, zlim=c(0,1))
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    if(!is.null(arealPreds)) {
      pdf(paste0("figures/ed/", plotNameRoot, "Quant", thisQuant, plotNameRootAreal, ".pdf"), width=5, height=3.8)
      par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
      plotMapDat(adm, arealQuants[,i], varAreas=orderedAreas, regionNames=orderedAreas, 
                 cols=quantCols, crosshatchNADensity=30, 
                 xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
                 ylab="Latitude", main=paste0(thisQuant, "th quantile"), legend.mar=4.7, zlim=range(arealQuants, na.rm=TRUE))
      if(length(arealMean) > 41) {
        plotMapDat(admFinal, new=FALSE)
      }
      dev.off()
    }
  }
  
  # CI widths
  lowsOrig = quantiles[quantiles < 0.5]
  lows = rev(lowsOrig)
  highs = quantiles[quantiles > 0.5]
  if(!all(highs == (1-lows))) {
    stop("quantiles must come in pairs or must be 0.5")
  }
  widths = quants[1:length(lows),]
  # for(i in 1:length(lows)) {
  #   thisLow = lows[i]
  #   thisHigh = highs[i]
  #   thisLowI = match(thisLow, lowsOrig)
  #   thishighI = which(quantiles > 0.5)[i]
  #   thisWidth = quants[thishighI,] - quants[thisLowI,]
  #   widths[i,] = thisWidth
  # }
  # zlim = c(0, max(widths, na.rm=TRUE))
  
  for(i in 1:length(lows)) {
    thisLow = lows[i]
    thisHigh = highs[i]
    thisLowI = match(thisLow, lowsOrig)
    thishighI = which(quantiles > 0.5)[i]
    thisWidth = quants[thishighI,] - quants[thisLowI,]
    thisWidthLevel = thisHigh - thisLow
    
    pdf(paste0("figures/ed/", plotNameRoot, "CI", thisWidthLevel*100, "Grid.pdf"), width=5, height=3.8)
    par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
    myQuiltPlot(popMat$lon, popMat$lat, thisWidth, colScale=quantCols, nx=259, 
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
                ylab="Latitude", main=paste0(100*thisWidthLevel, "% CI width"), 
                legend.mar=4.7)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    if(!is.null(arealPreds)) {
      thisArealWidth = arealQuants[thishighI,] - arealQuants[thisLowI,]
      
      pdf(paste0("figures/ed/", plotNameRoot, "CI", thisWidthLevel*100, plotNameRootAreal, ".pdf"), width=5, height=3.8)
      par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
      plotMapDat(adm, thisArealWidth, varAreas=orderedAreas, regionNames=orderedAreas, 
                 cols=quantCols, crosshatchNADensity=30, 
                 xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
                 ylab="Latitude", main=paste0(round(thisWidthLevel*100), "% CI Width"), legend.mar=4.7)
      if(length(arealMean) > 41) {
        plotMapDat(admFinal, new=FALSE)
      }
      dev.off()
    }
  }
}

summaryTabBYM2 = function(SD0, popMat=popMatNGAThresh, gridPreds=NULL, 
                          normalized=TRUE, extractMethod="bilinear", 
                          nsim=1000, quantiles=c(0.025, 0.975)) {
  
  # get grid level predictions if need be
  if(is.null(gridPreds)) {
    gridPreds = predGrid(SD0, popMat, normalized, 
                         extractMethod, nsim, quantiles)
  }
  
  # extract variables from gridPreds
  gridDraws = gridPreds$gridDraws
  epsDraws = gridPreds$epsDraws
  alphaDraws = gridPreds$alphaDraws
  betaDraws = gridPreds$betaDraws
  sigmaSqDraws = gridPreds$sigmaSqDraws
  phiDraws = gridPreds$phiDraws
  logitGridDrawsNoNug = gridPreds$logitGridDrawsNoNug
  sigmaEpsSqDraws = gridPreds$sigmaEpsSqDraws
  quants = gridPreds$quants
  pdHess = gridPreds$pdHess
  
  # print parameter summary table
  fixedMat = rbind(alphaDraws, 
                   betaDraws, 
                   sigmaSqDraws, 
                   phiDraws)
  betaNames = rep("beta", nrow(betaDraws))
  row.names(fixedMat) = c("(Int)", 
                          betaNames, 
                          "sigmaSq", 
                          "phi")
  
  hasNugget = !is.null(sigmaEpsSqDraws)
  if(hasNugget) {
    fixedMat = rbind(fixedMat, 
                     sigmaEpsSqDraws)
    row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
  }
  
  # Make parameter summary tables
  parMeans = rowMeans(fixedMat)
  parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
  parSummary = cbind(parMeans, parQuants)
  colnames(parSummary)[1] = "Est"
  colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
  print(xtable(parSummary, digits=2))
  
  invisible(NULL)
}

# normalized: whether covariates are normalized
# extractMethod: extraction method for covariates in terra:extract
# predAtArea: name of area to predict at, if only 1
predGrid = function(SD0, popMat=popMatNGAThresh, 
                    normalized=TRUE, extractMethod="bilinear", 
                    nsim=1000, quantiles=c(0.025, 0.1, 0.9, 0.975), 
                    splineApprox=TRUE, admLevel=c("stratMICS", "adm2"), 
                    predAtArea=NULL) {
  admLevel = match.arg(admLevel)
  
  # get parameters
  alpha = SD0$par.fixed[1]
  beta = SD0$par.fixed[which(names(SD0$par.fixed) == "beta")]
  beta = SD0$par.fixed[which(names(SD0$par.fixed) == "beta")]
  parnames = names(SD0$par.fixed)
  hasNugget = "log_tauEps" %in% parnames
  
  if(!is.null(predAtArea)) {
    # if(admLevel == "stratMICS") {
    #   popMat = popMat[popMat$stratumMICS == predAtArea,]
    # } else {
    #   popMat = popMat[popMat$subarea == predAtArea,]
    # }
    popMat = popMat[popMat$area == predAtArea,]
  }
  
  # load covariates at prediction locations
  LLcoords = cbind(popMat$lon, popMat$lat)
  Xmat = getDesignMat(LLcoords, normalized)
  
  out = load("savedOutput/global/popMeanSDCal.RData")
  popMean = popMeanCalThresh
  popSD = popSDCalThresh
  popValsNorm = (log1p(Xmat[,2]) - popMean) * (1/popSD)
  Xmat = cbind(Xmat[,3:6], popValsNorm)
  
  # get projection matrix at prediction locations
  if(admLevel == "stratMICS") {
    Amat = t(makeApointToArea(popMat$stratumMICS, admFinal$NAME_FINAL)) # nrow(popMat) x 41
  } else {
    Amat = t(makeApointToArea(popMat$subarea, adm2Full$NAME_2)) # nrow(popMat) x (# admin2 areas)
  }
  
  # generate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  
  rmvnorm_prec2 <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- solve(t(as(chol_prec, "sparseMatrix"))) # Sigma = L %*% L^T
    z <- as.matrix(L %*% z)
    mu + z
  }
  
  sigmaEpsSq_tmb_draws = NULL
  if(SD0$pdHess) {
    L <- Cholesky(SD0[['jointPrecision']], super = T)
    mu = summary(SD0)[,1]
    t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = nsim)
    
    # extract fixed effects and random effects from draws
    parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
    epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_bym2',]
    alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
    beta_tmb_draws    <- t.draws[parnames == 'beta',]
    sigmaSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tau',]), nrow = 1)
    phi_tmb_draws    <- matrix(expit(t.draws[parnames == 'logit_phi',]), nrow = 1)
    
    fixedMat = rbind(alpha_tmb_draws, 
                     beta_tmb_draws, 
                     sigmaSq_tmb_draws, 
                     phi_tmb_draws)
    betaNames = colnames(Xmat)
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi")
    
    if(hasNugget) {
      sigmaEpsSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEps',]), nrow = 1)
      fixedMat = rbind(fixedMat, 
                       sigmaEpsSq_tmb_draws)
      row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
    }
    
    # Make parameter summary tables
    parMeans = rowMeans(fixedMat)
    parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
    parSummary = cbind(parMeans, parQuants)
    colnames(parSummary)[1] = "Est"
    colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
    print(xtable(parSummary, digits=2))
    
    # add effects to predictions
    gridDraws_tmb <- as.matrix(Amat %*% epsilon_tmb_draws)
    gridDraws_tmb <- sweep(gridDraws_tmb, 2, alpha_tmb_draws, '+')
    gridDraws_tmb <- gridDraws_tmb + (Xmat%*% beta_tmb_draws)
    
    if(!hasNugget) {
      probDraws = expit(gridDraws_tmb)
    }
    else {
      browser()
      probDraws <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSq_tmb_draws), 
                                              gridDraws_tmb), logisticApprox=FALSE, 
                                        splineApprox=splineApprox)
      # logitNormMeanGrouped is muuuuuuch faster than:
      # system.time(probDraws <- matrix(logitNormMean(cbind(c(gridDraws_tmb[,1:500]), rep(sqrt(sigmaEpsSq_tmb_draws[1:500]), each=nrow(gridDraws_tmb))), logisticApprox=FALSE, splineApprox=splineApprox), nrow=nrow(gridDraws_tmb)))
      
      if(FALSE) {
        # test spline approximation timing and accuracy versus regular
        first2Sigmas = 1:(2*nrow(gridDraws_tmb))
        inMat = cbind(c(gridDraws_tmb), rep(sqrt(sigmaEpsSq_tmb_draws), each=nrow(gridDraws_tmb)))[first2Sigmas,]
        system.time(probDrawsSp <- matrix(logitNormMean(inMat, logisticApprox=FALSE, splineApprox=TRUE), nrow=nrow(gridDraws_tmb)))[3]
        # elapsed 
        # 0.038
        system.time(probDrawsReg <- matrix(logitNormMean(inMat, logisticApprox=FALSE, splineApprox=FALSE), nrow=nrow(gridDraws_tmb)))[3]
        # elapsed 
        # 3.966 
        # 0.038 / 3.966
        # 0.009581442 (less than 1% of the time!!!)
        
        mean(abs(probDrawsSp - probDrawsReg))
        # 8.469822e-08
        
        mean((probDrawsSp - probDrawsReg)^2)
        # 1.753772e-10
        
        orderI = order(inMat[,1])
        plot(inMat[orderI,1], probDrawsSp[orderI], type="l", col="purple")
        lines(inMat[orderI,1], probDrawsReg[orderI], col="green")
      }
    }
    
    preds = rowMeans(probDraws)
    quants = apply(probDraws, 1, quantile, probs=quantiles, na.rm=TRUE)
  }
  else {
    Eps = SD0$par.random[grepl("Epsilon", names(SD0$par.random))]
    alpha = SD0$par.fixed[grepl("alpha", names(SD0$par.fixed))]
    beta = SD0$par.fixed[grepl("beta", names(SD0$par.fixed))]
    
    # set "draws" to be just the fixed values
    epsilon_tmb_draws = Eps
    alpha_tmb_draws = alpha
    beta_tmb_draws = beta
    phi_tmb_draws = expit(SD0$par.fixed[grepl("logit_phi", names(SD0$par.fixed))])
    sigmaSq_tmb_draws = 1/exp(SD0$par.fixed[names(SD0$par.fixed) == "log_tau"])
    
    fixedMat = rbind(alpha_tmb_draws, 
                     beta_tmb_draws, 
                     sigmaSq_tmb_draws, 
                     phi_tmb_draws)
    betaNames = colnames(Xmat)
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi")
    
    # add effects to predictions
    gridDraws_tmb <- Amat %*% Eps
    gridDraws_tmb <- gridDraws_tmb + alpha
    gridDraws_tmb <- gridDraws_tmb + (Xmat %*% beta)
    
    if(!hasNugget) {
      probDraws = expit(gridDraws_tmb)
    }
    else {
      tauEps = exp(SD0$par.fixed[grepl("log_tauEps", names(SD0$par.fixed))])
      sigmaEps = 1/sqrt(tauEps)
      probDraws = logitNormMean(cbind(c(gridDraws_tmb), rep(sigmaEps, length(gridDraws_tmb))), logisticApprox=splineApprox)
      
      sigmaEpsSq_tmb_draws = sigmaEps^2
      
      fixedMat = rbind(fixedMat, 
                       sigmaEpsSq_tmb_draws)
      row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
    }
    
    preds = probDraws
    quants = NULL
  }
  
  list(popMat=popMat, 
       gridDraws=probDraws, 
       epsDraws=epsilon_tmb_draws, 
       fixedMat=fixedMat, 
       alphaDraws=alpha_tmb_draws, 
       betaDraws=beta_tmb_draws, 
       sigmaSqDraws=sigmaSq_tmb_draws, 
       phiDraws=phi_tmb_draws, 
       logitGridDrawsNoNug=gridDraws_tmb, 
       sigmaEpsSqDraws=sigmaEpsSq_tmb_draws, 
       quants=quants, 
       pdHess=SD0$pdHess)
}

# for now, just uses smooth latent aggregation model
predArea = function(gridPreds, areaVarName="stratumMICS", 
                    adm=NULL, orderedAreas=NULL, estZeroPopAreas=FALSE) {
  # set default shapefile associated with predictions
  if(is.null(adm)) {
    if(areaVarName == "stratumMICS") {
      adm = admFinal
    }
    else if(areaVarName == "subarea") {
      adm = adm2
    }
    else if(areaVarName == "area") {
      adm = adm1
    }
  }
  
  popMat = gridPreds$popMat
  areas = popMat[[areaVarName]]
  
  if(is.null(orderedAreas)) {
    orderedAreas = sort(unique(areas))
  }
  
  probDraws = gridPreds$gridDraws
  naPixels = apply(probDraws, 1, function(x) {any(is.na(x))})
  if(any(naPixels)) {
    warning("removing NA pixels from aggregation")
    probDraws[naPixels,] = 0
    popMat$pop[probDraws] = 0
  }
  
  popDenominators = matrix(rep(popMat$pop, ncol(probDraws)), ncol=ncol(probDraws))
  popNumerators = sweep(probDraws, 1, popMat$pop, FUN="*")
  
  urban=popMat$urban
  arealPreds = aggPreds(popNumerators=popNumerators, popDenominators=popDenominators, 
                        areas=areas, urban=urban, orderedAreas=orderedAreas, 
                        stratifyByUrban=TRUE, normalize=FALSE, 
                        estZeroPopAreas=estZeroPopAreas)
  
  c(arealPreds, list(adm=adm))
}









