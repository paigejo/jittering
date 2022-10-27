# BYM2 model

# Construct generalized precision matrix for Besag model
# See https://journals.sagepub.com/doi/pdf/10.1177/0962280216660421 and 
# ?inla.scale.model examples for details.
# constr: constraint for inla.scale.model
makeQBesag = function(graphObj, constr=TRUE, scale.model=TRUE) {
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
  
  Q = as.spam.dgCMatrix(as(Q, "dgCMatrix"))
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

prepareBYM2argumentsForTMB = function(graphObj, u=0.5, alpha=2/3, constr=TRUE, scale.model=TRUE, tol=1e-8) {
  # make Q
  Q = makeQBesag(graphObj, constr, scale.model)
  
  # take sparse eigendecomposition
  gammas = eigen.spam(Q, symmetric=TRUE, only.values = TRUE)$values
  # gammaTildes = Schur(Q, vectors=FALSE)$EValues
  gammaTildes = 1/gammas
  gammaTildes[abs(gammas) < tol] = 0
  gammaTildesm1 = gammaTildes - 1
  
  # calculate trace(Q^-1)
  tr = sum(gammaTildes)
  
  lambda = getLambdaPCphi(u, alpha, Q, gammaTildesm1, tr, tol=tol)
  
  list(lambda=lambda, Q=Q, tr=tr, gammaTildesm1=gammaTildesm1)
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
runModBYM2jitter = function(dat, popMat, adm1Map, stratByUrban=TRUE, 
                            numPointsUrbRur=c(11, 16), 
                            scalingFactor=1, 
                            JInner=3, JOuterUrbRur=c(0, 1), 
                            integrationPointType=c("mean", "midpoint"), 
                            verbose=TRUE) {
  
  ys = dat$ys
  ns = dat$ns
  
  
}






