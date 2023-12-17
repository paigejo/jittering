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

prepareBYM2argumentsForTMB = function(graphObj, u=0.5, alpha=2/3, constr=TRUE, 
                                      scale.model=TRUE, tol=1e-8, 
                                      matrixType=c("spam", "TsparseMatrix")) {
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

dBYM2phiPC = function(phi, lambda, logitPhi=FALSE, Q=NULL, gammaTildesm1=NULL, tr=NULL, doLog=FALSE, tol=1e-8, returnPts=FALSE) {
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
  
  if(!returnPts) {
    if(!doLog) {
      exp(ldensity)
    } else {
      ldensity
    }
  } else {
    if(!doLog) {
      c(density=exp(ldensity), jacobian=exp(ljacobian), d=d)
    } else {
      c(ldensity=ldensity, ljacobian=ljacobian, d=d)
    }
  }
  
}

debugINLABYM2phi = function(graph, Q, eigenvalues = NULL, marginal.variances = NULL, 
                            rankdef, alpha, u = 1/2, lambda, scale.model = TRUE, return.as.table = FALSE, 
                            adjust.for.con.comp = TRUE, eps = sqrt(.Machine$double.eps), 
                            debug = FALSE) {
  my.debug <- function(...) if (debug) 
    cat("*** debug *** inla.pc.bym.phi: ", ..., "\n")
  stopifnot(scale.model == TRUE)
  stopifnot(adjust.for.con.comp == TRUE)
  res <- NULL
  if (missing(eigenvalues) && missing(marginal.variances)) {
    if (missing(graph) && !missing(Q)) {
      Q <- inla.as.sparse(Q)
    }
    else if (!missing(graph) && missing(Q)) {
      Q <- inla.pc.bym.Q(graph)
    }
    else if (missing(graph) && missing(Q)) {
      stop("Either <Q> or <graph> must be given.")
    }
    else {
      stop("Only one of <Q> and <graph> can be given.")
    }
    res <- inla.bym.constr.internal(Q, adjust.for.con.comp = adjust.for.con.comp)
    if (missing(rankdef)) {
      rankdef <- res$rankdef
    }
    n <- dim(Q)[1]
    res <- inla.scale.model.bym.internal(Q, adjust.for.con.comp = adjust.for.con.comp)
    Q <- res$Q
    f <- mean(res$var, na.rm = TRUE) - 1
    use.eigenvalues <- FALSE
  }
  else {
    n <- length(eigenvalues)
    eigenvalues <- pmax(0, sort(eigenvalues, decreasing = TRUE))
    gamma.invm1 <- c(1/eigenvalues[1:(n - rankdef)], rep(0, 
                                                         rankdef)) - 1
    f <- mean(marginal.variances) - 1
    use.eigenvalues <- TRUE
  }
  if (use.eigenvalues) {
    phi.s <- 1/(1 + exp(-seq(-15, 12, len = 1000)))
    d <- numeric(length(phi.s))
    k <- 1
    for (phi in phi.s) {
      aa <- n * phi * f
      bb <- sum(log1p(phi * gamma.invm1))
      d[k] <- (if (aa >= bb) 
        sqrt(aa - bb)
        else NA)
      k <- k + 1
    }
  }
  else {
    phi.s <- 1/(1 + exp(-c(seq(-15, 0, len = 40), 1:12)))
    d <- numeric(length(phi.s))
    log.q1.det <- inla.sparse.det.bym(Q, adjust.for.con.comp = adjust.for.con.comp, 
                                      constr = res$constr, rankdef = rankdef)
    d <- c(unlist(inla.mclapply(phi.s, function(phi) {
      aa <- n * phi * f
      bb <- (n * log((1 - phi)/phi) + inla.sparse.det.bym(Q + 
                                                            phi/(1 - phi) * Diagonal(n), adjust.for.con.comp = adjust.for.con.comp, 
                                                          constr = res$constr, rankdef = rankdef) - (log.q1.det - 
                                                                                                       n * log(phi)))
      my.debug("aa=", aa, " bb=", bb, " sqrt(", aa - bb, 
               ")")
      return(if (aa >= bb) sqrt(aa - bb) else NA)
    })))
  }
  names(d) <- NULL
  remove <- is.na(d)
  my.debug(paste(sum(remove), "out of", length(remove), "removed"))
  d <- d[!remove]
  phi.s <- phi.s[!remove]
  remove <- 1:6
  d <- d[-remove]
  phi.s <- phi.s[-remove]
  phi.intern <- log(phi.s/(1 - phi.s))
  ff.d <- splinefun(phi.intern, log(d))
  phi.intern <- seq(min(phi.intern) - 0, max(phi.intern), 
                    len = 10000)
  phi.s <- 1/(1 + exp(-phi.intern))
  d <- exp(ff.d(phi.intern))
  ff.d.core <- splinefun(phi.intern, log(d))
  ff.d <- function(phi.intern, deriv = 0) {
    if (deriv == 0) {
      return(exp(ff.d.core(phi.intern)))
    }
    else if (deriv == 1) {
      return(exp(ff.d.core(phi.intern)) * ff.d.core(phi.intern, 
                                                    deriv = 1))
    }
    else {
      stop("ERROR")
    }
  }
  f.d <- function(phi.s) ff.d(log(phi.s/(1 - phi.s)))
  d <- f.d(phi.s)
  if (missing(lambda)) {
    stopifnot(alpha > 0 && alpha < 1 && u > 0)
    lambda <- -log(1 - alpha)/f.d(u)
  }
  log.jac <- log(abs(ff.d(phi.intern, deriv = 1)))
  log.prior <- log(lambda) - lambda * d + log.jac
  ssum <- sum(exp(log.prior) * (c(0, diff(phi.intern)) + c(diff(phi.intern), 
                                                           0)))/2
  my.debug("empirical integral of the prior = ", ssum, ". correct it.")
  log.prior <- log.prior - log(ssum)
  if (return.as.table) {
    my.debug("return prior as a table on phi.intern")
    theta.prior.table <- paste(c("table:", cbind(phi.intern, 
                                                 log.prior)), sep = "", collapse = " ")
    return(theta.prior.table)
  }
  else {
    my.debug("return prior as a function of phi")
    prior.theta.1 <- splinefun(phi.intern, log.prior + phi.intern + 
                                 2 * log(1 + exp(-phi.intern)))
    prior.phi <- function(phi) prior.theta.1(log(phi/(1 - 
                                                        phi)))
    return(prior.phi)
  }
}

dBYM2naive = function(x, phi, tau=1, Q=NULL, V=NULL, gammaTildesm1=NULL, 
                      logitPhi=FALSE, logTau=FALSE, precomputedValues=NULL, 
                      doLog=TRUE, returnComponents=FALSE) {
  if(logitPhi) {
    lPhi = phi
    phi = expit(lPhi)
  }
  if(logTau) {
    lTau = tau
    tau = exp(lTau)
  }
  
  if(!is.null(precomputedValues)) {
    Q = precomputedValues$Q
    V = precomputedValues$V
    gammaTildesm1 = precomputedValues$gammaTildesm1
  }
  
  # calculate Sigma
  # // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag^+]^(-1) Eps^T
  # // = Eps tau [I + phi (Q_besag^+ - I)]^(-1) Eps^T
  # // = Eps tau V [I + phi (GammaTilde - I)]^(-1) V^T Eps^T
  # // i.e. the sum of squares of tau^0.5 Eps V diag(1/sqrt(1 + phi*gammaTildesm1))
  Sigma = (V %*% diag(1 + phi*(gammaTildesm1)) %*% t(V))/tau
  Q = (V %*% diag(1/(1 + phi*(gammaTildesm1))) %*% t(V))*tau
  
  # get the log determinant of Sigma
  # logDet = log(det(Sigma))
  logDet = sum(log(1 + phi*(gammaTildesm1))) - length(x) * log(tau)
  
  # get the quadratic form:
  quadForm = t(x) %*% Q %*% x
  
  # combine to get the log likelihood
  k = nrow(Q)
  ldensity = -0.5 * (logDet + quadForm + k*log(2*pi))
  
  if(!doLog) {
    val = exp(ldensity)
  } else {
    val = ldensity
  }
  
  if(returnComponents) {
    list(val=val, logDet=logDet, quadForm=quadForm)
  } else {
    val
  }
}

dBYM2sepNaive = function(u, v, phi, tau=1, Q=NULL, V=NULL, gammaTildesm1=NULL, 
                         logitPhi=FALSE, logTau=FALSE, precomputedValues=NULL, 
                         doLog=TRUE, returnComponents=FALSE, includeDetQinv=FALSE) {
  if(logitPhi) {
    lPhi = phi
    phi = expit(lPhi)
  }
  if(logTau) {
    lTau = tau
    tau = exp(lTau)
  } else {
    lTau = log(tau)
  }
  
  if(!is.null(precomputedValues)) {
    Q = precomputedValues$Q
    V = precomputedValues$V
    gammaTildesm1 = precomputedValues$gammaTildesm1
  }
  
  nAreas = nrow(Q)
  
  # calculate Sigma for Besag part
  # // leave out the log determinant of the Q term, since it is constant, 
  # // but include the term for the parameters. In other words, calculate part of: 
  # // log |(1/tau) * phi * Q^+| + log|(1/tau) * (1-phi) * I| 
  # //   = nAreas * log(phi / tau) + log|Q^+| + nAreas * log((1-phi) / tau)
  # //   = nAreas * log(phi (1 - phi) / tau^2) + log|Q^+|
  # //   = nAreas * log(phi (1 - phi) / tau^2) + C
  # // Type logDetTau = Type(nAreas) * log(phi * (1.0 - phi) / pow(tau, 2));
  # Type logDetTau = Type(nAreas) * (log(phi) + log(Type(1.0) - phi) - Type(2.0) * log_tau);
  if(!includeDetQinv) {
    logDet = nAreas * (log(phi) + log(1 - phi) - 2 * lTau)
  } else {
    logDet = nAreas * (log(phi) + log(1 - phi) - 2 * lTau) + sum(log1p(gammaTildesm1))
  }
  
  # get the quadratic forms:
  quadFormV = (tau/phi) * (t(v) %*% Q %*% v)
  quadFormU = (tau/(1-phi)) * sum(u^2)
  quadForm = quadFormU + quadFormV
  
  # combine to get the log likelihood
  k = nAreas
  ldensity = -0.5 * (logDet + quadForm)
  
  if(includeDetQinv) {
    ldensity = ldensity + k*log(2*pi)
  }
  
  if(!doLog) {
    val = exp(ldensity)
  } else {
    val = ldensity
  }
  
  if(returnComponents) {
    list(val=as.numeric(as.matrix(val)), logDet=logDet, quadForm=quadForm)
  } else {
    as.numeric(as.matrix(val))
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
  margVars = diag(pcArgs$V %*% diag(pcArgs$gammaTildesm1+1) %*% t(pcArgs$V))
  out = INLA:::inla.pc.bym.phi(Q=Qtest, lambda=lambda, rankdef=1)
  out2 = INLA:::inla.pc.bym.phi(Q=Qtest, lambda=lambda, rankdef=1, eigenvalues=gammas, marginal.variances=margVars)
  # out = INLA:::inla.pc.bym.phi(Q=Qtest, lambda=lambda, eigenvalues=gammas, rankdef=1, marginal.variances=gammas)
  browser()
  
  logDensityVals2 = sapply(logit(phis), dBYM2phiPC, lambda=lambda, Q=Q, gammaTildesm1=gammaTildesm1, tr=tr, doLog=TRUE, logitPhi=TRUE)
  head(cbind(logDensityVals, logDensityVals2, out(phis), out2(phis)))
  
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
  
  pdf(file="~/git/jittering/figures/test/pcPhi2.pdf", width=5, height=5)
  ylim = range(c(logDensityVals, out(phis)))
  plot(phis, logDensityVals, type="l", col="blue", 
       main=expression(paste("Log PC prior ", phi)), 
       xlab=expression(phi), ylab="Log density", ylim=ylim)
  lines(phis, out(phis), col="black")
  legend("top", c("My version", "INLA version"), lty=1, col=c("blue", "black"))
  dev.off()
  
  pdf(file="~/git/jittering/figures/test/pcPhiDiff.pdf", width=5, height=5)
  plot(phis, logDensityVals - out(phis), type="l", col="blue", 
       main=expression(paste("Log PC prior ", phi, " difference")), 
       xlab=expression(phi), ylab="Log density difference")
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
# same as runBYM2, except fits a single data set (the ed global data frame)
# doPredsAtPostMean: if TRUE, fix all model hyperparameters at the posterior mean 
# getPosteriorDensity: EXPERIMENTAL: evaluate the posterior density of direct estimates using multivariate normal approximation
runBYM2simple = function(dat=ed, covMat=NULL, graph="savedOutput/global/adm2Graph.dat", 
                         areaOrder=sort(unique(dat$subarea)), previousResult=NULL, 
                         doValidation=FALSE, Nsim = 1000, 
                         strategy="eb", int.strategy="ccd", 
                         fast=TRUE, getJointDraws=FALSE) {
  
  # Define formula
  if(!is.null(covMat)) {
    formula = y ~ X +
      f(idx, model="bym2",
        graph=graph, scale.model=TRUE, constr=TRUE, 
        hyper=list(prec=list(param=c(1, 0.1), prior="pc.prec"), 
                   phi=list(param=c(0.5, 2/3), prior="pc"))) +
      f(idxEps, model = "iid",
        hyper = list(prec = list(prior = "pc.prec", param = c(1,0.1))))
  } else {
    formula = y ~ 1 +
      f(idx, model="bym2",
        graph=graph, scale.model=TRUE, constr=TRUE, 
        hyper=list(prec=list(param=c(1, 0.1), prior="pc.prec"), 
                   phi=list(param=c(0.5, 2/3), prior="pc"))) +
      f(idxEps, model = "iid",
        hyper = list(prec = list(prior = "pc.prec", param = c(1,0.1))))
  }
  
  # INLA data
  subareas = dat$subarea
  if(!is.null(covMat)) {
    dat = list(y = dat$y,
               Ntrials = dat$n,
               X = covMat, 
               idx = as.numeric(match(subareas, areaOrder)),
               idxEps = 1:length(dat$y))
  } else {
    dat = list(y = dat$y,
               Ntrials = dat$n,
               idx = as.numeric(match(subareas, areaOrder)),
               idxEps = 1:length(dat$y))
  }
  
  # set posterior approximation strategy
  control.inla = list(strategy=strategy, int.strategy=int.strategy, 
                      fast=fast)
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousResult)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousResult
    modeControl$theta = previousResult$mode$theta
    modeControl$x = previousResult$mode$x
    modeControl$restart = !fixedParameters
    modeControl$fixed = fixedParameters
  }
  
  # Run model
  print("fitting BYM2 model...")
  result = inla(formula = formula, 
                family="binomial",
                Ntrials = Ntrials,
                data=dat, 
                control.compute=list(config=TRUE), 
                quantiles=c(0.1, 0.5, 0.9), 
                control.mode=modeControl)
  
  if(getJointDraws) {
    samples = inla.posterior.sample(n = Nsim, result = result)
    
    browser()
    
    result = c(list(mod=result, samples=samples))
  }
  
  result
}

plotPreds = function(SD0=NULL, tmbObj=NULL, popMat=popMatNGAThresh, gridPreds=NULL, 
                     arealPreds=NULL, normalized=TRUE, extractMethod="bilinear", 
                     nsim=1000, quantiles=c(0.025, 0.1, 0.9, 0.975), 
                     plotNameRoot="edFusion", plotNameRootAreal="Strat", CIwidthLims=NULL) {
  
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
    if(!is.null(betaDraws)) {
      betaNames = rep("beta", nrow(betaDraws))
    } else {
      betaNames = NULL
    }
    
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
  
  print(paste0("mean predicted urban prob: ", weighted.mean(preds[popMat$urban], w=popMat$pop[popMat$urban])))
  print(paste0("mean predicted rural prob: ", weighted.mean(preds[!popMat$urban], w=popMat$pop[!popMat$urban])))
  
  if(!is.null(tmbObj)) {
    # temp = tryCatch({
    #   repOut = tmbObj$env$report(c(SD0$par.fixed, SD0$par.random))
    #   print(paste0("mean data urban prob: ", mean(expit(c(repOut$latentFieldUrbMICS, repOut$latentFieldUrbDHS)), na.rm=TRUE)))
    #   print(paste0("mean data rural prob: ", mean(expit(c(repOut$latentFieldRurMICS, repOut$latentFieldRurDHS)), na.rm=TRUE)))
    # }, error=function(e) {print(e)})
    print(paste0("mean data urban prob: ", (sum(tmbObj$env$data$y_iUrbanDHS)+sum(tmbObj$env$data$y_iUrbanMICS))/(sum(tmbObj$env$data$n_iUrbanDHS)+sum(tmbObj$env$data$n_iUrbanMICS))))
    print(paste0("mean data rural prob: ", (sum(tmbObj$env$data$y_iRuralDHS)+sum(tmbObj$env$data$y_iRuralMICS))/(sum(tmbObj$env$data$n_iRuralDHS)+sum(tmbObj$env$data$n_iRuralMICS))))
  }
  
  if(!is.null(arealPreds)) {
    # plot areal predictions
    adm = arealPreds$adm
    orderedAreas = arealPreds$aggregationResults$region
    arealDraws = arealPreds$aggregationResults$p
    doAdm2 = identical(adm, adm2)
    arealMean = rowMeans(arealDraws)
    arealQuants = apply(arealDraws, 1, quantile, prob=quantiles, na.rm=TRUE)
    arealQuants[!is.finite(arealQuants)] = NA
    
    lwd = ifelse(doAdm2, .6, 1)
    border = ifelse(doAdm2, rgb(.6, .6, .6), "black")
    
    pdf(paste0("figures/ed/", plotNameRoot, "Pred", plotNameRootAreal, ".pdf"), width=5, height=3.8)
    par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
    plotMapDat(adm, arealMean, varAreas=orderedAreas, regionNames=orderedAreas, 
               cols=predCols, crosshatchNADensity=30, lwd=lwd, border=border, 
               xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", 
               ylab="Latitude", main="Predictions", legend.mar=4.7)
    if(length(arealMean) > 41) {
      plotMapDat(admFinal, new=FALSE)
    }
    dev.off()
  }
  
  if(!is.null(SD0) && !SD0$pdHess) {
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
                 cols=quantCols, crosshatchNADensity=30, lwd=lwd, border=border, 
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
                xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", zlim=CIwidthLims, 
                ylab="Latitude", main=paste0(100*thisWidthLevel, "% CI width"), 
                legend.mar=4.7)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    if(!is.null(arealPreds)) {
      thisArealWidth = arealQuants[thishighI,] - arealQuants[thisLowI,]
      
      pdf(paste0("figures/ed/", plotNameRoot, "CI", thisWidthLevel*100, plotNameRootAreal, ".pdf"), width=5, height=3.8)
      par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
      plotMapDat(adm, thisArealWidth, varAreas=orderedAreas, regionNames=orderedAreas, 
                 cols=quantCols, crosshatchNADensity=30, lwd=lwd, border=border, 
                 xlim=lonLimNGA, ylim=latLimNGA, asp=1, xlab="Longitude", zlim=CIwidthLims, 
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
predGrid = function(SD0=NULL, popMat=popMatNGAThresh, 
                    normalized=TRUE, extractMethod="bilinear", 
                    nsim=1000, quantiles=c(0.025, 0.1, 0.9, 0.975), 
                    splineApprox=TRUE, admLevel=c("stratMICS", "adm2"), 
                    predAtArea=NULL, sep=FALSE, QinvSumsNorm=NULL, 
                    includedCovs=c("urb", "access", "elev", "distRiversLakes", "popValsNorm"), 
                    constrParameterization=FALSE, obj=NULL) {
  admLevel = match.arg(admLevel)
  
  # get parameters
  if(!is.null(SD0)) {
    alpha = SD0$par.fixed[1]
    beta = SD0$par.fixed[which(names(SD0$par.fixed) == "beta")]
    parnames = names(SD0$par.fixed)
    hasNugget = ("log_tauEps" %in% parnames) || ("log_tauEpsUrb" %in% parnames) || ("log_tauEpsUDHS" %in% parnames)
    URclust = "log_tauEpsUrb" %in% parnames
    varClust = "log_tauEpsUDHS" %in% parnames
    hasUrbDiffMICS = "diffUrbMICS" %in% names(SD0$par.random)
    hasRurDiffMICS = "diffRurMICS" %in% names(SD0$par.random)
  } else {
    allPar = obj$env$last.par
    alpha = allPar[grepl("alpha", names(allPar))]
    beta = allPar[grepl("beta", names(allPar))]
    log_tauEpsI = grep("log_tauEps", names(allPar))
    log_tauEps = allPar[log_tauEpsI]
    log_tau = allPar[-log_tauEpsI][grepl("log_tauEps", names(allPar))]
    parnames = names(SD0$par.fixed)
  }
  
  
  
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
  tempNames = colnames(Xmat)
  includeI = colnames(Xmat) %in% includedCovs
  Xmat = cbind(Xmat[,includeI])
  colnames(Xmat) = tempNames[includeI]
  
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
  sigmaEpsSqUrb_tmb_draws = NULL
  sigmaEpsSqRur_tmb_draws = NULL
  sigmaEpsSqUMICS_tmb_draws = NULL
  sigmaEpsSqRMICS_tmb_draws = NULL
  sigmaEpsSqUDHS_tmb_draws = NULL
  sigmaEpsSqRDHS_tmb_draws = NULL
  if(SD0$pdHess) {
    L <- Cholesky(SD0[['jointPrecision']], super = T)
    mu = summary(SD0)[,1]
    t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = nsim)
    
    # extract fixed effects and random effects from draws
    # parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
    parnames <- colnames(SD0$jointPrecision)
    
    alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
    beta_tmb_draws    <- t.draws[parnames == 'beta',]
    sigmaSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tau',]), nrow = 1)
    phi_tmb_draws    <- matrix(expit(t.draws[parnames == 'logit_phi',]), nrow = 1)
    
    if(hasUrbDiffMICS) {
      urbDiffDraws = t.draws[parnames == 'diffUrbMICS',]
    }
    if(hasRurDiffMICS) {
      rurDiffDraws = t.draws[parnames == 'diffRurMICS',]
    }
    
    # get the spatial effect
    if(!sep) {
      epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_bym2',]
    } else {
      wStar  <- t.draws[parnames == 'w_bym2Star',]
      uStar  <- t.draws[parnames == 'u_bym2Star',] # uStar is unit var scaled
      
      if(is.null(QinvSumsNorm)) {
        if(admLevel == "adm2") {
          out = load("savedOutput/global/adm2Mat.RData")
          admMat = adm2Mat
        } else if(admLevel == "stratMICS") {
          out = load("savedOutput/global/admFinalMat.RData")
          admMat = admFinalMat
        }
        
        bym2ArgsTMB = prepareBYM2argumentsForTMB(admMat, u=0.5, alpha=2/3, 
                                                 constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
        Qinv = bym2ArgsTMB$V %*% bym2ArgsTMB$Q %*% t(bym2ArgsTMB$V)
        QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
      }
      
      # get how much u reduced by sum to zero constraint to u, then scale
      uSums = colSums(uStar)
      uFacs = uSums * sqrt(phi_tmb_draws*sigmaSq_tmb_draws)
      reduceU = outer(QinvSumsNorm, c(uFacs))
      
      # adjust Epsilon = w for the constraint on u
      epsilon_tmb_draws = wStar - reduceU
    }
    
    includeBeta = any(parnames == "beta")
    
    if(includeBeta) {
      if(!hasUrbDiffMICS && !hasRurDiffMICS) {
        fixedMat = rbind(alpha_tmb_draws, 
                         beta_tmb_draws, 
                         sigmaSq_tmb_draws, 
                         phi_tmb_draws)
        betaNames = colnames(Xmat)
        row.names(fixedMat) = c("(Int)", 
                                betaNames, 
                                "sigmaSq", 
                                "phi")
      } else if(hasUrbDiffMICS) {
        fixedMat = rbind(alpha_tmb_draws, 
                         beta_tmb_draws[1:2,], 
                         urbDiffDraws, 
                         beta_tmb_draws[3:nrow(beta_tmb_draws),], 
                         sigmaSq_tmb_draws, 
                         phi_tmb_draws)
        betaNames = colnames(Xmat)
        row.names(fixedMat) = c("(Int)", 
                                betaNames, 
                                "urbDiffMICS", 
                                "sigmaSq", 
                                "phi")
      } else {
        fixedMat = rbind(alpha_tmb_draws, 
                         beta_tmb_draws[1:2,], 
                         urbDiffDraws, 
                         rurDiffDraws, 
                         beta_tmb_draws[3:nrow(beta_tmb_draws),], 
                         sigmaSq_tmb_draws, 
                         phi_tmb_draws)
        betaNames = colnames(Xmat)
        row.names(fixedMat) = c("(Int)", 
                                betaNames, 
                                "urbDiffMICS", 
                                "rurDiffMICS", 
                                "sigmaSq", 
                                "phi")
      }
      
    } else {
      fixedMat = rbind(alpha_tmb_draws, 
                       sigmaSq_tmb_draws, 
                       phi_tmb_draws)
      row.names(fixedMat) = c("(Int)", 
                              "sigmaSq", 
                              "phi")
    }
    
    if(hasNugget) {
      if((!URclust) && (!varClust)) {
        sigmaEpsSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEps',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSq_tmb_draws)
        row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
      } else if(URclust) {
        sigmaEpsSqUrb_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUrb',]), nrow = 1)
        sigmaEpsSqRur_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRur',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSqUrb_tmb_draws, 
                         sigmaEpsSqRur_tmb_draws)
        row.names(fixedMat)[(nrow(fixedMat)-1):nrow(fixedMat)] = c("sigmaEpsSqUrb", "sigmaEpsSqRur")
      } else if(varClust) {
        sigmaEpsSqUDHS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUDHS',]), nrow = 1)
        sigmaEpsSqRDHS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRDHS',]), nrow = 1)
        sigmaEpsSqUMICS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUMICS',]), nrow = 1)
        sigmaEpsSqRMICS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRMICS',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSqUMICS_tmb_draws, 
                         sigmaEpsSqRMICS_tmb_draws, 
                         sigmaEpsSqUDHS_tmb_draws, 
                         sigmaEpsSqRDHS_tmb_draws)
        row.names(fixedMat)[(nrow(fixedMat)-3):nrow(fixedMat)] = c("sigmaEpsSqUMICS", "sigmaEpsSqRMICS", 
                                                                   "sigmaEpsSqUDHS", "sigmaEpsSqRDHS")
      }
      
    }
    
    # Make parameter summary tables
    parMeans = rowMeans(fixedMat)
    parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
    parSummary = cbind(parMeans, parQuants)
    colnames(parSummary)[1] = "Est"
    colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
    print(xtable(parSummary, digits=2))
    
    # reparameterize epsilon if need be
    if(constrParameterization) {
      if(admLevel == "adm2") {
        out = load("savedOutput/global/adm2Mat.RData")
        admMat = adm2Mat
      } else if(admLevel == "stratMICS") {
        out = load("savedOutput/global/admFinalMat.RData")
        admMat = admFinalMat
      }
      
      bym2ArgsTMB = prepareBYM2argumentsForTMB(admMat, u=0.5, alpha=2/3, 
                                               constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
      Vtilde = bym2ArgsTMB$V[,-ncol(bym2ArgsTMB$V)]
      gammaTildesm1 = bym2ArgsTMB$gammaTildesm1
      lambdas = sweep((1.0 + outer(gammaTildesm1[1:(length(gammaTildesm1)-1)], c(phi_tmb_draws))), 2, sigmaSq_tmb_draws, "*")
      epsilon_tmb_draws = Vtilde %*% (sqrt(lambdas) * epsilon_tmb_draws)
    }
    
    # add effects to predictions
    gridDraws_tmb <- as.matrix(Amat %*% epsilon_tmb_draws)
    gridDraws_tmb <- sweep(gridDraws_tmb, 2, alpha_tmb_draws, '+')
    if(includeBeta) {
      gridDraws_tmb <- gridDraws_tmb + (Xmat%*% beta_tmb_draws)
    }
    gridDrawsMICS = NULL
    if(hasUrbDiffMICS) {
      gridDrawsMICS = gridDraws_tmb
      urbanGridI = popMat$urban
      gridDrawsMICS[urbanGridI,] = sweep(gridDrawsMICS[urbanGridI,], 2, urbDiffDraws, "+")
    }
    if(hasRurDiffMICS) {
      if(is.null(gridDrawsMICS)) {
        gridDrawsMICS = gridDraws_tmb
      }
      ruralGridI = !popMat$urban
      gridDrawsMICS[ruralGridI,] = sweep(gridDrawsMICS[ruralGridI,], 2, rurDiffDraws, "+")
    }
    
    probDrawsMICS = NULL
    if(!hasNugget) {
      probDraws = expit(gridDraws_tmb)
      if(!is.null(gridDrawsMICS)) {
        probDrawsMICS = expit(gridDrawsMICS)
      }
    }
    else {
      if((!URclust) && (!varClust)) {
        probDraws <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSq_tmb_draws), 
                                                gridDraws_tmb), logisticApprox=FALSE, 
                                          splineApprox=splineApprox)
        if(!is.null(gridDrawsMICS)) {
          probDrawsMICS <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSq_tmb_draws), 
                                                      gridDrawsMICS), logisticApprox=FALSE, 
                                            splineApprox=splineApprox)
        }
      } else {
        gridDraws_tmbUrb = gridDraws_tmb[popMat$urban,]
        gridDraws_tmbRur = gridDraws_tmb[!popMat$urban,]
        if(!is.null(gridDrawsMICS)) {
          gridDrawsMICSUrb = gridDrawsMICS[popMat$urban,]
          gridDrawsMICSRur = gridDrawsMICS[!popMat$urban,]
        }
        if(URclust) {
          probDrawsUrb <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqUrb_tmb_draws), 
                                                     gridDraws_tmbUrb), logisticApprox=FALSE, 
                                               splineApprox=splineApprox)
          probDrawsRur <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqRur_tmb_draws), 
                                                     gridDraws_tmbRur), logisticApprox=FALSE, 
                                               splineApprox=splineApprox)
          if(!is.null(gridDrawsMICS)) {
            probDrawsUrbMICS <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqUrb_tmb_draws), 
                                                       gridDrawsMICSUrb), logisticApprox=FALSE, 
                                                 splineApprox=splineApprox)
            probDrawsRurMICS <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqRur_tmb_draws), 
                                                       gridDrawsMICSRur), logisticApprox=FALSE, 
                                                 splineApprox=splineApprox)
          }
        } else if(varClust) {
          probDrawsUrb <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqUDHS_tmb_draws), 
                                                     gridDraws_tmbUrb), logisticApprox=FALSE, 
                                               splineApprox=splineApprox)
          probDrawsRur <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqRDHS_tmb_draws), 
                                                     gridDraws_tmbRur), logisticApprox=FALSE, 
                                               splineApprox=splineApprox)
          if(!is.null(gridDrawsMICS)) {
            probDrawsUrbMICS <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqUMICS_tmb_draws), 
                                                           gridDrawsMICSUrb), logisticApprox=FALSE, 
                                                     splineApprox=splineApprox)
            probDrawsRurMICS <- logitNormMeanGrouped(rbind(sqrt(sigmaEpsSqRMICS_tmb_draws), 
                                                           gridDrawsMICSRur), logisticApprox=FALSE, 
                                                     splineApprox=splineApprox)
          }
        }
        
        probDraws = gridDraws_tmb
        probDraws[popMat$urban,] = probDrawsUrb
        probDraws[!popMat$urban,] = probDrawsRur
        if(!is.null(gridDrawsMICS)) {
          probDrawsMICS = gridDrawsMICS
          probDrawsMICS[popMat$urban,] = probDrawsUrbMICS
          probDrawsMICS[!popMat$urban,] = probDrawsRurMICS
        }
      }
      
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
    predMICS = quantsMICS = NULL
    if(!is.null(probDrawsMICS)) {
      predsMICS = rowMeans(probDrawsMICS)
      quantsMICS = apply(probDrawsMICS, 1, quantile, probs=quantiles, na.rm=TRUE)
    }
  }
  else {
    # Prec not PD
    stop("Precision matrix should be PD")
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
       gridDrawsMICS=probDrawsMICS, 
       epsDraws=epsilon_tmb_draws, 
       fixedMat=fixedMat, 
       alphaDraws=alpha_tmb_draws, 
       betaDraws=beta_tmb_draws, 
       sigmaSqDraws=sigmaSq_tmb_draws, 
       phiDraws=phi_tmb_draws, 
       logitGridDrawsNoNug=gridDraws_tmb, 
       sigmaEpsSqDraws=sigmaEpsSq_tmb_draws, 
       sigmaEpsSqUrbDraws=sigmaEpsSqUrb_tmb_draws, 
       sigmaEpsSqRurDraws=sigmaEpsSqRur_tmb_draws, 
       sigmaEpsSqUMICSDraws = sigmaEpsSqUMICS_tmb_draws, 
       sigmaEpsSqRMICSDraws = sigmaEpsSqRMICS_tmb_draws, 
       sigmaEpsSqUDHSDraws = sigmaEpsSqUDHS_tmb_draws, 
       sigmaEpsSqRDHSDraws = sigmaEpsSqRDHS_tmb_draws, 
       preds=preds, 
       quants=quants, 
       predsMICS=predsMICS, 
       quantsMICS=quantsMICS, 
       pdHess=SD0$pdHess)
}

# normalized: whether covariates are normalized
# extractMethod: extraction method for covariates in terra:extract
# predAtArea: name of area to predict at, if only 1
predGridINLA = function(mod, popMat=popMatNGAThresh, 
                        normalized=TRUE, extractMethod="bilinear", 
                        nsim=1000, quantiles=c(0.025, 0.1, 0.9, 0.975), 
                        splineApprox=TRUE, admLevel=c("stratMICS", "adm2"), 
                        predAtArea=NULL, 
                        includedCovs=c("urb", "access", "elev", "distRiversLakes", "popValsNorm")) {
  admLevel = match.arg(admLevel)
  
  if(!is.null(predAtArea)) {
    # if(admLevel == "stratMICS") {
    #   popMat = popMat[popMat$stratumMICS == predAtArea,]
    # } else {
    #   popMat = popMat[popMat$subarea == predAtArea,]
    # }
    popMat = popMat[popMat$area == predAtArea,]
  }
  
  # get projection matrix at prediction locations
  if(admLevel == "stratMICS") {
    Amat = t(makeApointToArea(popMat$stratumMICS, admFinal$NAME_FINAL)) # nrow(popMat) x 41
  } else {
    Amat = t(makeApointToArea(popMat$subarea, adm2Full$NAME_2)) # nrow(popMat) x (# admin2 areas)
  }
  
  samples = inla.posterior.sample(n = nsim, result=mod)
  
  # extract fixed effects and random effects from draws
  # parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  # parnames <- colnames(SD0$jointPrecision)
  
  alphaI = which("(Intercept):1" == row.names(samples[[1]]$latent))
  tauI = which(names(samples[[1]]$hyperpar) == "Precision for idx")
  phiI = which(names(samples[[1]]$hyperpar) == "Phi for idx")
  tauEpsI = which(names(samples[[1]]$hyperpar) == "Precision for idxEps")
  alpha_tmb_draws    <- sapply(samples, function(x) {x$latent[alphaI]})
  sigmaSq_tmb_draws    <- 1/sapply(samples, function(x) {x$hyperpar[tauI]})
  phi_tmb_draws    <- sapply(samples, function(x) {x$hyperpar[phiI]})
  sigmaEpsSq_tmb_draws    <- 1/sapply(samples, function(x) {x$hyperpar[tauEpsI]})
  
  # get the spatial effect
  bym2I = grep("idx:", row.names(samples[[1]]$latent))
  wI = bym2I[1:ncol(Amat)]
  epsilon_tmb_draws  <- sapply(samples, function(x) {x$latent[wI]})
  
  includedCovNames = includedCovs
  includedCovNames[includedCovNames == "popValsNorm"] = "pop"
  includedCovNames = paste0(includedCovNames, ":1")
  betaI = which(row.names(samples[[1]]$latent) %in% includedCovNames)
  includeBeta = length(betaI) > 0
  if(includeBeta) {
    beta_tmb_draws    <- sapply(samples, function(x) {x$latent[betaI]})
  } else {
    beta_tmb_draws = NULL
  }
  
  if(includeBeta) {
    # load covariates at prediction locations
    LLcoords = cbind(popMat$lon, popMat$lat)
    Xmat = getDesignMat(LLcoords, normalized)
    
    out = load("savedOutput/global/popMeanSDCal.RData")
    popMean = popMeanCalThresh
    popSD = popSDCalThresh
    popValsNorm = (log1p(Xmat[,2]) - popMean) * (1/popSD)
    Xmat = cbind(Xmat[,3:6], popValsNorm)
    tempNames = colnames(Xmat)
    includeI = colnames(Xmat) %in% includedCovs
    Xmat = cbind(Xmat[,includeI])
    colnames(Xmat) = tempNames[includeI]
    
    fixedMat = rbind(alpha_tmb_draws, 
                     beta_tmb_draws, 
                     sigmaSq_tmb_draws, 
                     phi_tmb_draws, 
                     sigmaEpsSq_tmb_draws)
    betaNames = colnames(Xmat)
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi", 
                            "sigmaEpsSq")
  } else {
    fixedMat = rbind(alpha_tmb_draws, 
                     sigmaSq_tmb_draws, 
                     phi_tmb_draws, 
                     sigmaEpsSq_tmb_draws)
    row.names(fixedMat) = c("(Int)", 
                            "sigmaSq", 
                            "phi", 
                            "sigmaEpsSq")
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
  if(includeBeta) {
    gridDraws_tmb <- gridDraws_tmb + (Xmat%*% beta_tmb_draws)
  }
  
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
  
  preds = rowMeans(probDraws)
  quants = apply(probDraws, 1, quantile, probs=quantiles, na.rm=TRUE)
  
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

simBYM2fromPrior = function(nsim=1, bym2ArgsTMB=NULL, admLevel=c("adm2", "admFinal"), 
                            logTau=1, logitPhi=0) {
  admLevel = match.arg(admLevel)
  
  if(is.null(bym2ArgsTMB)) {
    if(admLevel == "admFinal") {
      out = load("savedOutput/global/admFinalMat.RData")
      bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                               constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
    } else {
      out = load("savedOutput/global/adm2Mat.RData")
      bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                               constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
    }
  }
  V = bym2ArgsTMB$V
  gammatildesm1 = bym2ArgsTMB$gammaTildesm1
  
  # out = eigen.spam(tempQ, symmetric=TRUE)
  # gammas = out$values
  # gammaTildes = 1/gammas
  # gammaTildes[abs(gammas) < tol] = 0
  # gammaTildesm1 = gammaTildes - 1
  # // Q_besag = V Lambda V^T is the eigendecomposition of Q_besag
  # // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag^+]^(-1) Eps^T
  # // = Eps tau [I + phi (Q_besag^+ - I)]^(-1) Eps^T
  # // = Eps tau V [I + phi (GammaTilde - I)]^(-1) V^T Eps^T
  # // i.e. the sum of squares of tau^0.5 Eps V diag(1/(1 + phi*gammaTildesm1))
  
  # We are interested in the Variance matrix:
  # tau^(-1) V [I + phi (GammaTilde - I)] V^T
  # = tau^(-1) V [I + phi (GammaTildesm1)]^(1/2) [I + phi (GammaTildesm1)]^(1/2) V^T
  zMat = matrix(rnorm(nsim*nrow(V)), ncol=nsim)
  sigma = 1/sqrt(exp(logTau))
  phi = expit(logitPhi)
  sims = sigma * V %*% sweep(zMat, 1, sqrt(1 + phi*gammatildesm1), "*")
  
  return(sims)
}

simBYM2fromPrior2 = function(nsim=1, bym2ArgsTMB=NULL, admLevel=c("adm2", "admFinal"), 
                             logTau=1, logitPhi=0, tol=1e-8) {
  admLevel = match.arg(admLevel)
  tau = exp(logTau)
  sigma = 1/sqrt(tau)
  phi = expit(logitPhi)
  
  if(admLevel == "admFinal") {
    out = load("savedOutput/global/admFinalMat.RData")
    Qbesag = makeQBesag(admFinalMat, constr=TRUE, scale.model=TRUE)
  } else {
    out = load("savedOutput/global/adm2Mat.RData")
    Qbesag = makeQBesag(adm2Mat, constr=TRUE, scale.model=TRUE)
  }
  
  eigenQ = eigen(Qbesag)
  V = eigenQ$vectors
  gammas = eigenQ$values
  
  nullVals = abs(gammas < tol)
  gammas[nullVals] = 0
  
  gammaTildes = 1/gammas
  gammaTildes[abs(gammas) < tol] = 0
  # gammaTildesm1 = gammaTildes - 1
  
  
  
  # out = eigen.spam(tempQ, symmetric=TRUE)
  # gammas = out$values
  # gammaTildes = 1/gammas
  # gammaTildes[abs(gammas) < tol] = 0
  # gammaTildesm1 = gammaTildes - 1
  # // Q_besag = V Lambda V^T is the eigendecomposition of Q_besag
  # // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag^+]^(-1) Eps^T
  # // = Eps tau [I + phi (Q_besag^+ - I)]^(-1) Eps^T
  # // = Eps tau V [I + phi (GammaTilde - I)]^(-1) V^T Eps^T
  # // i.e. the sum of squares of tau^0.5 Eps V diag(1/(1 + phi*gammaTildesm1))
  
  # We are interested in the Variance matrix:
  # tau^(-1) V [phi GammaTilde] V^T
  # tau^(-1) ( V [phi GammaTilde]^(-1/2) [phi GammaTilde]^(-1/2) V^T )
  zMat = matrix(rnorm(nsim*nrow(V)), ncol=nsim)
  sims = sigma * V %*% sweep(zMat, 1, sqrt(phi*gammaTildes), "*")
  sims = sims + sqrt(1-phi)*sigma * matrix(rnorm(nsim*nrow(V)), ncol=nsim)
  return(sims)
}

plotBYM2priorCIwidth80 = function(tau=1, phi=.5, nsim=1000, bym2ArgsTMB=NULL) {
  
  sims1 = simBYM2fromPrior(nsim=nsim, bym2ArgsTMB=bym2ArgsTMB, admLevel="adm2", 
                           logTau=log(tau), logitPhi=logit(phi))
  sims2 = simBYM2fromPrior2(nsim=nsim, bym2ArgsTMB=bym2ArgsTMB, admLevel="adm2", 
                            logTau=log(tau), logitPhi=logit(phi))
  
  upper80 = expit(apply(sims1, 1, quantile, prob=.9))
  lower80 = expit(apply(sims1, 1, quantile, prob=.1))
  width80 = upper80 - lower80
  
  upper80 = expit(apply(sims2, 1, quantile, prob=.9))
  lower80 = expit(apply(sims2, 1, quantile, prob=.1))
  width802 = upper80 - lower80
  
  # determine which is Lake Chad
  lakeChadI = which(adm2@data$NAME_2 == "Lake Chad")
  width80[lakeChadI] = NA
  width802[lakeChadI] = NA
  
  plotMapDat(adm2, width80, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             crosshatchNADensity=30)
  plotMapDat(adm2, width802, varAreas=adm2@data$NAME_2, regionNames=adm2@data$NAME_2, 
             crosshatchNADensity=30)
}





