
fitMM = function(datDHS=NULL, datMCIS, inputsMDM=NULL, 
                  intPtsMICS=NULL, intPtsDHS=NULL, 
                  KMICS=100,
                  KDHSurb = 11, # 3 rings of 5 each
                  JInnerUrban = 3,
                  KDHSrur = 16, # 3 inner + 1 outer rings of 5 each
                  JInnerRural = 3,
                  JOuterRural = 1, admMICS=admFinal, adm2DHS=adm2Full, 
                  alpha_pri = c(0, 100^2), 
                  beta_pri = c(0, sqrt(1000)), 
                  pc.bym2Phi=list(u=0.5, alpha=2/3), 
                  pc.bym2Prec=list(u=1, alpha=.1), 
                  pc.expPrec=list(u=1, alpha=.1), 
                  maxit=1000) {
  
  datMCIS = sortByCol(datMCIS, "Stratum", admStrat$NAME_FINAL)
  
  # first generate all necessary inputs if need be
  if(is.null(inputsMDM)) {
    inputsMDM = makeInputsMDM(datDHS, datMCIS, 
                              intPtsMICS=NULL, intPtsDHS=NULL, 
                              KMICS=100,
                              KDHSurb = 11, # 3 rings of 5 each
                              JInnerUrban = 3,
                              KDHSrur = 16, # 3 inner + 1 outer rings of 5 each
                              JInnerRural = 3,
                              JOuterRural = 1, 
                              admMICS=admMICS, adm2DHS=adm2DHS, 
                              tolSeq = 1e-06)
    
    thisEnv = environment()
    list2env(inputsMDM, envir=thisEnv)
  }
  
  # set priors ----
  
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Phi$u, alpha=pc.bym2Phi$alpha, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha) # get PC prior lambda for bym2 precision
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha) # get PC prior lambda for nugget precision
  
  # conditioning by Kriging from Eq (2.30) in Rue Held:
  # Ax = e (for A = (0^T 1^T), e = 0), x = (w^T u^T)^T
  # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x - e)
  # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x)
  # x* = x - (Q_x^-1 A^T A x) / sum(Q^+)
  # x* = x - (sqrt(phi/tau) Q_{+:}^+ \\ Q_{+:}^+) * sum(u) / sum(Q^+)
  # for Q_{+:}^+ = rowSums(Q^+), where * denotes the constrained version of the effect
  # Hence, we need Q_{+:}^+ / sum(Q^+):
  Qinv = bym2ArgsTMB$V %*% diag(bym2ArgsTMB$gammaTildesm1+1) %*% t(bym2ArgsTMB$V)
  QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('alpha', 'beta', 'w_bym2Star', 'u_bym2Star', 
                 'nuggetUrbMICS', 'nuggetRurMICS')
  
  # collect input data
  
  data_full = list(
    y_iUrbanMICS=ysUrbMICS, # observed binomial experiment at point i (clust)
    y_iRuralMICS=ysRurMICS, # 
    n_iUrbanMICS=nsUrbMICS, # number binomial trials
    n_iRuralMICS=nsRurMICS, # 
    # AprojUrbanMICS=AUrbMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    # AprojRuralMICS=ARurMICS, # 
    areaidxlocUrbanMICS=areaidxlocUrbanMICS, # [nIntegrationPointsUrban * nObsUrban] length vector of areal indices associated with each observation
    areaidxlocRuralMICS=areaidxlocRuralMICS, # [nIntegrationPointsRural * nObsRural] length vector of areal indices associated with each observation
    X_betaUrbanMICS=intPtsMICS$XUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralMICS=intPtsMICS$XRur, # 
    wUrbanMICS=intPtsMICS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralMICS=intPtsMICS$wRural, # 
    
    Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
    # V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
    alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
    beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
    tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
    gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
    QinvSumsNorm=QinvSumsNorm, 
    lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
    lambdaTau=lambdaTau, # determines PC prior for tau
    lambdaTauEps=lambdaTauEps, 
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  anyna = function(x) {any(is.na(x))}
  myDim = function(x) {
    dims = dim(x)
    if(is.null(dims)) {
      length(x)
    }
    else {
      dims
    }
  }
  
  if(FALSE) {
    # for testing purposes
    any(sapply(data_full, anyna))
    sapply(data_full, anyna)
    sapply(data_full, myDim)
    sapply(data_full, class)
    hist(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS, breaks=50)
    hist(data_full$y_iRuralDHS/data_full$n_iRuralDHS, breaks=50)
    hist(data_full$n_iUrbanDHS-data_full$y_iUrbanDHS, breaks=seq(-0.5, 17.5, by=1))
    hist(data_full$n_iRuralDHS-data_full$y_iRuralDHS, breaks=seq(-0.5, 24.5, by=1))
    hist(data_full$y_iRuralDHS, breaks=seq(-0.5, 16.5, by=1))
    mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1)
    mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0)
    mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1) + mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 0)
    mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0) + mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 1)
    
    hist(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, breaks=50)
    hist(data_full$y_iRuralMICS/data_full$n_iRuralMICS, breaks=50)
    hist(data_full$n_iUrbanMICS-data_full$y_iUrbanMICS, breaks=seq(-0.5, 17.5, by=1))
    hist(data_full$y_iRuralMICS, breaks=seq(-0.5, 16.5, by=1))
    mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1)
    mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0)
    mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1) + mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 0)
    mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0) + mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 1)
  }
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  tmb_params <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 0, # SPDE parameter related to the range
                     log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     alpha = initAlpha, # intercept
                     beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanMICS)-1)), 
                     w_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  # X.Int. & -1.80 & -1.96 & -1.91 & -1.68 & -1.61 \\ 
  # beta & 0.94 & 0.77 & 0.83 & 1.05 & 1.11 \\ 
  # beta.1 & -0.04 & -0.13 & -0.10 & 0.01 & 0.05 \\ 
  # beta.2 & 0.16 & 0.01 & 0.06 & 0.26 & 0.31 \\ 
  # beta.3 & -0.01 & -0.17 & -0.11 & 0.10 & 0.16 \\ 
  # beta.4 & 0.83 & 0.63 & 0.71 & 0.96 & 1.03 \\ 
  # sigmaSq & 1.36 & 1.00 & 1.10 & 1.63 & 1.81 \\ 
  # phi & 0.90 & 0.69 & 0.80 & 0.98 & 0.99 \\ 
  # sigmaEpsSq & 0.50 & 0.42 & 0.45 & 0.56 & 0.59 \\ 
  tmb_params <- list(log_tau = -.3, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 2.2, # SPDE parameter related to the range
                     log_tauEps = .7, # Log tau (i.e. log spatial precision, Epsilon)
                     alpha = -1.8, # intercept
                     beta = c(.94, -.04, .16, -.01, .83), 
                     w_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  # make TMB fun and grad ----
  # dyn.load( dynlib("code/modM_MSepsparse"))
  dyn.load( dynlib("code/modM_MSep"))
  TMB::config(tmbad.sparse_hessian_compress = 1)
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   hessian=TRUE,
                   DLL='modM_MSep')
  # objFull <- MakeADFun(data=data_full,
  #                      parameters=tmb_params,
  #                      hessian=TRUE,
  #                      DLL='modM_MSep')
  
  lower = rep(-10, length(obj[['par']]))
  upper = rep( 10, length(obj[['par']]))
  
  # make wrapper functions that print out parameters and function values
  funWrapper = function(par, badParVal=1e10) {
    if(any(par < lower) || any(par > upper)) {
      return(badParVal)
    }
    print(par)
    objVal = testObj[['fn']](par)
    parNames = names(par)
    parVals = par
    parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
    parStr = paste(parStrs, collapse=", ")
    print(paste0("objective: ", objVal, " for parameters, ", parStr))
    objVal
  }
  
  grWrapper = function(par, badParVal=1e10) {
    if(any(par < lower) || any(par > upper)) {
      return(rep(badParVal, length(par)))
    }
    print(par)
    grVal = testObj[['gr']](par)
    parNames = names(par)
    parVals = par
    parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
    parStr = paste(parStrs, collapse=", ")
    grStr = paste(grVal, collapse=", ")
    print(paste0("gradient: ", grStr, " for parameters, ", parStr))
    grVal
  }
  
  if(FALSE) {
    # testing before running optimization
    initAlphaBeta = rep(1, 6)
    initParFull = objFull$par
    initParFull[1:6] = initAlphaBeta
    objFull$fn(initParFull)
    testRep = obj$report(initParFull)
    
    system.time(test <- obj$fn(obj$par))
    # 363.236  for non-sparse
    # 355.211 for sparse
    
    system.time(test <- obj$gr(obj$par))
    # 54.5 for non-sparse
    # 55.3 for sparse
    
    testRep = obj$report()
    
    range(testRep$latentFieldUrbMICS)
    range(testRep$latentFieldRurMICS)
    range(testRep$latentFieldUrbDHS)
    range(testRep$latentFieldRurDHS)
    
    range(testRep$fe_iUrbanMICS)
    range(testRep$fe_iRuralMICS)
    range(testRep$fe_iUrbanDHS)
    range(testRep$fe_iRuralDHS)
    
    range(testRep$projepsilon_iUrbanMICS)
    range(testRep$projepsilon_iRuralMICS)
    range(testRep$projepsilon_iUrbanDHS)
    range(testRep$projepsilon_iRuralDHS)
    
    range(testRep$liksUrbMICS)
    range(testRep$liksRurMICS)
    range(testRep$liksUrbDHS)
    range(testRep$liksRurDHS)
    
    range(testRep$logDet)
    sapply(testRep, range)
    
    any(apply(data_full$wUrbanDHS, 1, function(x) {all(x == 0)}))
    any(apply(data_full$wRuralDHS, 1, function(x) {all(x == 0)}))
    any(apply(data_full$wUrbanMICS, 1, function(x) {all(x == 0)}))
    any(apply(data_full$wRuralMICS, 1, function(x) {all(x == 0)}))
  }
  
  # * Run TMB ----
  # set start of optimization
  # obj$par = optParINLA
  # obj$par = optPar
  
  {
    # tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
    tolSeq = 1e-06
    testObj = obj
    optPar = testObj$par
    startTime = proc.time()[3]
    for(thisTol in tolSeq) {
      testObj = obj
      testObj$env$inner.control = list(maxit=maxit, tol10=thisTol)
      testObj$env$tracepar = TRUE
      print(paste0("optimizing for tol = ", thisTol, "."))
      opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                    method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
      optPar = opt1$par
      if(!is.null(opt1$message)) {
        print(paste0("error for tol = ", thisTol, ". Message:"))
        print(opt1$message)
        next
      }
      else {
        print(paste0("completed optimization for tol = ", thisTol, ""))
        
        ## Get standard errors
        print("getting standard errors...")
        sdTime = system.time(
          SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                               bias.correct = TRUE,
                               bias.correct.control = list(sd = TRUE))
        )[3]
        # SD0
        print(sdTime/60)
        # 3.9447 minutes for intern=FALSE
        
        if(SD0$pdHess) {
          print("Optimization and PD hess calculation done!")
          break
        }
        else {
          print("Hessan not PD. Rerunning optimization with stricter tol...")
          
        }
      }
      
      
      
    }
    endTime = proc.time()[3]
    sdTime/60
    totalTime = endTime - startTime
    print(paste0("optimization took ", totalTime/60, " minutes"))
    # optimization took 21.7764833333333 minutes (for intern=FALSE)
  }
  
  if(FALSE) {
    badPar = c(SD0$par.fixed, SD0$par.random)
    
    tempGrad = obj$env$f(badPar,order=1)
    tempHess = obj$env$spHess(badPar,random=TRUE)
    range(tempHess)
    
    tempEig = eigen(tempHess)
    range(tempEig$values)
  }
  
  # opt0 <- nlminb(start       =    obj[['par']],
  #                objective   =    obj[['fn']],
  #                gradient    =    obj[['gr']],
  #                lower = rep(-10, length(obj[['par']])),
  #                upper = rep( 10, length(obj[['par']])),
  #                control     =    list(trace=1))
  
  # * TMB Posterior Sampling ----
  
  if(FALSE) {
    funWrapper(optPar)
    testRep = obj$report(obj$env$last.par)
    testPar = obj$env$last.par
    testPar[6:9]= 1
    testRep = obj$report(testPar)
    feUrb = testRep$fe_iUrbanMICS
    feRur = testRep$fe_iRuralMICS
    latentUrb = testRep$latentFieldUrbMICS
    latentRur = testRep$latentFieldRurMICS
    
    out = load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
    tempXUrb = intPtsMICS$XUrb[transformIUrb,]
    tempXRur = intPtsMICS$XRur[transformIRur,]
    
    pdf("figures/test/edM_Mlatent.pdf", width=6, height=6)
    zlim = range(c(latentUrb, latentRur))
    plotWithColor(tempXRur$lon, tempXRur$lat, latentRur, pch=19, cex=.3, zlim = zlim)
    plotWithColor(tempXUrb$lon, tempXUrb$lat, latentUrb, pch=19, cex=.3, add=TRUE)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_MpredUrban.pdf", width=6, height=6)
    plotWithColor(tempXRur$lon, tempXRur$lat, tempXRur$urban, pch=19, cex=.3, zlim=c(0,1))
    plotWithColor(tempXUrb$lon, tempXUrb$lat, tempXUrb$urban, pch=19, cex=.3, add=TRUE, zlim=c(0,1))
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_MempiricalProp.pdf", width=6, height=6)
    plotWithColor(jitter(tempXRur$lon, amount=.08), jitter(tempXRur$lat, amount=.08), rep(data_full$y_iRuralMICS/data_full$n_iRuralMICS, KMICS), pch=19, cex=c(data_full$wRuralMICS), zlim=c(0,1))
    plotWithColor(jitter(tempXUrb$lon, amount=.08), jitter(tempXUrb$lat, amount=.08), rep(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, KMICS), pch=19, cex=c(data_full$wUrbanMICS), add=TRUE, zlim=c(0,1))
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_MtrueUrban.pdf", width=6, height=6)
    quilt.plot(popMatNGAThresh$lon, popMatNGAThresh$lat, popMatNGAThresh$urban, nx=259, ny=212)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
  }
  
  
  
  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime)
}




