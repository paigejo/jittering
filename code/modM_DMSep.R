
fitMDM = function(datDHS, datMCIS, inputsMDM=NULL, 
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
                 'nuggetUrbMICS', 'nuggetRurMICS', 'nuggetUrbDHS', 'nuggetRurDHS')
  
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
    
    y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
    y_iRuralDHS=ysRurDHS, # 
    n_iUrbanDHS=nsUrbDHS, # number binomial trials
    n_iRuralDHS=nsRurDHS, # 
    # AprojUrbanDHS=AUrbDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    # AprojRuralDHS=ARurDHS, # 
    areaidxlocUrbanDHS=areaidxlocUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] length vector of areal indices associated with each observation
    areaidxlocRuralDHS=areaidxlocRuralDHS, # [nIntegrationPointsRural * nObsRural] length vector of areal indices associated with each observation
    X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralDHS=intPtsDHS$covsRur, # 
    wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralDHS=intPtsDHS$wRural, # 
    
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
                     beta = c(initBeta1, rep(0, ncol(intPtsDHS$covsUrb)-1)), 
                     w_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS)), 
                     nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                     nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
  )
  
  # make TMB fun and grad ----
  # dyn.load( dynlib("code/modM_DMSepsparse"))
  dyn.load( dynlib("code/modM_DMSep"))
  TMB::config(tmbad.sparse_hessian_compress = 1)
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   hessian=TRUE,
                   DLL='modM_DMSep')
  # objFull <- MakeADFun(data=data_full,
  #                      parameters=tmb_params,
  #                      hessian=TRUE,
  #                      DLL='modM_DMSep')
  
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
    # optimization took 22.9861833333334 minutes (for intern=FALSE)
  }
  
  if(FALSE) {
    badPar = c(SD0$par.fixed, SD0$par.random)
    
    tempGrad = obj$env$f(badPar,order=1)
    tempHess = obj$env$spHess(badPar,random=TRUE)
    range(tempHess)
    
    tempEig = eigen(tempHess)
    range(tempEig$values)
  }
  
  if(!SD0$pdHess) {
    thisTMBpar = tmb_params
    thisTMBpar$log_tauEps = optPar[names(optPar) == "log_tauEps"]
    
    map = as.list(factor(c(alpha=1, beta=2:6, log_tau=7, logit_phi=8, log_tauEps=NA)))
    temp = unlist(map[grepl("beta", names(map))])
    map[grepl("beta", names(map))] = NULL
    map$beta = temp
    map=list(factor(c(log_tauEps=NA)))
    names(map) = "log_tauEps"
    objFixed <- MakeADFun(data=data_full,
                          parameters=tmb_params,
                          random=rand_effs,
                          map=map, 
                          hessian=TRUE,
                          DLL='modBYM2JitterFusionNugget')
    testObj = objFixed
    thisOptPar = optPar[-which(names(optPar) == "log_tauEps")]
    lower = lower[-which(names(optPar) == "log_tauEps")]
    upper = upper[-which(names(optPar) == "log_tauEps")]
    
    newStartTime = proc.time()[3]
    for(thisTol in tolSeq) {
      testObj = objFixed
      testObj$env$inner.control = list(maxit=maxit, tol10=thisTol)
      testObj$env$tracepar = TRUE
      print(paste0("optimizing for tol = ", thisTol, "."))
      opt1 <- optim(par=thisOptPar, fn=funWrapper, gr=grWrapper,
                    method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
      # opt1 <- optim(par=optPar, fn=funWrapper, 
      #               hessian = FALSE, control=list(reltol=thisTol))
      thisOptPar = opt1$par
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
        # ~97? minutes minutes for intern=FALSE
        
        if(SD0$pdHess) {
          print("Optimization and PD hess calculation done!")
          break
        }
        else {
          print("Hessan not PD. Rerunning optimization with stricter tol...")
          
        }
      }
      
      
      
    }
    thisEndTime = proc.time()[3]
    sdTime/60
    thisTotalTime = thisEndTime - newStartTime
    totalTime = thisEndTime - startTime
    print(paste0("second optimization took ", thisTotalTime/60, " minutes"))
    print(paste0("all optimization took ", totalTime/60, " minutes"))
    # 102.33995 minutes
  }
  
  # opt0 <- nlminb(start       =    obj[['par']],
  #                objective   =    obj[['fn']],
  #                gradient    =    obj[['gr']],
  #                lower = rep(-10, length(obj[['par']])),
  #                upper = rep( 10, length(obj[['par']])),
  #                control     =    list(trace=1))
  
  # * TMB Posterior Sampling ----
  
  if(FALSE) {
    hessTime = system.time(testHess <- hessian(testObj$fn, optPar))[3]
    hessTime/60
    # 10.28792 minutes for intern=FALSE
    eig = eigen(testHess)
    eig
    
    for(i in 1:length(eig$values)) {
      thisVal = eig$values[i]
      thisVec = eig$vectors[,i]
      barplot(eig$vectors[,i], names.arg=names(SD0$par.fixed), 
              main=paste0("Eigenvector ", i, " with value ", thisVal))
    }
    
    testObj$fn(SD0$par.fixed)
    testObj$fn(SD0$par.fixed + eig$vectors[,9]*.05)
    testObj$fn(SD0$par.fixed - eig$vectors[,9]*.05)
  }
  
  
  
  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime)
}




