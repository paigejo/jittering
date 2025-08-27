
fitMD = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL, 
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
                 maxit=1000, repar=TRUE) {
  
  # make sure Stratum variable exists in MICS data
  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }
  
  # make sure all the names in our datasets we expect to be there are (convert over variables as need be)
  nameTab = rbind(c("N", "ns"), 
                  c("N", "n"), 
                  c("Z", "ys"), 
                  c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    theseNames = nameTab[i,]
    fromN = theseNames[1]
    toN = theseNames[2]
    if(!(toN %in% names(datMICS))) {
      datMICS[[toN]] = datMICS[[fromN]]
    }
    if(!(toN %in% names(datDHS))) {
      datDHS[[toN]] = datDHS[[fromN]]
    }
  }
  
  # first generate all necessary inputs if need be
  if(is.null(inputsMDM)) {
    print("Making M_DM inputs...")
    
    inputsMDM = makeInputsMDM(datDHS, datMICS, 
                              intPtsMICS=intPtsMICS, intPtsDHS=intPtsDHS, 
                              KMICS=KMICS,
                              KDHSurb = KDHSurb, # 3 rings of 5 each
                              JInnerUrban = JInnerUrban,
                              KDHSrur = KDHSrur, # 3 inner + 1 outer rings of 5 each
                              JInnerRural = JInnerRural,
                              JOuterRural = JOuterRural, 
                              admMICS=admMICS, adm2DHS=adm2DHS)
    
  }
  
  thisEnv = environment()
  list2env(inputsMDM, envir=thisEnv)
  
  # set priors ----
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Phi$u, alpha=pc.bym2Phi$alpha, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha) # get PC prior lambda for bym2 precision
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha) # get PC prior lambda for nugget precision
  
  if(!repar) {
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
  }
  
  # Specify inputs for TMB ----
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  if(!repar) {
    rand_effs <- c('alpha', 'beta', 'w_bym2Star', 'u_bym2Star', 
                   'nuggetUrbDHS', 'nuggetRurDHS')
  } else {
    rand_effs <- c('beta', 'w_bym2Star', 'u_bym2Star', 
                   'nuggetUrbDHS', 'nuggetRurDHS')
  }
  
  
  # collect input data
  
  if(!repar) {
    data_full = list(
      y_iUrbanDHS=ysUrbDHS, # observed binomial experiment at point i (clust)
      y_iRuralDHS=ysRurDHS, # 
      n_iUrbanDHS=nsUrbDHS, # number binomial trials
      n_iRuralDHS=nsRurDHS, # 
      # AprojUrbanMICS=AUrbMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
      # AprojRuralMICS=ARurMICS, # 
      areaidxlocUrban=areaidxlocUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] length vector of areal indices associated with each observation
      areaidxlocRural=areaidxlocRuralDHS, # [nIntegrationPointsRural * nObsRural] length vector of areal indices associated with each observation
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
  } else {
    # include everything from sep case except QinvSumsNorm and alpha_pri
    data_full = list(
      y_iUrbanDHS=ysUrbDHS, # observed binomial experiment at point i (clust)
      y_iRuralDHS=ysRurDHS, # 
      n_iUrbanDHS=nsUrbDHS, # number binomial trials
      n_iRuralDHS=nsRurDHS, # 
      # AprojUrbanMICS=AUrbMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
      # AprojRuralMICS=ARurMICS, # 
      areaidxlocUrban=areaidxlocUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] length vector of areal indices associated with each observation
      areaidxlocRural=areaidxlocRuralDHS, # [nIntegrationPointsRural * nObsRural] length vector of areal indices associated with each observation
      X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
      X_betaRuralDHS=intPtsDHS$covsRur, # 
      wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
      wRuralDHS=intPtsDHS$wRural, # 
      
      Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
      # V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
      beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
      tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
      gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
      lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
      lambdaTau=lambdaTau, # determines PC prior for tau
      lambdaTauEps=lambdaTauEps, 
      options=0 # 1 for adreport of log tau and logit phi
    )
  }
  
  
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
  initUrbP = sum(c(data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  if(!repar) {
    tmb_params <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                       logit_phi = 0, # SPDE parameter related to the range
                       log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                       alpha = initAlpha, # intercept
                       beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanDHS)-1)), 
                       w_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                       u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                       nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                       nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
    )
  } else {
    # don't include alpha, but incorporate it into w_bym2Star
    tmb_params <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                       logit_phi = 0, # SPDE parameter related to the range
                       log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                       beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanDHS)-1)), 
                       w_bym2Star = rep(initAlpha, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                       u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                       nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                       nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
    )
  }
  
  
  # make TMB fun and grad ----
  # dyn.load( dynlib("code/modM_DSepsparse"))
  
  # first make sure no TMB dynlib is loaded
  unloadDynlibs()
  
  if(!repar) {
    
    # compile dynlib if need be
    if(!file.exists("code/modM_DSep.so")) {
      print("compiling code/modM_DSep.cpp...")
      compile( "code/modM_DSep.cpp", 
               framework="TMBad", safebounds=FALSE)
    }
    
    # load dynlib, make tmb obj
    dyn.load( dynlib("code/modM_DSep"))
    TMB::config(tmbad.sparse_hessian_compress = 1)
    obj <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     random=rand_effs,
                     hessian=TRUE,
                     DLL='modM_DSep')
  } else {
    
    # compile dynlib if need be
    if(!file.exists("code/modM_DSepRepar.so")) {
      print("compiling code/modM_DSepRepar.cpp...")
      compile( "code/modM_DSepRepar.cpp", 
               framework="TMBad", safebounds=FALSE)
    }
    
    # load dynlib, make tmb obj
    dyn.load( dynlib("code/modM_DSepRepar"))
    TMB::config(tmbad.sparse_hessian_compress = 1)
    obj <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     random=rand_effs,
                     hessian=TRUE,
                     DLL='modM_DSepRepar')
  }
  
  # objFull <- MakeADFun(data=data_full,
  #                      parameters=tmb_params,
  #                      hessian=TRUE,
  #                      DLL='modM_DSepRepar')
  
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
    objFull <- MakeADFun(data=data_full,
                         parameters=tmb_params,
                         hessian=TRUE,
                         DLL='modM_DSepRepar')
    
    initBeta = rep(1, 5)
    initParFull = objFull$par
    initParFull[names(initParFull) == "beta"] = initBeta
    objFull$fn(initParFull)
    testRepMD = obj$report(initParFull)
    
    system.time(test <- obj$fn(obj$par))
    # 363.236  for non-sparse
    # 355.211 for sparse
    
    system.time(test <- obj$gr(obj$par))
    # 54.5 for non-sparse
    # 55.3 for sparse
    
    load("~/git/jittering/savedOutput/test/testRepMDM.RData")
    # range(testRepMDM$latentFieldUrbMICS)
    # range(testRepMD$latentFieldUrbMICS)
    # range(testRepMDM$latentFieldRurMICS)
    # range(testRepMD$latentFieldRurMICS)
    range(testRepMDM$latentFieldUrbDHS)
    range(testRepMD$latentFieldUrbDHS)
    range(testRepMDM$latentFieldRurDHS)
    range(testRepMD$latentFieldRurDHS)
    # range(testRepMDM$fe_iUrbanMICS)
    # range(testRepMD$fe_iUrbanMICS)
    # range(testRepMDM$fe_iRuralMICS)
    # range(testRepMD$fe_iRuralMICS)
    range(testRepMDM$fe_iUrbanDHS)
    range(testRepMD$fe_iUrbanDHS)
    range(testRepMDM$fe_iRuralDHS)
    range(testRepMD$fe_iRuralDHS)
    # range(testRepMDM$projepsilon_iUrbanMICS)
    # range(testRepMD$projepsilon_iUrbanMICS)
    # range(testRepMDM$projepsilon_iRuralMICS)
    # range(testRepMD$projepsilon_iRuralMICS)
    range(testRepMDM$projepsilon_iUrbanDHS)
    range(testRepMD$projepsilon_iUrbanDHS)
    range(testRepMDM$projepsilon_iRuralDHS)
    range(testRepMD$projepsilon_iRuralDHS)
    # range(testRepMDM$liksUrbMICS)
    # range(testRepMD$liksUrbMICS)
    # range(testRepMDM$liksRurMICS)
    # range(testRepMD$liksRurMICS)
    range(testRepMDM$liksUrbDHS)
    range(testRepMD$liksUrbDHS)
    range(testRepMDM$liksRurDHS)
    range(testRepMD$liksRurDHS)
    
    out=load("savedOutput/global/intPtsMICS_100.RData")
    intPtsMICS = straightenMICS(intPtsMICS)
    
    
    MICStoPtMat = function(matMICS, urb=TRUE) {
      
      # aggOut = aggregate(datMICS$Stratum, by=list(datMICS$Stratum), FUN=length)
      # # all.equal(aggOut$Group.1, intPtsMICS$strataMICS)
      # # TRUE
      
      # get number of integration points (urban and rural) per MICS stratum
      allNumPerStrat = aggregate(datMICS$Stratum, by=list(strat=datMICS$Stratum, urb=datMICS$urban), FUN=length, drop=FALSE)
      allNumPerStrat = straightenNumPerStrat(allNumPerStrat, admFinal$NAME_FINAL)
      numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
      numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
      numPerStratRur[is.na(numPerStratRur)] = 0
      numPerStratUrb[is.na(numPerStratUrb)] = 0
      
      if(urb) {
        thisNumPerStrat = numPerStratUrb
      } else {
        thisNumPerStrat = numPerStratRur
      }
      
      startStratInds = which(matMICS$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
      nAreas = nrow(matMICS)/KMICS
      areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=thisNumPerStrat[x])})) # length nUrb, range = 1:41. gives area index for each obs
      allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
      nObs = length(allAreaIs)/KMICS
      allIntIs = rep(1:KMICS, each=nObs) # length nObs*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
      transformI = allAreaIs + (allIntIs-1)*nAreas
      matMICS = matMICS[transformI,] # now XUrb is [K * nObs] x nVar
      
      matMICS
    }
    
    covsUrb = MICStoPtMat(intPtsMICS$XUrb, TRUE)
    covsRur = MICStoPtMat(intPtsMICS$XRur, FALSE)
    plotWithColor(covsUrb$east, covsUrb$north, data_full$X_betaUrbanMICS[,5], pch=19, cex=.5)
    plotWithColor(covsUrb$east, covsUrb$north, data_full$areaidxlocUrbanMICS, pch=19, cex=.5)
    plotWithColor(covsUrb$east, covsUrb$north, data_full$areaidxlocUrbanMICS, pch=19, cex=.5, 
                  zlim=c(0, 40))
    plotWithColor(covsRur$east, covsRur$north, data_full$areaidxlocRuralMICS, pch=19, cex=.5, 
                  zlim=c(0, 40), add=TRUE)
    
    plotWithColor(intPtsMICS$XUrb$east, intPtsMICS$XUrb$north, rep(1:41, KMICS), pch=19, cex=.5)
    plotWithColor(intPtsMICS$XRur$east, intPtsMICS$XRur$north, rep(1:41, KMICS), pch=19, cex=.5)
    
    areaI = 1
    areaSeq = seq(1, 4100, by=41) + areaI-1
    plotWithColor(intPtsMICS$XRur$east[areaSeq], intPtsMICS$XRur$north[areaSeq], intPtsMICS$XRur$pop[areaSeq], 
                  pch=19, cex=.5)
    plotWithColor(intPtsMICS$XUrb$east[areaSeq], intPtsMICS$XUrb$north[areaSeq], intPtsMICS$XUrb$pop[areaSeq], 
                  pch=19, cex=.5)
    
    plotWithColor(intPtsMICS$XRur$east[areaSeq], intPtsMICS$XRur$north[areaSeq], intPtsMICS$XRur$pop[areaSeq], 
                  pch=19, cex=.5)
    plotWithColor(intPtsMICS$XRur$east[areaSeq], intPtsMICS$XRur$north[areaSeq], intPtsMICS$wRural[areaI,], 
                  pch=19, cex=.5)
    plotWithColor(intPtsMICS$XUrb$east[areaSeq], intPtsMICS$XUrb$north[areaSeq], intPtsMICS$XUrb$pop[areaSeq], 
                  pch=19, cex=.5)
    plotWithColor(intPtsMICS$XUrb$east[areaSeq], intPtsMICS$XUrb$north[areaSeq], intPtsMICS$wUrban[areaI,], 
                  pch=19, cex=.5)
    
    
    plotWithColor(covsUrb$east, covsUrb$north, data_full$X_betaUrbanMICS %*% rep(1,5), pch=19, cex=.5)
    plotWithColor(covsUrb$east, covsUrb$north, testRepMDM$fe_iUrbanMICS, pch=19, cex=.5)
    plotWithColor(covsRur$east, covsRur$north, data_full$X_betaRuralMICS %*% rep(1,5), pch=19, cex=.5)
    plotWithColor(covsRur$east, covsRur$north, testRepMDM$fe_iRuralMICS, pch=19, cex=.5)
    
    all.equal(as.numeric(data_full$X_betaUrbanMICS %*% rep(1,5)), as.numeric(testRepMDM$fe_iUrbanMICS))
    all.equal(as.numeric(data_full$X_betaRuralMICS %*% rep(1,5)), as.numeric(testRepMDM$fe_iRuralMICS))
    all.equal(as.numeric(data_full$X_betaUrbanDHS %*% rep(1,5)), as.numeric(testRepMDM$fe_iUrbanDHS))
    all.equal(as.numeric(data_full$X_betaRuralDHS %*% rep(1,5)), as.numeric(testRepMDM$fe_iRuralDHS))
    
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), intPtsDHS$covsUrb %*% rep(1,5), pch=19, cex=.5)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), testRepMDM$fe_iUrbanDHS, pch=19, cex=.5)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), intPtsDHS$covsRur %*% rep(1,5), pch=19, cex=.5)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), testRepMDM$fe_iRuralDHS, pch=19, cex=.5)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), data_full$areaidxlocUrbanDHS, pch=19, cex=.5, 
                  zlim=c(0, 40))
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), data_full$areaidxlocRuralDHS, pch=19, cex=.5, 
                  zlim=c(0, 40), add=TRUE)
    
    load("~/git/jittering/savedOutput/test/data_fullMDM.RData")
    data_fullMD = data_full
    compNames = names(data_fullMD)
    for(i in 1:length(compNames)) {
      print(compNames[i])
      mdmName = compNames[i]
      if(grepl("areaidxloc", mdmName)) {
        mdmName = paste0(mdmName, "DHS")
      }
      print(all(c(as.matrix(data_fullMD[[compNames[i]]])) == c(as.matrix(data_fullMDM[[mdmName]]))))
    }
    
    compNames = names(testRepMD)
    for(i in 1:length(compNames)) {
      print(compNames[i])
      mdmName = compNames[i]
      # if(grepl("areaidxloc", mdmName)) {
      #   mdmName = paste0(mdmName, "DHS")
      # }
      print(all(c(testRepMD[[compNames[i]]]) == c(testRepMDM[[mdmName]])))
    }
    testRepMD$bym2LogLik
    testRepMDM$bym2LogLik
    
    head(testRepMD$fe_iRuralDHS)
    head(testRepMDM$fe_iRuralDHS)
    max(abs(testRepMD$fe_iRuralDHS - testRepMDM$fe_iRuralDHS))
    max(abs(testRepMD$fe_iUrbanDHS - testRepMDM$fe_iUrbanDHS))
    max(abs(testRepMD$latentFieldRurDHS - testRepMDM$latentFieldRurDHS))
    max(abs(testRepMD$latentFieldUrbDHS - testRepMDM$latentFieldUrbDHS))
    max(abs(testRepMD$w_bym2Star - testRepMDM$w_bym2Star))
    max(abs(testRepMD$liksRurDHS - testRepMDM$liksRurDHS))
    max(abs(testRepMD$liksUrbDHS - testRepMDM$liksUrbDHS))
    
    hist(intPtsMICS$XUrb$pop[seq(1, 4100, by=41)+39], breaks=30)
    hist(covsUrb$pop[data_full$areaidxlocUrbanMICS == 39], breaks=30)
    
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
  browser()
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
    # optimization took 2.36478333333313 minutes (for intern=FALSE)
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
    intPtsMICS = straightenMICS(intPtsMICS)
    tempXUrb = intPtsMICS$XUrb[transformIUrb,]
    tempXRur = intPtsMICS$XRur[transformIRur,]
    
    pdf("figures/test/edM_Dlatent.pdf", width=6, height=6)
    zlim = range(c(latentUrb, latentRur))
    plotWithColor(tempXRur$lon, tempXRur$lat, latentRur, pch=19, cex=.3, zlim = zlim)
    plotWithColor(tempXUrb$lon, tempXUrb$lat, latentUrb, pch=19, cex=.3, add=TRUE)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_DpredUrban.pdf", width=6, height=6)
    plotWithColor(tempXRur$lon, tempXRur$lat, tempXRur$urban, pch=19, cex=.3, zlim=c(0,1))
    plotWithColor(tempXUrb$lon, tempXUrb$lat, tempXUrb$urban, pch=19, cex=.3, add=TRUE, zlim=c(0,1))
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_DempiricalProp.pdf", width=6, height=6)
    plotWithColor(jitter(tempXRur$lon, amount=.08), jitter(tempXRur$lat, amount=.08), rep(data_full$y_iRuralMICS/data_full$n_iRuralMICS, KMICS), pch=19, cex=c(data_full$wRuralMICS), zlim=c(0,1))
    plotWithColor(jitter(tempXUrb$lon, amount=.08), jitter(tempXUrb$lat, amount=.08), rep(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, KMICS), pch=19, cex=c(data_full$wUrbanMICS), add=TRUE, zlim=c(0,1))
    plotMapDat(admFinal, new=FALSE)
    dev.off()
    
    pdf("figures/test/edM_DtrueUrban.pdf", width=6, height=6)
    quilt.plot(popMatNGAThresh$lon, popMatNGAThresh$lat, popMatNGAThresh$urban, nx=259, ny=212)
    plotMapDat(admFinal, new=FALSE)
    dev.off()
  }
  
  
  
  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime)
}




