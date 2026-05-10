# Marginalized BYM2 model for M_M (MICS-only) data
# Uses single n-dimensional Epsilon_bym2 with eigendecomposition
# Eliminates the (w,u) funnel from the separable parameterization

fitMMMarg = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL, 
                     intPtsMICS=NULL, intPtsDHS=NULL, 
                     KMICS=100,
                     KDHSurb = 11,
                     JInnerUrban = 3,
                     KDHSrur = 16,
                     JInnerRural = 3,
                     JOuterRural = 1, admMICS=admFinal, adm2DHS=adm2Full, 
                     alpha_pri = c(0, 100^2), 
                     beta_pri = c(0, sqrt(1000)), 
                     pc.bym2Phi=list(u=0.5, alpha=1/3), 
                     pc.bym2Prec=list(u=1, alpha=0.5), 
                     pc.expPrec=list(u=1, alpha=0.5), 
                     covariates=NULL,
                     fixedHyperparams=NULL,
                     maxit=1000, adm2AsCovariate=FALSE, 
                     tolSeq = c(1e-06), getSDs=TRUE, doMCMC=FALSE) {
  
  # make sure Stratum variable exists in MICS data
  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }
  
  # make sure all the names in our datasets we expect to be there are
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
                              KDHSurb = KDHSurb,
                              JInnerUrban = JInnerUrban,
                              KDHSrur = KDHSrur,
                              JInnerRural = JInnerRural,
                              JOuterRural = JOuterRural, 
                              admMICS=admMICS, adm2DHS=adm2DHS, 
                              adm2AsCovariate=adm2AsCovariate)
    
  }
  
  thisEnv = environment()
  list2env(inputsMDM, envir=thisEnv)
  
  datMICS = sortByCol(datMICS, "Stratum", admMICS$NAME_FINAL)
  
  # set priors ----
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Phi$u, alpha=pc.bym2Phi$alpha, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)
  
  nAreas = ncol(bym2ArgsTMB$Q)
  
  # Subset covariates if requested ----
  allCovNames = colnames(intPtsMICS$XUrb)
  if(!is.null(covariates)) {
    keepIdx = which(allCovNames %in% covariates)
    cat("  Subsetting covariates to:", allCovNames[keepIdx], "\n")
    intPtsMICS$XUrb = intPtsMICS$XUrb[, keepIdx, drop=FALSE]
    intPtsMICS$XRur = intPtsMICS$XRur[, keepIdx, drop=FALSE]
  } else {
    cat("  Using all covariates:", allCovNames, "\n")
  }
  
  # Specify inputs for TMB ----
  
  ## specify random effects (alpha is explicit, not absorbed)
  rand_effs <- c('alpha', 'beta', 'Epsilon_bym2', 
                 'nuggetUrbMICS', 'nuggetRurMICS')
  
  # collect input data: use V_bym2 (eigenvectors) instead of Q_bym2
  data_full = list(
    y_iUrbanMICS=ysUrbMICS,
    y_iRuralMICS=ysRurMICS,
    n_iUrbanMICS=nsUrbMICS,
    n_iRuralMICS=nsRurMICS,
    areaidxlocUrbanMICS=areaidxlocUrbanMICS,
    areaidxlocRuralMICS=areaidxlocRuralMICS,
    X_betaUrbanMICS=intPtsMICS$XUrb,
    X_betaRuralMICS=intPtsMICS$XRur,
    wUrbanMICS=intPtsMICS$wUrban,
    wRuralMICS=intPtsMICS$wRural,
    V_bym2=bym2ArgsTMB$V,       # eigenvectors of Q (NOT the sparse Q itself)
    alpha_pri=alpha_pri,
    beta_pri=beta_pri,
    tr=bym2ArgsTMB$tr,
    gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
    lambdaPhi=bym2ArgsTMB$lambda,
    lambdaTau=lambdaTau,
    lambdaTauEps=lambdaTauEps, 
    options=0
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
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS))/sum(c(data_full$n_iUrbanMICS))
  initRurP = sum(c(data_full$y_iRuralMICS))/sum(c(data_full$n_iRuralMICS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  tmb_params <- list(log_tau = 0,
                     logit_phi = 0,
                     log_tauEps = 0,
                     alpha = initAlpha,
                     beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanMICS)-1)), 
                     Epsilon_bym2 = rep(0, nAreas),
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  # Handle fixed hyperparameters via map ----
  mapList <- list()
  if(!is.null(fixedHyperparams)) {
    for(pName in names(fixedHyperparams)) {
      tmb_params[[pName]] <- fixedHyperparams[[pName]]
      mapList[[pName]] <- factor(NA)
      cat("  Fixing", pName, "=", fixedHyperparams[[pName]], "\n")
    }
  }
  
  # make TMB fun and grad ----
  
  # first make sure no TMB dynlib is loaded
  unloadDynlibs()
  
  # compile dynlib if need be
  if(!any(file.exists(paste0("code/modM_MSepMarg", c(".o", ".so", ".dll"))))) {
    print("compiling code/modM_MSepMarg.cpp...")
    compile( "code/modM_MSepMarg.cpp", 
             framework="TMBad", safebounds=FALSE)
  }
  
  # load dynlib, make tmb obj
  dyn.load( dynlib("code/modM_MSepMarg"))
  TMB::config(tmbad.sparse_hessian_compress = 1)
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   map=mapList,
                   hessian=TRUE,
                   DLL='modM_MSepMarg')
  
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
  
  # * Run TMB ----
  
  if(!doMCMC) {
    {
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
          
          if(getSDs) {
            ## Get standard errors
            print("getting standard errors...")
            sdTime = system.time(
              SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                   bias.correct = TRUE,
                                   bias.correct.control = list(sd = TRUE))
            )[3]
            
            print(sdTime/60)
            
            if(SD0$pdHess) {
              print("Optimization and PD hess calculation done!")
              break
            }
            else {
              print("Hessan not PD. Rerunning optimization with stricter tol...")
            }
          } else {
            ## Get report without computing standard errors
            print("getting report (no standard errors)...")
            sdTime = system.time(
              SD0 <- testObj$report()
            )[3]
            print(paste0("report time: ", sdTime/60, " minutes"))
          }
        }
      }
      endTime = proc.time()[3]
      sdTime/60
      totalTime = endTime - startTime
      print(paste0("optimization took ", totalTime/60, " minutes"))
    }
  } else {
    {
      testObj = obj
      optPar = testObj$par
      startTime = proc.time()[3]
      
      # Run MCMC with adequate iterations for proper inference
      fit <- tmbstan(obj=obj, silent=FALSE, laplace=FALSE, 
                     iter=4000, warmup=2000, chains=1)
      
      endTime = proc.time()[3]
      totalTime = endTime - startTime
      print(paste0("MCMC took ", totalTime/60, " minutes"))
    }
  }
  
  
  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime)
}
