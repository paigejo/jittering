# IID spatial effects only model for M_M (MICS-only) data
# Uses same integration point structure as the full models
# Replaces BYM2 with IID spatial random effects + cluster-level nugget

fitMM_IIDonly = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
                         intPtsMICS=NULL, intPtsDHS=NULL,
                         KMICS=100,
                         KDHSurb = 11,
                         JInnerUrban = 3,
                         KDHSrur = 16,
                         JInnerRural = 3,
                         JOuterRural = 1, admMICS=admFinal, adm2DHS=adm2Full,
                         alpha_pri = c(0, 100^2),
                         beta_pri = c(0, sqrt(1000)),
                         pc.iidPrec=list(u=1, alpha=0.5),
                         pc.expPrec=list(u=1, alpha=0.5),
                         covariates=NULL,
                         fixedEffectsOnly=FALSE,
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
  
  # generate all necessary inputs if need be
  if(is.null(inputsMDM)) {
    print("Making M_DM inputs for IID model...")
    
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
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)
  
  # number of areas from area indices (0-indexed)
  nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1
  
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
  if(fixedEffectsOnly) {
    rand_effs <- c('nuggetUrbMICS', 'nuggetRurMICS')
    mapList <- list(u_spatial = factor(rep(NA, nAreas)),
                    log_tau   = factor(NA))
  } else {
    rand_effs <- c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS')
    mapList <- list()
  }
  
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
    nAreas=nAreas,
    alpha_pri=alpha_pri,
    beta_pri=beta_pri,
    lambdaTau=lambdaTau,
    lambdaTauEps=lambdaTauEps
  )
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS))/sum(c(data_full$n_iUrbanMICS))
  initRurP = sum(c(data_full$y_iRuralMICS))/sum(c(data_full$n_iRuralMICS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  tmb_params <- list(log_tau = 0,
                     log_tauEps = 0,
                     alpha = initAlpha,
                     beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanMICS)-1)),
                     u_spatial = rep(0, nAreas),
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS)))
  
  # make TMB fun and grad ----
  unloadDynlibs()
  
  if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
    print("compiling code/modM_MIIDonly.cpp...")
    compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
  }
  
  dyn.load(dynlib("code/modM_MIIDonly"))
  TMB::config(tmbad.sparse_hessian_compress = 1)
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   map=mapList,
                   hessian=TRUE,
                   DLL='modM_MIIDonly')
  
  lower = rep(-10, length(obj[['par']]))
  upper = rep( 10, length(obj[['par']]))
  
  # wrapper functions that print out parameters and function values
  funWrapper = function(par, badParVal=1e10) {
    if(any(par < lower) || any(par > upper)) {
      return(badParVal)
    }
    objVal = obj[['fn']](par)
    print(paste0("objective: ", objVal, " for log_tau=", par))
    objVal
  }
  
  grWrapper = function(par, badParVal=1e10) {
    if(any(par < lower) || any(par > upper)) {
      return(rep(badParVal, length(par)))
    }
    obj[['gr']](par)
  }
  
  # Run TMB ----
  if(!doMCMC) {
    optPar = obj$par
    startTime = proc.time()[3]
    for(thisTol in tolSeq) {
      obj$env$inner.control = list(maxit=maxit, tol10=thisTol)
      obj$env$tracepar = TRUE
      print(paste0("optimizing for tol = ", thisTol, "."))
      opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                    method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
      optPar = opt1$par
      if(!is.null(opt1$message)) {
        print(paste0("error for tol = ", thisTol, ". Message:"))
        print(opt1$message)
        next
      } else {
        print(paste0("completed optimization for tol = ", thisTol, ""))
        
        if(getSDs) {
          print("getting standard errors...")
          sdTime = system.time(
            SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                 bias.correct = TRUE,
                                 bias.correct.control = list(sd = TRUE))
          )[3]
          print(sdTime/60)
          
          if(SD0$pdHess) {
            print("Optimization and PD hess calculation done!")
            break
          } else {
            print("Hessian not PD. Rerunning optimization with stricter tol...")
          }
        } else {
          print("getting report (no standard errors)...")
          sdTime = system.time(
            SD0 <- obj$report()
          )[3]
        }
      }
    }
    endTime = proc.time()[3]
    totalTime = endTime - startTime
    print(paste0("optimization took ", totalTime/60, " minutes"))
  } else {
    startTime = proc.time()[3]
    fit <- tmbstan(obj=obj, silent=FALSE, laplace=FALSE,
                   iter=4000, warmup=2000, chains=1)
    endTime = proc.time()[3]
    totalTime = endTime - startTime
    sdTime = 0
    SD0 = fit
    print(paste0("MCMC took ", totalTime/60, " minutes"))
  }
  
  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime)
}
