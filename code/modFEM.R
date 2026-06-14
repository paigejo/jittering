# FE + nugget model for MICS-only data
# fitFEM: GH quadrature for nuggets (no Laplace on nuggets)
# Uses modM_FEnug_GH.cpp TMB template

fitFEM = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
                  intPtsMICS=NULL, intPtsDHS=NULL,
                  KMICS=100,
                  KDHSurb=16, JInnerUrban=4,
                  KDHSrur=21, JInnerRural=4, JOuterRural=1,
                  admMICS=admFinal, adm2DHS=adm2Full,
                  alpha_pri=c(0, 100^2),
                  beta_pri=c(0, sqrt(1000)),
                  pc.iidPrec=list(u=1, alpha=0.5),
                  pc.expPrec=list(u=1, alpha=0.5),
                  covariates=NULL,
                  fixedEffectsOnly=TRUE,
                  maxit=1000, Qgh=25, getSDs=TRUE, verbose=FALSE,
                  doMCMC=FALSE, mcmcIter=4000, mcmcWarmup=2000, mcmcChains=1) {

  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }

  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]; toN = nameTab[i,2]
    if(!(toN %in% names(datMICS))) datMICS[[toN]] = datMICS[[fromN]]
    if(!(toN %in% names(datDHS))) datDHS[[toN]] = datDHS[[fromN]]
  }

  if(is.null(inputsMDM)) {
    inputsMDM = makeInputsMDM(datDHS, datMICS,
                              intPtsMICS=intPtsMICS, intPtsDHS=intPtsDHS,
                              KMICS=KMICS,
                              KDHSurb=KDHSurb, JInnerUrban=JInnerUrban,
                              KDHSrur=KDHSrur, JInnerRural=JInnerRural,
                              JOuterRural=JOuterRural,
                              admMICS=admMICS, adm2DHS=adm2DHS)
  }

  thisEnv = environment()
  list2env(inputsMDM, envir=thisEnv)

  allCovNames = colnames(intPtsMICS$XUrb)
  if(!is.null(covariates)) {
    keepIdx = which(allCovNames %in% covariates)
    intPtsMICS$XUrb = intPtsMICS$XUrb[, keepIdx, drop=FALSE]
    intPtsMICS$XRur = intPtsMICS$XRur[, keepIdx, drop=FALSE]
  }

  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

  nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1
  nBeta = ncol(intPtsMICS$XUrb)
  nUrb = length(ysUrbMICS)
  nRur = length(ysRurMICS)

  gh = fastGHQuad::gaussHermiteData(Qgh)

  data_gh = list(
    y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
    n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
    areaidxlocUrbanMICS=as.numeric(areaidxlocUrbanMICS),
    areaidxlocRuralMICS=as.numeric(areaidxlocRuralMICS),
    X_betaUrbanMICS=intPtsMICS$XUrb,
    X_betaRuralMICS=intPtsMICS$XRur,
    wUrbanMICS=intPtsMICS$wUrban,
    wRuralMICS=intPtsMICS$wRural,
    nAreas=nAreas,
    lchoose_urban=lchoose(nsUrbMICS, ysUrbMICS),
    lchoose_rural=lchoose(nsRurMICS, ysRurMICS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
  )

  initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
  initRurP = sum(ysRurMICS)/sum(nsRurMICS)
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha

  params_init = list(
    log_tau=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta-1)),
    u_spatial=rep(0, nAreas)
  )

  unloadDynlibs()
  dllName = "modM_FEnug_GH"
  if(!any(file.exists(paste0("code/", dllName, c(".o", ".so", ".dll"))))) {
    cat("Compiling code/modM_FEnug_GH.cpp...\n")
    compile(paste0("code/", dllName, ".cpp"), framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib(paste0("code/", dllName)))

  if(fixedEffectsOnly) {
    mapList = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))
    randList = NULL
  } else {
    mapList = list()
    randList = c("alpha", "beta", "u_spatial")
  }

  obj = MakeADFun(data=data_gh, parameters=params_init,
                  random=randList, map=mapList,
                  DLL=dllName, silent=!verbose)

  opt = NULL; SD0 = NULL; sdTime = 0
  if(!doMCMC) {
    startTime = proc.time()[3]
    opt = nlminb(obj$par, obj$fn, obj$gr,
                 control=list(eval.max=2000, iter.max=maxit, trace=as.integer(verbose)))
    totalTime = proc.time()[3] - startTime

    if(opt$convergence != 0) {
      cat("Warning: nlminb did not converge:", opt$message, "\n")
    }

    if(getSDs) {
      sdTime = system.time({
        SD0 <- try(TMB::sdreport(obj, getJointPrecision=TRUE,
                                 bias.correct = FALSE,
                                 bias.correct.control=list(sd=TRUE)), silent=TRUE)
      })[3]
      if(inherits(SD0, "try-error")) SD0 = NULL
    }
  } else {
    if(!requireNamespace("tmbstan", quietly=TRUE)) {
      stop("tmbstan package is required when doMCMC=TRUE")
    }
    startTime = proc.time()[3]
    SD0 = tmbstan::tmbstan(obj=obj, silent=!verbose, laplace=FALSE,
                           iter=mcmcIter, warmup=mcmcWarmup, chains=mcmcChains)
    totalTime = proc.time()[3] - startTime
    if(verbose) cat("MCMC took", totalTime/60, "minutes\n")
  }

  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime,
       opt=opt, inputsMDM=inputsMDM,
       covNames=colnames(intPtsMICS$XUrb),
       method=if(doMCMC) "MCMC" else "GH")
}

# Backward compatibility
fitMFEM = fitFEM
