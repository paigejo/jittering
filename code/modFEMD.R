# FE + nugget model for combined DHS+MICS data
# fitFEMD: GH quadrature for nuggets, BYM2 spatial structure (mapped out for FE-only)
# Uses modMDM_BYM2_GH_v2.cpp TMB template

fitFEMD = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
                   intPtsMICS=NULL, intPtsDHS=NULL,
                   KMICS=100,
                   KDHSurb=16, JInnerUrban=4,
                   KDHSrur=21, JInnerRural=4, JOuterRural=1,
                   admMICS=admFinal, adm2DHS=adm2Full,
                   alpha_pri=c(0, 100^2),
                   beta_pri=c(0, sqrt(1000)),
                   pc.bym2Prec=list(u=0.5, alpha=2/3),
                   pc.iidPrec=list(u=1, alpha=0.5),
                   pc.expPrec=list(u=1, alpha=0.5),
                   covariates=NULL,
                   fixedEffectsOnly=TRUE,
                   maxit=1000, Qgh=25, getSDs=TRUE, verbose=FALSE,
                   doMCMC=FALSE, mcmcIter=4000, mcmcWarmup=2000, mcmcChains=1) {

  # Ensure column names
  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]; toN = nameTab[i,2]
    if(!(toN %in% names(datDHS))) datDHS[[toN]] = datDHS[[fromN]]
    if(!(toN %in% names(datMICS))) datMICS[[toN]] = datMICS[[fromN]]
  }
  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }

  # Build inputsMDM if not provided
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
  intM = intPtsMICS
  intD = intPtsDHS

  # Covariate matrices
  XUrbM = as.matrix(intM$XUrb)
  XRurM = as.matrix(intM$XRur)
  XUrbD = as.matrix(intD$covsUrb)
  XRurD = as.matrix(intD$covsRur)

  # Subset covariates if requested
  if(!is.null(covariates)) {
    keepIdxM = which(colnames(XUrbM) %in% covariates)
    XUrbM = XUrbM[, keepIdxM, drop=FALSE]
    XRurM = XRurM[, keepIdxM, drop=FALSE]
    keepIdxD = which(colnames(XUrbD) %in% covariates)
    XUrbD = XUrbD[, keepIdxD, drop=FALSE]
    XRurD = XRurD[, keepIdxD, drop=FALSE]
  }
  covNames = colnames(XUrbM)  # MICS names as canonical

  # BYM2 scaffolding (needed by template even when mapped out)
  load("savedOutput/global/admFinalMat.RData")
  bym2Args = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha,
                                        constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  nFree = ncol(bym2Args$Q) - 1

  # Priors
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

  gh = fastGHQuad::gaussHermiteData(Qgh)

  data_gh = list(
    y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
    n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
    areaidxlocUrbanMICS=as.integer(areaidxlocUrbanMICS),
    areaidxlocRuralMICS=as.integer(areaidxlocRuralMICS),
    X_betaUrbanMICS=XUrbM, X_betaRuralMICS=XRurM,
    wUrbanMICS=intM$wUrban, wRuralMICS=intM$wRural,

    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    areaidxlocUrbanDHS=as.integer(areaidxlocUrbanDHS),
    areaidxlocRuralDHS=as.integer(areaidxlocRuralDHS),
    X_betaUrbanDHS=XUrbD, X_betaRuralDHS=XRurD,
    wUrbanDHS=intD$wUrban, wRuralDHS=intD$wRural,

    Q_bym2=bym2Args$Q,
    lchoose_urban_mics=lchoose(nsUrbMICS, ysUrbMICS),
    lchoose_rural_mics=lchoose(nsRurMICS, ysRurMICS),
    lchoose_urban_dhs=lchoose(nsUrbDHS, ysUrbDHS),
    lchoose_rural_dhs=lchoose(nsRurDHS, ysRurDHS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    tr=bym2Args$tr, gammaTildesm1=bym2Args$gammaTildesm1,
    lambdaPhi=bym2Args$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    uniformPhiPrior=as.integer(FALSE),
    options=0
  )

  # Initial values
  initUrbP = sum(c(ysUrbMICS, ysUrbDHS))/sum(c(nsUrbMICS, nsUrbDHS))
  initRurP = sum(c(ysRurMICS, ysRurDHS))/sum(c(nsRurMICS, nsRurDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  nBeta = ncol(XUrbM)

  params = list(
    log_tau=0, logit_phi=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta - 1)),
    w_bym2Free=rep(0, nFree), u_bym2Free=rep(0, nFree)
  )

  if(fixedEffectsOnly) {
    mapList = list(
      w_bym2Free = factor(rep(NA, nFree)),
      u_bym2Free = factor(rep(NA, nFree)),
      log_tau = factor(NA),
      logit_phi = factor(NA)
    )
  } else {
    mapList = list()
  }

  unloadDynlibs()
  dllName = "modMDM_BYM2_GH_v2"
  if(!any(file.exists(paste0("code/", dllName, c(".o", ".so", ".dll"))))) {
    cat("Compiling code/modMDM_BYM2_GH_v2.cpp...\n")
    compile(paste0("code/", dllName, ".cpp"), framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib(paste0("code/", dllName)))

  obj = MakeADFun(data=data_gh, parameters=params, map=mapList,
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
                                 bias.correct=TRUE,
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
       covNames=covNames,
       method=if(doMCMC) "MCMC" else "GH")
}
