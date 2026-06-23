# FE + nugget DHS+MICS fusion model with COMMENSURATE PRIORS on covariate effects.
# fitFEMD_comm: GH quadrature for nuggets, separate (alpha, beta) for MICS vs DHS,
# commensurate prior (beta_M - beta_D) ~ N(0, sigma_comm^2 I) with PC prior on sigma_comm.
# Predictions in the paper use alpha (DHS) and beta (DHS).
#
# Hyperparameter sigma_comm is selected by empirical Bayes: alpha, alpha_M, beta, beta_M
# are random effects (Laplace-marginalized) so the outer optimizer maximizes the marginal
# likelihood over (log_tauEps, log_sigma_comm) [+ BYM2 hyperparams when fixedEffectsOnly=FALSE].

fitFEMD_comm = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
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
                        pc.sigmaComm=list(u=1, alpha=0.1),
                        commIntercept=0L, commSlope=1L,  # default = original (indep intercept, comm slope)
                        initParams=NULL,                 # optional named overrides of initial params (multi-start)
                        fixSigmaComm=NULL,               # if set: pin sigma_comm at this value (forces alpha_M~alpha, beta_M~beta when both comm flags on) -> nesting/bug check
                        marginalizeFE=TRUE,              # TRUE: alpha/beta(/_M) Laplace-marginalized (random); FALSE: fixed outer params (sim-study style)
                        covariates=NULL,
                        fixedEffectsOnly=TRUE,
                        maxit=1000, Qgh=25, getSDs=TRUE, verbose=FALSE,
                        doMCMC=FALSE, mcmcIter=4000, mcmcWarmup=2000, mcmcChains=1) {

  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]; toN = nameTab[i,2]
    if(!(toN %in% names(datDHS))) datDHS[[toN]] = datDHS[[fromN]]
    if(!(toN %in% names(datMICS))) datMICS[[toN]] = datMICS[[fromN]]
  }
  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
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
  intM = intPtsMICS
  intD = intPtsDHS

  XUrbM = as.matrix(intM$XUrb)
  XRurM = as.matrix(intM$XRur)
  XUrbD = as.matrix(intD$covsUrb)
  XRurD = as.matrix(intD$covsRur)

  if(!is.null(covariates)) {
    keepIdxM = which(colnames(XUrbM) %in% covariates)
    XUrbM = XUrbM[, keepIdxM, drop=FALSE]
    XRurM = XRurM[, keepIdxM, drop=FALSE]
    keepIdxD = which(colnames(XUrbD) %in% covariates)
    XUrbD = XUrbD[, keepIdxD, drop=FALSE]
    XRurD = XRurD[, keepIdxD, drop=FALSE]
  }
  covNames = colnames(XUrbM)

  load("savedOutput/global/admFinalMat.RData")
  bym2Args = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha,
                                        constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  nFree = ncol(bym2Args$Q) - 1

  lambdaTau     = getLambdaPCprec(u=pc.iidPrec$u,  alpha=pc.iidPrec$alpha)
  lambdaTauEps  = getLambdaPCprec(u=pc.expPrec$u,  alpha=pc.expPrec$alpha)
  # PC prior directly on sigma_comm: P(sigma_comm > u) = alpha
  lambdaSigmaComm = -log(pc.sigmaComm$alpha) / pc.sigmaComm$u

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
    lambdaSigmaComm=lambdaSigmaComm,
    commIntercept=as.integer(commIntercept),
    commSlope=as.integer(commSlope),
    uniformPhiPrior=as.integer(FALSE),
    options=0
  )

  initUrbP = sum(c(ysUrbMICS, ysUrbDHS))/sum(c(nsUrbMICS, nsUrbDHS))
  initRurP = sum(c(ysRurMICS, ysRurDHS))/sum(c(nsRurMICS, nsRurDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  nBeta = ncol(XUrbM)

  params = list(
    log_tau=0, logit_phi=0, log_tauEps=0, log_sigma_comm=0,
    alpha=initAlpha, alpha_M=initAlpha,
    beta=c(initBeta1,   rep(0, nBeta - 1)),
    beta_M=c(initBeta1, rep(0, nBeta - 1)),
    w_bym2Free=rep(0, nFree), u_bym2Free=rep(0, nFree)
  )
  # optional multi-start: override any of the above init values by name
  if(!is.null(initParams)) for(nm in names(initParams)) params[[nm]] = initParams[[nm]]

  feRand = c("alpha", "alpha_M", "beta", "beta_M")  # the fixed effects
  if(fixedEffectsOnly) {
    mapList = list(
      w_bym2Free = factor(rep(NA, nFree)),
      u_bym2Free = factor(rep(NA, nFree)),
      log_tau    = factor(NA),
      logit_phi  = factor(NA)
    )
    randList = feRand
  } else {
    mapList = list()
    randList = c(feRand, "w_bym2Free", "u_bym2Free")
  }
  # marginalizeFE=FALSE: treat alpha/beta(/_M) as FIXED outer params (drop from
  # random=) instead of Laplace-marginalizing them -> matches the sim-study
  # standard fits. Lets us test fixed-vs-marginalized cleanly on the same template.
  if(!marginalizeFE) randList = setdiff(randList, feRand)
  # sigma_comm is only identified if at least one commensurate prior uses it.
  # Fully non-commensurate (both flags 0): fix log_sigma_comm out of the optimization.
  if(commIntercept == 0L && commSlope == 0L) {
    mapList$log_sigma_comm = factor(NA)
  }
  # nesting/bug check: pin sigma_comm at a tiny value so the commensurate prior
  # forces alpha_M~alpha and beta_M~beta -> should reproduce the standard model.
  if(!is.null(fixSigmaComm)) {
    params$log_sigma_comm = log(fixSigmaComm)
    mapList$log_sigma_comm = factor(NA)
  }

  unloadDynlibs()
  dllName = "modMDM_BYM2_GH_comm"
  if(!any(file.exists(paste0("code/", dllName, c(".o", ".so", ".dll"))))) {
    cat("Compiling code/", dllName, ".cpp...\n", sep="")
    compile(paste0("code/", dllName, ".cpp"), framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib(paste0("code/", dllName)))

  obj = MakeADFun(data=data_gh, parameters=params,
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
       method=if(doMCMC) "MCMC-comm" else "GH-comm")
}
