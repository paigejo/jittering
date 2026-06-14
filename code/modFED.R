# FE + nugget model for DHS-only data
# fitFED: GH quadrature for nuggets, BYM2 spatial structure (mapped out for FE-only)
# Uses modD_BYM2_GH_v2.cpp TMB template

fitFED = function(datDHS=ed, inputsMDM=NULL,
                  intPtsDHS=NULL,
                  KDHSurb=16, JInnerUrban=4,
                  KDHSrur=21, JInnerRural=4, JOuterRural=1,
                  adm2DHS=adm2Full,
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
  }

  # Get DHS integration points and responses
  if(!is.null(inputsMDM)) {
    intPtsDHS = inputsMDM$intPtsDHS
    ysUrbDHS = inputsMDM$ysUrbDHS
    ysRurDHS = inputsMDM$ysRurDHS
    nsUrbDHS = inputsMDM$nsUrbDHS
    nsRurDHS = inputsMDM$nsRurDHS
  } else {
    if(is.null(intPtsDHS)) {
      # Try loading saved canonical intPts
      label = paste0(KDHSurb, "_", KDHSrur)
      savedF = sprintf("savedOutput/global/intPtsDHS_%s.RData", label)
      if(file.exists(savedF)) {
        cat("Loading saved intPts:", savedF, "\n")
        load(savedF) # loads intPtsDHS
      } else {
        cat("Generating DHS intPts (", KDHSurb, ",", KDHSrur, ")...\n")
        intPtsDHS = makeAllIntegrationPointsDHS(
          cbind(datDHS$east, datDHS$north), datDHS$urban,
          areaNames=datDHS$subarea, popPrior=TRUE,
          numPointsUrban=KDHSurb, numPointsRural=KDHSrur,
          JInnerUrban=JInnerUrban, JInnerRural=JInnerRural,
          JOuterRural=JOuterRural, adminMap=adm2DHS, saveOutput=FALSE)
      }
    }
    # Remove intercept from DHS covariates if present
    if("int" %in% colnames(intPtsDHS$covsUrb))
      intPtsDHS$covsUrb = intPtsDHS$covsUrb[, colnames(intPtsDHS$covsUrb) != "int", drop=FALSE]
    if("int" %in% colnames(intPtsDHS$covsRur))
      intPtsDHS$covsRur = intPtsDHS$covsRur[, colnames(intPtsDHS$covsRur) != "int", drop=FALSE]
    if(!is.null(intPtsDHS$wUrban)) intPtsDHS$wUrban = as.matrix(intPtsDHS$wUrban)
    if(!is.null(intPtsDHS$wRural)) intPtsDHS$wRural = as.matrix(intPtsDHS$wRural)

    isU = datDHS$urban
    ysUrbDHS = datDHS$y[isU]; ysRurDHS = datDHS$y[!isU]
    nsUrbDHS = datDHS$n[isU]; nsRurDHS = datDHS$n[!isU]
  }

  # Covariate matrices
  XUrb_D = as.matrix(intPtsDHS$covsUrb)
  XRur_D = as.matrix(intPtsDHS$covsRur)

  # Subset covariates if requested
  allCovNames = colnames(XUrb_D)
  if(!is.null(covariates)) {
    keepIdx = which(allCovNames %in% covariates)
    XUrb_D = XUrb_D[, keepIdx, drop=FALSE]
    XRur_D = XRur_D[, keepIdx, drop=FALSE]
  }
  covNames = colnames(XUrb_D)

  # BYM2 scaffolding (needed by template even when mapped out)
  load("savedOutput/global/admFinalMat.RData")
  bym2Args = prepareBYM2argumentsForTMB(admFinalMat, u=pc.bym2Prec$u, alpha=pc.bym2Prec$alpha,
                                        constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  nAreas = ncol(bym2Args$Q)
  nFree = nAreas - 1

  # A-projection matrices (nObs x nAreas) for spatial effect projection
  # When BYM2 is mapped out (fixedEffectsOnly=TRUE), these multiply zero vectors
  nUrbDHS = length(ysUrbDHS)
  nRurDHS = length(ysRurDHS)
  if(fixedEffectsOnly) {
    AprojUrbanDHS = matrix(0, nUrbDHS, nAreas)
    AprojRuralDHS = matrix(0, nRurDHS, nAreas)
  } else {
    areasUrb = adm2ToStratumMICS(datDHS$subarea[datDHS$urban])
    areasRur = adm2ToStratumMICS(datDHS$subarea[!datDHS$urban])
    AprojUrbanDHS = t(makeApointToArea(areasUrb, admFinal$NAME_FINAL))
    AprojRuralDHS = t(makeApointToArea(areasRur, admFinal$NAME_FINAL))
    mode(AprojUrbanDHS) = "numeric"
    mode(AprojRuralDHS) = "numeric"
  }

  # Priors
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

  # Initial values
  initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
  initRurP = sum(ysRurDHS)/sum(nsRurDHS)
  initAlpha = log(initRurP/(1-initRurP))
  initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

  gh = fastGHQuad::gaussHermiteData(Qgh)
  nBeta = ncol(XUrb_D)

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

  data_gh = list(
    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    AprojUrbanDHS=AprojUrbanDHS, AprojRuralDHS=AprojRuralDHS,
    X_betaUrbanDHS=XUrb_D, X_betaRuralDHS=XRur_D,
    wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
    Q_bym2=bym2Args$Q,
    lchoose_urban=lchoose(nsUrbDHS, ysUrbDHS),
    lchoose_rural=lchoose(nsRurDHS, ysRurDHS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    tr=bym2Args$tr, gammaTildesm1=bym2Args$gammaTildesm1,
    lambdaPhi=bym2Args$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    uniformPhiPrior=as.integer(FALSE),
    options=0
  )

  unloadDynlibs()
  dllName = "modD_BYM2_GH_v2"
  if(!any(file.exists(paste0("code/", dllName, c(".o", ".so", ".dll"))))) {
    cat("Compiling code/modD_BYM2_GH_v2.cpp...\n")
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
       opt=opt, covNames=covNames,
       method=if(doMCMC) "MCMC" else "GH")
}
