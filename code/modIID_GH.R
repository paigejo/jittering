# IID spatial + GH nugget model wrappers (MICS-only, DHS-only, DHS+MICS fusion)
# Analogous to fitMM, fitMD, fitMDM in modM_MSep.R / modM_DSep.R / modM_DMSep.R
# but replacing BYM2 spatial effects with IID normal effects

fitMM_iid = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
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
                     maxit=1000, Qgh=10, getSDs=TRUE, verbose=FALSE, ...) {

  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }

  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]
    toN = nameTab[i,2]
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

  out = load("savedOutput/global/admFinalMat.RData")
  nAreas = nrow(admFinalMat)
  nFree = nAreas - 1
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

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
    lchoose_urban=lchoose(nsUrbMICS, ysUrbMICS),
    lchoose_rural=lchoose(nsRurMICS, ysRurMICS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    options=0
  )

  nBeta = ncol(data_gh$X_betaUrbanMICS)

  initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
  initRurP = sum(ysRurMICS)/sum(nsRurMICS)
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha

  params_init = list(
    log_tau=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta-1)),
    v_iidFree=rep(0, nFree)
  )

  unloadDynlibs()
  if(!any(file.exists(paste0("code/modM_IID_GH", c(".o", ".so", ".dll"))))) {
    compile("code/modM_IID_GH.cpp", framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib("code/modM_IID_GH"))

  map_fe = list(
    v_iidFree=factor(rep(NA, nFree)),
    log_tau=factor(NA)
  )
  obj_fe = MakeADFun(data=data_gh, parameters=params_init, map=map_fe,
                     DLL="modM_IID_GH", silent=!verbose)
  opt_fe = optim(obj_fe$par, obj_fe$fn, obj_fe$gr,
                 method="BFGS", control=list(maxit=maxit, reltol=1e-6))

  params_full = params_init
  params_full$alpha = as.numeric(opt_fe$par["alpha"])
  params_full$beta = as.numeric(opt_fe$par[grep("^beta", names(opt_fe$par))])
  params_full$log_tauEps = as.numeric(opt_fe$par["log_tauEps"])

  obj = MakeADFun(data=data_gh, parameters=params_full,
                  DLL="modM_IID_GH", silent=!verbose)
  startTime = proc.time()[3]
  opt = optim(obj$par, obj$fn, obj$gr,
              method="BFGS", control=list(maxit=maxit, reltol=1e-6))
  totalTime = proc.time()[3] - startTime
  obj$env$last.par = opt$par

  SD0 = NULL
  sdTime = 0
  if(getSDs) {
    sdTime = system.time({
      SD0 <- try(TMB::sdreport(obj, getJointPrecision=TRUE,
                               bias.correct=TRUE,
                               bias.correct.control=list(sd=TRUE)), silent=TRUE)
    })[3]
    if(inherits(SD0, "try-error")) SD0 = NULL
  }

  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime,
       initOpt=opt_fe, opt=opt, inputsMDM=inputsMDM)
}


fitMD_iid = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
                     intPtsMICS=NULL, intPtsDHS=NULL,
                     KMICS=100,
                     KDHSurb=16, JInnerUrban=4,
                     KDHSrur=21, JInnerRural=4, JOuterRural=1,
                     admMICS=admFinal, adm2DHS=adm2Full,
                     alpha_pri=c(0, 100^2),
                     beta_pri=c(0, sqrt(1000)),
                     pc.iidPrec=list(u=1, alpha=0.5),
                     pc.expPrec=list(u=1, alpha=0.5),
                     maxit=1000, Qgh=10, getSDs=TRUE, verbose=FALSE, ...) {

  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }

  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]
    toN = nameTab[i,2]
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

  out = load("savedOutput/global/admFinalMat.RData")
  nAreas = nrow(admFinalMat)
  nFree = nAreas - 1
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

  gh = fastGHQuad::gaussHermiteData(Qgh)
  XUrbDHS = as.matrix(intPtsDHS$covsUrb)
  XRurDHS = as.matrix(intPtsDHS$covsRur)
  if("int" %in% colnames(XUrbDHS)) XUrbDHS = XUrbDHS[, colnames(XUrbDHS) != "int", drop=FALSE]
  if("int" %in% colnames(XRurDHS)) XRurDHS = XRurDHS[, colnames(XRurDHS) != "int", drop=FALSE]

  data_gh = list(
    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    AprojUrbanDHS=t(makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)),
    AprojRuralDHS=t(makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)),
    X_betaUrbanDHS=XUrbDHS,
    X_betaRuralDHS=XRurDHS,
    wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
    lchoose_urban=lchoose(nsUrbDHS, ysUrbDHS),
    lchoose_rural=lchoose(nsRurDHS, ysRurDHS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    options=0
  )
  mode(data_gh$AprojUrbanDHS) = "numeric"
  mode(data_gh$AprojRuralDHS) = "numeric"

  nBeta = ncol(data_gh$X_betaUrbanDHS)

  initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
  initRurP = sum(ysRurDHS)/sum(nsRurDHS)
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha

  params_init = list(
    log_tau=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta-1)),
    v_iidFree=rep(0, nFree)
  )

  unloadDynlibs()
  if(!any(file.exists(paste0("code/modD_IID_GH", c(".o", ".so", ".dll"))))) {
    compile("code/modD_IID_GH.cpp", framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib("code/modD_IID_GH"))

  map_fe = list(
    v_iidFree=factor(rep(NA, nFree)),
    log_tau=factor(NA)
  )
  obj_fe = MakeADFun(data=data_gh, parameters=params_init, map=map_fe,
                     DLL="modD_IID_GH", silent=!verbose)
  opt_fe = optim(obj_fe$par, obj_fe$fn, obj_fe$gr,
                 method="BFGS", control=list(maxit=maxit, reltol=1e-6))

  params_full = params_init
  params_full$alpha = as.numeric(opt_fe$par["alpha"])
  params_full$beta = as.numeric(opt_fe$par[grep("^beta", names(opt_fe$par))])
  params_full$log_tauEps = as.numeric(opt_fe$par["log_tauEps"])

  obj = MakeADFun(data=data_gh, parameters=params_full,
                  DLL="modD_IID_GH", silent=!verbose)

  startTime = proc.time()[3]
  opt = optim(obj$par, obj$fn, obj$gr,
              method="BFGS", control=list(maxit=maxit, reltol=1e-6))
  totalTime = proc.time()[3] - startTime
  obj$env$last.par = opt$par

  SD0 = NULL
  sdTime = 0
  if(getSDs) {
    sdTime = system.time({
      SD0 <- try(TMB::sdreport(obj, getJointPrecision=TRUE,
                               bias.correct=TRUE,
                               bias.correct.control=list(sd=TRUE)), silent=TRUE)
    })[3]
    if(inherits(SD0, "try-error")) SD0 = NULL
  }

  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime,
       initOpt=opt_fe, opt=opt, inputsMDM=inputsMDM)
}


fitMDM_iid = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
                      intPtsMICS=NULL, intPtsDHS=NULL,
                      KMICS=100,
                      KDHSurb=16, JInnerUrban=4,
                      KDHSrur=21, JInnerRural=4, JOuterRural=1,
                      admMICS=admFinal, adm2DHS=adm2Full,
                      alpha_pri=c(0, 100^2),
                      beta_pri=c(0, sqrt(1000)),
                      pc.iidPrec=list(u=1, alpha=0.5),
                      pc.expPrec=list(u=1, alpha=0.5),
                      maxit=1000, Qgh=10, getSDs=TRUE, verbose=FALSE, ...) {

  if(!("Stratum" %in% names(datMICS))) {
    datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
  }
  nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
  for(i in 1:nrow(nameTab)) {
    fromN = nameTab[i,1]
    toN = nameTab[i,2]
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

  out = load("savedOutput/global/admFinalMat.RData")
  nAreas = nrow(admFinalMat)
  nFree = nAreas - 1
  lambdaTau = getLambdaPCprec(u=pc.iidPrec$u, alpha=pc.iidPrec$alpha)
  lambdaTauEps = getLambdaPCprec(u=pc.expPrec$u, alpha=pc.expPrec$alpha)

  gh = fastGHQuad::gaussHermiteData(Qgh)
  data_gh = list(
    y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
    n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
    areaidxlocUrbanMICS=as.integer(areaidxlocUrbanMICS),
    areaidxlocRuralMICS=as.integer(areaidxlocRuralMICS),
    X_betaUrbanMICS=intPtsMICS$XUrb,
    X_betaRuralMICS=intPtsMICS$XRur,
    wUrbanMICS=intPtsMICS$wUrban,
    wRuralMICS=intPtsMICS$wRural,

    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    areaidxlocUrbanDHS=as.integer(areaidxlocUrbanDHS),
    areaidxlocRuralDHS=as.integer(areaidxlocRuralDHS),
    X_betaUrbanDHS=intPtsDHS$covsUrb,
    X_betaRuralDHS=intPtsDHS$covsRur,
    wUrbanDHS=intPtsDHS$wUrban,
    wRuralDHS=intPtsDHS$wRural,

    lchoose_urban_mics=lchoose(nsUrbMICS, ysUrbMICS),
    lchoose_rural_mics=lchoose(nsRurMICS, ysRurMICS),
    lchoose_urban_dhs=lchoose(nsUrbDHS, ysUrbDHS),
    lchoose_rural_dhs=lchoose(nsRurDHS, ysRurDHS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    options=0
  )

  nBeta = ncol(data_gh$X_betaUrbanMICS)

  initUrbP = sum(c(ysUrbMICS, ysUrbDHS))/sum(c(nsUrbMICS, nsUrbDHS))
  initRurP = sum(c(ysRurMICS, ysRurDHS))/sum(c(nsRurMICS, nsRurDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha

  params_init = list(
    log_tau=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta-1)),
    v_iidFree=rep(0, nFree)
  )

  unloadDynlibs()
  if(!any(file.exists(paste0("code/modMDM_IID_GH", c(".o", ".so", ".dll"))))) {
    compile("code/modMDM_IID_GH.cpp", framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib("code/modMDM_IID_GH"))

  map_fe = list(
    v_iidFree=factor(rep(NA, nFree)),
    log_tau=factor(NA)
  )
  obj_fe = MakeADFun(data=data_gh, parameters=params_init, map=map_fe,
                     DLL="modMDM_IID_GH", silent=!verbose)
  opt_fe = optim(obj_fe$par, obj_fe$fn, obj_fe$gr,
                 method="BFGS", control=list(maxit=maxit, reltol=1e-6))

  params_full = params_init
  params_full$alpha = as.numeric(opt_fe$par["alpha"])
  params_full$beta = as.numeric(opt_fe$par[grep("^beta", names(opt_fe$par))])
  params_full$log_tauEps = as.numeric(opt_fe$par["log_tauEps"])

  obj = MakeADFun(data=data_gh, parameters=params_full,
                  DLL="modMDM_IID_GH", silent=!verbose)
  startTime = proc.time()[3]
  opt = optim(obj$par, obj$fn, obj$gr,
              method="BFGS", control=list(maxit=maxit, reltol=1e-6))
  totalTime = proc.time()[3] - startTime
  obj$env$last.par = opt$par

  SD0 = NULL
  sdTime = 0
  if(getSDs) {
    sdTime = system.time({
      SD0 <- try(TMB::sdreport(obj, getJointPrecision=TRUE,
                               bias.correct=TRUE,
                               bias.correct.control=list(sd=TRUE)), silent=TRUE)
    })[3]
    if(inherits(SD0, "try-error")) SD0 = NULL
  }

  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime,
       initOpt=opt_fe, opt=opt, inputsMDM=inputsMDM)
}
