# FE + nugget models
# fitFEM / fitMFEM: GH quadrature for nuggets — now in modFEM.R
# fitMFEM_laplace: Laplace approximation for nuggets (as in modM_MIIDonly.cpp)

source("code/modFEM.R")


fitMFEM_laplace = function(datDHS=ed, datMICS=edMICS, inputsMDM=NULL,
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
                           maxit=1000, getSDs=TRUE, verbose=FALSE) {

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

  data_full = list(
    y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
    n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
    areaidxlocUrbanMICS=areaidxlocUrbanMICS,
    areaidxlocRuralMICS=areaidxlocRuralMICS,
    X_betaUrbanMICS=intPtsMICS$XUrb,
    X_betaRuralMICS=intPtsMICS$XRur,
    wUrbanMICS=intPtsMICS$wUrban,
    wRuralMICS=intPtsMICS$wRural,
    nAreas=nAreas,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
  )

  initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
  initRurP = sum(ysRurMICS)/sum(nsRurMICS)
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha

  tmb_params = list(
    log_tau=0, log_tauEps=0,
    alpha=initAlpha,
    beta=c(initBeta1, rep(0, nBeta-1)),
    u_spatial=rep(0, nAreas),
    nuggetUrbMICS=rep(0, nUrb),
    nuggetRurMICS=rep(0, nRur)
  )

  if(fixedEffectsOnly) {
    randList = c("nuggetUrbMICS", "nuggetRurMICS")
    mapList = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))
  } else {
    randList = c("alpha", "beta", "u_spatial", "nuggetUrbMICS", "nuggetRurMICS")
    mapList = list()
  }

  unloadDynlibs()
  dllName = "modM_MIIDonly"
  if(!any(file.exists(paste0("code/", dllName, c(".o", ".so", ".dll"))))) {
    cat("Compiling code/modM_MIIDonly.cpp...\n")
    compile(paste0("code/", dllName, ".cpp"), framework="TMBad", safebounds=FALSE)
  }
  dyn.load(dynlib(paste0("code/", dllName)))
  TMB::config(tmbad.sparse_hessian_compress = 1)

  obj = MakeADFun(data=data_full, parameters=tmb_params,
                  random=randList, map=mapList,
                  hessian=TRUE,
                  DLL=dllName, silent=!verbose)

  lower = rep(-10, length(obj$par))
  upper = rep( 10, length(obj$par))

  startTime = proc.time()[3]
  opt = optim(par=obj$par, fn=obj$fn, gr=obj$gr,
              method="BFGS", hessian=FALSE,
              control=list(maxit=maxit, reltol=1e-6))
  totalTime = proc.time()[3] - startTime

  if(opt$convergence != 0) {
    cat("Warning: optim did not converge:", opt$message, "\n")
  }

  SD0 = NULL; sdTime = 0
  if(getSDs) {
    sdTime = system.time({
      SD0 <- try(TMB::sdreport(obj, getJointPrecision=TRUE,
                               bias.correct=TRUE,
                               bias.correct.control=list(sd=TRUE)), silent=TRUE)
    })[3]
    if(inherits(SD0, "try-error")) SD0 = NULL
  }

  list(TMBobj=obj, TMBsd=SD0, totalTime=totalTime, sdTime=sdTime,
       opt=opt, inputsMDM=inputsMDM,
       covNames=colnames(intPtsMICS$XUrb), method="Laplace")
}
