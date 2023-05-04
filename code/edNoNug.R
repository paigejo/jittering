# script for women's secondary education in Nigeria application

# load datasets ----
out = load("savedOutput/global/ed.RData")
out = load("savedOutput/global/edMICS.RData")

# set parameters ----
KMICS=25

if(FALSE) {
  # do some precomputation ----
  
  # make integration points if necessary
  intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart = 1, loadSavedIntPoints=FALSE, 
                                            numPtsRur=KMICS, numPtsUrb=KMICS)
  intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  # modify the integration points to be in the correct format for TMB
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K * nStrat] x nVar
  stratUrb = XUrb$strat
  AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  stratIndexUrb = unlist(mapply(rep, 1:(nrow(AUrbMICS)*KMICS), each=rep(numPerStratUrb, KMICS)))
  obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  intPtIndexUrb = rep(1:KMICS, each=sum(numPerStratUrb))
  XUrb = XUrb[stratIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [K * nStrat] x nVar
  stratRur = XRur$strat
  ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  numPerStratRur = rowSums(ARurMICS)
  stratIndexRur = unlist(mapply(rep, 1:(nrow(ARurMICS)*KMICS), each=rep(numPerStratRur, KMICS)))
  obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  intPtIndexRur = rep(1:KMICS, each=sum(numPerStratRur))
  XRur = XRur[stratIndexRur,] # now XRur is [K * nObsRur] x nVar
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  stratIDs = match(edMICS$Stratum, admFinal$NAME_FINAL)
  edMICS = edMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICS[edMICS$urban,]$ys
  nsUrbMICS = edMICS[edMICS$urban,]$ns
  ysRurMICS = edMICS[!edMICS$urban,]$ys
  nsRurMICS = edMICS[!edMICS$urban,]$ns
  
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbMICS = t(AUrbMICS)
  ARurMICS = t(ARurMICS)
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbMICS) = "numeric"
  mode(ARurMICS) = "numeric"
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRur = XRur[,-(2:3)]
  intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  save(AUrbMICS, ARurMICS, AUrbDHS, ARurDHS, intPtsDHS, intPtsMICS, 
       ysUrbMICS, nsUrbMICS, ysRurMICS, nsRurMICS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
       file="savedOutput/global/edInputs.RData")
  
  if(FALSE) {
    # do some quick checks to make sure indices are correct:
    X_betaUrbanDHS=intPtsDHS$covsUrb
    X_betaUrbanMICS=intPtsMICS$XUrb
    X_betaRuralDHS=intPtsDHS$covsRur
    X_betaRuralMICS=intPtsMICS$XRur
    
    nDHSurb = length(ysUrbDHS)
    KDHSurb = 11
    nDHSrur = length(ysRurDHS)
    KDHSrur = 16
    nMICSurb = length(ysUrbMICS)
    KMICSurb = KMICS
    nMICSrur = length(ysRurMICS)
    KMICSrur = KMICS
    
    # indices for first integration point
    iFirstIntUrbDHS = 1:nDHSurb
    iFirstIntRurDHS = 1:nDHSrur
    iFirstIntUrbMICS = 1:nMICSurb
    iFirstIntRurMICS = 1:nMICSrur
    
    # indices for first observation
    iFirstObsUrbDHS = seq(1, nrow(X_betaUrbanDHS), by=nDHSurb)
    iFirstObsRurDHS = seq(1, nrow(X_betaRuralDHS), by=nDHSrur)
    iFirstObsUrbMICS = seq(1, nrow(X_betaUrbanMICS), by=nMICSurb)
    iFirstObsRurMICS = seq(1, nrow(X_betaRuralMICS), by=nMICSrur)
    
    X_betaUrbanDHS[iFirstObsUrbDHS,]
    X_betaRuralDHS[iFirstObsRurDHS,]
    X_betaUrbanMICS[iFirstObsUrbMICS,]
    X_betaRuralMICS[iFirstObsRurMICS,]
    
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), X_betaUrbanDHS[,2], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), X_betaRuralDHS[,2], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), X_betaUrbanDHS[,3], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), X_betaRuralDHS[,3], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), X_betaUrbanDHS[,4], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), X_betaRuralDHS[,4], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), X_betaUrbanDHS[,5], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), X_betaRuralDHS[,5], pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(intPtsDHS$xUrban), c(intPtsDHS$yUrban), c(intPtsDHS$wUrban), pch=19, cex=.2, legend.mar=2, ordering="increasing")
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), c(intPtsDHS$wRural), pch=19, cex=.2, legend.mar=2, ordering="increasing")
    
    plotWithColor(c(XUrb$east), c(XUrb$north), XUrb$access, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XRur$east), c(XRur$north), XRur$access, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XUrb$east), c(XUrb$north), XUrb$elev, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XRur$east), c(XRur$north), XRur$elev, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XUrb$east), c(XUrb$north), XUrb$distRiversLakes, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XRur$east), c(XRur$north), XRur$distRiversLakes, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XUrb$east), c(XUrb$north), XUrb$normPop, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XRur$east), c(XRur$north), XRur$normPop, pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XUrb$east), c(XUrb$north), c(wUrban), pch=19, cex=.2, legend.mar=2)
    plotWithColor(c(XRur$east), c(XRur$north), c(wRural), pch=19, cex=.2, legend.mar=2)
    
    plotWithColor(c(intPtsDHS$xRural), c(intPtsDHS$yRural), X_betaRuralDHS[,2], pch=19, cex=.2, legend.mar=2)
    
    iIthObsUrbDHS = seq(7, nrow(X_betaUrbanDHS), by=nDHSurb)
    X_betaUrbanDHS[iIthObsUrbDHS,]
    plotWithColor(c(intPtsDHS$xUrban[iIthObsUrbDHS[1],]), c(intPtsDHS$yUrban[iIthObsUrbDHS[1],]), 
                  X_betaUrbanDHS[iIthObsUrbDHS,2], pch=19, legend.mar=2)
  }
  
  # compile model ----
  dyn.unload( dynlib("code/modBYM2JitterFusion"))
  compile( "code/modBYM2JitterFusion.cpp")
}

# load in TMB function inputs
out = load("savedOutput/global/edInputs.RData")

# set priors ----
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)

out = load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                         constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3)

# Specify inputs for TMB ----

## specify random effects
rand_effs <- c('Epsilon_bym2')

# collect input data

data_full = list(
  y_iUrbanMICS=ysUrbMICS, # observed binomial experiment at point i (clust)
  y_iRuralMICS=ysRurMICS, # 
  n_iUrbanMICS=nsUrbMICS, # number binomial trials
  n_iRuralMICS=nsRurMICS, # 
  AprojUrbanMICS=AUrbMICS, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  AprojRuralMICS=ARurMICS, # 
  X_betaUrbanMICS=intPtsMICS$XUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralMICS=intPtsMICS$XRur, # 
  wUrbanMICS=intPtsMICS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
  wRuralMICS=intPtsMICS$wRural, # 
  
  y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
  y_iRuralDHS=ysRurDHS, # 
  n_iUrbanDHS=nsUrbDHS, # number binomial trials
  n_iRuralDHS=nsRurDHS, # 
  AprojUrbanDHS=AUrbDHS, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  AprojRuralDHS=ARurDHS, # 
  X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralDHS=intPtsDHS$covsRur, # 
  wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
  wRuralDHS=intPtsDHS$wRural, # 
  
  Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
  V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
  alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
  beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
  tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
  gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
  lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
  lambdaTau=lambdaTau, # determines PC prior for tau
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
  sapply(data_full, anyna)
  sapply(data_full, myDim)
  hist(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS, breaks=50)
  hist(data_full$y_iRuralDHS/data_full$n_iRuralDHS, breaks=50)
  hist(data_full$n_iUrbanDHS-data_full$y_iUrbanDHS, breaks=seq(-0.5, 17.5, by=1))
  hist(data_full$y_iRuralDHS, breaks=seq(-0.5, 16.5, by=1))
  mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1)
  mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0)
}

# initial parameters
initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

tmb_params <- list(alpha = initAlpha, # intercept
                   beta = c(initBeta1, rep(0, ncol(intPtsMICS$XUrb)-1)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)) # RE on mesh vertices
)

# make TMB fun and grad ----
dyn.load( dynlib("code/modBYM2JitterFusion"))
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterFusion')
objFull <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     hessian=TRUE,
                     DLL='modBYM2JitterFusion')

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

{
  tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
  testObj = obj
  optPar = testObj$par
  startTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = obj
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
    testObj$env$tracepar = TRUE
    print(paste0("optimizing for tol = ", thisTol, "."))
    opt1 <- optim(par=testObj$par, fn=funWrapper, gr=grWrapper,
                  method = c("BFGS"), hessian = FALSE)
    optPar = opt1$par
    if(!is.null(opt1$message)) {
      print(paste0("error for tol = ", thisTol, ". Message:"))
      print(opt1$message)
      next
    }
    else {
      print(paste0("completed optimization for tol = ", thisTol, ". Running one final time..."))
      funWrapper(optPar)
      break
    }
  }
  endTime = proc.time()[3]
  totalTime = endTime - startTime
  totalTime # 104.406  seconds
}

# opt0 <- nlminb(start       =    obj[['par']],
#                objective   =    obj[['fn']],
#                gradient    =    obj[['gr']],
#                lower = rep(-10, length(obj[['par']])),
#                upper = rep( 10, length(obj[['par']])),
#                control     =    list(trace=1))

# * TMB Posterior Sampling ----
## Get standard errors
SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                     bias.correct = TRUE,
                     bias.correct.control = list(sd = TRUE))
save(SD0, obj, objFull, file="savedOutput/ed/fitNoNug.RData")
load("savedOutput/ed/fitNoNug.RData")
## summary(SD0, 'report')
## summary(SD0, 'fixed')

gridPredsNoNug = predGrid(SD0, obj)
save(gridPredsNoNug, file="savedOutput/ed/gridPredsNoNug.RData")
out = load("savedOutput/ed/gridPredsNoNug.RData")

stratPredsNoNug = predArea(gridPredsNoNug, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin2PredsNoNug = predArea(gridPredsNoNug, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPredsNoNug, file="savedOutput/ed/stratPredsNoNug.RData")
save(admin2PredsNoNug, file="savedOutput/ed/admin2PredsNoNug.RData")
out = load("savedOutput/ed/stratPredsNoNug.RData")
out = load("savedOutput/ed/admin2PredsNoNug.RData")

plotPreds(SD0, obj, 
          gridPreds=gridPredsNoNug, arealPreds=stratPredsNoNug, 
          plotNameRoot="edFusionNoNug", plotNameRootAreal="Strat")
plotPreds(SD0, obj, 
          gridPreds=gridPredsNoNug, arealPreds=admin2PredsNoNug, 
          plotNameRoot="edFusionNoNug", plotNameRootAreal="Admin2")
