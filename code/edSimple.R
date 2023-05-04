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
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [nStrat * K] x nVar
  stratUrb = XUrb$strat
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  XUrb = XUrb[stratIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  numPerStratRur = rowSums(ARurMICS)
  stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
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
  
  # compile model ----
  dyn.unload( dynlib("code/modBYM2JitterFusionSimple"))
  compile( "code/modBYM2JitterFusionSimple.cpp")
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
  X_betaUrbanMICS=matrix(intPtsMICS$XUrb[,1], ncol=1), # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralMICS=matrix(intPtsMICS$XRur[,1], ncol=1), # 
  wUrbanMICS=intPtsMICS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
  wRuralMICS=intPtsMICS$wRural, # 
  
  y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
  y_iRuralDHS=ysRurDHS, # 
  n_iUrbanDHS=nsUrbDHS, # number binomial trials
  n_iRuralDHS=nsRurDHS, # 
  AprojUrbanDHS=AUrbDHS, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  AprojRuralDHS=ARurDHS, # 
  X_betaUrbanDHS=matrix(intPtsDHS$covsUrb[,1], ncol=1), # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralDHS=matrix(intPtsDHS$covsRur[,1], ncol=1), # 
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
                   beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanMICS)-1)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)) # RE on mesh vertices
)

# make TMB fun and grad ----
dyn.load( dynlib("code/modBYM2JitterFusionSimple"))
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterFusionSimple')
objFull <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     hessian=TRUE,
                     DLL='modBYM2JitterFusionSimple')

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
  totalTime
}

# compare with svyglm
allYs = c(data_full$y_iUrbanMICS, data_full$y_iRuralMICS, 
          data_full$y_iUrbanDHS, data_full$y_iRuralDHS)
allNs = c(data_full$n_iUrbanMICS, data_full$n_iRuralMICS, 
          data_full$n_iUrbanDHS, data_full$n_iRuralDHS)
responseMat = cbind(successes=allYs, failures=allNs-allYs)
allUrbs = c(rep(1, length(data_full$y_iUrbanMICS)), rep(0, length(data_full$y_iRuralMICS)), 
            rep(1, length(data_full$y_iUrbanDHS)), rep(0, length(data_full$y_iRuralDHS)))
out = glm(responseMat ~ allUrbs, family="binomial")
summary(out)
optPar
outPar = out$coefficients
outAlpha = outPar[1]
outBeta = outPar[2]
ps = expit(outAlpha + outBeta * allUrbs)
sum(sapply(1:length(allNs), function(ind) {dbinom(allYs[ind], allNs[ind], ps[ind], log=TRUE)}))

lastPar = testObj$env$last.par
test = testObj$report(lastPar)
head(test$fe_iUrbanMICS)
head(test$testUrbanMICS)
tail(test$fe_iUrbanMICS)
tail(test$testUrbanMICS)
head(test$fe_iRuralMICS)
head(test$testRuralMICS)
tail(test$fe_iRuralMICS)
tail(test$testRuralMICS)

head(test$fe_iUrbanMICS)
head(test$latentFieldUrbMICS)
head(test$fe_iUrbanDHS)
head(test$latentFieldUrbDHS)
head(test$fe_iRuralMICS)
head(test$latentFieldRurMICS)
head(test$fe_iRuralDHS)
head(test$latentFieldRurDHS)

tail(test$fe_iUrbanMICS)
tail(test$latentFieldUrbMICS)
tail(test$fe_iUrbanDHS)
tail(test$latentFieldUrbDHS)
tail(test$fe_iRuralMICS)
tail(test$latentFieldRurMICS)
tail(test$fe_iRuralDHS)
tail(test$latentFieldRurDHS)

lastAlphaBeta = lastPar[1:2]
psUrbMICS=rep(sum(lastAlphaBeta), length(data_full$n_iUrbanMICS))
test$liksUrbMICS[1:5,1:5]
psUrbMICS = rep(expit(sum(lastAlphaBeta)), length(data_full$n_iUrbanMICS))
head(sapply(1:length(psUrbMICS), function(i) {dbinom(data_full$y_iUrbanMICS[i], data_full$n_iUrbanMICS[i], psUrbMICS[i])}))
test$liksUrbMICS[1:6,1]
sapply(1:6, function(i) {dbinom(data_full$y_iUrbanMICS[i], data_full$n_iUrbanMICS[i], expit(test$latentFieldUrbMICS[i]))})
range(abs(sapply(1:length(psUrbMICS), function(i) {dbinom(data_full$y_iUrbanMICS[i], data_full$n_iUrbanMICS[i], psUrbMICS[i])}) -
            test$liksUrbMICS[,1]))

psRurMICS = rep(expit(lastAlphaBeta[1]), length(data_full$n_iRuralMICS))
head(sapply(1:length(psRurMICS), function(i) {dbinom(data_full$y_iRuralMICS[i], data_full$n_iRuralMICS[i], psRurMICS[i])}))
test$liksRurMICS[1:6,1]
sapply(1:6, function(i) {dbinom(data_full$y_iRuralMICS[i], data_full$n_iRuralMICS[i], expit(test$latentFieldRurMICS[i]))})
range(abs(sapply(1:length(psRurMICS), function(i) {dbinom(data_full$y_iRuralMICS[i], data_full$n_iRuralMICS[i], psRurMICS[i])}) -
            test$liksRurMICS[,1]))

psUrbDHS=rep(sum(lastAlphaBeta), length(data_full$n_iUrbanDHS))
test$liksUrbDHS[1:5,1:5]
psUrbDHS = rep(expit(sum(lastAlphaBeta)), length(data_full$n_iUrbanDHS))
head(sapply(1:length(psUrbDHS), function(i) {dbinom(data_full$y_iUrbanDHS[i], data_full$n_iUrbanDHS[i], psUrbDHS[i])}))
test$liksUrbDHS[1:6,1]
sapply(1:6, function(i) {dbinom(data_full$y_iUrbanDHS[i], data_full$n_iUrbanDHS[i], expit(test$latentFieldUrbDHS[i]))})
range(abs(sapply(1:length(psUrbDHS), function(i) {dbinom(data_full$y_iUrbanDHS[i], data_full$n_iUrbanDHS[i], psUrbDHS[i])}) -
            test$liksUrbDHS[,1]))

psRurDHS = rep(expit(lastAlphaBeta[1]), length(data_full$n_iRuralDHS))
head(sapply(1:length(psRurDHS), function(i) {dbinom(data_full$y_iRuralDHS[i], data_full$n_iRuralDHS[i], psRurDHS[i])}))
test$liksRurDHS[1:6,1]
sapply(1:6, function(i) {dbinom(data_full$y_iRuralDHS[i], data_full$n_iRuralDHS[i], expit(test$latentFieldRurDHS[i]))})
mean(abs(sapply(1:length(psRurDHS), function(i) {dbinom(data_full$y_iRuralDHS[i], data_full$n_iRuralDHS[i], psRurDHS[i])}) -
      test$liksRurDHS[,1]))

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
## summary(SD0, 'report')
## summary(SD0, 'fixed')

## take samples from fitted model
mu <- c(SD0$par.fixed,SD0$par.random)
test = hessian(funWrapper, SD0$par.fixed)
eig = eigen(test)
for(i in 1:length(eig$values)) {
  thisVal = eig$values[i]
  thisVec = eig$vectors[,i]
  barplot(eig$vectors[,i], names.arg=names(SD0$par.fixed), 
          main=paste0("Eigenvector ", i, " with value ", thisVal))
}


sum(SD0$par.random)
sumTo0 = rnorm(length(SD0$par.random))
sumTo0 = sumTo0 - mean(sumTo0)

objFull$fn(c(SD0$par.fixed, SD0$par.random))
test1 = objFull$report()

objFull$fn(c(SD0$par.fixed, SD0$par.random + 1))
objFull$fn(c(SD0$par.fixed, SD0$par.random + sumTo0))
objFull$fn(c(SD0$par.fixed + c(-1, rep(0, 7)), SD0$par.random + 1))
objFull$fn(c(SD0$par.fixed, SD0$par.random)) - objFull$fn(c(SD0$par.fixed, SD0$par.random + sumTo0))
objFull$fn(c(SD0$par.fixed + c(1, rep(0, 7)), SD0$par.random + 0))
objFull$fn(c(SD0$par.fixed + c(-1, 10, rep(0, 6)), SD0$par.random + 0))
objFull$fn(c(SD0$par.fixed + c(-1, 10, rep(0, 6)), SD0$par.random + 0))
objFull$fn(c(SD0$par.fixed + c(initAlpha-SD0$par.fixed[1], initBeta1-SD0$par.fixed[2], 
                               rep(0, 6)), SD0$par.random + 0))
test2 = objFull$report()

library(weights)
test1$quad
test2$quad
test1$transformedEpsilon
test2$transformedEpsilon
head(test1$latentFieldRurDHS)
head(test2$latentFieldRurDHS)
head(data_full$y_iRuralDHS)
head(data_full$n_iRuralDHS)
latMat1 = matrix(test1$latentFieldRurDHS, nrow=nrow(data_full$wRuralDHS))
latMat2 = matrix(test2$latentFieldRurDHS, nrow=nrow(data_full$wRuralDHS))
latMat3 = matrix(test1$latentFieldRurMICS, nrow=nrow(data_full$wRuralMICS))
latMat4 = matrix(test2$latentFieldRurMICS, nrow=nrow(data_full$wRuralMICS))
hist(expit(latMat1))
hist(expit(latMat1[data_full$wUrbanDHS!=0]))
wtd.hist(expit(c(latMat1)), weight=data_full$wRuralDHS)
hist(expit(latMat2))
wtd.hist(expit(c(latMat2)), weight=data_full$wRuralDHS)
hist(expit(latMat3))
wtd.hist(expit(c(latMat3)), weight=data_full$wRuralMICS)
hist(expit(latMat4))
wtd.hist(expit(c(latMat4)), weight=data_full$wRuralMICS)
hist(data_full$y_iRuralDHS/data_full$n_iRuralDHS)
hist(data_full$y_iRuralMICS/data_full$n_iRuralMICS)
wtd.hist(c(expit(latMat1), expit(latMat3)), 
         weight=c(data_full$wRuralDHS, data_full$wRuralMICS))
wtd.hist(c(expit(latMat2), expit(latMat4)), 
         weight=c(data_full$wRuralDHS, data_full$wRuralMICS))
hist(c(data_full$y_iRuralDHS/data_full$n_iRuralDHS, data_full$y_iRuralMICS/data_full$n_iRuralMICS), breaks=50)
mean(rowSums(expit(latMat1)*data_full$wRuralDHS))
mean(rowSums(expit(latMat2)*data_full$wRuralDHS))
mean(rowSums(expit(latMat3)*data_full$wRuralMICS))
mean(rowSums(expit(latMat4)*data_full$wRuralMICS))
sum(test1$liksRurDHS * data_full$wRuralDHS)
sum(test2$liksRurDHS * data_full$wRuralDHS)
sum(mapply(i=rep(1:nrow(data_full$wRuralDHS), ncol(data_full$wRuralDHS)), 
           j=rep(1:ncol(data_full$wRuralDHS), each=nrow(data_full$wRuralDHS)), 
           FUN=function(i, j) {
             data_full$wRuralDHS[i,j]*dbinom(data_full$y_iRuralDHS[i], data_full$n_iRuralDHS[i], prob=expit(latMat1[i,j]))
           }))
sum(mapply(i=rep(1:nrow(data_full$wRuralDHS), ncol(data_full$wRuralDHS)), 
           j=rep(1:ncol(data_full$wRuralDHS), each=nrow(data_full$wRuralDHS)), 
           FUN=function(i, j) {
             data_full$wRuralDHS[i,j]*dbinom(data_full$y_iRuralDHS[i], data_full$n_iRuralDHS[i], prob=expit(latMat2[i,j]))
           }))

head(test1$latentFieldUrbDHS)
head(test2$latentFieldUrbDHS)
head(data_full$y_iUrbanDHS)
head(data_full$n_iUrbanDHS)
latMat1 = matrix(test1$latentFieldUrbDHS, nrow=nrow(data_full$wUrbanDHS))
latMat2 = matrix(test2$latentFieldUrbDHS, nrow=nrow(data_full$wUrbanDHS))
latMat3 = matrix(test1$latentFieldUrbMICS, nrow=nrow(data_full$wUrbanMICS))
latMat4 = matrix(test2$latentFieldUrbMICS, nrow=nrow(data_full$wUrbanMICS))
hist(expit(latMat1))
hist(expit(latMat1[data_full$wUrbanDHS!=0]))
wtd.hist(expit(c(latMat1)), weight=data_full$wUrbanDHS)
hist(expit(latMat2))
wtd.hist(expit(c(latMat2)), weight=data_full$wUrbanDHS)
hist(expit(latMat3))
wtd.hist(expit(c(latMat3)), weight=data_full$wUrbanMICS)
hist(expit(latMat4))
wtd.hist(expit(c(latMat4)), weight=data_full$wUrbanMICS)
hist(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS)
hist(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS)
wtd.hist(c(expit(latMat1), expit(latMat3)), 
         weight=c(data_full$wUrbanDHS, data_full$wUrbanMICS))
wtd.hist(c(expit(latMat2), expit(latMat4)), 
         weight=c(data_full$wUrbanDHS, data_full$wUrbanMICS))
hist(c(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS, data_full$y_iUrbanMICS/data_full$n_iUrbanMICS), breaks=50)
mean(rowSums(expit(latMat1)*data_full$wUrbanDHS))
mean(rowSums(expit(latMat2)*data_full$wUrbanDHS))
mean(rowSums(expit(latMat3)*data_full$wUrbanMICS))
mean(rowSums(expit(latMat4)*data_full$wUrbanMICS))
sum(test1$liksUrbDHS * data_full$wUrbanDHS)
sum(test2$liksUrbDHS * data_full$wUrbanDHS)
sum(test1$liksUrbMICS * data_full$wUrbanMICS)
sum(test2$liksUrbMICS * data_full$wUrbanMICS)
sum(mapply(i=rep(1:nrow(data_full$wUrbanDHS), ncol(data_full$wUrbanDHS)), 
           j=rep(1:ncol(data_full$wUrbanDHS), each=nrow(data_full$wUrbanDHS)), 
           FUN=function(i, j) {
             data_full$wUrbanDHS[i,j]*dbinom(data_full$y_iUrbanDHS[i], data_full$n_iUrbanDHS[i], prob=expit(latMat1[i,j]))
           }))
sum(mapply(i=rep(1:nrow(data_full$wUrbanDHS), ncol(data_full$wUrbanDHS)), 
           j=rep(1:ncol(data_full$wUrbanDHS), each=nrow(data_full$wUrbanDHS)), 
           FUN=function(i, j) {
             data_full$wUrbanDHS[i,j]*dbinom(data_full$y_iUrbanDHS[i], data_full$n_iUrbanDHS[i], prob=expit(latMat2[i,j]))
           }))

## simulate draws
rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- chol_prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

L <- Cholesky(SD0[['jointPrecision']], super = T)
t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 500)

## summarize the draws
parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_s',]
alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)

# project from mesh to raster, add intercept
data(kenyaPopulationData)
predCoords = cbind(popMatKenya$east, popMatKenya$north)
A.pred = inla.spde.make.A(mesh = mesh.s, loc = predCoords)
pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)
pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')

# convert to probability scale
pred_tmb = expit(pred_tmb)

## find the median and sd across draws, as well as 90% intervals
summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                  sd     = (apply(pred_tmb, 1, sd)),
                  lower = (apply(pred_tmb, 1, quantile, .05)),
                  upper = (apply(pred_tmb, 1, quantile, .95)))

# * Plot TMB results ----
## make summary rasters
# kenyaExtent = extent(cbind(popMatKenya$east, popMatKenya$north))
# pop.rast <- raster(nrows=length(unique(popMatKenya$north)), ncols=length(length(unique(popMatKenya$north))),
#                   xmn=min(popMatKenya$east), xmx=max(length(unique(popMatKenya$east))), 
#                   ymn=min(popMatKenya$north), ymx=max(popMatKenya$north),
#                   vals=popMatKenya$pop, ext = kenyaExtent)
pop.rast = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=popMatKenya$pop), 
                         res=c(5, 5))
ras_med_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 1]), 
                            res=c(5, 5))
ras_sdv_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 2]), 
                            res=c(5, 5))
ras_lower_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 3]), 
                              res=c(5, 5))
ras_upper_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 4]), 
                              res=c(5, 5))

## plot truth, pixels falling within/without the 90% interval,
##  post. median, and post sd

# set the range for the truth and median
rast.zrange <- range(c(dat$y/dat$n, values(ras_med_tmb)), na.rm = T)

# plot tmb
par(mfrow = c(2, 2))
quilt.plot(dat$east, dat$north, dat$y/dat$n, main = 'Observations', col = (viridis(100)), 
           nx=30, ny=30)
# plot(ras_inInt_tmb, main = 'Pixels where 90% CIs did not cover Truth')
# points(dat$east, dat$north, pch=".")
plot(ras_med_tmb, main = 'TMB Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main='TMB Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

## plot TMB meds and stdevs
med.zrange <- range(values(ras_med_tmb), na.rm = T)
sdv.zrange <- range(values(ras_sdv_tmb), na.rm = T)

par(mfrow = c(1, 2))
plot(ras_med_tmb, main = 'TMB Posterior Median',
     zlim = med.zrange, col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main = 'TMB Posterior Standard Deviation',
     zlim = sdv.zrange)
points(dat$east, dat$north, pch=".")
