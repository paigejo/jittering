# script for women's secondary education in Nigeria application

# load datasets ----
out = load("savedOutput/global/ed.RData")

if(FALSE) {
  # do some precomputation ----
  
  # make integration points if necessary
  intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  load("savedOutput/global/intPtsDHS.RData")
  
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  # extract cluster information (in the correct order)
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  save(AUrbDHS, ARurDHS, intPtsDHS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
       file="savedOutput/global/edInputsDHS.RData")
  
  # compile model ----
  dyn.unload( dynlib("code/modBYM2JitterDHS"))
  compile( "code/modBYM2JitterDHS.cpp")
  # compile("code/modBYM2JitterDHS.cpp","-O0 -g") # for using gbdsource
}

# load in TMB function inputs
out = load("savedOutput/global/edInputsDHS.RData")

# set priors ----
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)

out = load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                         constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
lambdaTau = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for bym2 precision
lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision

# collect input data ----

data_full = list(
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
  lambdaTauEps=lambdaTauEps, # determines PC prior for tauEps, the nugget precision
  options=0 # 1 for adreport of log tau and logit phi
)

# Specify starting values for TMB params ----
# initial parameters
initUrbP = sum(c(data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanDHS))
initRurP = sum(c(data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

tmb_params <- list(alpha = initAlpha, # intercept
                   beta = c(initBeta1, rep(0, ncol(intPtsDHS$covsUrb)-1)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                   nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                   nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
)

## specify random effects
rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS')

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

# * Run TMB ----
dyn.load( dynlib("code/modBYM2JitterDHS"))
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterDHS')
objFull <- MakeADFun(data=data_full,
                     parameters=tmb_params,
                     hessian=TRUE,
                     DLL='modBYM2JitterDHS')

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

{
  # tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
  tolSeq = 1e-06
  testObj = obj
  optPar = testObj$par
  startTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = obj
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
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
      
      ## Get standard errors
      print("getting standard errors...")
      sdTime = system.time({
        SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                             bias.correct = TRUE,
                             bias.correct.control = list(sd = TRUE))
        
        if(!SD0$pdHess) {
          # try recalculating for fixed parameters numerically
          warning("Using alternative method for hessian calculation...")
          Hess = numDeriv::hessian( func=testObj$fn, x=optPar )
          SD0 <- sdreport( testObj, hessian.fixed=Hess,
                           getJointPrecision=TRUE,
                           bias.correct = TRUE,
                           bias.correct.control = list(sd = TRUE) )
          
          # system.time(HessTest <- numDeriv::jacobian(func=objFull$gr, x=testObj$env$last.par.best))
          # 15.82288 minutes, not PD
          # system.time(HessTest2 <- numDeriv::hessian(func=objFull$fn, x=testObj$env$last.par.best))
        }
      }
      )[3]
      # SD0
      print(sdTime/60)
      # 3.9447 minutes for intern=FALSE
      
      if(SD0$pdHess) {
        print("Optimization and PD hess calculation done!")
        break
      }
      else {
        print("Hessan not PD. Rerunning optimization with stricter tol...")
      }
    }
    
    
    
  }
  endTime = proc.time()[3]
  sdTime/60
  totalTime = endTime - startTime
  print(paste0("optimization took ", totalTime/60, " minutes"))
  # optimization took  minutes (for intern=FALSE)
}

if(FALSE) {
  badPar = c(SD0$par.fixed, SD0$par.random)
  
  tempGrad = obj$env$f(badPar,order=1)
  tempHess = obj$env$spHess(badPar,random=TRUE)
  range(tempHess)
  
  tempEig = eigen(tempHess)
  range(tempEig$values)
}

if(!SD0$pdHess) {
  thisTMBpar = tmb_params
  thisTMBpar$log_tauEps = optPar[names(optPar) == "log_tauEps"]
  
  map = as.list(factor(c(alpha=1, beta=2:6, log_tau=7, logit_phi=8, log_tauEps=NA)))
  temp = unlist(map[grepl("beta", names(map))])
  map[grepl("beta", names(map))] = NULL
  map$beta = temp
  map=list(factor(c(log_tauEps=NA)))
  names(map) = "log_tauEps"
  objFixed <- MakeADFun(data=data_full,
                        parameters=tmb_params,
                        random=rand_effs,
                        map=map, 
                        hessian=TRUE,
                        DLL='modBYM2JitterDHS')
  testObj = objFixed
  thisOptPar = optPar[-which(names(optPar) == "log_tauEps")]
  lower = lower[-which(names(optPar) == "log_tauEps")]
  upper = upper[-which(names(optPar) == "log_tauEps")]
  
  newStartTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = objFixed
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
    testObj$env$tracepar = TRUE
    print(paste0("optimizing for tol = ", thisTol, "."))
    opt1 <- optim(par=thisOptPar, fn=funWrapper, gr=grWrapper,
                  method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
    # opt1 <- optim(par=optPar, fn=funWrapper, 
    #               hessian = FALSE, control=list(reltol=thisTol))
    thisOptPar = opt1$par
    if(!is.null(opt1$message)) {
      print(paste0("error for tol = ", thisTol, ". Message:"))
      print(opt1$message)
      next
    }
    else {
      print(paste0("completed optimization for tol = ", thisTol, ""))
      
      ## Get standard errors
      print("getting standard errors...")
      sdTime = system.time(
        SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                             bias.correct = TRUE,
                             bias.correct.control = list(sd = TRUE))
      )[3]
      # SD0
      print(sdTime/60)
      #  minutes for intern=FALSE
      
      if(SD0$pdHess) {
        print("Optimization and PD hess calculation done!")
        break
      }
      else {
        print("Hessan not PD. Rerunning optimization with stricter tol...")
        
      }
    }
    
    
    
  }
  thisEndTime = proc.time()[3]
  sdTime/60
  thisTotalTime = thisEndTime - newStartTime
  totalTime = thisEndTime - startTime
  print(paste0("second optimization took ", thisTotalTime/60, " minutes"))
  print(paste0("all optimization took ", totalTime/60, " minutes"))
}

# opt0 <- nlminb(start       =    obj[['par']],
#                objective   =    obj[['fn']],
#                gradient    =    obj[['gr']],
#                lower = rep(-10, length(obj[['par']])),
#                upper = rep( 10, length(obj[['par']])),
#                control     =    list(trace=1))

# * TMB Posterior Sampling ----

if(FALSE) {
  hessTime = system.time(testHess <- hessian(testObj$fn, optPar))[3]
  hessTime/60
  # 10.28792 minutes for intern=FALSE
  eig = eigen(testHess)
  eig
  
  for(i in 1:length(eig$values)) {
    thisVal = eig$values[i]
    thisVec = eig$vectors[,i]
    barplot(eig$vectors[,i], names.arg=names(SD0$par.fixed), 
            main=paste0("Eigenvector ", i, " with value ", thisVal))
  }
  
  testObj$fn(SD0$par.fixed)
  testObj$fn(SD0$par.fixed + eig$vectors[,9]*.05)
  testObj$fn(SD0$par.fixed - eig$vectors[,9]*.05)
}



## summary(SD0, 'report')
## summary(SD0, 'fixed')

save(SD0, obj, objFull, totalTime, sdTime, file="savedOutput/ed/fitDHS.RData")
out = load("savedOutput/ed/fitDHS.RData")

gridPreds = predGrid(SD0)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\
# \hline
# (Int) & -2.51 & -2.73 & -2.66 & -2.37 & -2.29 \\
# urb & 0.47 & 0.20 & 0.30 & 0.65 & 0.72 \\
# access & -0.34 & -0.46 & -0.41 & -0.26 & -0.22 \\
# elev & -0.41 & -0.51 & -0.47 & -0.34 & -0.31 \\
# distRiversLakes & 0.41 & 0.31 & 0.35 & 0.47 & 0.50 \\
# popValsNorm & 1.08 & 0.82 & 0.91 & 1.24 & 1.33 \\
# sigmaSq & 1.83 & 0.20 & 0.35 & 3.71 & 7.86 \\
# phi & 0.23 & 0.03 & 0.06 & 0.46 & 0.65 \\
# sigmaEpsSq & 2.36 & 2.07 & 2.17 & 2.57 & 2.67 \\
# \hline
# \end{tabular}
# \end{table}
save(gridPreds, file="savedOutput/ed/gridPredsDHS.RData")
out = load("savedOutput/ed/gridPredsDHS.RData")

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPreds, file="savedOutput/ed/stratPredsDHS.RData")
save(admin1Preds, file="savedOutput/ed/admin1PredsDHS.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsDHS.RData")
out = load("savedOutput/ed/stratPredsDHS.RData")
out = load("savedOutput/ed/admin1PredsDHS.RData")
out = load("savedOutput/ed/admin2PredsDHS.RData")

summaryTabBYM2(SD0, popMat=popMatNGAThresh, 
               gridPreds=gridPreds)
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edDHS", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edDHS", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edDHS", plotNameRootAreal="Admin2")