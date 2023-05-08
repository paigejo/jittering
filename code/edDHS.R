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
lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3)
lambdaTauEps = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for nugget precision

# Specify starting values for TMB params ----
# initial parameters
initUrbP = sum(c(data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanDHS))
initRurP = sum(c(data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

tmb_params <- list(alpha = initAlpha, # intercept
                   beta = c(initBeta1, rep(0, ncol(intPtsMICS$XUrb)-1)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                   nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                   nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
)

## specify random effects
rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS')

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

## make the autodiff generated likelihood func & gradient
dyn.load( dynlib("code/modBYM2JitterDHS"))
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterDHS')

# * Run TMB ----
lower = rep(-10, length(obj[['par']]))
upper = rep( 10, length(obj[['par']]))
funWrapper = function(pars) {
  print(pars)
  obj[['fn']](pars)
}
grWrapper = function(pars) {
  print(pars)
  obj[['gr']](pars)
}

badPar = c(-1.46335890,  1.73814507, -0.10405578,  0.08469777, -0.02948654,  0.16331678, -1.60453975, -0.22204565 )
obj[['fn']](badPar)
obj[['gr']](badPar)
obj[['he']](badPar)
library(numDeriv)
grad(obj[['fn']], badPar)

test = obj$report()
test$quad
test$logDet
test$logDetTau
test$logPriPhi
test$KLD
test$d
test$logPriPhi
test$lliksUrb

opt0 <- nlminb(start       =    obj[['par']],
               objective   =    funWrapper,
               gradient    =    grWrapper,
               lower       =    lower,
               upper       =    upper,
               control     =    list(trace=1))

testOut = optim(obj[['par']], funWrapper, control=list(trace=4))
testOut$value
# [1] -10.89829
optPar1 = testOut$par
optPar1
#      alpha       beta       beta       beta       beta       beta    log_tau  logit_phi 
# -0.3507440  0.6236230  0.7150636  1.2023273 -0.1570947 -0.2156714 -2.5024082 -0.3730605 

opt0 <- nlminb(start       =    optPar1,
               objective   =    funWrapper,
               gradient    =    grWrapper,
               lower       =    lower,
               upper       =    upper,
               control     =    list(trace=1))
opt0
# $par
# alpha       beta       beta       beta       beta       beta    log_tau  logit_phi 
# -0.3507549  0.6236147  0.7150841  1.2023305 -0.1570978 -0.2156813 -2.5023804 -0.3730772 
# $objective
# [1] -7.839766
# $convergence
# [1] 1
# $message
# [1] "false convergence (8)"
obj[["gr"]](opt0$par)
# outer mgc:  2287113 
# [1]  1020303.31   744758.58 -2287112.89 -1011361.45    38339.58  1826511.89 -1695171.90  1564203.59
grad(funWrapper, opt0$par)
load("savedOutput/h.RData")
funWrapper(opt0$par + h)
funWrapper(opt0$par - h)
test = obj$report()
test$quad
test$logDet
test$logDetTau
test$logPriPhi
test$KLD
test$d
test$logPriPhi
test$lliksUrb

testObjFunWrapper = function(par, badParVal=1e10) {
  if(any(par < lower) || any(par > upper)) {
    return(badParVal)
  }
  print(par)
  objVal = testObj[['fn']](par)
  parNames = names(par)
  parVals = par
  parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
  parStr = paste(parStrs, sep=", ")
  print(paste0("objective: ", objVal, " for parameters, ", parStr))
  objVal
}

tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
for(thisTol in tolSeq) {
  testObj = obj
  testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
  testObj$env$tracepar = TRUE
  print(paste0("optimizing for tol = ", thisTol, "."))
  tryCatch({
    opt1 <- optim(par=testObj$par, fn = testObjFunWrapper, gr = testObj[['gr']],
                 method = c("BFGS"), hessian = FALSE)
  }, error={
    print(paste0("error for tol = ", thisTol, "."))
  }, finally={
    print(paste0("completed optimization for tol = ", thisTol, "."))
  })
}



# * TMB Posterior Sampling ----
## Get standard errors
SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                     bias.correct = TRUE,
                     bias.correct.control = list(sd = TRUE))
## summary(SD0, 'report')
## summary(SD0, 'fixed')

# plot predictions


## take samples from fitted model
mu <- c(SD0$par.fixed,SD0$par.random)

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
