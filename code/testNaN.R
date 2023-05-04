# script for testing women's secondary education in Nigeria application
# see: https://www.ics.uci.edu/~pattis/common/handouts/macmingweclipse/allexperimental/mac-gdb-install.html
# for installing GBD
# run TMB::gdbsource("code/testNaN.R", interactive=TRUE)

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
  compile("code/modBYM2JitterDHS.cpp","-O0 -g") # for using gbdsource
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

# Specify starting values for TMB params ----
tmb_params <- list(alpha = 0, # intercept
                   beta = rep(0, ncol(intPtsDHS$covsUrb)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)) # RE on mesh vertices
)

## specify random effects
rand_effs <- c('Epsilon_bym2')

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

badPar = c(-0.0437403001,  0.0061653050,  0.0148973164, 0.0081539719, 
           -0.0006482087, -0.0229490176, -0.1992381562, 0.0117170980)
print(obj[['fn']](badPar))
print(obj[['gr']](badPar))