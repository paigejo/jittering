# script for validation

# generates K sets of indices in the MICS dataset corresponding to K folds of 
# stratified validation. The fold index is added to a new column of the edMICS 
# dataset named 'fold'. The new dataset is names edMICSval.
# Inputs:
# K: number of folds
# seed: random number seed for reproducibility
stratKmics = function(K=10, seed=123) {
  set.seed(seed)
  
  out = load("savedOutput/global/edMICS.RData")
  
  allInds = 1:nrow(edMICS)
  allStrata = edMICS$StratumID
  allUrb = edMICS$urban
  uniqueStrata = sort(unique(edMICS$StratumID))
  urbVals = c(FALSE, TRUE)
  
  # make a new dataset with 'fold' column
  foldVec = numeric(nrow(edMICS))
  
  thisIndList = function(thisStrat, thisUrb) {
    thisInds = allInds[(allStrata == thisStrat) & (allUrb == thisUrb)]
    
    thisN = length(thisInds)
    
    if(thisN == 0) {
      return(foldVec)
    } else if(thisN < K) {
      warning(paste0("Some folds will have no data for strat ", thisStrat, ", urb=", thisUrb, " with ", thisN, " pts."))
    }
    
    nBigFolds = ((thisN-1) %% K) + 1
    nSmallFolds = K - nBigFolds
    bigSize = ceiling(thisN / K)
    smallSize = bigSize - 1
    
    bigInds = sample(1:K, nBigFolds)
    
    foldIndsList = list()
    for(i in 1:K) {
      thisFoldSize = ifelse(i %in% bigInds, bigSize, smallSize)
      
      if(thisFoldSize == 0) {
        # do nothing, don't update foldVec
        next
      } else if(length(thisInds) == 1) {
        # sample() behaves weirdly if length(thisInds) == 1. I think it calls sample.int in that case
        if(thisFoldSize == 1) {
          tempInds = thisInds
        } else {
          stop(paste0("there is 1 index left in the stratum but thisFoldSize=", thisFoldSize))
        }
      } else {
        tempInds = sample(thisInds, thisFoldSize)
      }
      
      thisInds = thisInds[-match(tempInds, thisInds)]
      foldVec[tempInds] = i
    }
    
    foldVec
  }
  
  for(i in 1:length(uniqueStrata)) {
    thisStrat = uniqueStrata[i]
    
    foldVec = thisIndList(thisStrat, TRUE)
    foldVec = thisIndList(thisStrat, FALSE)
  }
  
  edMICSval = edMICS
  edMICSval$fold = foldVec
  
  save(edMICSval, file="savedOutput/validation/edMICSval.RData")
  
  invisible(edMICSval)
}




# generates K sets of indices in the DHS dataset corresponding to K folds of 
# stratified validation. The fold index is added to a new column of the edDHS 
# dataset named 'fold'. The new dataset is names edMICSval.
# Inputs:
# K: number of folds
# seed: random number seed for reproducibility
stratKdhs = function(K=10, seed=12) {
  set.seed(seed)
  
  out = load("savedOutput/global/ed.RData")
  
  allInds = 1:nrow(ed)
  allStrata = ed$admin1
  allUrb = ed$urban
  uniqueStrata = sort(unique(ed$admin1))
  urbVals = c(FALSE, TRUE)
  
  # make a new dataset with 'fold' column
  foldVec = numeric(nrow(ed))
  
  thisIndList = function(thisStrat, thisUrb) {
    thisInds = allInds[(allStrata == thisStrat) & (allUrb == thisUrb)]
    
    thisN = length(thisInds)
    
    if(thisN == 0) {
      return(foldVec)
    } else if(thisN < K) {
      warning(paste0("Some folds will have no data for strat ", thisStrat, ", urb=", thisUrb, " with ", thisN, " pts."))
    }
    
    nBigFolds = ((thisN-1) %% K) + 1
    nSmallFolds = K - nBigFolds
    bigSize = ceiling(thisN / K)
    smallSize = bigSize - 1
    
    bigInds = sample(1:K, nBigFolds)
    
    foldIndsList = list()
    for(i in 1:K) {
      thisFoldSize = ifelse(i %in% bigInds, bigSize, smallSize)
      
      if(thisFoldSize == 0) {
        # do nothing, don't update foldVec
        next
      } else if(length(thisInds) == 1) {
        # sample() behaves weirdly if length(thisInds) == 1. I think it calls sample.int in that case
        if(thisFoldSize == 1) {
          tempInds = thisInds
        } else {
          stop(paste0("there is 1 index left in the stratum but thisFoldSize=", thisFoldSize))
        }
      } else {
        tempInds = sample(thisInds, thisFoldSize)
      }
      
      thisInds = thisInds[-match(tempInds, thisInds)]
      foldVec[tempInds] = i
    }
    
    foldVec
  }
  
  for(i in 1:length(uniqueStrata)) {
    thisStrat = uniqueStrata[i]
    
    foldVec = thisIndList(thisStrat, TRUE)
    foldVec = thisIndList(thisStrat, FALSE)
  }
  
  edVal = ed
  edVal$fold = foldVec
  
  save(edVal, file="savedOutput/validation/edVal.RData")
  
  invisible(edVal)
}

# submodels of the BYM2 model

getValidationDataM_d = function(fold, inSample=TRUE) {
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  # load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLInds = edVal$fold != fold
  inSampleLIndsUrb = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRur = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
  inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
  outOfSampleLInds = edVal$fold == fold
  edInSample = edVal[inSampleLInds,]
  edOutOfSample = edVal[outOfSampleLInds,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrb)
  nPtsRurDHS = sum(inSampleLIndsRur)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now construct the standard inputs to TMB (we will modify them later based on 
  # in sample indices)
  
  # do some precomputation ----
  
  # use integration points
  out = load("savedOutput/global/intPtsDHS.RData")
  
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2]
  ARurDHS = ARurDHS[,inSampleLIndsRur2]
  
  # extract cluster information (in the correct order)
  ysUrbDHS = edInSample$y[edInSample$urban]
  ysRurDHS = edInSample$y[!edInSample$urban]
  nsUrbDHS = edInSample$n[edInSample$urban]
  nsRurDHS = edInSample$n[!edInSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2, times=KrurDHS),-1]
  
  # save(AUrbDHS, ARurDHS, intPtsDHS, 
  #      ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
  #      file="savedOutput/global/edInputsDHS.RData")
  
  # load in TMB function inputs
  # out = load("savedOutput/global/edInputsDHS.RData")
  
  # modify weights based on fold indices ----
  intPtsDHS$wRural = intPtsDHS$wRural[inSampleLIndsRur2,]
  intPtsDHS$wUrban = intPtsDHS$wUrban[inSampleLIndsUrb2,]
  intPtsDHS$wUrban[,1] = 1
  intPtsDHS$wUrban[,-1] = 0
  intPtsDHS$wRural[,1] = 1
  intPtsDHS$wRural[,-1] = 0
  
  
  # set priors ----
  alpha_pri = c(0, 100^2)
  beta_pri = c(0, 10^2)
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3)
  lambdaTauEps = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for nugget precision
  
  # Specify starting values for TMB params ----
  tmb_params <- list(alpha = 0, # intercept
                     beta = rep(0, ncol(intPtsDHS$covsUrb)), 
                     log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 0, # SPDE parameter related to the range
                     Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)) # RE on mesh vertices
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
    lambdaTauEps=lambdaTauEps, 
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
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
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterDHS'))
}

getValidationDataM_D = function(fold, inSample=TRUE) {
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  # load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLInds = edVal$fold != fold
  inSampleLIndsUrb = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRur = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
  inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
  outOfSampleLInds = edVal$fold == fold
  edInSample = edVal[inSampleLInds,]
  edOutOfSample = edVal[outOfSampleLInds,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrb)
  nPtsRurDHS = sum(inSampleLIndsRur)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now construct the standard inputs to TMB (we will modify them later based on 
  # in sample indices)
  
  # do some precomputation ----
  
  # use integration points
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2]
  ARurDHS = ARurDHS[,inSampleLIndsRur2]
  
  # extract cluster information (in the correct order)
  ysUrbDHS = edInSample$y[edInSample$urban]
  ysRurDHS = edInSample$y[!edInSample$urban]
  nsUrbDHS = edInSample$n[edInSample$urban]
  nsRurDHS = edInSample$n[!edInSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2, times=KrurDHS),-1]
  # save(AUrbDHS, ARurDHS, intPtsDHS, 
  #      ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
  #      file="savedOutput/global/edInputsDHS.RData")
  
  # load in TMB function inputs
  # out = load("savedOutput/global/edInputsDHS.RData")
  
  # set priors ----
  alpha_pri = c(0, 100^2)
  beta_pri = c(0, 10^2)
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3)
  lambdaTauEps = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for nugget precision
  
  # Specify starting values for TMB params ----
  tmb_params <- list(alpha = 0, # intercept
                     beta = rep(0, ncol(intPtsDHS$covsUrb)), 
                     log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 0, # SPDE parameter related to the range
                     Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)) # RE on mesh vertices
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
    lambdaTau=lambdaTau, # determines PC prior for tau, BYM2 precision
    lambdaTauEps=lambdaTauEps, # determines PC prior for tauEps, nugget precision
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
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
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterDHS'))
}

# fold is 1-20, with 1-10 removing part of DHS data and 11-20 removing part of MICS data
getValidationDataM_dm = function(fold) {
  foldMICS = fold - 10
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  # load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLIndsDHS = edVal$fold != fold
  inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
  inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
  outOfSampleLIndsDHS = edVal$fold == fold
  outOfSampleLIndsUrbDHS = (edVal$fold == fold) & (edVal$urban)
  outOfSampleLIndsRurDHS = (edVal$fold == fold) & (!edVal$urban)
  outOfSampleLIndsUrb2DHS = outOfSampleLIndsDHS[edVal$urban]
  outOfSampleLIndsRur2DHS = outOfSampleLIndsDHS[!edVal$urban]
  edInSample = edVal[inSampleLIndsDHS,]
  edOutOfSample = edVal[outOfSampleLIndsDHS,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrbDHS)
  nPtsRurDHS = sum(inSampleLIndsRurDHS)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now load MICS data
  # out = load("savedOutput/validation/edMICSval.RData")
  out = load("savedOutput/validation/simEdMICS.RData")
  edMICSval = simEdMICS[[fold]]
  
  inSampleLIndsMICS = edMICSval$fold != foldMICS
  inSampleLIndsUrbMICS = (edMICSval$fold != foldMICS) & (edMICSval$urban)
  inSampleLIndsRurMICS = (edMICSval$fold != foldMICS) & (!edMICSval$urban)
  inSampleLIndsUrb2MICS = inSampleLIndsMICS[edMICSval$urban]
  inSampleLIndsRur2MICS = inSampleLIndsMICS[!edMICSval$urban]
  outOfSampleLIndsMICS = edMICSval$fold == foldMICS
  outOfSampleLIndsUrbMICS = (edMICSval$fold == foldMICS) & (edMICSval$urban)
  outOfSampleLIndsRurMICS = (edMICSval$fold == foldMICS) & (!edMICSval$urban)
  outOfSampleLIndsUrb2MICS = outOfSampleLIndsMICS[edMICSval$urban]
  outOfSampleLIndsRur2MICS = outOfSampleLIndsMICS[!edMICSval$urban]
  edMICSInSample = edMICSval[inSampleLIndsMICS,]
  edMICSOutOfSample = edMICSval[outOfSampleLIndsMICS,]
  
  nPtsUrbMICS = sum(inSampleLIndsUrbMICS)
  nPtsRurMICS = sum(inSampleLIndsRurMICS)
  # KMICS = nrow(intPtsMICS$covsUrb)/sum(edMICSval$urban)
  KMICS=25
  
  # do some precomputation ----
  
  # make integration points if necessary
  # intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
  #                                           numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  out = load("savedOutput/global/intPtsDHS.RData")
  out = load("savedOutput/global/intPtsMICS.RData")
  
  
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHSoutOfSample = t(AUrbDHS[,outOfSampleLIndsUrb2DHS])
  ARurDHSoutOfSample = t(ARurDHS[,outOfSampleLIndsRur2DHS])
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2DHS]
  ARurDHS = ARurDHS[,inSampleLIndsRur2DHS]
  
  # extract cluster information (in the correct order)
  ysUrbDHS = edInSample$y[edInSample$urban]
  ysRurDHS = edInSample$y[!edInSample$urban]
  nsUrbDHS = edInSample$n[edInSample$urban]
  nsRurDHS = edInSample$n[!edInSample$urban]
  
  ysUrbDHSoutOfSample = edOutOfSample$y[edOutOfSample$urban]
  ysRurDHSoutOfSample = edOutOfSample$y[!edOutOfSample$urban]
  nsUrbDHSoutOfSample = edOutOfSample$n[edOutOfSample$urban]
  nsRurDHSoutOfSample = edOutOfSample$n[!edOutOfSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # remove rows of out of sample covariates
  intPtsDHS$covsUrbOutOfSample = intPtsDHS$covsUrb[rep(outOfSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRurOutOfSample = intPtsDHS$covsRur[rep(outOfSampleLIndsRur2DHS, times=KrurDHS),-1]
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratRur = rowSums(ARurMICS)
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  
  # now do the same for the out of sample data if need be
  if(fold > 10) {
    # first extract only the relevant covariates for the in sample data
    XUrbOutOfSample = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
    stratUrb = XUrbOutOfSample$strat
    XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    AUrbMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    numPerStratUrbOutOfSample = rowSums(AUrbMICSOutOfSample)
    # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample * KMICS))
    # obsIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), KMICS)
    # intPtIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), each=KMICS)
    actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    XUrbOutOfSample = XUrbOutOfSample[actualIndexUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    
    XRurOutOfSample = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
    stratRur = XRurOutOfSample$strat
    XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    ARurMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[!edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    numPerStratRurOutOfSample = rowSums(ARurMICSOutOfSample)
    # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample * KMICS))
    # obsIndexRur = rep(1:sum(numPerStratRurOutOfSample), KMICS)
    # intPtIndexRur = rep(1:sum(numPerStratRurOutOfSample), each=KMICS)
    actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    XRurOutOfSample = XRurOutOfSample[actualIndexRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
  }
  else {
    AUrbMICSOutOfSample = NULL
    XUrbOutOfSample = NULL
    ARurMICSOutOfSample = NULL
    XRurOutOfSample = NULL
  }
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  if(fold > 10) {
    # Do the same with the out of sample observations
    wUrbanOutOfSample = intPtsMICS$wUrban
    stratIndexUrbW = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample))
    wUrbanOutOfSample = wUrban[stratIndexUrbW,]
    
    wRuralOutOfSample = intPtsMICS$wRural
    stratIndexRurW = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample))
    wRuralOutOfSample = wRural[stratIndexRurW,]
  } else {
    wUrbanOutOfSample = NULL
    wRuralOutOfSample = NULL
  }
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  edMICSInSample = edMICSInSample[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICSInSample[edMICSInSample$urban,]$ys
  nsUrbMICS = edMICSInSample[edMICSInSample$urban,]$ns
  ysRurMICS = edMICSInSample[!edMICSInSample$urban,]$ys
  nsRurMICS = edMICSInSample[!edMICSInSample$urban,]$ns
  
  ysUrbMICSoutOfSample = edMICSOutOfSample$ys[edMICSOutOfSample$urban]
  ysRurMICSoutOfSample = edMICSOutOfSample$ys[!edMICSOutOfSample$urban]
  nsUrbMICSoutOfSample = edMICSOutOfSample$ns[edMICSOutOfSample$urban]
  nsRurMICSoutOfSample = edMICSOutOfSample$ns[!edMICSOutOfSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  if(fold > 10) {
    AUrbMICSOutOfSample = t(AUrbMICSOutOfSample)
    ARurMICSOutOfSample = t(ARurMICSOutOfSample)
  }
  AUrbMICS = t(AUrbMICS)
  ARurMICS = t(ARurMICS)
  mode(AUrbMICS) = "numeric"
  mode(ARurMICS) = "numeric"
  
  # create final list of objects for analysis to save
  intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRur = XRur[,-(2:3)]
  intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  
  intPtsMICS$XUrbOutOfSample = XUrbOutOfSample[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRurOutOfSample = XRurOutOfSample[,-(2:3)]
  if(fold > 10) {
    intPtsMICS$XUrbOutOfSample = as.matrix(intPtsMICS$XUrbOutOfSample)
    intPtsMICS$XRurOutOfSample = as.matrix(intPtsMICS$XRurOutOfSample)
  }
  intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
  intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
  
  # put weight only on the simulated locations for both DHS and MICS data
  intPtsDHS$wUrbanOutOfSample = intPtsDHS$wUrban[outOfSampleLIndsUrb2DHS,]
  intPtsDHS$wRuralOutOfSample = intPtsDHS$wRural[outOfSampleLIndsRur2DHS,]
  if(fold < 11) {
    intPtsDHS$wUrbanOutOfSample[,1] = 1
    intPtsDHS$wUrbanOutOfSample[,-1] = 0
    intPtsDHS$wRuralOutOfSample[,1] = 1
    intPtsDHS$wRuralOutOfSample[,-1] = 0
  }
  
  intPtsDHS$wUrban = intPtsDHS$wUrban[inSampleLIndsUrb2DHS,]
  intPtsDHS$wRural = intPtsDHS$wRural[inSampleLIndsRur2DHS,]
  intPtsDHS$wUrban[,-1] = 0
  intPtsDHS$wUrban[,1] = 1
  intPtsDHS$wRural[,-1] = 0
  intPtsDHS$wRural[,1] = 1
  
  intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
  intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
  
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  intPtsMICS$wUrban[,-1] = 0
  intPtsMICS$wUrban[,1] = 1
  intPtsMICS$wRural[,-1] = 0
  intPtsMICS$wRural[,1] = 1
  
  if(fold > 10) {
    intPtsMICS$wUrbanOutOfSample[,-1] = 0
    intPtsMICS$wUrbanOutOfSample[,1] = 1
    intPtsMICS$wRuralOutOfSample[,-1] = 0
    intPtsMICS$wRuralOutOfSample[,1] = 1
  }
  
  # update the MICS covariates to the ones from the simulated locations for the first 
  # column
  urbCovsMICS = edMICSInSample[edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XUrb[1:sum(edMICSInSample$urban),] = matrix(unlist(urbCovsMICS), ncol=ncol(urbCovsMICS))
  rurCovsMICS = edMICSInSample[!edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XRur[1:sum(!edMICSInSample$urban),] = matrix(unlist(rurCovsMICS), ncol=ncol(rurCovsMICS))
  
  if(fold > 10) {
    # same for the out of sample simulated locations
    urbCovsMICSOutOfSample = edMICSOutOfSample[edMICSOutOfSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
    intPtsMICS$XUrbOutOfSample[1:sum(edMICSOutOfSample$urban),] = matrix(unlist(urbCovsMICSOutOfSample), ncol=ncol(urbCovsMICSOutOfSample))
    rurCovsMICSOutOfSample = edMICSOutOfSample[!edMICSOutOfSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
    intPtsMICS$XRurOutOfSample[1:sum(!edMICSOutOfSample$urban),] = matrix(unlist(rurCovsMICSOutOfSample), ncol=ncol(rurCovsMICSOutOfSample))
  }
  
  # set priors ----
  alpha_pri = c(0, 100^2)
  beta_pri = c(0, 10^2)
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for bym2 precision
  lambdaTauEps = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for nugget precision
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('Epsilon_bym2', 'nuggetUrbMICS', 'nuggetRurMICS', 
                 'nuggetUrbDHS', 'nuggetRurDHS')
  
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
    lambdaTauEps=lambdaTauEps, 
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  dataOutOfSample = list(
    y_iUrbanDHS=ysUrbDHSoutOfSample, # same as above but for DHS survey
    y_iRuralDHS=ysRurDHSoutOfSample, # 
    n_iUrbanDHS=nsUrbDHSoutOfSample, # number binomial trials
    n_iRuralDHS=nsRurDHSoutOfSample, # 
    AprojUrbanDHS=AUrbDHSoutOfSample, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralDHS=ARurDHSoutOfSample, # 
    X_betaUrbanDHS=intPtsDHS$covsUrbOutOfSample, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralDHS=intPtsDHS$covsRurOutOfSample, # 
    wUrbanDHS=intPtsDHS$wUrbanOutOfSample, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralDHS=intPtsDHS$wRuralOutOfSample, 
    
    y_iUrbanMICS=ysUrbMICSoutOfSample, # same as above but for MICS survey
    y_iRuralMICS=ysRurMICSoutOfSample, # 
    n_iUrbanMICS=nsUrbMICSoutOfSample, # number binomial trials
    n_iRuralMICS=nsRurMICSoutOfSample, # 
    AprojUrbanMICS=AUrbMICSOutOfSample, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralMICS=ARurMICSOutOfSample, # 
    X_betaUrbanMICS=intPtsMICS$XUrbOutOfSample, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralMICS=intPtsMICS$XRurOutOfSample, # 
    wUrbanMICS=intPtsMICS$wUrbanOutOfSample, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralMICS=intPtsMICS$wRuralOutOfSample
  )
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  tmb_params <- list(alpha = initAlpha, # intercept
                     beta = c(initBeta1, rep(0, ncol(intPtsMICS$XUrb)-1)), 
                     log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 0, # SPDE parameter related to the range
                     log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS)), 
                     nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                     nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
  )
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterFusionNugget'), 
       dataOutOfSample=dataOutOfSample)
}

# fold is 1-20, with 1-10 removing part of DHS data and 11-20 removing part of MICS data
getValidationDataM_DM = function(fold) {
  foldMICS = fold - 10
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  # load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLIndsDHS = edVal$fold != fold
  inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
  inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
  outOfSampleLIndsDHS = edVal$fold == fold
  outOfSampleLIndsUrbDHS = (edVal$fold == fold) & (edVal$urban)
  outOfSampleLIndsRurDHS = (edVal$fold == fold) & (!edVal$urban)
  outOfSampleLIndsUrb2DHS = outOfSampleLIndsDHS[edVal$urban]
  outOfSampleLIndsRur2DHS = outOfSampleLIndsDHS[!edVal$urban]
  edInSample = edVal[inSampleLIndsDHS,]
  edOutOfSample = edVal[outOfSampleLIndsDHS,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrbDHS)
  nPtsRurDHS = sum(inSampleLIndsRurDHS)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now load MICS data
  out = load("savedOutput/validation/edMICSval.RData")
  
  inSampleLIndsMICS = edMICSval$fold != foldMICS
  inSampleLIndsUrbMICS = (edMICSval$fold != foldMICS) & (edMICSval$urban)
  inSampleLIndsRurMICS = (edMICSval$fold != foldMICS) & (!edMICSval$urban)
  inSampleLIndsUrb2MICS = inSampleLIndsMICS[edMICSval$urban]
  inSampleLIndsRur2MICS = inSampleLIndsMICS[!edMICSval$urban]
  outOfSampleLIndsMICS = edMICSval$fold == foldMICS
  outOfSampleLIndsUrbMICS = (edMICSval$fold == foldMICS) & (edMICSval$urban)
  outOfSampleLIndsRurMICS = (edMICSval$fold == foldMICS) & (!edMICSval$urban)
  outOfSampleLIndsUrb2MICS = outOfSampleLIndsMICS[edMICSval$urban]
  outOfSampleLIndsRur2MICS = outOfSampleLIndsMICS[!edMICSval$urban]
  edMICSInSample = edMICSval[inSampleLIndsMICS,]
  edMICSOutOfSample = edMICSval[outOfSampleLIndsMICS,]
  
  nPtsUrbMICS = sum(inSampleLIndsUrbMICS)
  nPtsRurMICS = sum(inSampleLIndsRurMICS)
  # KMICS = nrow(intPtsMICS$covsUrb)/sum(edMICSval$urban)
  KMICS=25
  
  # do some precomputation ----
  
  # make integration points if necessary
  # intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
  #                                           numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  out = load("savedOutput/global/intPtsDHS.RData")
  out = load("savedOutput/global/intPtsMICS.RData")
  
  
  AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHSoutOfSample = t(AUrbDHS[,outOfSampleLIndsUrb2DHS])
  ARurDHSoutOfSample = t(ARurDHS[,outOfSampleLIndsRur2DHS])
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2DHS]
  ARurDHS = ARurDHS[,inSampleLIndsRur2DHS]
  
  # extract cluster information (in the correct order)
  ysUrbDHS = edInSample$y[edInSample$urban]
  ysRurDHS = edInSample$y[!edInSample$urban]
  nsUrbDHS = edInSample$n[edInSample$urban]
  nsRurDHS = edInSample$n[!edInSample$urban]
  
  ysUrbDHSoutOfSample = edOutOfSample$y[edOutOfSample$urban]
  ysRurDHSoutOfSample = edOutOfSample$y[!edOutOfSample$urban]
  nsUrbDHSoutOfSample = edOutOfSample$n[edOutOfSample$urban]
  nsRurDHSoutOfSample = edOutOfSample$n[!edOutOfSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # remove rows of out of sample covariates
  intPtsDHS$covsUrbOutOfSample = intPtsDHS$covsUrb[rep(outOfSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRurOutOfSample = intPtsDHS$covsRur[rep(outOfSampleLIndsRur2DHS, times=KrurDHS),-1]
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  
  # first extract only the relevant covariates for the in sample data
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratRur = rowSums(ARurMICS)
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  
  # now do the same for the out of sample data if need be
  if(fold > 10) {
    # first extract only the relevant covariates for the in sample data
    XUrbOutOfSample = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
    stratUrb = XUrbOutOfSample$strat
    XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    AUrbMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    numPerStratUrbOutOfSample = rowSums(AUrbMICSOutOfSample)
    # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample * KMICS))
    # obsIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), KMICS)
    # intPtIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), each=KMICS)
    actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    XUrbOutOfSample = XUrbOutOfSample[actualIndexUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    
    XRurOutOfSample = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
    stratRur = XRurOutOfSample$strat
    XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    ARurMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[!edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    numPerStratRurOutOfSample = rowSums(ARurMICSOutOfSample)
    # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample * KMICS))
    # obsIndexRur = rep(1:sum(numPerStratRurOutOfSample), KMICS)
    # intPtIndexRur = rep(1:sum(numPerStratRurOutOfSample), each=KMICS)
    actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    XRurOutOfSample = XRurOutOfSample[actualIndexRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
  }
  else {
    AUrbMICSOutOfSample = NULL
    XUrbOutOfSample = NULL
    ARurMICSOutOfSample = NULL
    XRurOutOfSample = NULL
  }
  
  # w matrices are nStrata x K. They should be nObs x K. Start with the in sample observations
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  if(fold > 10) {
    # Do the same with the out of sample observations
    wUrbanOutOfSample = intPtsMICS$wUrban
    stratIndexUrbW = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample))
    wUrbanOutOfSample = wUrban[stratIndexUrbW,]
    
    wRuralOutOfSample = intPtsMICS$wRural
    stratIndexRurW = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample))
    wRuralOutOfSample = wRural[stratIndexRurW,]
  } else {
    wUrbanOutOfSample = NULL
    wRuralOutOfSample = NULL
  }
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  edMICSInSample = edMICSInSample[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICSInSample[edMICSInSample$urban,]$ys
  nsUrbMICS = edMICSInSample[edMICSInSample$urban,]$ns
  ysRurMICS = edMICSInSample[!edMICSInSample$urban,]$ys
  nsRurMICS = edMICSInSample[!edMICSInSample$urban,]$ns
  
  ysUrbMICSoutOfSample = edMICSOutOfSample$y[edMICSOutOfSample$urban]
  ysRurMICSoutOfSample = edMICSOutOfSample$y[!edMICSOutOfSample$urban]
  nsUrbMICSoutOfSample = edMICSOutOfSample$ns[edMICSOutOfSample$urban]
  nsRurMICSoutOfSample = edMICSOutOfSample$ns[!edMICSOutOfSample$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  if(fold > 10) {
    AUrbMICSOutOfSample = t(AUrbMICSOutOfSample)
    ARurMICSOutOfSample = t(ARurMICSOutOfSample)
  }
  AUrbMICS = t(AUrbMICS)
  ARurMICS = t(ARurMICS)
  mode(AUrbMICS) = "numeric"
  mode(ARurMICS) = "numeric"
  
  # create final list of objects for analysis to save
  intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRur = XRur[,-(2:3)]
  intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  
  intPtsMICS$XUrbOutOfSample = XUrbOutOfSample[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRurOutOfSample = XRurOutOfSample[,-(2:3)]
  if(fold > 10) {
    intPtsMICS$XUrbOutOfSample = as.matrix(intPtsMICS$XUrbOutOfSample)
    intPtsMICS$XRurOutOfSample = as.matrix(intPtsMICS$XRurOutOfSample)
  }
  intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
  intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
  
  # put weight only on the simulated locations for both DHS and MICS data
  intPtsDHS$wUrbanOutOfSample = intPtsDHS$wUrban[outOfSampleLIndsUrb2DHS,]
  intPtsDHS$wRuralOutOfSample = intPtsDHS$wRural[outOfSampleLIndsRur2DHS,]
  # intPtsDHS$wUrbanOutOfSample[,1] = 1
  # intPtsDHS$wUrbanOutOfSample[,-1] = 0
  # intPtsDHS$wRuralOutOfSample[,1] = 1
  # intPtsDHS$wRuralOutOfSample[,-1] = 0
  
  intPtsDHS$wUrban = intPtsDHS$wUrban[inSampleLIndsUrb2DHS,]
  intPtsDHS$wRural = intPtsDHS$wRural[inSampleLIndsRur2DHS,]
  # intPtsDHS$wUrban[,-1] = 0
  # intPtsDHS$wUrban[,1] = 1
  # intPtsDHS$wRural[,-1] = 0
  # intPtsDHS$wRural[,1] = 1
  
  intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
  intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
  # intPtsMICS$wUrbanOutOfSample[,1] = 1
  # intPtsMICS$wUrbanOutOfSample[,-1] = 0
  # intPtsMICS$wRuralOutOfSample[,1] = 1
  # intPtsMICS$wRuralOutOfSample[,-1] = 0
  
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  # intPtsMICS$wUrban[,-1] = 0
  # intPtsMICS$wUrban[,1] = 1
  # intPtsMICS$wRural[,-1] = 0
  # intPtsMICS$wRural[,1] = 1
  
  # set priors ----
  alpha_pri = c(0, 100^2)
  beta_pri = c(0, 10^2)
  
  out = load("savedOutput/global/admFinalMat.RData")
  bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                           constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  lambdaTau = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for bym2 precision
  lambdaTauEps = getLambdaPCprec(u=0.5, alpha=2/3) # get PC prior lambda for nugget precision
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('Epsilon_bym2', 'nuggetUrbMICS', 'nuggetRurMICS', 
                 'nuggetUrbDHS', 'nuggetRurDHS')
  
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
    lambdaTauEps=lambdaTauEps, 
    options=0 # 1 for adreport of log tau and logit phi
  )
  
  dataOutOfSample = list(
    y_iUrbanDHS=ysUrbDHSoutOfSample, # same as above but for DHS survey
    y_iRuralDHS=ysRurDHSoutOfSample, # 
    n_iUrbanDHS=nsUrbDHSoutOfSample, # number binomial trials
    n_iRuralDHS=nsRurDHSoutOfSample, # 
    AprojUrbanDHS=AUrbDHSoutOfSample, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralDHS=ARurDHSoutOfSample, # 
    X_betaUrbanDHS=intPtsDHS$covsUrbOutOfSample, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralDHS=intPtsDHS$covsRurOutOfSample, # 
    wUrbanDHS=intPtsDHS$wUrbanOutOfSample, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralDHS=intPtsDHS$wRuralOutOfSample, 
    
    y_iUrbanMICS=ysUrbMICSoutOfSample, # same as above but for MICS survey
    y_iRuralMICS=ysRurMICSoutOfSample, # 
    n_iUrbanMICS=nsUrbMICSoutOfSample, # number binomial trials
    n_iRuralMICS=nsRurMICSoutOfSample, # 
    AprojUrbanMICS=AUrbMICSOutOfSample, # nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
    AprojRuralMICS=ARurMICSOutOfSample, # 
    X_betaUrbanMICS=intPtsMICS$XUrbOutOfSample, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
    X_betaRuralMICS=intPtsMICS$XRurOutOfSample, # 
    wUrbanMICS=intPtsMICS$wUrbanOutOfSample, # nObsUrban x nIntegrationPointsUrban weight matrix
    wRuralMICS=intPtsMICS$wRuralOutOfSample
  )
  
  # initial parameters
  initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
  initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
  initAlpha = logit(initRurP)
  initBeta1 = logit(initUrbP) - initAlpha
  
  tmb_params <- list(alpha = initAlpha, # intercept
                     beta = c(initBeta1, rep(0, ncol(intPtsMICS$XUrb)-1)), 
                     log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = 0, # SPDE parameter related to the range
                     log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                     Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)), 
                     nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS)), 
                     nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                     nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
  )
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterFusionNugget'), 
       dataOutOfSample=dataOutOfSample)
}

# This function generates and saves all the validation datasets and input 
# parameters for all models
getAllValidationData = function(folds=1:20) {
  
  # first generate M_d data
  print("generating data for M_d...")
  time1 = system.time(datMd <- lapply(folds, getValidationDataM_d))[3]
  print(paste0("Took ", time1, " seconds. Now generating data for M_D..."))
  time2 = system.time(datMD <- lapply(folds, getValidationDataM_D))[3]
  print(paste0("Took ", time2, " seconds. Now generating data for M_dm..."))
  time3 = system.time(datMdm <- lapply(folds, getValidationDataM_dm))[3]
  print(paste0("Took ", time3, " seconds. Now generating data for M_DM..."))
  time4 = system.time(datMDM <- lapply(folds, getValidationDataM_DM))[3]
  print(paste0("Took ", time4, " seconds. Now saving results..."))
  
  save(datMd, file="savedOutput/validation/datMd.RData")
  save(datMD, file="savedOutput/validation/datM_D.RData")
  save(datMdm, file="savedOutput/validation/datMdm.RData")
  save(datMDM, file="savedOutput/validation/datM_DM.RData")
  
  invisible(list(datMd, datMD, datMdm, datMDM))
}

getValidationFit = function(fold, model=c("Md", "MD", "Mdm", "MDM"), regenModFit=FALSE, randomBeta=FALSE) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  # load in the data for the appropriate model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  }
  fname = paste0("savedOutput/validation/dat", fnameRoot, ".RData")
  out = load(fname)
  
  # get the data from the appropriate variable. The data is in the following format:
  # list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
  #      edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
  #      MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
  #                           hessian=TRUE, DLL='modBYM2JitterFusionNugget'))
  varname = paste0("dat", model, collapse="")
  dat = get(varname)[[fold]]
  
  edInSample = dat$edInSample
  edMICSInSample = dat$edMICSInSample
  edOutOfSample = dat$edOutOfSample
  edMICSOutOfSample = dat$edMICSOutOfSample
  
  # set starting values to full M_DM optimum
  dat$MakeADFunInputs$parameters$alpha = -1.7905689 # intercept
  dat$MakeADFunInputs$parameters$beta = c(.6, -.62, .19, .13, .12)
  dat$MakeADFunInputs$parameters$log_tau = 0.6306069 # Log tau (i.e. log spatial precision, Epsilon)
  dat$MakeADFunInputs$parameters$logit_phi = 2.4493192 # SPDE parameter related to the range
  dat$MakeADFunInputs$parameters$log_tauEps = -0.4203604 # Log tau (i.e. log spatial precision, Epsilon)
  if(randomBeta) {
    dat$MakeADFunInputs$random = c("beta", dat$MakeADFunInputs$random)
  }
  
  MakeADFunInputs = dat$MakeADFunInputs
  MakeADFunInputsFull = MakeADFunInputs
  MakeADFunInputsFull$random = NULL
  
  # make sure we set the out of sample data correctly for 'DHS data only' models
  if(is.null(edMICSOutOfSample)) {
    out = load("savedOutput/validation/edMICSval.RData")
    edMICSOutOfSample = edMICSval
  }
  
  if(regenModFit || !file.exists(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))) {
    # now fit the model. First we load DLLs and build the functions then we optimize
    dyn.load(dynlib(paste0("code/", MakeADFunInputs$DLL)))
    # browser()
    if(FALSE) {
      dimLen = function(x) {
        out = dim(x)
        if(is.null(out)) {
          length(x)
        } else {
          out
        }
      }
      
      tempDataFull = tempMakeADFunInputs$data
      tempPar = tempMakeADFunInputs$parameters
      
      sapply(MakeADFunInputs$data, dimLen)
      sapply(tempDataFull, dimLen)
      sapply(MakeADFunInputs$data, class)
      sapply(tempDataFull, class)
      lapply(MakeADFunInputs$data, str)
      sapply(tempDataFull, str)
      sapply(MakeADFunInputs$parameters, dimLen)
      sapply(tempPar, dimLen)
      # Browse[1]>       sapply(MakeADFunInputs$parameters, dimLen)
      # alpha          beta       log_tau     logit_phi    log_tauEps  Epsilon_bym2 
      # 1             5             1             1             1            41 
      # nuggetUrbMICS nuggetRurMICS  nuggetUrbDHS  nuggetRurDHS 
      # 662          1307           569           810 
      # Browse[1]>       sapply(tempPar, dimLen)
      # alpha          beta       log_tau     logit_phi    log_tauEps  Epsilon_bym2 
      # 1             5             1             1             1            41 
      # nuggetUrbMICS nuggetRurMICS  nuggetUrbDHS  nuggetRurDHS 
      # 734          1449           569           810 
      sapply(MakeADFunInputs$data, anyNA)
      sapply(MakeADFunInputs$parameters, anyNA)
      
      head(tempDataFull$X_betaUrbanMICS)
    }
    obj <- do.call("MakeADFun", MakeADFunInputs)
    objFull <- do.call("MakeADFun", MakeADFunInputsFull)
    
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
      tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
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
          sdTime = system.time(
            SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                 bias.correct = TRUE,
                                 bias.correct.control = list(sd = TRUE))
          )[3]
          # SD0
          print(paste0("SE calculations took ", sdTime/60, " minutes"))
          
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
    }
    
    if(!SD0$pdHess) {
      # stop("Hessian not PD")
      warning("Hessian not PD...")
      hessPD = FALSE
    } else {
      hessPD = TRUE
    }
    
    # save fit model before generating predictions
    save(SD0, obj, objFull, totalTime, sdTime, hessPD, file=paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  } else {
    out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  }
  
  # predict at the left out clusters
  if(hessPD) {
    preds = predClusters(nsim=1000, fold, SD0, obj, 
                         model=model, 
                         quantiles=c(0.025, 0.1, 0.9, 0.975))
  } else {
    preds = NULL
  }
  
  save(SD0, obj, objFull, totalTime, sdTime, hessPD, preds, file=paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  
  allScores = scoreValidationPreds(fold, model=model, regenScores=TRUE)
  
  list(SD0, obj, objFull, totalTime, sdTime, hessPD, allScores)
}

# make predictions for a set of clusters of one type (MICS or DHS)
# SD0: TMB sdreport object
# obj: TMB fun object
# fold: cross validation fold from 1-20
# ptType: The type of integration points. either MICS or DHS
# model: which model we're making predictions with
predClusters = function(nsim=1000, fold, SD0, obj, 
                        model=c("Md", "MD", "Mdm", "MDM"), 
                        quantiles=c(0.025, 0.1, 0.9, 0.975), 
                        addBinVar=TRUE) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  # set the file name root depending on the model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  }
  
  # for loading the info for predicting at the left out data locations, if the 
  # model doesn't use MICS data use the M_DM info to make the predictions for 
  # curiosity's sake.
  fnameRootLeftOut = ifelse(fnameRoot %in% c("Md", "M_D"), "M_DM", fnameRoot)
  modelLeftOut = ifelse(fnameRoot %in% c("Md", "M_D"), "MDM", model)
  
  # load relevant data and model fit
  out = load(paste0("savedOutput/validation/dat", fnameRootLeftOut, ".RData", collapse=""))
  leftOutDat = get(paste0("dat", modelLeftOut))[[fold]]$dataOutOfSample
  out = load(paste0("savedOutput/validation/dat", fnameRoot, ".RData", collapse=""))
  load("~/git/jittering/savedOutput/validation/edMICSval.RData")
  load("~/git/jittering/savedOutput/validation/edVal.RData")
  
  varname = paste0("dat", model)
  dat = get(varname)[[fold]]
  
  foldMod = ifelse((fold > 11) && (model %in% c("Md", "MD")), 11, fold)
  out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", foldMod, ".RData"))
  
  # generate predictions at the left out clusters
  
  # generate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  
  sigmaEpsSq_tmb_draws = NULL
  if(SD0$pdHess) {
    L <- Cholesky(SD0[['jointPrecision']], super = T)
    mu = summary(SD0)[,1]
    t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = nsim)
    
    # extract fixed effects and random effects from draws
    parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
    epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_bym2',]
    alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
    beta_tmb_draws    <- t.draws[parnames == 'beta',]
    sigmaSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tau',]), nrow = 1)
    phi_tmb_draws    <- matrix(expit(t.draws[parnames == 'logit_phi',]), nrow = 1)
    
    fixedMat = rbind(alpha_tmb_draws, 
                     beta_tmb_draws, 
                     sigmaSq_tmb_draws, 
                     phi_tmb_draws)
    betaNames = colnames(obj$env$.data$X_betaUrbanDHS)
    row.names(fixedMat) = c("(Int)", 
                            betaNames, 
                            "sigmaSq", 
                            "phi")
    
    hasNugget = "log_tauEps" %in% row.names(summary(SD0))
    if(hasNugget) {
      sigmaEpsSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEps',]), nrow = 1)
      fixedMat = rbind(fixedMat, 
                       sigmaEpsSq_tmb_draws)
      row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
    }
    
    # Make parameter summary tables
    parMeans = rowMeans(fixedMat)
    parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
    parSummary = cbind(parMeans, parQuants)
    colnames(parSummary)[1] = "Est"
    colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
    print(xtable(parSummary, digits=2))
    
    # add effects to predictions
    if(fold <= 10) {
      Aurb = leftOutDat$AprojUrbanDHS
      Arur = leftOutDat$AprojRuralDHS
      Xurb = leftOutDat$X_betaUrbanDHS
      Xrur = leftOutDat$X_betaRuralDHS
      wUrb = leftOutDat$wUrbanDHS
      wRur = leftOutDat$wRuralDHS
      yUrb = leftOutDat$y_iUrbanDHS
      yRur = leftOutDat$y_iRuralDHS
      nUrb = leftOutDat$n_iUrbanDHS
      nRur = leftOutDat$n_iRuralDHS
    } else {
      Aurb = leftOutDat$AprojUrbanMICS
      Arur = leftOutDat$AprojRuralMICS
      Xurb = leftOutDat$X_betaUrbanMICS
      Xrur = leftOutDat$X_betaRuralMICS
      wUrb = leftOutDat$wUrbanMICS
      wRur = leftOutDat$wRuralMICS
      yUrb = leftOutDat$y_iUrbanMICS
      yRur = leftOutDat$y_iRuralMICS
      nUrb = leftOutDat$n_iUrbanMICS
      nRur = leftOutDat$n_iRuralMICS
    }
    
    # expand A matrices to size of Xmat (project to integration points instead of clusters)
    Kurb = nrow(Xurb) / nrow(Aurb)
    Krur = nrow(Xrur) / nrow(Arur)
    bigAurb = matrix(rep(Aurb, times=Kurb), ncol=ncol(Aurb))
    bigArur = matrix(rep(Arur, times=Krur), ncol=ncol(Arur))
    
    # get latent preds at cluster integration points
    clustIntDrawsUrb <- as.matrix(bigAurb %*% epsilon_tmb_draws)
    clustIntDrawsUrb <- sweep(clustIntDrawsUrb, 2, alpha_tmb_draws, '+')
    clustIntDrawsUrb <- clustIntDrawsUrb + (Xurb %*% beta_tmb_draws)
    
    clustIntDrawsRur <- as.matrix(bigArur %*% epsilon_tmb_draws)
    clustIntDrawsRur <- sweep(clustIntDrawsRur, 2, alpha_tmb_draws, '+')
    clustIntDrawsRur <- clustIntDrawsRur + (Xrur %*% beta_tmb_draws)
    
    # convert predictions to probability scale
    if(!hasNugget) {
      probIntDrawsUrb = expit(clustIntDrawsUrb)
      probIntDrawsRur = expit(clustIntDrawsRur)
    }
    else {
      nClustUrb = nrow(clustIntDrawsUrb)/Kurb
      nClustRur = nrow(clustIntDrawsRur)/Krur
      clustIDUrb = rep(1:nClustUrb, Kurb)
      clustIDRur = rep(1:nClustRur, Krur)
      # probIntDrawsUrb = matrix(logitNormMean(cbind(c(clustIntDrawsUrb), rep(sqrt(sigmaEpsSq_tmb_draws), each=nrow(clustIntDrawsUrb))), logisticApprox=FALSE, splineApprox=TRUE), nrow=nrow(clustIntDrawsUrb))
      # probIntDrawsRur = matrix(logitNormMean(cbind(c(clustIntDrawsRur), rep(sqrt(sigmaEpsSq_tmb_draws), each=nrow(clustIntDrawsRur))), logisticApprox=FALSE, splineApprox=TRUE), nrow=nrow(clustIntDrawsRur))
      logitIntDrawsUrb = sapply(1:ncol(clustIntDrawsUrb), function(colI) {
        thisNuggetSD = sqrt(sigmaEpsSq_tmb_draws[colI])
        nugsUrb = rnorm(nClustUrb, sd=thisNuggetSD)
        clustIntDrawsUrb[,colI] + rep(nugsUrb, times=Kurb)
      })
      probIntDrawsUrb = expit(logitIntDrawsUrb)
      
      logitIntDrawsRur = sapply(1:ncol(clustIntDrawsRur), function(colI) {
        thisNuggetSD = sqrt(sigmaEpsSq_tmb_draws[colI])
        nugsRur = rnorm(nClustRur, sd=thisNuggetSD)
        clustIntDrawsRur[,colI] + rep(nugsRur, times=Krur)
      })
      probIntDrawsRur = expit(logitIntDrawsRur)
    }
    
    # take weighted average of predictions at integration points (i.e. evaluate integral of predictions for each cluster numerically)
    # We will make block diagonal Wurb and Wrur matrices, where element ij is the integration weight for cluster i associated with integration point j
    # buildRowUrb = c(rep(1, Kurb), rep(0, nrow(Xurb)))
    # Wurb = matrix(c(rep(buildRowUrb, times=length(yUrb)-1), rep(1, Kurb)), byrow=TRUE, ncol=nrow(Xurb))
    # Wurb = sweep(Wurb, 2, c(t(wUrb)), FUN="*")
    
    buildMatUrb = rbind(c(wUrb), 
                        matrix(0, ncol=Kurb*nrow(wUrb), nrow=nrow(wUrb)))
    leaveOutInds = (length(buildMatUrb)-nrow(wUrb) + 1):length(buildMatUrb)
    Wurb = matrix(c(buildMatUrb)[-leaveOutInds], nrow=nrow(wUrb))
    zeroCols = seq(nrow(Wurb)+1, ncol(Wurb), by=nrow(Wurb)+1)
    Wurb = Wurb[,-zeroCols]
    probDrawsUrb = Wurb %*% probIntDrawsUrb
    
    # buildRowRur = c(rep(1, Krur), rep(0, nrow(Xrur)))
    # Wrur = matrix(c(rep(buildRowRur, times=length(yRur)-1), rep(1, Krur)), byrow=TRUE, ncol=nrow(Xrur))
    # Wrur = sweep(Wrur, 2, c(t(wRur)), FUN="*")
    
    buildMatRur = rbind(c(wRur), 
                        matrix(0, ncol=Krur*nrow(wRur), nrow=nrow(wRur)))
    leaveOutInds = (length(buildMatRur)-nrow(wRur) + 1):length(buildMatRur)
    Wrur = matrix(c(buildMatRur)[-leaveOutInds], nrow=nrow(wRur))
    zeroCols = seq(nrow(Wrur)+1, ncol(Wrur), by=nrow(Wrur)+1)
    Wrur = Wrur[,-zeroCols]
    probDrawsRur = Wrur %*% probIntDrawsRur
    
    # calculate central prediction before binomial variation is added in
    predsUrb = rowMeans(probDrawsUrb)
    predsRur = rowMeans(probDrawsRur)
    
    # add in binomial variation
    probDrawsUrb = addBinomialVar(probDrawsUrb, nUrb)
    probDrawsRur = addBinomialVar(probDrawsRur, nRur)
    
    quantsUrb = apply(probDrawsUrb, 1, quantile, probs=quantiles, na.rm=TRUE)
    quantsRur = apply(probDrawsRur, 1, quantile, probs=quantiles, na.rm=TRUE)
  }
  else {
    # in the case that the hessian is not PD
    stop("Hessian not PD")
    browser()
    Eps = SD0$par.random[grepl("Epsilon", names(SD0$par.random))]
    alpha = SD0$par.fixed[grepl("alpha", names(SD0$par.fixed))]
    beta = SD0$par.fixed[grepl("beta", names(SD0$par.fixed))]
    
    # set "draws" to be just the fixed values
    epsilon_tmb_draws = Eps
    alpha_tmb_draws = alpha
    beta_tmb_draws = beta
    phi_tmb_draws = expit(SD0$par.fixed[grepl("logit_phi", names(SD0$par.fixed))])
    sigmaSq_tmb_draws = 1/exp(SD0$par.fixed[grepl("log_tau", names(SD0$par.fixed))])
    
    # add effects to predictions
    gridDraws_tmb <- Amat %*% Eps
    gridDraws_tmb <- gridDraws_tmb + alpha
    gridDraws_tmb <- gridDraws_tmb + (Xmat %*% beta)
    
    if(!hasNugget) {
      probDraws = expit(gridDraws_tmb)
    }
    else {
      tauEps = exp(SD0$par.fixed[grepl("log_tauEps", names(SD0$par.fixed))])
      sigmaEps = 1/sqrt(tauEps)
      probDraws = logitNormMean(cbind(c(gridDraws_tmb), rep(sigmaEps, length(gridDraws_tmb))), logisticApprox=FALSE)
      
      sigmaEpsSq_tmb_draws = sigmaEps^2
    }
    
    preds = probDraws
    quants = NULL
  }
  
  list(probDrawsUrb=probDrawsUrb, probDrawsRur=probDrawsRur, 
       predsUrb=predsUrb, predsRur=predsRur, 
       parSummary=parSummary, fixedMat=fixedMat, 
       quantsUrb=quantsUrb, quantsRur=quantsRur, 
       yUrb=yUrb, yRur=yRur, nUrb=nUrb, nRur=nRur)
}

scoreValidationPreds = function(fold, model=c("Md", "MD", "Mdm", "MDM"), regenScores=FALSE) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  # load in the data for the appropriate model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  }
  
  # load predictions
  out = load(paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  # list(probDrawsUrb, probDrawsRur, predsUrb, predsRur, 
  #      quantsUrb, quantsRur, yUrb, yRur, nUrb, nRur)
  
  if(!is.null(preds)) {
    probDrawsUrb = preds$probDrawsUrb
    predsUrb = preds$predsUrb
    quantsUrb = preds$quantsUrb
    yUrb = preds$yUrb
    nUrb = preds$nUrb
    probDrawsRur = preds$probDrawsRur
    predsRur = preds$predsRur
    quantsRur = preds$quantsRur
    yRur = preds$yRur
    nRur = preds$nRur
    
    # combine predictions from urban and rural areas
    probDraws = rbind(probDrawsUrb, probDrawsRur)
    preds = c(predsUrb, predsRur)
    quants = c(quantsUrb, quantsRur)
    ys = c(yUrb, yRur)
    ns = c(nUrb, nRur)
    
    # calculate score
    scoresUrb = getScores(truth=yUrb/nUrb, estMat=probDrawsUrb, weights=nUrb, 
                          significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, getAverage=TRUE, na.rm=TRUE)
    scoresRur = getScores(truth=yRur/nRur, estMat=probDrawsRur, weights=nRur, 
                          significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, getAverage=TRUE, na.rm=TRUE)
    scores = getScores(truth=ys/ns, estMat=probDraws, weights=ns, 
                       significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, getAverage=TRUE, na.rm=TRUE)
    
    # add computation time
    scoresUrb = cbind(scoresUrb, Time=totalTime/60)
    scoresRur = cbind(scoresRur, Time=totalTime/60)
    scores = cbind(scores, Time=totalTime/60)
  } else {
    scoresUrb = NULL
    scoresRur = NULL
    scores = NULL
  }
  
  save(scoresUrb, scoresRur, scores, preds, file=paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData"))
  
  invisible(NULL)
}


# function for collecting validation results for each model
validationTable = function(quantiles=c(0.025, 0.1, 0.9, 0.975)) {
  models=c("Md", "MD", "Mdm", "MDM")
  folds = 1:20
  
  thisSummaryFun = function(x) {
    c(mean(x), sd(x), quantile(x, probs=quantiles))
  }
  
  # calculate the fold weights (overall, urban, and rural)
  # Note that the scores within each fold average have already been weighted by 
  # the cluster ns so we just have to weight by the fold total ns here
  out = load("savedOutput/validation/edVal.RData")
  out = load("savedOutput/validation/edMICSval.RData")
  weightsDHS = aggregate(edVal$n, by=list(fold=edVal$fold), FUN=sum)$x
  weightsMICS = aggregate(edMICSval$ns, by=list(fold=edMICSval$fold), FUN=sum)$x
  weightsUrbDHS = aggregate(edVal$n[edVal$urban], by=list(fold=edVal$fold[edVal$urban]), FUN=sum)$x
  weightsUrbMICS = aggregate(edMICSval$ns[edMICSval$urban], by=list(fold=edMICSval$fold[edMICSval$urban]), FUN=sum)$x
  weightsRurDHS = aggregate(edVal$n[!edVal$urban], by=list(fold=edVal$fold[!edVal$urban]), FUN=sum)$x
  weightsRurMICS = aggregate(edMICSval$ns[!edMICSval$urban], by=list(fold=edMICSval$fold[!edMICSval$urban]), FUN=sum)$x
  weightsDHS = weightsDHS/sum(weightsDHS)
  weightsUrbDHS = weightsUrbDHS/sum(weightsUrbDHS)
  weightsRurDHS = weightsRurDHS/sum(weightsRurDHS)
  weightsMICS = weightsMICS/sum(weightsMICS)
  weightsUrbMICS = weightsUrbMICS/sum(weightsUrbMICS)
  weightsRurMICS = weightsRurMICS/sum(weightsRurMICS)
  
  # aggregate the results of all folds for each model
  scoresTabsDHS = list()
  scoresTabsUrbDHS = list()
  scoresTabsRurDHS = list()
  parTabsDHS = list()
  scoresTabsMICS = list()
  scoresTabsUrbMICS = list()
  scoresTabsRurMICS = list()
  parTabsMICS = list()
  for(i in 1:length(models)) {
    model = models[i]
    
    # get the file name root for the results of this model
    fnameRoot = model
    if(fnameRoot == "MD") {
      fnameRoot = "M_D"
    } else if(fnameRoot == "MDM") {
      fnameRoot = "M_DM"
    }
    
    # collect the results over all folds
    thisScoresTabDHS = list()
    thisScoresTabUrbDHS = list()
    thisScoresTabRurDHS = list()
    thisParTabDHS = list()
    thisScoresTabMICS = list()
    thisScoresTabUrbMICS = list()
    thisScoresTabRurMICS = list()
    thisParTabMICS = list()
    for(j in 1:length(folds)) {
      fold = folds[j]
      isMICS = j <= 10
      
      # load the scores
      out = load(paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData"))
      out = load(paste0("~/git/jittering/savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
      
      # calculate parameter summary statistics
      foldParTab = t(apply(t(preds$fixedMat), 2, thisSummaryFun))
      colnames(foldParTab) = c("Est", "SD", paste0("Q", quantiles*100))
      
      if(!isMICS) {
        thisScoresTabDHS = rbind(thisScoresTabDHS, scores)
        thisScoresTabUrbDHS = rbind(thisScoresTabUrbDHS, scoresUrb)
        thisScoresTabRurDHS = rbind(thisScoresTabRurDHS, scoresRur)
        
        thisParTabDHS = c(thisParTabDHS, list(foldParTab))
      } else {
        thisScoresTabMICS = rbind(thisScoresTabMICS, scores)
        thisScoresTabUrbMICS = rbind(thisScoresTabUrbMICS, scoresUrb)
        thisScoresTabRurMICS = rbind(thisScoresTabRurMICS, scoresRur)
        
        thisParTabMICS = c(thisParTabMICS, list(foldParTab))
      }
    }
    scoresTabsDHS = c(scoresTabsDHS, list(thisScoresTabDHS))
    scoresTabsUrbDHS = c(scoresTabsUrbDHS, list(thisScoresTabUrbDHS))
    scoresTabsRurDHS = c(scoresTabsRurDHS, list(thisScoresTabRurDHS))
    parTabsDHS = c(parTabsDHS, list(thisParTabDHS))
    scoresTabsMICS = c(scoresTabsMICS, list(thisScoresTabMICS))
    scoresTabsUrbMICS = c(scoresTabsUrbMICS, list(thisScoresTabUrbMICS))
    scoresTabsRurMICS = c(scoresTabsRurMICS, list(thisScoresTabRurMICS))
    parTabsMICS = c(parTabsMICS, list(thisParTabMICS))
  }
  
  # calculate averages
  scoresTabsAvgDHS = do.call("rbind", lapply(scoresTabsDHS, function(x) {colSums(sweep(x, 1, weightsDHS, "*"))}))
  scoresTabsUrbAvgDHS = do.call("rbind", lapply(scoresTabsUrbDHS, function(x) {colSums(sweep(x, 1, weightsUrbDHS, "*"))}))
  scoresTabsRurAvgDHS = do.call("rbind", lapply(scoresTabsRurDHS, function(x) {colSums(sweep(x, 1, weightsRurDHS, "*"))}))
  scoresTabsAvgMICS = do.call("rbind", lapply(scoresTabsMICS, function(x) {colSums(sweep(x, 1, weightsUrbMICS, "*"))}))
  scoresTabsUrbAvgMICS = do.call("rbind", lapply(scoresTabsUrbMICS, function(x) {colSums(sweep(x, 1, weightsUrbMICS, "*"))}))
  scoresTabsRurAvgMICS = do.call("rbind", lapply(scoresTabsRurMICS, function(x) {colSums(sweep(x, 1, weightsRurMICS, "*"))}))
  
  parTabsAvgDHS = lapply(parTabsDHS, function(x) {
    rNames = row.names(x[[1]])
    cNames = colnames(x[[1]])
    x = lapply(x, function(y) {array(y, dim=c(dim(y), 1))})
    scoreArray = do.call("abind", list(x, along=3))
    out = apply(scoreArray, 1:2, mean)
    row.names(out) = rNames
    colnames(out) = cNames
    out
  })
  names(parTabsAvgDHS) = models
  parTabsAvgMICS = lapply(parTabsMICS, function(x) {
    rNames = row.names(x[[1]])
    cNames = colnames(x[[1]])
    x = lapply(x, function(y) {array(y, dim=c(dim(y), 1))})
    scoreArray = do.call("abind", list(x, along=3))
    out = apply(scoreArray, 1:2, mean)
    row.names(out) = rNames
    colnames(out) = cNames
    out
  })
  names(parTabsAvgMICS) = models
  
  parTabsAvg = lapply(1:length(models), function (i) {
    rNames = row.names(parTabsAvgDHS[[1]])
    cNames = colnames(parTabsAvgDHS[[1]])
    thisParTab = abind(parTabsAvgDHS[[i]], parTabsAvgMICS[[i]], along=3)
    thisParTab = apply(thisParTab, 1:2, mean)
    
    row.names(thisParTab) = rNames
    colnames(thisParTab) = cNames
    thisParTab
  })
  names(parTabsAvg) = models
  
  browser()
  
}

