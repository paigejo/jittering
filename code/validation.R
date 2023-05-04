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
      warning(paste("Some folds will have no data for strat ", thisStrat, ", urb=", thisUrb, " with ", thisN, " pts."))
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
stratKdhs = function(K=10, seed=123) {
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
      warning(paste("Some folds will have no data for strat ", thisStrat, ", urb=", thisUrb, " with ", thisN, " pts."))
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

getValidationDataM_d = function(fold) {
  # first load the DHS data
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
  
  # make integration points if necessary
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
    options=0 # 1 for adreport of log tau and logit phi
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

getValidationDataM_D = function(fold) {
  # first load the DHS data
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
  
  # make integration points if necessary
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
    options=0 # 1 for adreport of log tau and logit phi
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
  
  # first load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLIndsDHS = edVal$fold != fold
  inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2DHS = inSampleLInds[edVal$urban]
  inSampleLIndsRur2DHS = inSampleLInds[!edVal$urban]
  outOfSampleLIndsDHS = edVal$fold == fold
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
  inSampleLIndsUrb2MICS = inSampleLInds[edMICSval$urban]
  inSampleLIndsRur2MICS = inSampleLInds[!edMICSval$urban]
  outOfSampleLIndsMICS = edMICSval$fold == foldMICS
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
  
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2DHS]
  ARurDHS = ARurDHS[,inSampleLIndsRur2DHS]
  
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
  
  # remove rows of out of sample covariates
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  XUrb = XUrb[stratIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
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
  stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  edMICSInSample = edMICSInSample[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICSInSample[edMICSInSample$urban,]$ys
  nsUrbMICS = edMICSInSample[edMICSInSample$urban,]$ns
  ysRurMICS = edMICSInSample[!edMICSInSample$urban,]$ys
  nsRurMICS = edMICSInSample[!edMICSInSample$urban,]$ns
  
  # make sure A matrices are nArea x nObs, as TMB expects
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
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  
  # put weight only on the simulated locations for both DHS and MICS data
  intPtsMICS$wUrban[,-1] = 0
  intPtsMICS$wUrban[,1] = 1
  intPtsMICS$wRural[,-1] = 0
  intPtsMICS$wRural[,1] = 1
  
  intPtsDHS$wUrban[,-1] = 0
  intPtsDHS$wUrban[,1] = 1
  intPtsDHS$wRural[,-1] = 0
  intPtsDHS$wRural[,1] = 1
  
  # update the MICS covariates to the ones from the simulated locations for the first 
  # column
  urbCovsMICS = edMICSInSample[edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XUrb[1:sum(edMICSInSample$urban),] = matrix(unlist(urbCovsMICS), ncol=ncol(urbCovsMICS))
  rurCovsMICS = edMICSInSample[!edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XRur[1:sum(!edMICSInSample$urban),] = matrix(unlist(rurCovsMICS), ncol=ncol(rurCovsMICS))
  
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
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterFusionNugget'))
}

# fold is 1-20, with 1-10 removing part of DHS data and 11-20 removing part of MICS data
getValidationDataM_DM = function(fold) {
  foldMICS = fold - 10
  
  # first load the DHS data
  out = load("savedOutput/validation/edVal.RData")
  
  inSampleLIndsDHS = edVal$fold != fold
  inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
  inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
  inSampleLIndsUrb2DHS = inSampleLInds[edVal$urban]
  inSampleLIndsRur2DHS = inSampleLInds[!edVal$urban]
  outOfSampleLIndsDHS = edVal$fold == fold
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
  inSampleLIndsUrb2MICS = inSampleLInds[edMICSval$urban]
  inSampleLIndsRur2MICS = inSampleLInds[!edMICSval$urban]
  outOfSampleLIndsMICS = edMICSval$fold == foldMICS
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
  
  AUrbDHS = AUrbDHS[,inSampleLIndsUrb2DHS]
  ARurDHS = ARurDHS[,inSampleLIndsRur2DHS]
  
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
  
  # remove rows of out of sample covariates
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  numPerStratUrb = rowSums(AUrbMICS)
  stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  XUrb = XUrb[stratIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
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
  stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  edMICSInSample = edMICSInSample[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICSInSample[edMICSInSample$urban,]$ys
  nsUrbMICS = edMICSInSample[edMICSInSample$urban,]$ns
  ysRurMICS = edMICSInSample[!edMICSInSample$urban,]$ys
  nsRurMICS = edMICSInSample[!edMICSInSample$urban,]$ns
  
  # make sure A matrices are nArea x nObs, as TMB expects
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
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  
  # put weight only on the simulated locations for both DHS and MICS data
  intPtsMICS$wUrban[,-1] = 0
  intPtsMICS$wUrban[,1] = 1
  intPtsMICS$wRural[,-1] = 0
  intPtsMICS$wRural[,1] = 1
  
  intPtsDHS$wUrban[,-1] = 0
  intPtsDHS$wUrban[,1] = 1
  intPtsDHS$wRural[,-1] = 0
  intPtsDHS$wRural[,1] = 1
  
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
  
  # return TMB inputs and in/out of sample datasets. Here are the TMB inputs:
  # obj <- MakeADFun(data=data_full,
  #                  parameters=tmb_params,
  #                  random=rand_effs,
  #                  hessian=TRUE,
  #                  DLL='modBYM2JitterDHS')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL='modBYM2JitterFusionNugget'))
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

getValidationFit = function(fold, model=c("Md", "MD", "Mdm", "MDM"), regenModFit=FALSE) {
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
  fname = paste("savedOutput/validation/dat", fnameRoot, ".RData", collapse="")
  out = load(fname)
  
  # get the data from the appropriate variable. The data is in the following format:
  # list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
  #      edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
  #      MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
  #                           hessian=TRUE, DLL='modBYM2JitterFusionNugget'))
  varname = paste0("dat", model, collapse="")
  dat = get(varname)
  
  edInSample = dat$edInSample
  edMICSInSample = dat$edMICSInSample
  edOutOfSample = dat$edOutOfSample
  edMICSOutOfSample = dat$edMICSOutOfSample
  MakeADFunInputs = dat$MakeADFunInputs
  MakeADFunInputsFull = MakeADFunInputs
  MakeADFunInputsFull$random = NULL
  
  # make sure we set the out of sample data correctly for 'DHS data only' models
  if(is.null(edMICSOutOfSample)) {
    out = load("savedOutput/validation/edMICSval.RData")
    edMICSOutOfSample = edMICSval
  }
  
  if(regenModFit || !file.exists(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))) {
    # now fit the model. First we build the functions then we optimize
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
      stop("Hessian not PD")
    }
    
    # save fit model before generating predictions
    save(SD0, obj, objFull, totalTime, sdTime, file=paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  } else {
    out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  }
  
  list(SD0, obj, objFull, totalTime, sdTime)
}

# make predictions for a set of clusters of one type (MICS or DHS)
# SD0: TMB sdreport object
# obj: TMB fun object
# fold: cross validation fold from 1-20
# ptType: The type of integration points. either MICS or DHS
# model: which model we're making predictions with
predClusters = function(nsim=1000, SD0, obj, fold, 
                        model=c("Md", "MD", "Mdm", "MDM"), 
                        quantiles=c(0.025, 0.1, 0.9, 0.975)) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  }
  
  # load relevant data and model fit
  out = load(paste("savedOutput/validation/dat", fnameRoot, ".RData", collapse=""))
  out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  
  # generate predictions
}


