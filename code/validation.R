# script for validation

# generates K sets of indices in the MICS dataset corresponding to K folds of 
# stratified validation. The fold index is added to a new column of the edMICS 
# dataset named 'fold'. The new dataset is names edMICSval.
# Inputs:
# K: number of folds
# seed: random number seed for reproducibility
stratKmics = function(K=10, seed=1234) {
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
  
  # sort it by stratum
  # stratIDs = match(edMICSval$Stratum, admFinal$NAME_FINAL)
  # edMICSval = edMICSval[order(stratIDs),]
  edMICSval = edMICSval[order(edMICSval$Stratum),]
  
  save(edMICSval, file="savedOutput/validation/edMICSval.RData")
  
  invisible(edMICSval)
}




# generates K sets of indices in the DHS dataset corresponding to K folds of 
# stratified validation. The fold index is added to a new column of the edDHS 
# dataset named 'fold'. The new dataset is names edMICSval.
# Inputs:
# K: number of folds
# seed: random number seed for reproducibility
stratKdhs = function(K=10, seed=1234) {
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
  
  # sort it by stratum
  # stratIDs = match(edVal$Stratum, admFinal$NAME_FINAL)
  # edVal = edVal[order(stratIDs),]
  
  save(edVal, file="savedOutput/validation/edVal.RData")
  
  invisible(edVal)
}

# submodels of the BYM2 model

getValidationDataM_d = function(fold, admLevel=c("admFinal", "adm2"), areal=FALSE, res=300, adm2AsCovariate=TRUE) {
  KDHSurb = 11 # 3 rings of 5 each
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  admLevel = match.arg(admLevel)
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  # load("savedOutput/global/intPtsMICS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  
  # load the DHS data
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/validation/edVal.RData")
  out = load("savedOutput/validation/edMICSval.RData")
  
  # order edVal so it matches with ordering of ed
  sortI = match(ed$clusterID, edVal$clusterID)
  temp = edVal[sortI,]
  # all.equal(temp$subarea, ed$subarea)
  edVal = edVal[sortI,]
  
  areas = sort(unique(edVal$area))
  foldArea = areas[fold]
  
  if(!areal) {
    inSampleLInds = edVal$fold != fold
    inSampleLIndsUrb = (edVal$fold != fold) & (edVal$urban)
    inSampleLIndsRur = (edVal$fold != fold) & (!edVal$urban)
    inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
    inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
    outOfSampleLInds = edVal$fold == fold
  } else {
    inSampleLInds = (edVal$area != foldArea) | (edVal$fold <= 5)
    inSampleLIndsUrb = ((edVal$area != foldArea) & (edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsRur = ((edVal$area != foldArea) & (!edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
    inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
    outOfSampleLInds = (edVal$area == foldArea) & (edVal$fold > 5)
  }
  
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
  
  if(admLevel == "admFinal") {
    AUrbDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasUrban), admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
    ARurDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasRural), admFinal$NAME_FINAL) # 41 x 810
  } else {
    AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
    ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  }
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
  alpha_pri = c(0, 10^2)
  beta_pri = c(0, 5^2)
  
  if(admLevel == "admFinal") {
    out = load("savedOutput/global/admFinalMat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  } else {
    out = load("savedOutput/global/adm2Mat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  }
  lambdaTau = getLambdaPCprec(u=1, alpha=.1)
  lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision
  
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
  
  if(admLevel == "adm2") {
    AUrbDHS=as(AUrbDHS, "sparseMatrix")
    ARurDHS=as(ARurDHS, "sparseMatrix")
  }
  
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
  
  DLL = ifelse(admLevel == "admFinal", 'modBYM2JitterDHS', 'modBYM2JitterDHS2')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL=DLL))
}

getValidationDataM_D = function(fold, admLevel=c("admFinal", "adm2"), areal=FALSE, res=300, adm2AsCovariate=TRUE) {
  KDHSurb = 11 # 3 rings of 5 each
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  admLevel = match.arg(admLevel)
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  # load("savedOutput/global/intPtsMICS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  
  # load the DHS data
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/validation/edVal.RData")
  
  # order edVal so it matches with ordering of ed
  sortI = match(ed$clusterID, edVal$clusterID)
  temp = edVal[sortI,]
  # all.equal(temp$subarea, ed$subarea)
  edVal = edVal[sortI,]
  
  strata = sort(unique(edVal$Stratum))
  foldArea = strata[fold]
  
  if(!areal) {
    inSampleLInds = edVal$fold != fold
    inSampleLIndsUrb = (edVal$fold != fold) & (edVal$urban)
    inSampleLIndsRur = (edVal$fold != fold) & (!edVal$urban)
    inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
    inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
    outOfSampleLInds = edVal$fold == fold
  } else {
    inSampleLInds = (edVal$area != foldArea) | (edVal$fold <= 5)
    inSampleLIndsUrb = ((edVal$area != foldArea) & (edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsRur = ((edVal$area != foldArea) & (!edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsUrb2 = inSampleLInds[edVal$urban]
    inSampleLIndsRur2 = inSampleLInds[!edVal$urban]
    outOfSampleLInds = (edVal$area == foldArea) & (edVal$fold > 5)
  }
  
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
  if(admLevel == "admFinal") {
    AUrbDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasUrban), admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
    ARurDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasRural), admFinal$NAME_FINAL) # 41 x 810
  } else {
    AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
    ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  }
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
  alpha_pri = c(0, 10^2)
  beta_pri = c(0, 5^2)
  
  if(admLevel == "admFinal") {
    out = load("savedOutput/global/admFinalMat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  } else {
    out = load("savedOutput/global/adm2Mat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  }
  lambdaTau = getLambdaPCprec(u=1, alpha=.1)
  lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision
  
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
  if(admLevel == "adm2") {
    AUrbDHS=as(AUrbDHS, "sparseMatrix")
    ARurDHS=as(ARurDHS, "sparseMatrix")
  }
  
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
  
  DLL = ifelse(admLevel == "admFinal", 'modBYM2JitterDHS', 'modBYM2JitterDHS2')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL=DLL))
}

# fold is 1-20, with 1-10 removing part of DHS data and 11-20 removing part of MICS data
getValidationDataM_dm = function(fold, admLevel=c("admFinal", "adm2"), areal=FALSE, res=300, adm2AsCovariate=TRUE) {
  KDHSurb = 11 # 3 rings of 5 each
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  admLevel = match.arg(admLevel)
  foldMICS = fold - 10
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  # load("savedOutput/global/intPtsMICS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  
  # load the DHS data
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/validation/edVal.RData")
  
  # order edVal so it matches with ordering of ed
  sortI = match(ed$clusterID, edVal$clusterID)
  temp = edVal[sortI,]
  # all.equal(temp$subarea, ed$subarea)
  edVal = edVal[sortI,]
  
  areas = sort(unique(edVal$area))
  foldArea = areas[fold]
  
  if(!areal) {
    inSampleLIndsDHS = edVal$fold != fold
    inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
    inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
    inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
    inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
    outOfSampleLIndsDHS = edVal$fold == fold
    outOfSampleLIndsUrb2DHS = outOfSampleLIndsDHS[edVal$urban]
    outOfSampleLIndsRur2DHS = outOfSampleLIndsDHS[!edVal$urban]
  } else {
    inSampleLIndsDHS = (edVal$area != foldArea) | (edVal$fold <= 5)
    inSampleLIndsUrbDHS = ((edVal$area != foldArea) & (edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsRurDHS = ((edVal$area != foldArea) & (!edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
    inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
    outOfSampleLIndsDHS = (edVal$area == foldArea) & (edVal$fold > 5)
  }
  
  edInSample = edVal[inSampleLIndsDHS,]
  edOutOfSample = edVal[outOfSampleLIndsDHS,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrbDHS)
  nPtsRurDHS = sum(inSampleLIndsRurDHS)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now load MICS data
  # out = load("savedOutput/validation/edMICSval.RData")
  out = load("savedOutput/validation/simEdMICS.RData")
  if(!areal) {
    edMICSval = simEdMICS[[fold]]
  } else {
    edMICSval = simEdMICS[[1]]
  }
  
  if(!areal) {
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
  } else {
    inSampleLIndsMICS = (edMICSval$Area != foldArea) | (edMICSval$fold <= 5)
    inSampleLIndsUrbMICS = ((edMICSval$Area != foldArea) & (edMICSval$urban)) | (edMICSval$fold <= 5)
    inSampleLIndsRurMICS = ((edMICSval$Area != foldArea) & (!edMICSval$urban)) | (edMICSval$fold <= 5)
    inSampleLIndsUrb2MICS = inSampleLIndsMICS[edMICSval$urban]
    inSampleLIndsRur2MICS = inSampleLIndsMICS[!edMICSval$urban]
    outOfSampleLIndsMICS = (edMICSval$Area == foldArea) & (edMICSval$fold > 5)
  }
  
  edMICSInSample = edMICSval[inSampleLIndsMICS,]
  edMICSOutOfSample = edMICSval[outOfSampleLIndsMICS,]
  
  nPtsUrbMICS = sum(inSampleLIndsUrbMICS)
  nPtsRurMICS = sum(inSampleLIndsRurMICS)
  KMICS = nrow(intPtsMICS$XUrb)/41
  # KMICS=25
  
  # do some precomputation ----
  
  # make integration points if necessary
  # intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
  #                                           numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  out = load("savedOutput/global/intPtsDHS.RData")
  # out = load("savedOutput/global/intPtsMICS.RData")
  out = load(paste0("savedOutput/global/intPtsMICS", "_", res, ifelse(adm2AsCovariate, "_adm2Cov", ""), ".RData"))
  
  
  if(admLevel == "admFinal") {
    AUrbDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasUrban), admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
    ARurDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasRural), admFinal$NAME_FINAL) # 41 x 810
  } else {
    AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
    ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  }
  
  # store projection matrices to convenient prediction at out of sample clusters 
  # if necessary
  if(!areal) {
    AUrbDHSoutOfSample = t(AUrbDHS[,outOfSampleLIndsUrb2DHS])
    ARurDHSoutOfSample = t(ARurDHS[,outOfSampleLIndsRur2DHS])
  } else {
    AUrbDHSoutOfSample = NULL
    ARurDHSoutOfSample = NULL
  }
  
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
  
  if(!areal) {
    # remove rows of out of sample covariates
    intPtsDHS$covsUrbOutOfSample = intPtsDHS$covsUrb[rep(outOfSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
    intPtsDHS$covsRurOutOfSample = intPtsDHS$covsRur[rep(outOfSampleLIndsRur2DHS, times=KrurDHS),-1]
  }
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICSInSample$Stratum, by=list(strat=edMICSInSample$Stratum, urb=edMICSInSample$urban), FUN=length, drop=FALSE)
  if(areal && (nrow(allNumPerStrat) != 82)) {
    stop("check this...")
  }
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratUrb[is.na(numPerStratUrb)] = 0
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  if(admLevel == "adm2") {
    AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  } else {
    AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  }
  # numPerStratUrb = rowSums(AUrbMICS)
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # startInds = seq(1, KMICS*length(admFinal@data$NAME_FINAL), by=KMICS)
  # getInds = function(intPtI = 1, numPerStrat) {
  #   unlist(mapply(rep, startInds+intPtI-1, each=numPerStrat))
  # }
  # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrb))
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  # AUrbMICS = AUrbMICS[,actualIndexUrb]
  startStratInds = which(XUrb$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XUrb)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrb[x])})) # length nUrb, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
  nUrb = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIUrb = allAreaIs + (allIntIs-1)*nAreas
  XUrb = XUrb[transformIUrb,] # now XUrb is [K * nObsUrb] x nVar
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  if(admLevel == "adm2") {
    AUrbMICS = AUrbMICS[,transformIUrb]
  }
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  if(admLevel == "adm2") {
    ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  } else {
    ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
  }
  # numPerStratRur = rowSums(ARurMICS)
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRur))
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  # XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  # ARurMICS = ARurMICS[,actualIndexRur]
  startStratInds = which(XRur$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XRur)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRur[x])})) # length nRur, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
  nRur = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIRur = allAreaIs + (allIntIs-1)*nAreas
  XRur = XRur[transformIRur,]
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  if(admLevel == "adm2") {
    ARurMICS = ARurMICS[,transformIRur]
  }
  
  # now do the same for the out of sample data if need be
  if((fold > 10) && !areal) {
    allNumPerStrat = aggregate(edMICSOutOfSample$Stratum, by=list(strat=edMICSOutOfSample$Stratum, urb=edMICSOutOfSample$urban), FUN=length, drop=FALSE)
    numPerStratUrbOutOfSample = allNumPerStrat[allNumPerStrat[,2], 3]
    numPerStratUrbOutOfSample[is.na(numPerStratUrbOutOfSample)] = 0
    numPerStratRurOutOfSample = allNumPerStrat[!allNumPerStrat[,2], 3]
    numPerStratRurOutOfSample[is.na(numPerStratRurOutOfSample)] = 0
    
    # first extract only the relevant covariates for the in sample data
    XUrbOutOfSample = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
    stratUrb = XUrbOutOfSample$strat
    if(admLevel == "adm2") {
      AUrbMICSOutOfSample = makeApointToArea(XUrbOutOfSample$subarea, adm2$NAME_2)
    } else {
      AUrbMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    }
    # numPerStratUrbOutOfSample = rowSums(AUrbMICSOutOfSample)
    # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample * KMICS))
    # obsIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), KMICS)
    # intPtIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), each=KMICS)
    # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrbOutOfSample))
    # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    # XUrbOutOfSample = XUrbOutOfSample[actualIndexUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    # XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    # AUrbMICSOutOfSample = AUrbMICSOutOfSample[,actualIndexUrb]
    startStratInds = which(XUrbOutOfSample$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
    nAreas = nrow(XUrbOutOfSample)/KMICS
    areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrbOutOfSample[x])})) # length nUrb, range = 1:41. gives area index for each obs
    allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
    nUrb = length(allAreaIs)/KMICS
    allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
    transformIUrb = allAreaIs + (allIntIs-1)*nAreas
    XUrbOutOfSample = XUrbOutOfSample[transformIUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    if(admLevel == "adm2") {
      AUrbMICSOutOfSample = AUrbMICSOutOfSample[,transformIUrb]
    }
    
    XRurOutOfSample = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
    stratRur = XRurOutOfSample$strat
    if(admLevel == "adm2") {
      ARurMICSOutOfSample = makeApointToArea(XRurOutOfSample$subarea, adm2$NAME_2)
    } else {
      ARurMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[!edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    }
    # numPerStratRurOutOfSample = rowSums(ARurMICSOutOfSample)
    # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample * KMICS))
    # obsIndexRur = rep(1:sum(numPerStratRurOutOfSample), KMICS)
    # intPtIndexRur = rep(1:sum(numPerStratRurOutOfSample), each=KMICS)
    # actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRurOutOfSample))
    # actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    # XRurOutOfSample = XRurOutOfSample[actualIndexRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
    # XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    # ARurMICSOutOfSample = ARurMICSOutOfSample[,actualIndexRur]
    
    startStratInds = which(XRurOutOfSample$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
    nAreas = nrow(XRurOutOfSample)/KMICS
    areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRurOutOfSample[x])})) # length nRur, range = 1:41. gives area index for each obs
    allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
    nRur = length(allAreaIs)/KMICS
    allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
    transformIRur = allAreaIs + (allIntIs-1)*nAreas
    XRurOutOfSample = XRurOutOfSample[transformIRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
    XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    if(admLevel == "adm2") {
      ARurMICSOutOfSample = ARurMICSOutOfSample[,transformIRur]
    }
  }
  else {
    AUrbMICSOutOfSample = NULL
    XUrbOutOfSample = NULL
    ARurMICSOutOfSample = NULL
    XRurOutOfSample = NULL
  }
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  if((fold > 10) && !areal) {
    # Do the same with the out of sample MICS observations
    wUrbanOutOfSample = intPtsMICS$wUrban
    stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrbanOutOfSample), each=numPerStratUrbOutOfSample))
    wUrbanOutOfSample = wUrbanOutOfSample[stratIndexUrbW,]
    
    wRuralOutOfSample = intPtsMICS$wRural
    stratIndexRurW = unlist(mapply(rep, 1:nrow(wRuralOutOfSample), each=numPerStratRurOutOfSample))
    wRuralOutOfSample = wRuralOutOfSample[stratIndexRurW,]
  } else {
    wUrbanOutOfSample = NULL
    wRuralOutOfSample = NULL
  }
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  # edMICSInSample = edMICSInSample[order(stratIDs),]
  
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
  if((fold > 10) && !areal) {
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
  
  if(!areal) {
    # set up out of sample covariates and weights
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
    
    intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
    intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
    
    if(fold > 10) {
      intPtsMICS$wUrbanOutOfSample[,-1] = 0
      intPtsMICS$wUrbanOutOfSample[,1] = 1
      intPtsMICS$wRuralOutOfSample[,-1] = 0
      intPtsMICS$wRuralOutOfSample[,1] = 1
    }
  } else {
    intPtsMICS$XUrbOutOfSample = NULL
    intPtsMICS$XRurOutOfSample = NULL
    intPtsMICS$wUrbanOutOfSample = NULL
    intPtsMICS$wRuralOutOfSample = NULL
  }
  
  intPtsDHS$wUrban = intPtsDHS$wUrban[inSampleLIndsUrb2DHS,]
  intPtsDHS$wRural = intPtsDHS$wRural[inSampleLIndsRur2DHS,]
  intPtsDHS$wUrban[,-1] = 0
  intPtsDHS$wUrban[,1] = 1
  intPtsDHS$wRural[,-1] = 0
  intPtsDHS$wRural[,1] = 1
  
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  intPtsMICS$wUrban[,-1] = 0
  intPtsMICS$wUrban[,1] = 1
  intPtsMICS$wRural[,-1] = 0
  intPtsMICS$wRural[,1] = 1
  
  # update the MICS covariates to the ones from the simulated locations for the first 
  # column
  subareasUrb = edMICSInSample[edMICSInSample$urban,]$subarea
  urbCovsMICS = edMICSInSample[edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XUrb[1:sum(edMICSInSample$urban),] = matrix(unlist(urbCovsMICS), ncol=ncol(urbCovsMICS))
  subareasRur = edMICSInSample[!edMICSInSample$urban,]$subarea
  rurCovsMICS = edMICSInSample[!edMICSInSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
  intPtsMICS$XRur[1:sum(!edMICSInSample$urban),] = matrix(unlist(rurCovsMICS), ncol=ncol(rurCovsMICS))
  
  if(admLevel == "adm2") {
    AUrbMICS[1:sum(edMICSInSample$urban),] = t(makeApointToArea(subareasUrb,  adm2$NAME_2))
    ARurMICS[1:sum(!edMICSInSample$urban),] = t(makeApointToArea(subareasRur,  adm2$NAME_2))
  }
  
  if(!areal) {
    if(fold > 10) {
      # same for the out of sample simulated locations
      subareasUrb = edMICSOutOfSample[edMICSOutOfSample$urban,]$subarea
      urbCovsMICSOutOfSample = edMICSOutOfSample[edMICSOutOfSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
      intPtsMICS$XUrbOutOfSample[1:sum(edMICSOutOfSample$urban),] = matrix(unlist(urbCovsMICSOutOfSample), ncol=ncol(urbCovsMICSOutOfSample))
      subareasRur = edMICSOutOfSample[!edMICSOutOfSample$urban,]$subarea
      rurCovsMICSOutOfSample = edMICSOutOfSample[!edMICSOutOfSample$urban,][c("urb", "access", "elev", "distRiversLakes", "pop")]
      intPtsMICS$XRurOutOfSample[1:sum(!edMICSOutOfSample$urban),] = matrix(unlist(rurCovsMICSOutOfSample), ncol=ncol(rurCovsMICSOutOfSample))
      
      if(admLevel == "adm2") {
        AUrbMICSOutOfSample[1:sum(edMICSOutOfSample$urban),] = t(makeApointToArea(subareasUrb,  adm2$NAME_2))
        ARurMICSOutOfSample[1:sum(!edMICSOutOfSample$urban),] = t(makeApointToArea(subareasRur, adm2$NAME_2))
      }
    }
  }
  
  # set priors ----
  alpha_pri = c(0, 10^2)
  beta_pri = c(0, 5^2)
  
  if(admLevel == "admFinal") {
    out = load("savedOutput/global/admFinalMat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  } else {
    out = load("savedOutput/global/adm2Mat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  }
  lambdaTau = getLambdaPCprec(u=1, alpha=.1)
  lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('Epsilon_bym2', 'nuggetUrbMICS', 'nuggetRurMICS', 
                 'nuggetUrbDHS', 'nuggetRurDHS')
  
  # make sure A matrices are sparse if BYM2 is at the admin2 level
  if(admLevel == "adm2") {
    # in sample
    AUrbDHS=as(AUrbDHS, "sparseMatrix")
    ARurDHS=as(ARurDHS, "sparseMatrix")
    AUrbMICS=as(AUrbMICS, "sparseMatrix")
    ARurMICS=as(ARurMICS, "sparseMatrix")
    
    # out of sample
    if(!areal) {
      AUrbDHSoutOfSample=as(AUrbDHSoutOfSample, "sparseMatrix")
      ARurDHSoutOfSample=as(ARurDHSoutOfSample, "sparseMatrix")
      if(fold > 10) {
        AUrbMICSOutOfSample=as(AUrbMICSOutOfSample, "sparseMatrix")
        ARurMICSOutOfSample=as(ARurMICSOutOfSample, "sparseMatrix")
      }
    }
  }
  
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
  
  DLL = ifelse(admLevel == "admFinal", 'modBYM2JitterFusionNugget', 'modBYM2JitterFusionNugget2')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL=DLL), 
       dataOutOfSample=dataOutOfSample)
}

# fold is 1-20, with 1-10 removing part of DHS data and 11-20 removing part of MICS data
getValidationDataM_DM = function(fold, admLevel=c("admFinal", "adm2"), areal=FALSE, res=300, adm2AsCovariate=TRUE) {
  KDHSurb = 11 # 3 rings of 5 each
  KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
  admLevel = match.arg(admLevel)
  foldMICS = fold - 10
  
  # load in the precomputed integration points
  load("savedOutput/global/intPtsDHS.RData")
  # load("savedOutput/global/intPtsMICS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  
  # load the DHS data
  out = load("savedOutput/global/ed.RData")
  out = load("savedOutput/validation/edVal.RData")
  
  # order edVal so it matches with ordering of ed
  sortI = match(ed$clusterID, edVal$clusterID)
  temp = edVal[sortI,]
  # all.equal(temp$subarea, ed$subarea)
  edVal = edVal[sortI,]
  
  areas = sort(unique(edVal$area))
  foldArea = areas[fold]
  
  if(!areal) {
    inSampleLIndsDHS = edVal$fold != fold
    inSampleLIndsUrbDHS = (edVal$fold != fold) & (edVal$urban)
    inSampleLIndsRurDHS = (edVal$fold != fold) & (!edVal$urban)
    inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
    inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
    outOfSampleLIndsDHS = edVal$fold == fold
    outOfSampleLIndsUrb2DHS = outOfSampleLIndsDHS[edVal$urban]
    outOfSampleLIndsRur2DHS = outOfSampleLIndsDHS[!edVal$urban]
  } else {
    inSampleLIndsDHS = (edVal$area != foldArea) | (edVal$fold <= 5)
    inSampleLIndsUrbDHS = ((edVal$area != foldArea) & (edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsRurDHS = ((edVal$area != foldArea) & (!edVal$urban)) | (edVal$fold <= 5)
    inSampleLIndsUrb2DHS = inSampleLIndsDHS[edVal$urban]
    inSampleLIndsRur2DHS = inSampleLIndsDHS[!edVal$urban]
    outOfSampleLIndsDHS = (edVal$area == foldArea) & (edVal$fold > 5)
  }
  
  edInSample = edVal[inSampleLIndsDHS,]
  edOutOfSample = edVal[outOfSampleLIndsDHS,]
  
  nPtsUrbDHS = sum(inSampleLIndsUrbDHS)
  nPtsRurDHS = sum(inSampleLIndsRurDHS)
  KurbDHS = nrow(intPtsDHS$covsUrb)/sum(edVal$urban)
  KrurDHS = nrow(intPtsDHS$covsRur)/sum(!edVal$urban)
  
  # now load MICS data
  out = load("savedOutput/validation/edMICSval.RData")
  
  if(!areal) {
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
  } else {
    inSampleLIndsMICS = (edMICSval$Area != foldArea) | (edMICSval$fold <= 5)
    inSampleLIndsUrbMICS = ((edMICSval$Area != foldArea) & (edMICSval$urban)) | (edMICSval$fold <= 5)
    inSampleLIndsRurMICS = ((edMICSval$Area != foldArea) & (!edMICSval$urban)) | (edMICSval$fold <= 5)
    inSampleLIndsUrb2MICS = inSampleLIndsMICS[edMICSval$urban]
    inSampleLIndsRur2MICS = inSampleLIndsMICS[!edMICSval$urban]
    outOfSampleLIndsMICS = (edMICSval$Area == foldArea) & (edMICSval$fold > 5)
  }
  
  edMICSInSample = edMICSval[inSampleLIndsMICS,]
  edMICSOutOfSample = edMICSval[outOfSampleLIndsMICS,]
  
  nPtsUrbMICS = sum(inSampleLIndsUrbMICS)
  nPtsRurMICS = sum(inSampleLIndsRurMICS)
  KMICS = nrow(intPtsMICS$XUrb)/41
  # KMICS=25
  
  # do some precomputation ----
  
  # make integration points if necessary
  # intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
  #                                           numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  
  out = load("savedOutput/global/intPtsDHS.RData")
  out = load(paste0("savedOutput/global/intPtsMICS", "_", res, ifelse(adm2AsCovariate, "_adm2Cov", ""), ".RData"))
  
  
  if(admLevel == "admFinal") {
    AUrbDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasUrban), admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
    ARurDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasRural), admFinal$NAME_FINAL) # 41 x 810
  } else {
    AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
    ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  }
  
  # store projection matrices to convenient prediction at out of sample clusters 
  # if necessary
  if(!areal) {
    AUrbDHSoutOfSample = t(AUrbDHS[,outOfSampleLIndsUrb2DHS])
    ARurDHSoutOfSample = t(ARurDHS[,outOfSampleLIndsRur2DHS])
  } else {
    AUrbDHSoutOfSample = NULL
    ARurDHSoutOfSample = NULL
  }
  
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
  
  if(!areal) {
    # remove rows of out of sample covariates
    intPtsDHS$covsUrbOutOfSample = intPtsDHS$covsUrb[rep(outOfSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
    intPtsDHS$covsRurOutOfSample = intPtsDHS$covsRur[rep(outOfSampleLIndsRur2DHS, times=KrurDHS),-1]
  }
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[rep(inSampleLIndsUrb2DHS, times=KurbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[rep(inSampleLIndsRur2DHS, times=KrurDHS),-1]
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICSInSample$Stratum, by=list(strat=edMICSInSample$Stratum, urb=edMICSInSample$urban), FUN=length, drop=FALSE)
  if(areal && (nrow(allNumPerStrat) != 82)) {
    stop("check this...")
  }
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratUrb[is.na(numPerStratUrb)] = 0
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  stratUrb = XUrb$strat
  if(admLevel == "adm2") {
    AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  } else {
    AUrbMICS = makeApointToArea(edMICSInSample$Stratum[edMICSInSample$urban], admFinal$NAME_FINAL)
  }
  # numPerStratUrb = rowSums(AUrbMICS)
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # startInds = seq(1, KMICS*length(admFinal@data$NAME_FINAL), by=KMICS)
  # getInds = function(intPtI = 1, numPerStrat) {
  #   unlist(mapply(rep, startInds+intPtI-1, each=numPerStrat))
  # }
  # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrb))
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  # AUrbMICS = AUrbMICS[,actualIndexUrb]
  startStratInds = which(XUrb$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XUrb)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrb[x])})) # length nUrb, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
  nUrb = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIUrb = allAreaIs + (allIntIs-1)*nAreas
  XUrb = XUrb[transformIUrb,] # now XUrb is [K * nObsUrb] x nVar
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  if(admLevel == "adm2") {
    AUrbMICS = AUrbMICS[,transformIUrb]
  }
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  stratRur = XRur$strat
  if(admLevel == "adm2") {
    ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  } else {
    ARurMICS = makeApointToArea(edMICSInSample$Stratum[!edMICSInSample$urban], admFinal$NAME_FINAL)
  }
  # numPerStratRur = rowSums(ARurMICS)
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRur))
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  # XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  # ARurMICS = ARurMICS[,actualIndexRur]
  startStratInds = which(XRur$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XRur)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRur[x])})) # length nRur, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
  nRur = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIRur = allAreaIs + (allIntIs-1)*nAreas
  XRur = XRur[transformIRur,]
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  if(admLevel == "adm2") {
    ARurMICS = ARurMICS[,transformIRur] 
  }
  
  # now do the same for the out of sample data if need be
  if((fold > 10) && !areal) {
    allNumPerStrat = aggregate(edMICSOutOfSample$Stratum, by=list(strat=edMICSOutOfSample$Stratum, urb=edMICSOutOfSample$urban), FUN=length, drop=FALSE)
    if(nrow(allNumPerStrat) != 82) {
      stop("bad number of rows in out of sample nPerStrat")
    }
    numPerStratUrbOutOfSample = allNumPerStrat[allNumPerStrat[,2], 3]
    numPerStratUrbOutOfSample[is.na(numPerStratUrbOutOfSample)] = 0
    numPerStratRurOutOfSample = allNumPerStrat[!allNumPerStrat[,2], 3]
    numPerStratRurOutOfSample[is.na(numPerStratRurOutOfSample)] = 0
    
    # first extract only the relevant covariates for the in sample data
    XUrbOutOfSample = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
    stratUrb = XUrbOutOfSample$strat
    if(admLevel == "adm2") {
      AUrbMICSOutOfSample = makeApointToArea(XUrbOutOfSample$subarea, adm2$NAME_2)
    } else {
      AUrbMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    }
    # numPerStratUrbOutOfSample = rowSums(AUrbMICSOutOfSample)
    # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICSOutOfSample), each=numPerStratUrbOutOfSample * KMICS))
    # obsIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), KMICS)
    # intPtIndexUrb = rep(1:sum(numPerStratUrbOutOfSample), each=KMICS)
    # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrbOutOfSample))
    # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrbOutOfSample), each=rep(numPerStratUrbOutOfSample, times=KMICS)))
    # XUrbOutOfSample = XUrbOutOfSample[actualIndexUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    # XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    # AUrbMICSOutOfSample = AUrbMICSOutOfSample[,actualIndexUrb]
    
    startStratInds = which(XUrbOutOfSample$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
    nAreas = nrow(XUrbOutOfSample)/KMICS
    areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrbOutOfSample[x])})) # length nUrb, range = 1:41. gives area index for each obs
    allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
    nUrb = length(allAreaIs)/KMICS
    allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
    transformIUrb = allAreaIs + (allIntIs-1)*nAreas
    XUrbOutOfSample = XUrbOutOfSample[transformIUrb,] # now XUrbOutOfSample is [K * nObsUrb] x nVar
    XUrbOutOfSample = XUrbOutOfSample[,names(XUrbOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    if(admLevel == "adm2") {
      AUrbMICSOutOfSample = AUrbMICSOutOfSample[,transformIUrb]
    }
    
    XRurOutOfSample = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
    stratRur = XRurOutOfSample$strat
    if(admLevel == "adm2") {
      ARurMICSOutOfSample = makeApointToArea(XRurOutOfSample$subarea, adm2$NAME_2)
    } else {
      ARurMICSOutOfSample = makeApointToArea(edMICSOutOfSample$Stratum[!edMICSOutOfSample$urban], admFinal$NAME_FINAL)
    }
    # numPerStratRurOutOfSample = rowSums(ARurMICSOutOfSample)
    # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICSOutOfSample), each=numPerStratRurOutOfSample * KMICS))
    # obsIndexRur = rep(1:sum(numPerStratRurOutOfSample), KMICS)
    # intPtIndexRur = rep(1:sum(numPerStratRurOutOfSample), each=KMICS)
    # actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRurOutOfSample))
    # actualIndexRur = unlist(mapply(rep, 1:nrow(XRurOutOfSample), each=rep(numPerStratRurOutOfSample, times=KMICS)))
    # XRurOutOfSample = XRurOutOfSample[actualIndexRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
    # XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    # ARurMICSOutOfSample = ARurMICSOutOfSample[,actualIndexRur]
    
    startStratInds = which(XRurOutOfSample$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
    nAreas = nrow(XRurOutOfSample)/KMICS
    areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRurOutOfSample[x])})) # length nRur, range = 1:41. gives area index for each obs
    allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
    nRur = length(allAreaIs)/KMICS
    allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
    transformIRur = allAreaIs + (allIntIs-1)*nAreas
    XRurOutOfSample = XRurOutOfSample[transformIRur,] # now XRurOutOfSample is [K * nObsRur] x nVar
    XRurOutOfSample = XRurOutOfSample[,names(XRurOutOfSample) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
    if(admLevel == "adm2") {
      ARurMICSOutOfSample = ARurMICSOutOfSample[,transformIRur]
    }
  }
  else {
    AUrbMICSOutOfSample = NULL
    XUrbOutOfSample = NULL
    ARurMICSOutOfSample = NULL
    XRurOutOfSample = NULL
  }
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  if((fold > 10) && !areal) {
    # Do the same with the out of sample MICS observations
    wUrbanOutOfSample = intPtsMICS$wUrban
    stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrbanOutOfSample), each=numPerStratUrbOutOfSample))
    wUrbanOutOfSample = wUrbanOutOfSample[stratIndexUrbW,]
    
    wRuralOutOfSample = intPtsMICS$wRural
    stratIndexRurW = unlist(mapply(rep, 1:nrow(wRuralOutOfSample), each=numPerStratRurOutOfSample))
    wRuralOutOfSample = wRuralOutOfSample[stratIndexRurW,]
  } else {
    wUrbanOutOfSample = NULL
    wRuralOutOfSample = NULL
  }
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICSInSample$Stratum, admFinal$NAME_FINAL)
  # edMICSInSample = edMICSInSample[order(stratIDs),]
  
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
  if((fold > 10) && !areal) {
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
  
  if(!areal) {
    # set up out of sample covariates and weights
    intPtsMICS$XUrbOutOfSample = XUrbOutOfSample[,-(2:3)] # don't include strata or intercept
    intPtsMICS$XRurOutOfSample = XRurOutOfSample[,-(2:3)]
    if(fold > 10) {
      intPtsMICS$XUrbOutOfSample = as.matrix(intPtsMICS$XUrbOutOfSample)
      intPtsMICS$XRurOutOfSample = as.matrix(intPtsMICS$XRurOutOfSample)
    }
    
    intPtsMICS$wUrbanOutOfSample = wUrbanOutOfSample
    intPtsMICS$wRuralOutOfSample = wRuralOutOfSample
    
    intPtsDHS$wUrbanOutOfSample = intPtsDHS$wUrban[outOfSampleLIndsUrb2DHS,]
    intPtsDHS$wRuralOutOfSample = intPtsDHS$wRural[outOfSampleLIndsRur2DHS,]
  } else {
    intPtsMICS$XUrbOutOfSample = NULL
    intPtsMICS$XRurOutOfSample = NULL
    intPtsMICS$wUrbanOutOfSample = NULL
    intPtsMICS$wRuralOutOfSample = NULL
  }
  
  intPtsDHS$wUrban = intPtsDHS$wUrban[inSampleLIndsUrb2DHS,]
  intPtsDHS$wRural = intPtsDHS$wRural[inSampleLIndsRur2DHS,]
  
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  
  # set priors ----
  alpha_pri = c(0, 10^2)
  beta_pri = c(0, 5^2)
  
  if(admLevel == "admFinal") {
    out = load("savedOutput/global/admFinalMat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  } else {
    out = load("savedOutput/global/adm2Mat.RData")
    bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
  }
  lambdaTau = getLambdaPCprec(u=1, alpha=.1)
  lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision
  
  # Specify inputs for TMB ----
  
  ## specify random effects
  rand_effs <- c('Epsilon_bym2', 'nuggetUrbMICS', 'nuggetRurMICS', 
                 'nuggetUrbDHS', 'nuggetRurDHS')
  
  # make sure A matrices are sparse if BYM2 is at the admin2 level
  if(admLevel == "adm2") {
    # in sample
    AUrbDHS=as(AUrbDHS, "sparseMatrix")
    ARurDHS=as(ARurDHS, "sparseMatrix")
    AUrbMICS=as(AUrbMICS, "sparseMatrix")
    ARurMICS=as(ARurMICS, "sparseMatrix")
    
    # out of sample
    if(!areal) {
      AUrbDHSoutOfSample=as(AUrbDHSoutOfSample, "sparseMatrix")
      ARurDHSoutOfSample=as(ARurDHSoutOfSample, "sparseMatrix")
      if(fold > 10) {
        AUrbMICSOutOfSample=as(AUrbMICSOutOfSample, "sparseMatrix")
        ARurMICSOutOfSample=as(ARurMICSOutOfSample, "sparseMatrix")
      }
    }
  }
  
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
  
  DLL = ifelse(admLevel == "admFinal", 'modBYM2JitterFusionNugget', 'modBYM2JitterFusionNugget2')
  
  list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
       edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
       MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
                            hessian=TRUE, DLL=DLL), 
       dataOutOfSample=dataOutOfSample)
}

# This function generates and saves all the validation datasets and input 
# parameters for all models
getAllValidationData = function(folds=1:20, res=100, adm2AsCovariate=TRUE) {
  
  # first generate M_d data
  print("generating data for M_d...")
  time1 = system.time(datMd <- lapply(folds, getValidationDataM_d, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time1, " seconds. Now generating data for M_D..."))
  time2 = system.time(datMD <- lapply(folds, getValidationDataM_D, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time2, " seconds. Now generating data for M_dm..."))
  time3 = system.time(datMdm <- lapply(folds, getValidationDataM_dm, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time3, " seconds. Now generating data for M_DM..."))
  time4 = system.time(datMDM <- lapply(folds, getValidationDataM_DM, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time4, " seconds. Now saving results..."))
  
  save(datMd, file="savedOutput/validation/datMd.RData")
  save(datMD, file="savedOutput/validation/datM_D.RData")
  save(datMdm, file="savedOutput/validation/datMdm.RData")
  save(datMDM, file="savedOutput/validation/datM_DM.RData")
  
  invisible(list(datMd, datMD, datMdm, datMDM))
}

# This function generates and saves all the validation datasets and input 
# parameters for all models
getAllValidationDataAreal = function(folds=1:37, res=100, adm2AsCovariate=TRUE) {
  
  # first generate M_d data
  print("generating data for M_d...")
  time1 = system.time(datMdareal <- lapply(folds, getValidationDataM_d, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate, areal=TRUE))[3]
  print(paste0("Took ", time1, " seconds. Now generating data for M_D..."))
  time2 = system.time(datMDareal <- lapply(folds, getValidationDataM_D, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate, areal=TRUE))[3]
  print(paste0("Took ", time2, " seconds. Now generating data for M_dm..."))
  time3 = system.time(datMdmareal <- lapply(folds, getValidationDataM_dm, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate, areal=TRUE))[3]
  print(paste0("Took ", time3, " seconds. Now generating data for M_DM..."))
  time4 = system.time(datMDMareal <- lapply(folds, getValidationDataM_DM, admLevel="admFinal", res=res, adm2AsCovariate=adm2AsCovariate, areal=TRUE))[3]
  print(paste0("Took ", time4, " seconds. Now saving results..."))
  
  save(datMdareal, file="savedOutput/validation/datMdareal.RData")
  save(datMDareal, file="savedOutput/validation/datM_Dareal.RData")
  save(datMdmareal, file="savedOutput/validation/datMdmareal.RData")
  save(datMDMareal, file="savedOutput/validation/datM_DMareal.RData")
  
  invisible(list(datMdareal, datMDareal, datMdmareal, datMDMareal))
}

# This function generates and saves all the validation datasets and input 
# parameters for all models
getAllValidationData2 = function(folds=1:20, res=300, adm2AsCovariate=TRUE) {
  
  # first generate M_d data
  print("generating data for M_d2...")
  time1 = system.time(datMd2 <- lapply(folds, getValidationDataM_d, admLevel="adm2", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time1, " seconds. Now generating data for M_D2..."))
  time2 = system.time(datMD2 <- lapply(folds, getValidationDataM_D, admLevel="adm2", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time2, " seconds. Now generating data for M_dm2..."))
  time3 = system.time(datMdm2 <- lapply(folds, getValidationDataM_dm, admLevel="adm2", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time3, " seconds. Now generating data for M_DM2..."))
  time4 = system.time(datMDM2 <- lapply(folds, getValidationDataM_DM, admLevel="adm2", res=res, adm2AsCovariate=adm2AsCovariate))[3]
  print(paste0("Took ", time4, " seconds. Now saving results..."))
  
  # intParText = ifelse(res == 25, "", paste0("_", res, ""))
  # intParText = ifelse(adm2AsCovariate, paste0(intParText, "_adm2Cov"), intParText)
  save(datMd2, file=paste0("savedOutput/validation/datMd2.RData"))
  save(datMD2, file=paste0("savedOutput/validation/datM_D2.RData"))
  save(datMdm2, file=paste0("savedOutput/validation/datMdm2.RData"))
  save(datMDM2, file=paste0("savedOutput/validation/datM_DM2.RData"))
  
  invisible(list(datMd2, datMD2, datMdm2, datMDM2))
}

# This function generates and saves all the areal validation datasets and input 
# parameters for all models, 1 for each of the 37 admin1 States in Nigeria
getAllValidationData2Areal = function(folds=1:37, res=300) {
  
  # first generate M_d data
  print("generating data for M_d2...")
  time1 = system.time(datMd2areal <- lapply(folds, getValidationDataM_d, res=res, admLevel="adm2", areal=TRUE))[3]
  print(paste0("Took ", time1, " seconds. Now generating data for M_D2..."))
  time2 = system.time(datMD2areal <- lapply(folds, getValidationDataM_D, res=res, admLevel="adm2", areal=TRUE))[3]
  print(paste0("Took ", time2, " seconds. Now generating data for M_dm2..."))
  time3 = system.time(datMdm2areal <- lapply(folds, getValidationDataM_dm, res=res, admLevel="adm2", areal=TRUE))[3]
  print(paste0("Took ", time3, " seconds. Now generating data for M_DM2..."))
  time4 = system.time(datMDM2areal <- lapply(folds, getValidationDataM_DM, res=res, admLevel="adm2", areal=TRUE))[3]
  print(paste0("Took ", time4, " seconds. Now saving results..."))
  
  save(datMd2areal, file="savedOutput/validation/datMd2areal.RData")
  save(datMD2areal, file="savedOutput/validation/datM_D2areal.RData")
  save(datMdm2areal, file="savedOutput/validation/datMdm2areal.RData")
  save(datMDM2areal, file="savedOutput/validation/datM_DM2areal.RData")
  
  invisible(list(datMd2areal, datMD2areal, datMdm2areal, datMDM2areal))
}

getValidationFit = function(fold, 
                            model=c("Md", "MD", "Mdm", "MDM", "Md2", "MD2", "Mdm2", "MDM2"), 
                            regenModFit=FALSE, regenPreds=TRUE, randomBeta=FALSE, randomAlpha=FALSE, 
                            fromOptPar=FALSE, areal=FALSE, nsim=10000, sep=TRUE, forceRegenIfMissing=TRUE, 
                            varClust=FALSE) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  admLevel = ifelse(model %in% c("Md", "MD", "Mdm", "MDM"), 1, 2)
  
  out = load("savedOutput/validation/edVal.RData")
  areas = sort(unique(edVal$area))
  foldArea = areas[fold]
  
  # load in the data for the appropriate model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  } else if(fnameRoot == "MD2") {
    fnameRoot = "M_D2"
  } else if(fnameRoot == "MDM2") {
    fnameRoot = "M_DM2"
  }
  
  if(areal) {
    fnameRoot = paste0(fnameRoot, "areal", collapse="")
  }
  # fnameRoot = paste0(fnameRoot, "_", res, collapse="")
  
  fname = paste0("savedOutput/validation/dat", fnameRoot, ".RData")
  out = load(fname)
  
  if(varClust) {
    fnameRoot = paste0(fnameRoot, "VarClust", collapse="")
  }
  
  # get the data from the appropriate variable. The data is in the following format:
  # list(edInSample=edInSample, edOutOfSample=edOutOfSample, 
  #      edMICSInSample=edMICSInSample, edMICSOutOfSample=edMICSOutOfSample, 
  #      MakeADFunInputs=list(data=data_full, parameters=tmb_params, random=rand_effs, 
  #                           hessian=TRUE, DLL='modBYM2JitterFusionNugget'))
  varname = paste0("dat", model, collapse="")
  if(areal) {
    varname = paste0(varname, "areal", collapse = "")
  }
  dat = get(varname)[[fold]]
  
  edInSample = dat$edInSample
  edMICSInSample = dat$edMICSInSample
  edOutOfSample = dat$edOutOfSample
  edMICSOutOfSample = dat$edMICSOutOfSample
  
  # specify random effects
  if(!sep) {
    rand_effsStart <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS', 'beta', 'alpha')
    rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS', 'beta', 'alpha')
    if(randomBeta) {
      rand_effsStart = c("beta", rand_effsStart)
    }
    if(randomAlpha) {
      rand_effsStart = c("alpha", rand_effsStart)
    }
  } else {
    rand_effsStart <- c('w_bym2Star', 'u_bym2Star', 'nuggetUrbDHS', 'nuggetRurDHS')
    if(randomBeta) {
      rand_effsStart = c("beta", rand_effsStart)
    }
    if(randomAlpha) {
      rand_effsStart = c("alpha", rand_effsStart)
    }
    
    # in this case, must also change random effects input into main TMB call...
    if(model %in% c("Md", "MD", "Md2", "MD2")) {
      dat$MakeADFunInputs$random = c('w_bym2Star', 'u_bym2Star', 'nuggetUrbDHS', 'nuggetRurDHS')
    } else {
      dat$MakeADFunInputs$random = c('w_bym2Star', 'u_bym2Star', 'nuggetUrbMICS', 'nuggetRurMICS', 'nuggetUrbDHS', 'nuggetRurDHS')
    }
    if(randomBeta) {
      dat$MakeADFunInputs$random = c("beta", dat$MakeADFunInputs$random)
    }
    if(randomAlpha) {
      dat$MakeADFunInputs$random = c("alpha", dat$MakeADFunInputs$random)
    }
    
    if(varClust) {
      # we have separate cluster level variances. Adjust parameters accordingly
      if(model %in% c("Md", "MD")) {
        tauEpsPar = list(log_tauEpsUrb=0, log_tauEpsRur=0)
      } else if(model %in% c("Mdm", "MDM")) {
        tauEpsPar = list(log_tauEpsUMICS=0, log_tauEpsRMICS=0, 
                         log_tauEpsUDHS=0, log_tauEpsRDHS=0)
      }
      whichI = which(names(dat$MakeADFunInputs$parameters) == "log_tauEps")
      newPar = c(dat$MakeADFunInputs$parameters[1:(whichI-1)], 
                 tauEpsPar, 
                 dat$MakeADFunInputs$parameters[(whichI+1):length(dat$MakeADFunInputs$parameters)])
      dat$MakeADFunInputs$parameters = newPar
    }
  }
  
  # if the file doesn't exist, recreate it.
  if(forceRegenIfMissing) {
    regenModFit = regenModFit || !file.exists(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
    regenPreds = regenPreds || !file.exists(paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  }
  
  # initialize with simple/unadjusted model ----
  if(!varClust) {
    optParStart = c(0, 0, 0)
  } else if(model %in% c("Md", "MD")) {
    optParStart = c(0, 0, 0, 0)
  } else if(model %in% c("Mdm", "MDM")) {
    optParStart = c(0, 0, 0, 0)
  }
  
  if(((model == "Md2") || (model == "Md")) && regenModFit) {
    # now set the initial parameters
    print("Initializing optimization for the unadjusted DHS model")
    initUrbP = sum(c(edInSample$y[edInSample$urban]))/sum(c(edInSample$n[edInSample$urban]))
    initRurP = sum(c(edInSample$y[!edInSample$urban]))/sum(c(edInSample$n[!edInSample$urban]))
    initAlpha = logit(initRurP)
    initBeta1 = logit(initUrbP) - initAlpha
    
    if(!sep) {
      tmb_params <- list(alpha = initAlpha, # intercept
                         beta = c(initBeta1, rep(0, ncol(dat$MakeADFunInputs$data$X_betaUrbanDHS)-1)), 
                         log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                         logit_phi = 0, # SPDE parameter related to the range
                         log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                         Epsilon_bym2 = rep(0, ncol(dat$MakeADFunInputs$data$Q_bym2)), # RE on mesh vertices
                         nuggetUrbDHS = rep(0, sum(edInSample$urban)), 
                         nuggetRurDHS = rep(0, sum(!edInSample$urban))
      )
    } else {
      tmb_params <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                         logit_phi = 0, # SPDE parameter related to the range
                         log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                         alpha = initAlpha, # intercept
                         beta = c(initBeta1, rep(0, ncol(dat$MakeADFunInputs$data$X_betaUrbanDHS)-1)), 
                         w_bym2Star = rep(0, ncol(dat$MakeADFunInputs$data$Q)), # RE on mesh vertices
                         u_bym2Star = rep(0, ncol(dat$MakeADFunInputs$data$Q)), # RE on mesh vertices
                         nuggetUrbDHS = rep(0, sum(edInSample$urban)), 
                         nuggetRurDHS = rep(0, sum(!edInSample$urban))
      )
      
      areaidxlocUrban = apply(dat$MakeADFunInputs$data$AprojUrbanDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
      areaidxlocRural = apply(dat$MakeADFunInputs$data$AprojRuralDHS, 1, function(x) {match(1, x)}) - 1
      areaidxlocUrban = as.integer(areaidxlocUrban)
      areaidxlocRural = as.integer(areaidxlocRural)
      
      # conditioning by Kriging from Eq (2.30) in Rue Held:
      # Ax = e (for A = (0^T 1^T), e = 0), x = (w^T u^T)^T
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x - e)
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x)
      # x* = x - (Q_x^-1 A^T A x) / sum(Q^+)
      # x* = x - (sqrt(phi/tau) Q_{+:}^+ \\ Q_{+:}^+) * sum(u) / sum(Q^+)
      # for Q_{+:}^+ = rowSums(Q^+), where * denotes the constrained version of the effect
      # Hence, we need Q_{+:}^+ / sum(Q^+):
      Qinv = dat$MakeADFunInputs$data$V %*% dat$MakeADFunInputs$data$Q %*% t(dat$MakeADFunInputs$data$V)
      QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
      
      # make sure prior agrees with INLA
      beta_pri = c(0, sqrt(1000))
      
      data_start = list(
        y_iUrbanDHS=dat$MakeADFunInputs$data$y_iUrbanDHS, # same as above but for DHS survey
        y_iRuralDHS=dat$MakeADFunInputs$data$y_iRuralDHS, # 
        n_iUrbanDHS=dat$MakeADFunInputs$data$n_iUrbanDHS, # number binomial trials
        n_iRuralDHS=dat$MakeADFunInputs$data$n_iRuralDHS, # 
        # AprojUrbanDHS=dat$MakeADFunInputs$data$AprojUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
        # AprojRuralDHS=dat$MakeADFunInputs$data$AprojRuralDHS, # 
        areaidxlocUrban = areaidxlocUrban, 
        areaidxlocRural = areaidxlocRural, 
        X_betaUrbanDHS=dat$MakeADFunInputs$data$X_betaUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
        X_betaRuralDHS=dat$MakeADFunInputs$data$X_betaRuralDHS, # 
        wUrbanDHS=dat$MakeADFunInputs$data$wUrbanDHS, # nObsUrban x nIntegrationPointsUrban weight matrix
        wRuralDHS=dat$MakeADFunInputs$data$wRuralDHS, # 
        
        # V_bym2=dat$MakeADFunInputs$data$V_bym2, # eigenvectors of Q (i.e. Q = V Lambda V^T)
        Q_bym2=dat$MakeADFunInputs$data$Q_bym2, # BYM2 unit scaled structure matrix
        alpha_pri=dat$MakeADFunInputs$data$alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
        beta_pri=dat$MakeADFunInputs$data$beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
        tr=dat$MakeADFunInputs$data$tr, # precomputed for Q_bym2
        gammaTildesm1=dat$MakeADFunInputs$data$gammaTildesm1, # precomputed for Q_bym2
        QinvSumsNorm=QinvSumsNorm, 
        lambdaPhi=dat$MakeADFunInputs$data$lambdaPhi, # precomputed for Q_bym2
        lambdaTau=dat$MakeADFunInputs$data$lambdaTau, # determines PC prior for tau
        lambdaTauEps=dat$MakeADFunInputs$data$lambdaTauEps, 
        options=0 # 1 for adreport of log tau and logit phi
      )
    }
    
    if(varClust) {
      # we have separate cluster level variances. Adjust parameters accordingly
      tauEpsPar = list(log_tauEpsUrb=0, log_tauEpsRur=0)
      whichI = which(names(tmb_params) == "log_tauEps")
      newPar = c(tmb_params[1:(whichI-1)], 
                 tauEpsPar, 
                 tmb_params[(whichI+1):length(tmb_params)])
      tmb_params = newPar
    }
  } else if((model %in% c("MD", "Mdm", "MDM", "MD2", "Mdm2", "MDM2")) && regenModFit && !fromOptPar) {
    print("Initializing optimization via the unadjusted DHS model")
    initUrbP = sum(c(edInSample$y[edInSample$urban]))/sum(c(edInSample$n[edInSample$urban]))
    initRurP = sum(c(edInSample$y[!edInSample$urban]))/sum(c(edInSample$n[!edInSample$urban]))
    initAlpha = logit(initRurP)
    initBeta1 = logit(initUrbP) - initAlpha
    
    if(!sep) {
      tmb_paramsStart <- list(alpha = initAlpha, # intercept
                              beta = c(initBeta1, rep(0, ncol(dat$MakeADFunInputs$data$X_betaUrbanDHS)-1)), 
                              log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              logit_phi = 0, # SPDE parameter related to the range
                              log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              Epsilon_bym2 = rep(0, ncol(dat$MakeADFunInputs$data$Q_bym2)), # RE on mesh vertices
                              nuggetUrbDHS = rep(0, sum(edInSample$urban)), 
                              nuggetRurDHS = rep(0, sum(!edInSample$urban))
      )
    } else {
      tmb_paramsStart <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              logit_phi = 0, # SPDE parameter related to the range
                              log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                              alpha = initAlpha, # intercept
                              beta = c(initBeta1, rep(0, ncol(dat$MakeADFunInputs$data$X_betaUrbanDHS)-1)), 
                              w_bym2Star = rep(0, ncol(dat$MakeADFunInputs$data$Q)), # RE on mesh vertices
                              u_bym2Star = rep(0, ncol(dat$MakeADFunInputs$data$Q)), # RE on mesh vertices
                              nuggetUrbDHS = rep(0, sum(edInSample$urban)), 
                              nuggetRurDHS = rep(0, sum(!edInSample$urban))
      )
    }
    
    if(varClust) {
      # we have separate cluster level variances. Adjust parameters accordingly
      tauEpsPar = list(log_tauEpsUrb=0, log_tauEpsRur=0)
      whichI = which(names(tmb_paramsStart) == "log_tauEps")
      newPar = c(tmb_paramsStart[1:(whichI-1)], 
                 tauEpsPar, 
                 tmb_paramsStart[(whichI+1):length(tmb_paramsStart)])
      tmb_paramsStart = newPar
    }
    
    # collect input data, setting only first weights as nonzero (to 1)
    wUrbanDHStemp=dat$MakeADFunInputs$data$wUrbanDHS
    wRuralDHStemp=dat$MakeADFunInputs$data$wRuralDHS
    wUrbanDHStemp[,1] = 1
    wRuralDHStemp[,1] = 1
    wUrbanDHStemp[,-1] = 0
    wRuralDHStemp[,-1] = 0
    
    if(!sep) {
      data_start = list(
        y_iUrbanDHS=dat$MakeADFunInputs$data$y_iUrbanDHS, # same as above but for DHS survey
        y_iRuralDHS=dat$MakeADFunInputs$data$y_iRuralDHS, # 
        n_iUrbanDHS=dat$MakeADFunInputs$data$n_iUrbanDHS, # number binomial trials
        n_iRuralDHS=dat$MakeADFunInputs$data$n_iRuralDHS, # 
        AprojUrbanDHS=dat$MakeADFunInputs$data$AprojUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
        AprojRuralDHS=dat$MakeADFunInputs$data$AprojRuralDHS, # 
        X_betaUrbanDHS=dat$MakeADFunInputs$data$X_betaUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
        X_betaRuralDHS=dat$MakeADFunInputs$data$X_betaRuralDHS, # 
        wUrbanDHS=dat$MakeADFunInputs$data$wUrbanDHS, # nObsUrban x nIntegrationPointsUrban weight matrix
        wRuralDHS=dat$MakeADFunInputs$data$wRuralDHS, # 
        
        V_bym2=dat$MakeADFunInputs$data$V_bym2, # eigenvectors of Q (i.e. Q = V Lambda V^T)
        Q_bym2=dat$MakeADFunInputs$data$Q_bym2, # BYM2 unit scaled structure matrix
        alpha_pri=dat$MakeADFunInputs$data$alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
        beta_pri=dat$MakeADFunInputs$data$beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
        tr=dat$MakeADFunInputs$data$tr, # precomputed for Q_bym2
        gammaTildesm1=dat$MakeADFunInputs$data$gammaTildesm1, # precomputed for Q_bym2
        lambdaPhi=dat$MakeADFunInputs$data$lambdaPhi, # precomputed for Q_bym2
        lambdaTau=dat$MakeADFunInputs$data$lambdaTau, # determines PC prior for tau
        lambdaTauEps=dat$MakeADFunInputs$data$lambdaTauEps, 
        options=0 # 1 for adreport of log tau and logit phi
      )
    } else {
      areaidxlocUrban = apply(dat$MakeADFunInputs$data$AprojUrbanDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
      areaidxlocRural = apply(dat$MakeADFunInputs$data$AprojRuralDHS, 1, function(x) {match(1, x)}) - 1
      areaidxlocUrban = as.integer(areaidxlocUrban)
      areaidxlocRural = as.integer(areaidxlocRural)
      
      # conditioning by Kriging from Eq (2.30) in Rue Held:
      # Ax = e (for A = (0^T 1^T), e = 0), x = (w^T u^T)^T
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x - e)
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x)
      # x* = x - (Q_x^-1 A^T A x) / sum(Q^+)
      # x* = x - (sqrt(phi/tau) Q_{+:}^+ \\ Q_{+:}^+) * sum(u) / sum(Q^+)
      # for Q_{+:}^+ = rowSums(Q^+), where * denotes the constrained version of the effect
      # Hence, we need Q_{+:}^+ / sum(Q^+):
      Qinv = dat$MakeADFunInputs$data$V %*% dat$MakeADFunInputs$data$Q %*% t(dat$MakeADFunInputs$data$V)
      QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
      
      # make sure prior agrees with INLA
      beta_pri = c(0, sqrt(1000))
      
      data_start = list(
        y_iUrbanDHS=dat$MakeADFunInputs$data$y_iUrbanDHS, # same as above but for DHS survey
        y_iRuralDHS=dat$MakeADFunInputs$data$y_iRuralDHS, # 
        n_iUrbanDHS=dat$MakeADFunInputs$data$n_iUrbanDHS, # number binomial trials
        n_iRuralDHS=dat$MakeADFunInputs$data$n_iRuralDHS, # 
        # AprojUrbanDHS=dat$MakeADFunInputs$data$AprojUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
        # AprojRuralDHS=dat$MakeADFunInputs$data$AprojRuralDHS, # 
        areaidxlocUrban = areaidxlocUrban, 
        areaidxlocRural = areaidxlocRural, 
        X_betaUrbanDHS=dat$MakeADFunInputs$data$X_betaUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
        X_betaRuralDHS=dat$MakeADFunInputs$data$X_betaRuralDHS, # 
        wUrbanDHS=wUrbanDHStemp, # nObsUrban x nIntegrationPointsUrban weight matrix
        wRuralDHS=wRuralDHStemp, # 
        
        # V_bym2=dat$MakeADFunInputs$data$V_bym2, # eigenvectors of Q (i.e. Q = V Lambda V^T)
        Q_bym2=dat$MakeADFunInputs$data$Q_bym2, # BYM2 unit scaled structure matrix
        alpha_pri=dat$MakeADFunInputs$data$alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
        beta_pri=dat$MakeADFunInputs$data$beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
        tr=dat$MakeADFunInputs$data$tr, # precomputed for Q_bym2
        gammaTildesm1=dat$MakeADFunInputs$data$gammaTildesm1, # precomputed for Q_bym2
        QinvSumsNorm=QinvSumsNorm, 
        lambdaPhi=dat$MakeADFunInputs$data$lambdaPhi, # precomputed for Q_bym2
        lambdaTau=dat$MakeADFunInputs$data$lambdaTau, # determines PC prior for tau
        lambdaTauEps=dat$MakeADFunInputs$data$lambdaTauEps, 
        options=0 # 1 for adreport of log tau and logit phi
      )
    }
    
    if(!sep) {
      if(admLevel == 2) {
        dyn.load( dynlib("code/modBYM2JitterDHS2"))
        TMB::config(tmbad.sparse_hessian_compress = 1)
        objStart <- MakeADFun(data=data_start,
                              parameters=tmb_paramsStart,
                              random=rand_effsStart,
                              hessian=TRUE,
                              DLL='modBYM2JitterDHS2')
      } else {
        stop("admLevel == 1 and !sep not currently supported")
      }
    } else {
      if(admLevel == 2) {
        dyn.load( dynlib("code/modM_D2Sep"))
        TMB::config(tmbad.sparse_hessian_compress = 1)
        objStart <- MakeADFun(data=data_start,
                              parameters=tmb_paramsStart,
                              random=rand_effsStart,
                              hessian=TRUE,
                              DLL='modM_D2Sep')
      } else {
        if(!varClust) {
          dyn.load(dynlib("code/modM_DSep"))
          TMB::config(tmbad.sparse_hessian_compress = 1)
          objStart <- MakeADFun(data=data_start,
                                parameters=tmb_paramsStart,
                                random=rand_effsStart,
                                hessian=TRUE,
                                DLL='modM_DSep')
        } else {
          dyn.load(dynlib("code/modM_DSepURClust"))
          TMB::config(tmbad.sparse_hessian_compress = 1)
          objStart <- MakeADFun(data=data_start,
                                parameters=tmb_paramsStart,
                                random=rand_effsStart,
                                hessian=TRUE,
                                DLL='modM_DSepURClust')
        }
        
      }
    }
    
    
    lower = rep(-10, length(objStart[['par']]))
    upper = rep( 10, length(objStart[['par']]))
    
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
    
    testObj = objStart
    optPar = testObj$par
    testObj = objStart
    testObj$env$inner.control = list(maxit=1000, tol10=1e-06)
    testObj$env$tracepar = TRUE
    optStart <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                      method = c("BFGS"), hessian = FALSE, control=list(reltol=1e-06))
    optParStart = optStart$par
    print("Optimization/model fitting for the unadjusted DHS model complete")
    
    # make sure last.par is at the optimum
    if(!varClust) {
      parI = match(c("log_tau", "logit_phi", "log_tauEps"), names(testObj$env$last.par))
    } else {
      parI = match(c("log_tau", "logit_phi", "log_tauEpsUrb", "log_tauEpsRur"), names(testObj$env$last.par))
    }
    
    if(!all(testObj$env$last.par[parI] == optParStart)) {
      # last.par is not the optimum. Check last.par.best
      
      if(!all(testObj$env$last.par.best[parI] == optParStart)) {
        # last.par.best is not the optimum.
        # must run function one last time at the optimum so last.par is correct
        invisible(funWrapper(optParStart))
      } else {
        # last.par.best is at the optimum. set last.par to be last.par.best
        testObj$env$last.par = testObj$env$last.par.best
      }
    }
    
    # now set the initial parameters based on the previous optimization
    if(model %in% c("MD", "MD2")) {
      if(!sep) {
        tmb_params <- list(alpha = testObj$env$last.par[grepl("alpha", names(testObj$env$last.par))], # intercept
                           beta = testObj$env$last.par[grepl("beta", names(testObj$env$last.par))], 
                           log_tau = testObj$env$last.par[names(testObj$env$last.par) == "log_tau"], # Log tau (i.e. log spatial precision, Epsilon)
                           logit_phi = testObj$env$last.par[grepl("logit_phi", names(testObj$env$last.par))], # SPDE parameter related to the range
                           log_tauEps = testObj$env$last.par[grepl("log_tauEps", names(testObj$env$last.par))], # Log tau (i.e. log spatial precision, Epsilon)
                           Epsilon_bym2 = testObj$env$last.par[grepl("Epsilon_bym2", names(testObj$env$last.par))], # RE on mesh vertices
                           nuggetUrbDHS = testObj$env$last.par[grepl("nuggetUrbDHS", names(testObj$env$last.par))], 
                           nuggetRurDHS = testObj$env$last.par[grepl("nuggetRurDHS", names(testObj$env$last.par))]
        )
      } else {
        tmb_params <- list(log_tau = testObj$env$last.par[names(testObj$env$last.par) == "log_tau"], # Log tau (i.e. log spatial precision, Epsilon)
                           logit_phi = testObj$env$last.par[grepl("logit_phi", names(testObj$env$last.par))], # SPDE parameter related to the range
                           log_tauEps = testObj$env$last.par[grepl("log_tauEps", names(testObj$env$last.par))], # Log tau (i.e. log spatial precision, Epsilon)
                           alpha = testObj$env$last.par[grepl("alpha", names(testObj$env$last.par))], # intercept
                           beta = testObj$env$last.par[grepl("beta", names(testObj$env$last.par))], 
                           w_bym2Star = testObj$env$last.par[grepl("w_bym2Star", names(testObj$env$last.par))], # RE on mesh vertices
                           u_bym2Star = testObj$env$last.par[grepl("u_bym2Star", names(testObj$env$last.par))], # RE on mesh vertices
                           nuggetUrbDHS = testObj$env$last.par[grepl("nuggetUrbDHS", names(testObj$env$last.par))], 
                           nuggetRurDHS = testObj$env$last.par[grepl("nuggetRurDHS", names(testObj$env$last.par))]
        )
      }
      
    } else if(model %in% c("Mdm", "MDM", "Mdm2", "MDM2")) {
      if(!sep) {
        tmb_params <- list(alpha = testObj$env$last.par[grepl("alpha", names(testObj$env$last.par))], # intercept
                           beta = testObj$env$last.par[grepl("beta", names(testObj$env$last.par))], 
                           log_tau = testObj$env$last.par[names(testObj$env$last.par) == "log_tau"], # Log tau (i.e. log spatial precision, Epsilon)
                           logit_phi = testObj$env$last.par[grepl("logit_phi", names(testObj$env$last.par))], # SPDE parameter related to the range
                           log_tauEps = testObj$env$last.par[grepl("log_tauEps", names(testObj$env$last.par))], # Log tau (i.e. log spatial precision, Epsilon)
                           Epsilon_bym2 = testObj$env$last.par[grepl("Epsilon_bym2", names(testObj$env$last.par))], # RE on mesh vertices
                           nuggetUrbMICS = rep(0, length(dat$MakeADFunInputs$data$y_iUrbanMICS)), 
                           nuggetRurMICS = rep(0, length(dat$MakeADFunInputs$data$y_iRuralMICS)), 
                           nuggetUrbDHS = testObj$env$last.par[grepl("nuggetUrbDHS", names(testObj$env$last.par))], 
                           nuggetRurDHS = testObj$env$last.par[grepl("nuggetRurDHS", names(testObj$env$last.par))]
        )
      } else {
        tmb_params <- list(log_tau = testObj$env$last.par[names(testObj$env$last.par) == "log_tau"], # Log tau (i.e. log spatial precision, Epsilon)
                           logit_phi = testObj$env$last.par[grepl("logit_phi", names(testObj$env$last.par))], # SPDE parameter related to the range
                           log_tauEps = testObj$env$last.par[grepl("log_tauEps", names(testObj$env$last.par))], # Log tau (i.e. log spatial precision, Epsilon)
                           alpha = testObj$env$last.par[grepl("alpha", names(testObj$env$last.par))], # intercept
                           beta = testObj$env$last.par[grepl("beta", names(testObj$env$last.par))], 
                           w_bym2Star = testObj$env$last.par[grepl("w_bym2Star", names(testObj$env$last.par))], # RE on mesh vertices
                           u_bym2Star = testObj$env$last.par[grepl("u_bym2Star", names(testObj$env$last.par))], # RE on mesh vertices
                           nuggetUrbMICS = rep(0, length(dat$MakeADFunInputs$data$y_iUrbanMICS)), 
                           nuggetRurMICS = rep(0, length(dat$MakeADFunInputs$data$y_iRuralMICS)), 
                           nuggetUrbDHS = testObj$env$last.par[grepl("nuggetUrbDHS", names(testObj$env$last.par))], 
                           nuggetRurDHS = testObj$env$last.par[grepl("nuggetRurDHS", names(testObj$env$last.par))]
        )
      }
    }
    
    if(varClust) {
      # we have separate cluster level variances. Adjust parameters accordingly
      if(model %in% c("Md", "MD")) {
        tauEpsPar = list(log_tauEpsUrb=0, log_tauEpsRur=0)
      } else if(model %in% c("Mdm", "MDM")) {
        tauEpsPar = list(log_tauEpsUMICS=0, log_tauEpsRMICS=0, 
                         log_tauEpsUDHS=0, log_tauEpsRDHS=0)
      }
      whichI = which(names(tmb_params) == "log_tauEps")
      newPar = c(tmb_params[1:(whichI-1)], 
                 tauEpsPar, 
                 tmb_params[(whichI+1):length(tmb_params)])
      tmb_params = newPar
    }
    
    if(admLevel == 2) {
      if(!sep) {
        dyn.unload( dynlib("code/modBYM2JitterDHS2"))
      } else {
        dyn.unload( dynlib("code/modM_D2Sep"))
      }
    } else {
      if(!sep) {
        stop("adm1 level and !sep not supported")
      } else {
        if(!varClust) {
          dyn.unload( dynlib("code/modM_DSep"))
        } else {
          dyn.unload( dynlib("code/modM_DSepURClust"))
        }
      }
    }
    
  } else {
    # we still need inputs to make the AD function
    tmb_params = dat$MakeADFunInputs$parameters
  }
  
  # 
  # dat$MakeADFunInputs$parameters$alpha = -2.19427268 # intercept
  # dat$MakeADFunInputs$parameters$beta = c(0.56250562, 0.01287842, 0.10090683, 0.09207191, 1.25325119)
  # dat$MakeADFunInputs$parameters$log_tau = -0.08209775 # Log tau (i.e. log spatial precision, Epsilon)
  # dat$MakeADFunInputs$parameters$logit_phi = -1.78428461 # SPDE parameter related to the range
  # dat$MakeADFunInputs$parameters$log_tauEps = 0.64672006 # Log tau (i.e. log spatial precision, Epsilon)
  
  dat$MakeADFunInputs$parameters = tmb_params
  MakeADFunInputs = dat$MakeADFunInputs
  # MakeADFunInputsFull = MakeADFunInputs
  # MakeADFunInputsFull$random = NULL
  
  # make sure we set the out of sample data correctly for 'DHS data only' models
  if(is.null(edMICSOutOfSample)) {
    out = load("savedOutput/validation/edMICSval.RData")
    edMICSOutOfSample = edMICSval
  }
  
  # adjust the DLL if we are using the seperated parameterization
  if(sep) {
    if(admLevel == 2) {
      if(model %in% c("Mdm2", "MDM2")) {
        MakeADFunInputs$DLL = "modM_DM2Sep"
      } else {
        MakeADFunInputs$DLL = "modM_D2Sep"
      }
    } else {
      if(!varClust) {
        if(model %in% c("Mdm", "MDM")) {
          MakeADFunInputs$DLL = "modM_DMSep"
        } else {
          MakeADFunInputs$DLL = "modM_DSep"
        }
      } else {
        if(model %in% c("Mdm", "MDM")) {
          MakeADFunInputs$DLL = "modM_DMSepVarClust"
        } else {
          MakeADFunInputs$DLL = "modM_DSepURClust"
        }
      }
    }
    
    dat$MakeADFunInputs$DLL = MakeADFunInputs$DLL
  } else {
    stop("!sep Deprecated")
  }
  
  # regen model fit ----
  if(regenModFit || !file.exists(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))) {
    # now fit the model. First we load DLLs and build the functions then we optimize
    dyn.load(dynlib(paste0("code/", MakeADFunInputs$DLL)))
    
    TMB::config(tmbad.sparse_hessian_compress = 1)
    MakeADFunInputs$parameters = tmb_params
    
    if(sep) {
      # reconstruct data to include indices rather than A matrices. Also add in 
      # QinvRowSums:
      if(model %in% c("Mdm", "MDM", "Mdm2", "MDM2")) {
        areaidxlocUrbanMICS = apply(dat$MakeADFunInputs$data$AprojUrbanMICS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
        areaidxlocRuralMICS = apply(dat$MakeADFunInputs$data$AprojRuralMICS, 1, function(x) {match(1, x)}) - 1
        areaidxlocUrbanMICS = as.integer(areaidxlocUrbanMICS)
        areaidxlocRuralMICS = as.integer(areaidxlocRuralMICS)
        
        areaidxlocUrbanDHS = apply(dat$MakeADFunInputs$data$AprojUrbanDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
        areaidxlocRuralDHS = apply(dat$MakeADFunInputs$data$AprojRuralDHS, 1, function(x) {match(1, x)}) - 1
        areaidxlocUrbanDHS = as.integer(areaidxlocUrbanDHS)
        areaidxlocRuralDHS = as.integer(areaidxlocRuralDHS)
      } else {
        areaidxlocUrban = apply(dat$MakeADFunInputs$data$AprojUrbanDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
        areaidxlocRural = apply(dat$MakeADFunInputs$data$AprojRuralDHS, 1, function(x) {match(1, x)}) - 1
        areaidxlocUrban = as.integer(areaidxlocUrban)
        areaidxlocRural = as.integer(areaidxlocRural)
      }
      
      # conditioning by Kriging from Eq (2.30) in Rue Held:
      # Ax = e (for A = (0^T 1^T), e = 0), x = (w^T u^T)^T
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x - e)
      # x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x)
      # x* = x - (Q_x^-1 A^T A x) / sum(Q^+)
      # x* = x - (sqrt(phi/tau) Q_{+:}^+ \\ Q_{+:}^+) * sum(u) / sum(Q^+)
      # for Q_{+:}^+ = rowSums(Q^+), where * denotes the constrained version of the effect
      # Hence, we need Q_{+:}^+ / sum(Q^+):
      Qinv = dat$MakeADFunInputs$data$V %*% dat$MakeADFunInputs$data$Q %*% t(dat$MakeADFunInputs$data$V)
      QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
      
      # make sure prior agrees with INLA
      beta_pri = c(0, sqrt(1000))
      
      if(model %in% c("Mdm", "MDM", "Mdm2", "MDM2")) {
        newDat = list(
          y_iUrbanMICS=dat$MakeADFunInputs$data$y_iUrbanMICS, # same as above but for MICS survey
          y_iRuralMICS=dat$MakeADFunInputs$data$y_iRuralMICS, # 
          n_iUrbanMICS=dat$MakeADFunInputs$data$n_iUrbanMICS, # number binomial trials
          n_iRuralMICS=dat$MakeADFunInputs$data$n_iRuralMICS, # 
          # AprojUrbanMICS=dat$MakeADFunInputs$data$AprojUrbanMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
          # AprojRuralMICS=dat$MakeADFunInputs$data$AprojRuralMICS, # 
          areaidxlocUrbanMICS = areaidxlocUrbanMICS, 
          areaidxlocRuralMICS = areaidxlocRuralMICS, 
          X_betaUrbanMICS=dat$MakeADFunInputs$data$X_betaUrbanMICS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
          X_betaRuralMICS=dat$MakeADFunInputs$data$X_betaRuralMICS, # 
          wUrbanMICS=dat$MakeADFunInputs$data$wUrbanMICS, # nObsUrban x nIntegrationPointsUrban weight matrix
          wRuralMICS=dat$MakeADFunInputs$data$wRuralMICS, # 
          
          y_iUrbanDHS=dat$MakeADFunInputs$data$y_iUrbanDHS, # same as above but for DHS survey
          y_iRuralDHS=dat$MakeADFunInputs$data$y_iRuralDHS, # 
          n_iUrbanDHS=dat$MakeADFunInputs$data$n_iUrbanDHS, # number binomial trials
          n_iRuralDHS=dat$MakeADFunInputs$data$n_iRuralDHS, # 
          # AprojUrbanDHS=dat$MakeADFunInputs$data$AprojUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
          # AprojRuralDHS=dat$MakeADFunInputs$data$AprojRuralDHS, # 
          areaidxlocUrbanDHS = areaidxlocUrbanDHS, 
          areaidxlocRuralDHS = areaidxlocRuralDHS, 
          X_betaUrbanDHS=dat$MakeADFunInputs$data$X_betaUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
          X_betaRuralDHS=dat$MakeADFunInputs$data$X_betaRuralDHS, # 
          wUrbanDHS=dat$MakeADFunInputs$data$wUrbanDHS, # nObsUrban x nIntegrationPointsUrban weight matrix
          wRuralDHS=dat$MakeADFunInputs$data$wRuralDHS, # 
          
          # V_bym2=dat$MakeADFunInputs$data$V_bym2, # eigenvectors of Q (i.e. Q = V Lambda V^T)
          Q_bym2=dat$MakeADFunInputs$data$Q_bym2, # BYM2 unit scaled structure matrix
          alpha_pri=dat$MakeADFunInputs$data$alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
          beta_pri=dat$MakeADFunInputs$data$beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
          tr=dat$MakeADFunInputs$data$tr, # precomputed for Q_bym2
          gammaTildesm1=dat$MakeADFunInputs$data$gammaTildesm1, # precomputed for Q_bym2
          QinvSumsNorm=QinvSumsNorm, 
          lambdaPhi=dat$MakeADFunInputs$data$lambdaPhi, # precomputed for Q_bym2
          lambdaTau=dat$MakeADFunInputs$data$lambdaTau, # determines PC prior for tau
          lambdaTauEps=dat$MakeADFunInputs$data$lambdaTauEps, 
          options=0 # 1 for adreport of log tau and logit phi
        )
      } else {
        newDat = list(
          y_iUrbanDHS=dat$MakeADFunInputs$data$y_iUrbanDHS, # same as above but for DHS survey
          y_iRuralDHS=dat$MakeADFunInputs$data$y_iRuralDHS, # 
          n_iUrbanDHS=dat$MakeADFunInputs$data$n_iUrbanDHS, # number binomial trials
          n_iRuralDHS=dat$MakeADFunInputs$data$n_iRuralDHS, # 
          # AprojUrbanDHS=dat$MakeADFunInputs$data$AprojUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
          # AprojRuralDHS=dat$MakeADFunInputs$data$AprojRuralDHS, # 
          areaidxlocUrban = areaidxlocUrban, 
          areaidxlocRural = areaidxlocRural, 
          X_betaUrbanDHS=dat$MakeADFunInputs$data$X_betaUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
          X_betaRuralDHS=dat$MakeADFunInputs$data$X_betaRuralDHS, # 
          wUrbanDHS=dat$MakeADFunInputs$data$wUrbanDHS, # nObsUrban x nIntegrationPointsUrban weight matrix
          wRuralDHS=dat$MakeADFunInputs$data$wRuralDHS, # 
          
          # V_bym2=dat$MakeADFunInputs$data$V_bym2, # eigenvectors of Q (i.e. Q = V Lambda V^T)
          Q_bym2=dat$MakeADFunInputs$data$Q_bym2, # BYM2 unit scaled structure matrix
          alpha_pri=dat$MakeADFunInputs$data$alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
          beta_pri=dat$MakeADFunInputs$data$beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
          tr=dat$MakeADFunInputs$data$tr, # precomputed for Q_bym2
          gammaTildesm1=dat$MakeADFunInputs$data$gammaTildesm1, # precomputed for Q_bym2
          QinvSumsNorm=QinvSumsNorm, 
          lambdaPhi=dat$MakeADFunInputs$data$lambdaPhi, # precomputed for Q_bym2
          lambdaTau=dat$MakeADFunInputs$data$lambdaTau, # determines PC prior for tau
          lambdaTauEps=dat$MakeADFunInputs$data$lambdaTauEps, 
          options=0 # 1 for adreport of log tau and logit phi
        )
      }
      MakeADFunInputs$data = newDat
    }
    
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
    
    # Error in getParameterOrder(data, parameters, new.env(), DLL = DLL) : 
    #   Error when reading the variable: 'AprojUrbanDHS'. Please check data and parameters.
    obj <- do.call("MakeADFun", MakeADFunInputs)
    # objFull <- do.call("MakeADFun", MakeADFunInputsFull)
    
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
      if(!varClust) {
        optPar = optParStart
      } else if(model %in% c("Md", "MD")) {
        optPar = optParStart
      } else if(model %in% c("Mdm", "MDM")) {
        # initialize MICS cluster variance at DHS cluster variance
        optPar = c(optParStart[1:2], optParStart[3:4], optParStart[3:4])
      }
      
      if(fromOptPar) {
        if(file.exists(paste0("savedOutput/validation/folds/optPar", fnameRoot, "_fold", fold, ".RData"))) {
          # load optPar (replacing the default value) from file if requested
          load(paste0("savedOutput/validation/folds/optPar", fnameRoot, "_fold", fold, ".RData"))
        }
      }
      
      startTime = proc.time()[3]
      for(thisTol in tolSeq) {
        testObj = obj
        testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
        testObj$env$tracepar = TRUE
        print(paste0("optimizing for tol = ", thisTol, "."))
        opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                      method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
        optPar = opt1$par
        
        save(optPar, opt1, file=paste0("savedOutput/validation/folds/optPar", fnameRoot, "_fold", fold, ".RData"))
        if(!is.null(opt1$message)) {
          print(paste0("error for tol = ", thisTol, ". Message:"))
          print(opt1$message)
          next
        } else {
          # make sure last.par is at the optimum
          parI = match(c("log_tau", "logit_phi", "log_tauEps"), names(obj$env$last.par))
          if(!all(obj$env$last.par[parI] == optPar)) {
            # last.par is not the optimum. Check last.par.best
            
            if(!all(obj$env$last.par.best[parI] == optPar)) {
              # last.par.best is not the optimum.
              # must run function one last time at the optimum so last.par is correct
              invisible(funWrapper(optPar))
            } else {
              # last.par.best is at the optimum. set last.par to be last.par.best
              obj$env$last.par = obj$env$last.par.best
            }
          }
          
          print(paste0("completed optimization for tol = ", thisTol, ""))
          
          ## Get standard errors
          print("getting standard errors...")
          sdTime = system.time({
            SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                 bias.correct = TRUE,
                                 bias.correct.control = list(sd = TRUE))
            
            if(!SD0$pdHess) {
              # try recalculating for fixed parameters numerically
              warning("initial hessian non-PD. Trying another way...")
              Hess = numDeriv::hessian( func=testObj$fn, x=optPar )
              SD0 <- sdreport( testObj, hessian.fixed=Hess,
                               getJointPrecision=TRUE,
                               bias.correct = TRUE,
                               bias.correct.control = list(sd = TRUE) )
            }
            
            if(!SD0$pdHess) {
              # try fixing parameter at estimated value
              stop("non-PD hessian...")
            }
          }
          )[3]
          # SD0
          print(paste0("SE calculations took ", sdTime/60, " minutes"))
          
          if(SD0$pdHess) {
            print("Optimization and PD hess calculation done!")
            break
          }
          else {
            # try some other ways of calculating the hessian:
            print("Hessian not PD. Using empirical Bayes inference.")
            
            
            print("Hessian not PD. Testing other SD calculations...")
            print(SD0)
            SD0testSkipDelta <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                 bias.correct = TRUE,
                                 bias.correct.control = list(sd = TRUE), 
                                 skip.delta.method=TRUE)
            print(SD0testSkipDelta)
            print(paste0("SD0testSkipDelta PD: ", SD0testSkipDelta$pdHess))
            
            save(SD0, SD0testSkipDelta, file=paste0("savedOutput/validation/folds/testSkipDelta", fnameRoot, "_fold", fold, ".RData"))
            SD0testSkipBiasCorrect <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                              bias.correct = FALSE, 
                                              skip.delta.method=TRUE)
            save(SD0testSkipBiasCorrect, file=paste0("savedOutput/validation/folds/testSkipBias", fnameRoot, "_fold", fold, ".RData"))
            print(SD0testSkipBiasCorrect)
            print(paste0("SD0testSkipBiasCorrect PD: ", SD0testSkipBiasCorrect$pdHess))
            SD0testNoSDCorrect <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                                              bias.correct = TRUE,
                                              bias.correct.control = list(sd = FALSE), 
                                              skip.delta.method=TRUE)
            save(SD0testNoSDCorrect, file=paste0("savedOutput/validation/folds/testNoSDCorrect", fnameRoot, "_fold", fold, ".RData"))
            print(SD0testNoSDCorrect)
            print(paste0("SD0testNoSDCorrect PD: ", SD0testNoSDCorrect$pdHess))
            
            hessTest <- numDeriv::jacobian(objFull$gr, )
            save(SD0testNoSDCorrect, file=paste0("savedOutput/validation/folds/testNoSDCorrect", fnameRoot, "_fold", fold, ".RData"))
            print(SD0testNoSDCorrect)
            print(paste0("SD0testNoSDCorrect PD: ", SD0testNoSDCorrect$pdHess))
            
            browser()
            # print("Hessian not PD. Rerunning optimization with stricter tol...")
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
    save(SD0, obj, totalTime, sdTime, hessPD, file=paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  } else {
    out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", fold, ".RData"))
  }
  
  # predict at the left out clusters/areas
  admLevelString = ifelse(admLevel == 1, "stratMICS", "adm2")
  if(regenPreds) {
    if(hessPD) {
      if(!areal) {
        preds = predClusters(nsim=nsim, fold, SD0, obj, 
                             model=model, sep=sep, varClust=varClust, 
                             quantiles=c(0.025, 0.1, 0.9, 0.975))
      } else {
        gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=nsim, admLevel=admLevelString, 
                             predAtArea=foldArea,
                             quantiles=c(0.025, 0.1, 0.9, 0.975), sep=sep)
        preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
        preds$fixedMat = gridPreds$fixedMat
      }
    } else {
      preds = NULL
    }
    
    save(SD0, obj, totalTime, sdTime, hessPD, preds, file=paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  } else {
    out = load(paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  }
  
  
  allScores = scoreValidationPreds(fold, model=model, regenScores=TRUE, areal=areal, varClust=varClust)
  
  dyn.unload( dynlib(paste0("code/", dat$MakeADFunInputs$DLL)))
  
  list(SD0, obj, totalTime, sdTime, hessPD, allScores)
}

# make predictions for a set of clusters of one type (MICS or DHS)
# SD0: TMB sdreport object
# obj: TMB fun object
# fold: cross validation fold from 1-20
# ptType: The type of integration points. either MICS or DHS
# model: which model we're making predictions with
predClusters = function(nsim=1000, fold, SD0, obj, 
                        model=c("Md", "MD", "Mdm", "MDM", "Md2", "MD2", "Mdm2", "MDM2"), 
                        quantiles=c(0.025, 0.1, 0.9, 0.975), 
                        addBinVar=TRUE, maxIterChunk=1000, 
                        sep=TRUE, QinvSumsNorm=NULL, verbose=TRUE, varClust=FALSE) {
  
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  admLevel = ifelse(grepl("2", model), "adm2", "stratMICS")
  
  # set the file name root depending on the model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  } else if(fnameRoot == "MD2") {
    fnameRoot = "M_D2"
  } else if(fnameRoot == "MDM2") {
    fnameRoot = "M_DM2"
  }
  
  # for loading the info for predicting at the left out data locations, if the 
  # model doesn't use MICS data use the M_DM info to make the predictions for 
  # curiosity's sake.
  
  if((fold <= 10) && (fnameRoot == "Md")) {
    # for Md model, predict normally for DHS data, but use M_DM prediction info for MICS
    fnameRootLeftOut = "Mdm"
    modelLeftOut = "Mdm"
  } else if(fnameRoot == "Md") {
    fnameRootLeftOut = "M_DM"
    modelLeftOut = "MDM"
  } else if(fnameRoot == "M_D") {
    # for MD model, use M_DM prediction info for both DHS and MICS
    fnameRootLeftOut = "M_DM"
    modelLeftOut = "MDM"
  } else if((fold <= 10) && (fnameRoot == "Md2")) {
    # for Md model, predict normally for DHS data, but use M_DM prediction info for MICS
    fnameRootLeftOut = "Mdm2"
    modelLeftOut = "Mdm2"
  } else if(fnameRoot == "Md2") {
    fnameRootLeftOut = "M_DM2"
    modelLeftOut = "MDM2"
  } else if(fnameRoot == "M_D2") {
    # for MD model, use M_DM prediction info for both DHS and MICS
    fnameRootLeftOut = "M_DM2"
    modelLeftOut = "MDM2"
  } else {
    fnameRootLeftOut = fnameRoot
    modelLeftOut = model
  }
  
  # load relevant data and model fit
  out = load(paste0("savedOutput/validation/dat", fnameRootLeftOut, ".RData", collapse=""))
  leftOutDat = get(paste0("dat", modelLeftOut))[[fold]]$dataOutOfSample
  out = load(paste0("savedOutput/validation/dat", fnameRoot, ".RData", collapse=""))
  load("~/git/jittering/savedOutput/validation/edMICSval.RData")
  load("~/git/jittering/savedOutput/validation/edVal.RData")
  
  varname = paste0("dat", model)
  dat = get(varname)[[fold]]
  
  if(varClust) {
    fnameRoot = paste0(fnameRoot, "VarClust", collapse="")
  }
  
  foldMod = ifelse((fold > 11) && (model %in% c("Md", "MD", "Md2", "MD2")), 11, fold)
  out = load(paste0("savedOutput/validation/folds/fit", fnameRoot, "_fold", foldMod, ".RData"))
  
  # compute QinvSumsNorm if necessary
  if(sep && is.null(QinvSumsNorm)) {
    if(admLevel == "adm2") {
      out = load("savedOutput/global/adm2Mat.RData")
      admMat = adm2Mat
    } else if(admLevel == "stratMICS") {
      out = load("savedOutput/global/admFinalMat.RData")
      admMat = admFinalMat
    }
    
    bym2ArgsTMB = prepareBYM2argumentsForTMB(admMat, u=0.5, alpha=2/3, 
                                             constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
    Qinv = bym2ArgsTMB$V %*% bym2ArgsTMB$Q %*% t(bym2ArgsTMB$V)
    QinvSumsNorm = rowSums(Qinv)/sum(Qinv)
  }
  
  # get predictions in chunks for MICS folds because they are memory intensive to create
  if((fold > 10) && maxIterChunk < nsim) {
    
    startI = 1
    endI = maxIterChunk
    
    while(startI <= nsim) {
      
      # calculate this chunk of predictions
      thisNsim = endI - startI + 1
      out = predClusters(nsim=thisNsim, fold=fold, SD0=SD0, obj=obj, 
                         model=model, 
                         quantiles=quantiles, 
                         addBinVar=addBinVar, maxIterChunk=maxIterChunk, 
                         sep=sep, QinvSumsNorm=QinvSumsNorm, verbose=FALSE,
                         varClust=varClust)
      
      # concatenate results
      if(startI == 1) {
        probDrawsUrb = matrix(nrow=nrow(out$probDrawsUrb), ncol=nsim)
        probDrawsRur = matrix(nrow=nrow(out$probDrawsRur), ncol=nsim)
        probDrawsUrb[,startI:endI] = out$probDrawsUrb
        probDrawsRur[,startI:endI] = out$probDrawsRur
        predsUrb = out$predsUrb
        predsRur = out$predsRur
        fixedMatNames = row.names(out$fixedMat)
        fixedMat = matrix(nrow=nrow(out$fixedMat), ncol=nsim)
        fixedMat[,startI:endI] = out$fixedMat
      } else {
        probDrawsUrb[,startI:endI] = out$probDrawsUrb
        probDrawsRur[,startI:endI] = out$probDrawsRur
        thisWt = thisNsim/endI
        predsUrb = thisWt * out$predsUrb + (1-thisWt) * predsUrb
        predsRur = thisWt * out$predsRur + (1-thisWt) * predsRur
        fixedMat[,startI:endI] = out$fixedMat
      }
      
      # update indices
      startI = endI + 1
      endI = min(c(endI + maxIterChunk, nsim))
    }
    
    # do final postprocessing/get summary statistics
    
    # Make parameter summary tables
    row.names(fixedMat) = fixedMatNames
    parMeans = rowMeans(fixedMat)
    parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
    parSummary = cbind(parMeans, parQuants)
    colnames(parSummary)[1] = "Est"
    colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
    print(xtable(parSummary, digits=2))
    
    # calculate predictive quantiles
    quantsUrb = apply(probDrawsUrb, 1, quantile, probs=quantiles, na.rm=TRUE)
    quantsRur = apply(probDrawsRur, 1, quantile, probs=quantiles, na.rm=TRUE)
    
    # get responses and ns for each cluster
    yUrb = leftOutDat$y_iUrbanMICS
    yRur = leftOutDat$y_iRuralMICS
    nUrb = leftOutDat$n_iUrbanMICS
    nRur = leftOutDat$n_iRuralMICS
    
    return(list(probDrawsUrb=probDrawsUrb, probDrawsRur=probDrawsRur,
                predsUrb=predsUrb, predsRur=predsRur,
                parSummary=parSummary, fixedMat=fixedMat,
                quantsUrb=quantsUrb, quantsRur=quantsRur,
                yUrb=yUrb, yRur=yRur, nUrb=nUrb, nRur=nRur))
  }
  
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
    # parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
    parnames <- colnames(SD0$jointPrecision)
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
    
    # get the spatial effect
    if(!sep) {
      epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_bym2',]
    } else {
      wStar  <- t.draws[parnames == 'w_bym2Star',]
      uStar  <- t.draws[parnames == 'u_bym2Star',] # uStar is unit var scaled
      
      # get how much u reduced by sum to zero constraint to u, then scale
      uSums = colSums(uStar)
      uFacs = uSums * sqrt(phi_tmb_draws*sigmaSq_tmb_draws)
      reduceU = outer(QinvSumsNorm, c(uFacs))
      
      # adjust Epsilon = w for the constraint on u
      epsilon_tmb_draws = wStar - reduceU
    }
    
    hasNugget = any(grepl("log_tauEps", row.names(summary(SD0))))
    if(hasNugget) {
      URclust = "log_tauEpsUrb" %in% parnames
      if(!varClust) {
        sigmaEpsSq_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEps',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSq_tmb_draws)
        row.names(fixedMat)[nrow(fixedMat)] = "sigmaEpsSq"
      } else if(URclust) {
        sigmaEpsSqUrb_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUrb',]), nrow = 1)
        sigmaEpsSqRur_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRur',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSqUrb_tmb_draws, 
                         sigmaEpsSqRur_tmb_draws)
        row.names(fixedMat)[(nrow(fixedMat)-1):nrow(fixedMat)] = c("sigmaEpsSqUrb", "sigmaEpsSqRur")
      } else {
        sigmaEpsSqUDHS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUDHS',]), nrow = 1)
        sigmaEpsSqRDHS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRDHS',]), nrow = 1)
        sigmaEpsSqUMICS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsUMICS',]), nrow = 1)
        sigmaEpsSqRMICS_tmb_draws    <- matrix(1/exp(t.draws[parnames == 'log_tauEpsRMICS',]), nrow = 1)
        fixedMat = rbind(fixedMat, 
                         sigmaEpsSqUMICS_tmb_draws, 
                         sigmaEpsSqRMICS_tmb_draws, 
                         sigmaEpsSqUDHS_tmb_draws, 
                         sigmaEpsSqRDHS_tmb_draws)
        row.names(fixedMat)[(nrow(fixedMat)-3):nrow(fixedMat)] = c("sigmaEpsSqUMICS", "sigmaEpsSqRMICS", 
                                                                   "sigmaEpsSqUDHS", "sigmaEpsSqRDHS")
      }
    }
    
    # Make parameter summary tables
    parMeans = rowMeans(fixedMat)
    parQuants = t(apply(fixedMat, 1, quantile, probs=quantiles))
    parSummary = cbind(parMeans, parQuants)
    colnames(parSummary)[1] = "Est"
    colnames(parSummary)[2:ncol(parSummary)] = paste0("Q", quantiles)
    if(verbose) {
      print(xtable(parSummary, digits=2))
    }
    
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
    Kurb = ncol(wUrb)
    Krur = ncol(wRur)
    if(model %in% c("Md", "MD", "Mdm", "MDM")) {
      Kurb = nrow(Xurb) / nrow(Aurb)
      Krur = nrow(Xrur) / nrow(Arur)
      # bigAurb = matrix(rep(Aurb, times=Kurb), ncol=ncol(Aurb))
      # bigArur = matrix(rep(Arur, times=Krur), ncol=ncol(Arur))
      bigAurb = sapply(1:ncol(Aurb), function(x) {rep(Aurb[,x], times=Kurb)})
      bigArur = sapply(1:ncol(Arur), function(x) {rep(Arur[,x], times=Krur)})
    } else {
      bigAurb = Aurb
      bigArur = Arur
    }
    
    # get latent preds at cluster integration points
    # matMultChunk = function(X1, X2, nColMax=maxIterChunk) {
    #   matMultChunkHelper = function(colStart=1) {
    #     inds = colStart:min(c(colStart + nColMax-1, ncol(X2)))
    #     X1 %*% X2[,inds]
    #   }
    #   
    #   startIs = seq(1, ncol(X2), by=nColMax)
    #   do.call("cbind", lapply(startIs, matMultChunkHelper))
    # }
    # browser()
    clustIntDrawsUrb <- as.matrix(bigAurb %*% epsilon_tmb_draws)
    # clustIntDrawsUrb <- as.matrix(matMultChunk(bigAurb, epsilon_tmb_draws))
    rm(bigAurb)
    clustIntDrawsUrb <- sweep(clustIntDrawsUrb, 2, alpha_tmb_draws, '+')
    clustIntDrawsUrb <- clustIntDrawsUrb + (Xurb %*% beta_tmb_draws)
    
    clustIntDrawsRur <- as.matrix(bigArur %*% epsilon_tmb_draws)
    # clustIntDrawsRur <- as.matrix(matMultChunk(bigArur, epsilon_tmb_draws))
    rm(bigArur)
    clustIntDrawsRur <- sweep(clustIntDrawsRur, 2, alpha_tmb_draws, '+')
    clustIntDrawsRur <- clustIntDrawsRur + (Xrur %*% beta_tmb_draws)
    
    # convert predictions to probability scale
    # browser()
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
      # browser()
      if(!varClust) {
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
      } else if(URclust) {
        logitIntDrawsUrb = sapply(1:ncol(clustIntDrawsUrb), function(colI) {
          thisNuggetSD = sqrt(sigmaEpsSqUrb_tmb_draws[colI])
          nugsUrb = rnorm(nClustUrb, sd=thisNuggetSD)
          clustIntDrawsUrb[,colI] + rep(nugsUrb, times=Kurb)
        })
        probIntDrawsUrb = expit(logitIntDrawsUrb)
        
        logitIntDrawsRur = sapply(1:ncol(clustIntDrawsRur), function(colI) {
          thisNuggetSD = sqrt(sigmaEpsSqRur_tmb_draws[colI])
          nugsRur = rnorm(nClustRur, sd=thisNuggetSD)
          clustIntDrawsRur[,colI] + rep(nugsRur, times=Krur)
        })
        probIntDrawsRur = expit(logitIntDrawsRur)
      } else {
        # set U/R cluster effect variances depending of if predicting DHS or 
        # MICS clusters
        if(fold <= 10) {
          # DHS fold
          sigmaEpsSqUrb_tmb_draws = sigmaEpsSqUDHS_tmb_draws
          sigmaEpsSqRur_tmb_draws = sigmaEpsSqRDHS_tmb_draws
        } else {
          # MICS fold
          sigmaEpsSqUrb_tmb_draws = sigmaEpsSqUMICS_tmb_draws
          sigmaEpsSqRur_tmb_draws = sigmaEpsSqRMICS_tmb_draws
        }
        
        logitIntDrawsUrb = sapply(1:ncol(clustIntDrawsUrb), function(colI) {
          thisNuggetSD = sqrt(sigmaEpsSqUrb_tmb_draws[colI])
          nugsUrb = rnorm(nClustUrb, sd=thisNuggetSD)
          clustIntDrawsUrb[,colI] + rep(nugsUrb, times=Kurb)
        })
        probIntDrawsUrb = expit(logitIntDrawsUrb)
        
        logitIntDrawsRur = sapply(1:ncol(clustIntDrawsRur), function(colI) {
          thisNuggetSD = sqrt(sigmaEpsSqRur_tmb_draws[colI])
          nugsRur = rnorm(nClustRur, sd=thisNuggetSD)
          clustIntDrawsRur[,colI] + rep(nugsRur, times=Krur)
        })
        probIntDrawsRur = expit(logitIntDrawsRur)
      }
      # browser()
    }
    
    # take weighted average of predictions at integration points (i.e. evaluate integral of predictions for each cluster numerically)
    # We will make block diagonal Wurb and Wrur matrices, where element ij is the integration weight for cluster i associated with integration point j
    # buildRowUrb = c(rep(1, Kurb), rep(0, nrow(Xurb)))
    # Wurb = matrix(c(rep(buildRowUrb, times=length(yUrb)-1), rep(1, Kurb)), byrow=TRUE, ncol=nrow(Xurb))
    # Wurb = sweep(Wurb, 2, c(t(wUrb)), FUN="*")
    # browser()
    buildMatUrb = rbind(c(wUrb), 
                        matrix(0, ncol=Kurb*nrow(wUrb), nrow=nrow(wUrb)))
    leaveOutInds = (length(buildMatUrb)-nrow(wUrb) + 1):length(buildMatUrb)
    Wurb = matrix(c(buildMatUrb)[-leaveOutInds], nrow=nrow(wUrb))
    zeroCols = seq(nrow(Wurb)+1, ncol(Wurb), by=nrow(Wurb)+1)
    Wurb = Wurb[,-zeroCols]
    probDrawsUrb = Wurb %*% probIntDrawsUrb
    browser()
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
    
    # sometimes numerical roundoff makes the above sums just barely outside of [0,1]
    probDrawsUrb[probDrawsUrb < 0] = 0
    probDrawsUrb[probDrawsUrb > 1] = 1
    probDrawsRur[probDrawsRur < 0] = 0
    probDrawsRur[probDrawsRur > 1] = 1
    
    # calculate central prediction before binomial variation is added in
    predsUrb = rowMeans(probDrawsUrb)
    predsRur = rowMeans(probDrawsRur)
    browser()
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
  
  # browser()
  if(FALSE) {
    thisNs = c(nUrb, nRur)
    thisTruth = c(yUrb, yRur)/thisNs
    thisEstMat = rbind(probDrawsUrb, probDrawsRur)
    thisWeights = thisNs/sum(thisNs)
    theseScores = getScores(thisTruth, estMat=thisEstMat, weights=thisWeights, 
                            significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, 
                            getAverage=TRUE, na.rm=TRUE)
    #          Bias       Var       MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95
    # 1 -0.04969721 0.1036997 0.1061695 0.3258366 0.1793202        0.819808       0.9952007        1.035166        1.047229
    #   Coverage50 Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95
    # 1  0.3596456  0.6744459  0.8208786  0.8961908 0.4018325 0.6749564 0.7925829 0.8727094
    
    # i4, j11, current
    #          Bias       Var       MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95
    # 1 0.005108033 0.1014577 0.1014838 0.3185652 0.1800558       0.8058633       0.8776547       0.8773777       0.9649815
    #   Coverage50 Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95
    # 1  0.2946404  0.6378262  0.8295828  0.9208956 0.4337488 0.7400739 0.8589104 0.9280471
    
    # i3, j11, current
    #          Bias       Var       MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95
    # 1 0.002739965 0.1344116 0.1344191 0.3666321 0.2086547       0.9349913        1.161867        1.237347         1.16586
    # Coverage50 Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95
    # 1  0.3727854  0.6535637  0.7704153  0.8565431 0.4157941 0.7011344 0.8184991 0.8866274
  }
  
  list(probDrawsUrb=probDrawsUrb, probDrawsRur=probDrawsRur, 
       predsUrb=predsUrb, predsRur=predsRur, 
       parSummary=parSummary, fixedMat=fixedMat, 
       quantsUrb=quantsUrb, quantsRur=quantsRur, 
       yUrb=yUrb, yRur=yRur, nUrb=nUrb, nRur=nRur)
}

# make predictions for a set of clusters of one type (MICS or DHS)
# SD0: TMB sdreport object
# obj: TMB fun object
# fold: cross validation fold from 1-41 represent the MICS stratum in alphabetical order
# model: which model we're making predictions with
predStratum = function(nsim=1000, fold, SD0, obj, 
                       model=c("Md", "MD", "Mdm", "MDM", "Md2", "MD2", "Mdm2", "MDM2"), 
                       quantiles=c(0.025, 0.1, 0.9, 0.975)) {
  
  # clean input arguments
  model = match.arg(model)
  
  # load in the data for the appropriate model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  } else if(fnameRoot == "MD2") {
    fnameRoot = "M_D2"
  } else if(fnameRoot == "MDM2") {
    fnameRoot = "M_DM2"
  }
  
  # get what stratum we're predicting at
  out = load("savedOutput/validation/edMICSval.RData")
  strata = sort(unique(edVal$Stratum))
  foldArea = strata[fold]
  
  gridPreds = predGrid(SD0=SD0, tmbObj=obj, 
                       normalized=TRUE, extractMethod="bilinear", 
                       nsim=nsim, quantiles=quantiles, 
                       splineApprox=TRUE, admLevel="stratMICS", 
                       predAtArea=foldArea)
  
  stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
  
  stratPreds
}

scoreValidationPreds = function(fold, 
                                model=c("Md", "MD", "Mdm", "MDM", "Md2", "MD2", "Mdm2", "MDM2"), 
                                regenScores=FALSE, areal=FALSE, varClust=FALSE) {
  # clean input arguments
  model = match.arg(model)
  foldMICS = fold - 10
  
  # load in the data for the appropriate model
  fnameRoot = model
  if(fnameRoot == "MD") {
    fnameRoot = "M_D"
  } else if(fnameRoot == "MDM") {
    fnameRoot = "M_DM"
  } else if(fnameRoot == "MD2") {
    fnameRoot = "M_D2"
  } else if(fnameRoot == "MDM2") {
    fnameRoot = "M_DM2"
  }
  
  if(areal) {
    fnameRoot = paste0(fnameRoot, "areal", collapse="")
  }
  
  if(varClust) {
    fnameRoot = paste0(fnameRoot, "VarClust", collapse="")
  }
  
  # load predictions
  out = load(paste0("savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
  # list(probDrawsUrb, probDrawsRur, predsUrb, predsRur, 
  #      quantsUrb, quantsRur, yUrb, yRur, nUrb, nRur)
  
  if(!is.null(preds)) {
    
    if(!areal) {
      probDrawsUrb = preds$probDrawsUrb
      predsUrb = preds$predsUrb
      yUrb = preds$yUrb
      nUrb = preds$nUrb
      probDrawsRur = preds$probDrawsRur
      predsRur = preds$predsRur
      yRur = preds$yRur
      nRur = preds$nRur
      ys = c(yUrb, yRur)
      ns = c(nUrb, nRur)
      
      # combine predictions from urban and rural areas
      probDraws = rbind(probDrawsUrb, probDrawsRur)
      preds = c(predsUrb, predsRur)
      
      # calculate score
      scoresUrb = getScores(truth=yUrb/nUrb, estMat=probDrawsUrb, weights=nUrb, 
                            significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, 
                            getAverage=TRUE, na.rm=TRUE, ns=nUrb)
      scoresRur = getScores(truth=yRur/nRur, estMat=probDrawsRur, weights=nRur, 
                            significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, 
                            getAverage=TRUE, na.rm=TRUE, ns=nRur)
      scores = getScores(truth=ys/ns, estMat=probDraws, weights=ns, 
                         significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE, 
                         getAverage=TRUE, na.rm=TRUE, ns=ns)
      
      # add computation time
      scoresUrb = cbind(scoresUrb, Time=totalTime/60)
      scoresRur = cbind(scoresRur, Time=totalTime/60)
      scores = cbind(scores, Time=totalTime/60)
      
      save(scoresUrb, scoresRur, scores, preds, file=paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData"))
    } else {
      out = load("savedOutput/validation/directEsts.RData")
      
      probDraws = preds$aggregationResults$p
      
      # calculate score
      scoresDHS = getScoresDirectEstimates(logitDirectEsts=estsDHS$logit.est[fold], 
                                           logitDirectEstsVar=estsDHS$logit.var[fold], 
                                           estMat=matrix(probDraws[fold,], nrow=1), 
                                           significance=c(.5, .8, .9, .95), 
                                           getAverage=TRUE, na.rm=TRUE)
      scoresMICS = getScoresDirectEstimates(logitDirectEsts=estsMICS$logit.est[fold], 
                                           logitDirectEstsVar=estsMICS$logit.var[fold], 
                                           estMat=matrix(probDraws[fold,], nrow=1), 
                                           significance=c(.5, .8, .9, .95), 
                                           getAverage=TRUE, na.rm=TRUE)
      scores = getScoresDirectEstimates(logitDirectEsts=ests$logit.est[fold], 
                                        logitDirectEstsVar=ests$logit.var[fold], 
                                        estMat=matrix(probDraws[fold,], nrow=1), 
                                        significance=c(.5, .8, .9, .95), 
                                        getAverage=TRUE, na.rm=TRUE)
      
      # add computation time
      scoresDHS = cbind(scoresDHS, Time=totalTime/60)
      scoresMICS = cbind(scoresMICS, Time=totalTime/60)
      scores = cbind(scores, Time=totalTime/60)
      
      save(scoresDHS, scoresMICS, scores, preds, file=paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData"))
    }
  }
  
  invisible(NULL)
}


# function for collecting validation results for each model
validationTable = function(quantiles=c(0.025, 0.1, 0.9, 0.975), areal=FALSE, 
                           admLevel=c("adm2", "admFinal", "all"), varClust=FALSE) {
  admLevel = match.arg(admLevel)
  
  if(admLevel == "adm2") {
    models=c("Md2", "MD2", "Mdm2", "MDM2")
  } else if(admLevel == "admFinal") {
    models=c("Md", "MD", "Mdm", "MDM")
  } else if(admLevel == "all") {
    models=c("Md", "MD", "Mdm", "MDM", "Md2", "MD2", "Mdm2", "MDM2")
  }
  
  if(!areal) {
    folds = 1:20
  } else {
    folds = 1:37
  }
  
  thisSummaryFun = function(x) {
    c(mean(x), sd(x), quantile(x, probs=quantiles))
  }
  
  # calculate the fold weights (overall, urban, and rural)
  # Note that the scores within each fold average have already been weighted by 
  # the cluster ns so we just have to weight by the fold total ns here
  out = load("savedOutput/validation/edVal.RData")
  out = load("savedOutput/validation/edMICSval.RData")
  if(!areal) {
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
  }
  
  # aggregate the results of all folds for each model
  if(areal) {
    scoresTabsFull = list()
    parTabsFull = list()
  }
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
    } else if(fnameRoot == "MD2") {
      fnameRoot = "M_D2"
    } else if(fnameRoot == "MDM2") {
      fnameRoot = "M_DM2"
    }
    
    if(areal) {
      fnameRoot = paste0(fnameRoot, "areal", collapse="")
    }
    
    if(varClust) {
      fnameRoot = paste0(fnameRoot, "VarClust", collapse="")
    }
    
    # collect the results over all folds
    if(areal) {
      thisScoresTabFull = list() # only for areal case
      thisParTabFull = list()
    }
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
      
      if(!areal) {
        isMICS = j > 10
      }
      
      # load the scores
      if(file.exists(paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData")) &&
         file.exists(paste0("~/git/jittering/savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))) {
        out = load(paste0("savedOutput/validation/folds/scores", fnameRoot, "_fold", fold, ".RData"))
        out = load(paste0("~/git/jittering/savedOutput/validation/folds/preds", fnameRoot, "_fold", fold, ".RData"))
        
        if(is.null(preds) || is.null(scores)) {
          warning(paste0("NULL scores or preds for model ", fnameRoot, " fold ", fold))
          
          # add filler
          if(!areal && !isMICS) {
            thisScoresTabDHS = rbind(thisScoresTabDHS, NA)
            thisScoresTabUrbDHS = rbind(thisScoresTabUrbDHS, NA)
            thisScoresTabRurDHS = rbind(thisScoresTabRurDHS, NA)
            
            thisParTabDHS = c(thisParTabDHS, list(NA))
          } else if(!areal) {
            thisScoresTabMICS = rbind(thisScoresTabMICS, NA)
            thisScoresTabUrbMICS = rbind(thisScoresTabUrbMICS, NA)
            thisScoresTabRurMICS = rbind(thisScoresTabRurMICS, NA)
            
            thisParTabMICS = c(thisParTabMICS, list(NA))
          } else {
            thisScoresTabDHS = rbind(thisScoresTabDHS, NA)
            thisScoresTabMICS = rbind(thisScoresTabDHS, NA)
            thisScoresTabFull = rbind(thisScoresTabFull, NA)
          }
          
          next
        }
      } else {
        warning(paste0("no scores or preds file for model ", fnameRoot, " fold ", fold))
        
        # add filler
        if(!areal && !isMICS) {
          thisScoresTabDHS = rbind(thisScoresTabDHS, NA)
          thisScoresTabUrbDHS = rbind(thisScoresTabUrbDHS, NA)
          thisScoresTabRurDHS = rbind(thisScoresTabRurDHS, NA)
          
          thisParTabDHS = c(thisParTabDHS, list(NA))
        } else if(!areal) {
          thisScoresTabMICS = rbind(thisScoresTabMICS, NA)
          thisScoresTabUrbMICS = rbind(thisScoresTabUrbMICS, NA)
          thisScoresTabRurMICS = rbind(thisScoresTabRurMICS, NA)
          
          thisParTabMICS = c(thisParTabMICS, list(NA))
        } else {
          thisScoresTabDHS = rbind(thisScoresTabDHS, NA)
          thisScoresTabMICS = rbind(thisScoresTabDHS, NA)
          thisScoresTabFull = rbind(thisScoresTabFull, NA)
          
          thisParTabFull = rbind(thisParTabFull, NA)
        }
        
        next
      }
      
      # calculate parameter summary statistics
      foldParTab = t(apply(t(preds$fixedMat), 2, thisSummaryFun))
      colnames(foldParTab) = c("Est", "SD", paste0("Q", quantiles*100))
      
      if(!areal && !isMICS) {
        thisScoresTabDHS = rbind(thisScoresTabDHS, scores)
        thisScoresTabUrbDHS = rbind(thisScoresTabUrbDHS, scoresUrb)
        thisScoresTabRurDHS = rbind(thisScoresTabRurDHS, scoresRur)
        
        thisParTabDHS = c(thisParTabDHS, list(foldParTab))
      } else if(!areal) {
        thisScoresTabMICS = rbind(thisScoresTabMICS, scores)
        thisScoresTabUrbMICS = rbind(thisScoresTabUrbMICS, scoresUrb)
        thisScoresTabRurMICS = rbind(thisScoresTabRurMICS, scoresRur)
        
        thisParTabMICS = c(thisParTabMICS, list(foldParTab))
      } else {
        # if((model == "Md2") || (model == "MDM2")) {
        #   browser()
        # }
        
        thisScoresTabDHS = rbind(thisScoresTabDHS, scoresDHS)
        thisScoresTabMICS = rbind(thisScoresTabMICS, scoresMICS)
        thisScoresTabFull = rbind(thisScoresTabFull, scores)
        
        thisParTabFull = c(thisParTabFull, list(foldParTab))
      }
    }
    
    scoresTabsDHS = c(scoresTabsDHS, list(thisScoresTabDHS))
    scoresTabsMICS = c(scoresTabsMICS, list(thisScoresTabMICS))
    if(!areal) {
      parTabsDHS = c(parTabsDHS, list(thisParTabDHS))
      parTabsMICS = c(parTabsMICS, list(thisParTabMICS))
      scoresTabsUrbDHS = c(scoresTabsUrbDHS, list(thisScoresTabUrbDHS))
      scoresTabsRurDHS = c(scoresTabsRurDHS, list(thisScoresTabRurDHS))
      scoresTabsUrbMICS = c(scoresTabsUrbMICS, list(thisScoresTabUrbMICS))
      scoresTabsRurMICS = c(scoresTabsRurMICS, list(thisScoresTabRurMICS))
    } else {
      parTabsFull = c(parTabsFull, list(thisParTabFull))
      scoresTabsFull = c(scoresTabsFull, list(thisScoresTabFull))
    }
  }
  
  # calculate averages
  if(!areal) {
    scoresTabsAvgDHS = do.call("rbind", lapply(scoresTabsDHS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsDHS[!is.na(x[,1])]/sum(weightsDHS[!is.na(x[,1])]), "*"))}))
    scoresTabsUrbAvgDHS = do.call("rbind", lapply(scoresTabsUrbDHS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsUrbDHS[!is.na(x[,1])]/sum(weightsUrbDHS[!is.na(x[,1])]), "*"))}))
    scoresTabsRurAvgDHS = do.call("rbind", lapply(scoresTabsRurDHS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsRurDHS[!is.na(x[,1])]/sum(weightsRurDHS[!is.na(x[,1])]), "*"))}))
    scoresTabsAvgMICS = do.call("rbind", lapply(scoresTabsMICS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsMICS[!is.na(x[,1])]/sum(weightsMICS[!is.na(x[,1])]), "*"))}))
    scoresTabsUrbAvgMICS = do.call("rbind", lapply(scoresTabsUrbMICS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsUrbMICS[!is.na(x[,1])]/sum(weightsUrbMICS[!is.na(x[,1])]), "*"))}))
    scoresTabsRurAvgMICS = do.call("rbind", lapply(scoresTabsRurMICS, function(x) {colSums(sweep(x[!is.na(x[,1]),], 1, weightsRurMICS[!is.na(x[,1])]/sum(weightsRurMICS[!is.na(x[,1])]), "*"))}))
    
    colMeds = function(x) {apply(x, 2, median, na.rm=TRUE)}
    scoresTabsMedDHS = do.call("rbind", lapply(scoresTabsDHS, colMeds))
    scoresTabsMedMICS = do.call("rbind", lapply(scoresTabsMICS, colMeds))
    scoresTabsUrbMedDHS = do.call("rbind", lapply(scoresTabsUrbDHS, colMeds))
    scoresTabsRurMedDHS = do.call("rbind", lapply(scoresTabsRurDHS, colMeds))
    scoresTabsUrbMedMICS = do.call("rbind", lapply(scoresTabsUrbMICS, colMeds))
    scoresTabsRurMedMICS = do.call("rbind", lapply(scoresTabsRurMICS, colMeds))
  } else {
    # no need for weighted averaging in this case
    scoresTabsAvgDHS = do.call("rbind", lapply(scoresTabsDHS, colMeans, na.rm=TRUE))
    scoresTabsAvgMICS = do.call("rbind", lapply(scoresTabsMICS, colMeans, na.rm=TRUE))
    scoresTabsAvgFull = do.call("rbind", lapply(scoresTabsFull, colMeans, na.rm=TRUE))
    
    colMeds = function(x) {apply(x, 2, median, na.rm=TRUE)}
    scoresTabsMedDHS = do.call("rbind", lapply(scoresTabsDHS, colMeds))
    scoresTabsMedMICS = do.call("rbind", lapply(scoresTabsMICS, colMeds))
    scoresTabsMedFull = do.call("rbind", lapply(scoresTabsFull, colMeds))
  }
  
  if(!areal) {
    parTabsAvgDHS = lapply(parTabsDHS, function(x) {
      badFits = sapply(x, function(l) {any(is.na(unlist(l)))})
      x[badFits] = NULL
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
      badFits = sapply(x, function(l) {any(is.na(unlist(l)))})
      x[badFits] = NULL
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
      rNames = row.names(parTabsAvgDHS[[i]])
      cNames = colnames(parTabsAvgDHS[[i]])
      thisParTab = abind(parTabsAvgDHS[[i]], parTabsAvgMICS[[i]], along=3)
      thisParTab = apply(thisParTab, 1:2, mean)
      
      row.names(thisParTab) = rNames
      colnames(thisParTab) = cNames
      thisParTab
    })
    names(parTabsAvg) = models
  } else {
    parTabsAvg = lapply(parTabsFull, function(x) {
      badFits = sapply(x, function(l) {any(is.na(unlist(l)))})
      x[badFits] = NULL
      rNames = row.names(x[[1]])
      cNames = colnames(x[[1]])
      x = lapply(x, function(y) {array(y, dim=c(dim(y), 1))})
      scoreArray = do.call("abind", list(x, along=3))
      out = apply(scoreArray, 1:2, mean)
      row.names(out) = rNames
      colnames(out) = cNames
      out
    })
    names(parTabsAvg) = models
  }
  
  # take a weighted average of MICS and DHS scores by total people in datasets 
  # in case of cluster validation
  nMICS = sum(edMICSval$ns)
  nUrbMICS = sum(edMICSval$ns[edMICSval$urban])
  nRurMICS = sum(edMICSval$ns[!edMICSval$urban])
  nDHS = sum(edVal$n)
  nUrbDHS = sum(edVal$n[edVal$urban])
  nRurDHS = sum(edVal$n[!edVal$urban])
  nTotal = nMICS + nDHS
  nUrbTotal = nUrbMICS + nUrbDHS
  nRurTotal = nRurMICS + nRurDHS
  
  wMICS = nMICS/nTotal
  wDHS = nDHS/nTotal
  wUrbMICS = nUrbMICS/nUrbTotal
  wUrbDHS = nUrbDHS/nUrbTotal
  wRurMICS = nRurMICS/nRurTotal
  wRurDHS = nRurDHS/nRurTotal
  
  if(!areal) {
    finalTabAvg = scoresTabsAvgDHS * wDHS + scoresTabsAvgMICS * wMICS
    finalTabUrbAvg = scoresTabsUrbAvgDHS * wUrbDHS + scoresTabsUrbAvgMICS * wUrbMICS
    finalTabRurAvg = scoresTabsRurAvgDHS * wRurDHS + scoresTabsRurAvgMICS * wRurMICS
    
    finalTabDHSAvg = scoresTabsAvgDHS
    finalTabMICSAvg = scoresTabsAvgMICS
    finalTabDHSMed = scoresTabsMedDHS
    finalTabMICSMed = scoresTabsMedMICS
    finalTabUrbDHSMed = scoresTabsUrbMedDHS
    finalTabRurDHSMed = scoresTabsRurMedDHS
    finalTabUrbMICSMed = scoresTabsUrbMedMICS
    finalTabRurMICSMed = scoresTabsRurMedMICS
  } else {
    finalTabAvg = scoresTabsAvgFull
    finalTabDHSAvg = scoresTabsAvgDHS
    finalTabMICSAvg = scoresTabsAvgMICS
    
    finalTabMed = scoresTabsMedFull
    finalTabDHSMed = scoresTabsMedDHS
    finalTabMICSMed = scoresTabsMedMICS
    
    finalTabUrbAvg = NULL
    finalTabRurAvg = NULL
  }
  
  # plot scores spatially ----
  predCols = makeBlueGreenYellowSequentialColors(64)
  quantCols = makePurpleYellowSequentialColors(64)
  arealText = ifelse(areal, "areal", "")
  admLevelText = admLevel
  if(areal) {
    pdf(paste0("figures/validation/", admLevelText, arealText, "MSEdiff.pdf"), width=6, height=6)
    theseCols = makeRedBlueDivergingColors(64, valRange=range(scoresTabsFull[[4]]$MSE-scoresTabsFull[[2]]$MSE), 
                                           center=0, rev=TRUE)
    plotMapDat(adm1, scoresTabsFull[[4]]$MSE-scoresTabsFull[[2]]$MSE, 
               varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, cols=theseCols, 
               main="MSE M_DM2 - MSE M_D2")
    dev.off()
  }
  
  # save/print out tables ----
  
  browser()
  
  
  # old results: ----
  
  #          Bias        Var       MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95
  # 1 0.003537896 0.09405488 0.0940674 0.3067041 0.1726749       0.7729948       0.9801325        1.110449        1.239974
  #   Coverage50 Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95
  # 1  0.3782958  0.6722866  0.8125554  0.8765654 0.3623988 0.6269316 0.7425313 0.8279065
  
  # finalTabAvg
  #              Bias        Var        MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,] -0.048923416 0.10258375 0.10724543 0.3269813 0.1848153       0.8074385       1.0257699       1.1268245       1.1880048  0.3220954
  # [2,] -0.042409437 0.08963965 0.09597390 0.3078246 0.1739119       0.7602655       0.9553505       1.0343832       1.0760787  0.3216234
  # [3,]  0.008536307 0.12283088 0.12347722 0.3498838 0.2006246       0.8801756       1.1053651       1.1743620       1.2067275  0.3412146
  # [4,]  0.016442123 0.09399670 0.09507257 0.3076165 0.1742858       0.7629908       0.9098012       0.9606157       0.9930743  0.3339502
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95       Time
  # [1,]  0.6241462  0.7698937  0.8659569 0.3304475 0.6010558 0.7258421 0.8153675 0.09974754
  # [2,]  0.6114158  0.7642869  0.8640286 0.3305149 0.6013947 0.7313789 0.8207659 0.84239516
  # [3,]  0.6279874  0.7664817  0.8696364 0.3663356 0.6523856 0.7778282 0.8572857 0.19906408
  # [4,]  0.6212378  0.7842337  0.8822892 0.3775339 0.6754577 0.8055434 0.8871717 2.95259512
  
  # scoresTabsAvgDHS
  #              Bias        Var        MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,] -0.010947783 0.09962442 0.10035338 0.3165307 0.1792139       0.7960109       1.0360634       1.1507908       1.2423856  0.3457359
  # [2,]  0.015759156 0.07685016 0.07760338 0.2783876 0.1575662       0.7033438       0.8709920       0.9333493       0.9734245  0.3400876
  # [3,] -0.001068097 0.10353460 0.10405081 0.3223983 0.1824954       0.8136387       1.0223140       1.0884226       1.1213663  0.3582090
  # [4,]  0.031055866 0.08425341 0.08576563 0.2925230 0.1664032       0.7415413       0.9233491       0.9609172       1.0130063  0.3475848
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95       Time
  # [1,]  0.6383757  0.7765516  0.8688758 0.3307033 0.5951478 0.7165991 0.8049035 0.09359287
  # [2,]  0.6199160  0.7745838  0.8749643 0.3436871 0.6169147 0.7476554 0.8363716 0.86390443
  # [3,]  0.6494250  0.7872871  0.8881316 0.3627088 0.6400537 0.7670350 0.8499071 0.20021130
  # [4,]  0.6179197  0.7768636  0.8741421 0.3727287 0.6630486 0.7923907 0.8760249 3.00214026
  
  # scoresTabsAvgMICS
  #              Bias       Var       MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,] -0.094478535 0.1061337 0.1155131 0.3395178 0.1915348       0.8211469       1.0134220        1.098075       1.1227702  0.2937365
  # [2,] -0.112187791 0.1049818 0.1180110 0.3431368 0.1935200       0.8285482       1.0565461        1.155582       1.1992213  0.2994740
  # [3,]  0.020057637 0.1459785 0.1467809 0.3828551 0.2223721       0.9599925       1.2049921        1.277454       1.3091258  0.3208282
  # [4,] -0.001088348 0.1056846 0.1062371 0.3257226 0.1837417       0.7887214       0.8935493        0.960254       0.9691641  0.3175942
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95      Time
  # [1,]  0.6070766  0.7619071  0.8624555 0.3301406 0.6081430 0.7369298 0.8279200 0.1071306
  # [2,]  0.6012190  0.7519349  0.8509102 0.3147137 0.5827772 0.7118537 0.8020455 0.8165929
  # [3,]  0.6022712  0.7415239  0.8474497 0.3706861 0.6671788 0.7907756 0.8661370 0.1976879
  # [4,]  0.6252181  0.7930748  0.8920624 0.3832982 0.6903434 0.8213213 0.9005434 2.8931613
  
  # parTabsAvg
  # $Md
  # Est         SD       Q2.5        Q10         Q90       Q97.5
  # (Int)           -1.94129596 0.11905241 -2.1730412 -2.0938546 -1.78924275 -1.70940508
  # urban            0.29642889 0.06230221  0.1749010  0.2174880  0.37667100  0.41814940
  # access          -0.18153894 0.11248702 -0.4004782 -0.3237846 -0.03714779  0.03855007
  # elev             0.14067058 0.27463587 -0.3977896 -0.2064539  0.49131870  0.67851991
  # distRiversLakes  0.08114496 0.44009669 -0.7794160 -0.4836371  0.64039444  0.93638650
  # pop              0.93009097 0.07523636  0.7832854  0.8341408  1.02681785  1.07645384
  # sigmaSq          0.43019316 0.05426750  0.3333176  0.3635974  0.50129899  0.54509118
  # phi              0.48659818 0.01390103  0.4591779  0.4687317  0.50425673  0.51378021
  # sigmaEpsSq       1.44432079 0.09746750  1.2657981  1.3209650  1.57211626  1.64474648
  # 
  # $MD
  # Est         SD       Q2.5        Q10         Q90       Q97.5
  # (Int)           -2.02993926 0.12300228 -2.2678335 -2.1858132 -1.87239691 -1.78734017
  # urban            0.43797318 0.06108454  0.3190559  0.3594872  0.51567260  0.55730990
  # access          -0.11460163 0.11120885 -0.3351878 -0.2584252  0.02797801  0.09872751
  # elev             0.12492345 0.27461989 -0.4066872 -0.2289754  0.47616777  0.66681025
  # distRiversLakes  0.08358633 0.44389183 -0.7796613 -0.4895262  0.65152610  0.94650479
  # pop              0.99104510 0.07654125  0.8419479  0.8937704  1.08978246  1.13990138
  # sigmaSq          0.40426090 0.04945289  0.3155787  0.3426645  0.46865263  0.50715050
  # phi              0.51621109 0.01197842  0.4930573  0.5008839  0.53150162  0.53950840
  # sigmaEpsSq       1.36000680 0.08821376  1.1972001  1.2486907  1.47400694  1.53953688
  # 
  # $Mdm
  # Est         SD       Q2.5        Q10         Q90       Q97.5
  # (Int)           -1.52283279 0.08935393 -1.6982142 -1.6364734 -1.40916937 -1.34622667
  # urban            0.87915197 0.04746879  0.7871418  0.8182994  0.93922191  0.97252775
  # access          -0.09121571 0.07689023 -0.2405250 -0.1890602  0.00754624  0.06050912
  # elev             0.04916828 0.25289755 -0.4377771 -0.2771477  0.37470865  0.54548374
  # distRiversLakes  0.03758685 0.42766601 -0.7975181 -0.5094472  0.58472839  0.86837598
  # pop              0.52891551 0.05170684  0.4283231  0.4622907  0.59536389  0.62983590
  # sigmaSq          0.39978964 0.03721313  0.3324792  0.3533561  0.44823183  0.47805137
  # phi              0.51885521 0.00987176  0.4994156  0.5061621  0.53142987  0.53798595
  # sigmaEpsSq       1.85094506 0.09167776  1.6772827  1.7343914  1.96950545  2.03388824
  # 
  # $MDM
  # Est         SD       Q2.5        Q10        Q90        Q97.5
  # (Int)           -1.38323725 0.07151436 -1.5216664 -1.4741821 -1.2912558 -1.244226439
  # urban            1.20249733 0.04191640  1.1200387  1.1488269  1.2560754  1.283440691
  # access          -0.12546053 0.06130758 -0.2444885 -0.2036496 -0.0470810 -0.005727111
  # elev             0.01906404 0.28396156 -0.5329500 -0.3444793  0.3815009  0.569318587
  # distRiversLakes  0.02815097 0.55128121 -1.0430663 -0.6833002  0.7348523  1.104860023
  # pop              0.29653954 0.06048748  0.1791659  0.2182861  0.3734627  0.414436811
  # sigmaSq          0.35601707 0.03147675  0.2992667  0.3165689  0.3970594  0.421319398
  # phi              0.69022354 0.00670695  0.6770004  0.6816572  0.6987463  0.703217893
  # sigmaEpsSq       1.63425590 0.07128903  1.5005885  1.5434381  1.7262691  1.777046946
  
  # newer old results: ----
  # finalTabAvg
  #              Bias        Var        MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,] -0.001195913 0.10194118 0.10255866 0.3199920 0.1815174       0.7936561       0.9894571       1.0665475       1.1303076  0.3305173
  # [2,]  0.009899270 0.08876116 0.08935446 0.2979548 0.1686452       0.7382367       0.8907166       0.9524627       0.9946624  0.3207881
  # [3,]  0.008636098 0.12281674 0.12344696 0.3498537 0.2005877       0.8813509       1.1052730       1.1646444       1.2105253  0.3429644
  # [4,] -0.002387789 0.08462564 0.08506660 0.2906171 0.1619909       0.7023527       0.8655062       0.9284315       0.9907109  0.3399711
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95       Time
  # [1,]  0.6201958  0.7744234  0.8703859 0.3473067 0.6296035 0.7605387 0.8435679 0.06937045
  # [2,]  0.6148707  0.7732745  0.8735706 0.3478748 0.6357671 0.7704163 0.8541860 0.52524633
  # [3,]  0.6274462  0.7690926  0.8685646 0.3664714 0.6526301 0.7780663 0.8575582 0.13546164
  # [4,]  0.6449592  0.8004512  0.8881724 0.3439799 0.6258084 0.7573977 0.8426754 1.77649791
  
  # scoresTabsAvgDHS
  #               Bias        Var        MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,] -0.0109389127 0.09954326 0.10023341 0.3163467 0.1791701       0.7958612       1.0330904       1.1513947       1.2422438  0.3504919
  # [2,]  0.0166653791 0.07693741 0.07767767 0.2785133 0.1576946       0.7044831       0.8668311       0.9413000       0.9799486  0.3356106
  # [3,] -0.0005162375 0.10361923 0.10412636 0.3225182 0.1825279       0.8140095       1.0200468       1.0841962       1.1337893  0.3607627
  # [4,]  0.0002172389 0.07342022 0.07385371 0.2714786 0.1513177       0.6699246       0.8480120       0.9113472       0.9746350  0.3573118
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95       Time
  # [1,]  0.6372166  0.7799045  0.8672647 0.3315550 0.5933909 0.7197066 0.8050478 0.06626508
  # [2,]  0.6213675  0.7746574  0.8727310 0.3418937 0.6162588 0.7483100 0.8340258 0.52347791
  # [3,]  0.6499489  0.7901513  0.8862632 0.3627219 0.6401469 0.7685309 0.8508136 0.13644521
  # [4,]  0.6513264  0.8027665  0.8890711 0.3376643 0.6072395 0.7370268 0.8229898 1.79488796
  
  # scoresTabsAvgMICS
  #              Bias        Var        MSE      RMSE      CRPS IntervalScore50 IntervalScore80 IntervalScore90 IntervalScore95 Coverage50
  # [1,]  0.010491674 0.10481769 0.10534801 0.3243648 0.1843333       0.7910109       0.9371151       0.9647657       0.9960302  0.3065559
  # [2,]  0.001782726 0.10294479 0.10336179 0.3212767 0.1817814       0.7787271       0.9193693       0.9658533       1.0123129  0.3030072
  # [3,]  0.019615132 0.14584585 0.14662372 0.3826451 0.2222521       0.9621329       1.2075094       1.2611491       1.3025769  0.3216137
  # [4,] -0.005512750 0.09806752 0.09851745 0.3135756 0.1747943       0.7412530       0.8864920       0.9489256       1.0099954  0.3191694
  #      Coverage80 Coverage90 Coverage95   Width50   Width80   Width90   Width95       Time
  # [1,]  0.5997779  0.7678484  0.8741301 0.3662022 0.6730438 0.8095204 0.8897761 0.07309562
  # [2,]  0.6070771  0.7716157  0.8745777 0.3550497 0.6591690 0.7969346 0.8783698 0.52736771
  # [3,]  0.6004523  0.7438308  0.8473336 0.3709693 0.6676047 0.7895048 0.8656490 0.13428177
  # [4,]  0.6373212  0.7976738  0.8870943 0.3515560 0.6480835 0.7818343 0.8662901 1.75443743
  
  # parTabsAvg
  # $Md
  # Est         SD       Q2.5        Q10         Q90       Q97.5
  # (Int)           -1.94278262 0.11792990 -2.1723209 -2.0939462 -1.79110063 -1.71240072
  # urb              0.29569057 0.06170367  0.1761480  0.2173431  0.37487710  0.41664975
  # access          -0.17989018 0.11187726 -0.3968374 -0.3246100 -0.03650993  0.03955828
  # elev             0.13729947 0.27852403 -0.4043243 -0.2194224  0.49458244  0.68639357
  # distRiversLakes  0.08095797 0.43777081 -0.7787264 -0.4854279  0.63581016  0.92753396
  # pop              0.92950334 0.07549213  0.7829236  0.8334573  1.02580797  1.07724773
  # sigmaSq          0.42995568 0.05379267  0.3344933  0.3639561  0.50001503  0.54460481
  # phi              0.48681625 0.01390506  0.4596718  0.4689324  0.50462285  0.51387525
  # sigmaEpsSq       1.44409397 0.09571045  1.2666877  1.3232987  1.56821778  1.63796377
  # 
  # $MD
  # Est         SD       Q2.5        Q10         Q90      Q97.5
  # (Int)           -2.02872730 0.12268758 -2.2674606 -2.1855941 -1.87043313 -1.7860693
  # urb              0.43780777 0.06066964  0.3209290  0.3603642  0.51604613  0.5573030
  # access          -0.11512487 0.11077485 -0.3301644 -0.2569739  0.02633578  0.1008438
  # elev             0.12559561 0.27334361 -0.4094415 -0.2247518  0.47629437  0.6576237
  # distRiversLakes  0.08648196 0.44859381 -0.7818074 -0.4933770  0.65877585  0.9646978
  # pop              0.99162538 0.07737600  0.8407588  0.8935913  1.09028762  1.1433731
  # sigmaSq          0.40442384 0.04960489  0.3163684  0.3434661  0.46921096  0.5093302
  # phi              0.51620363 0.01204519  0.4924113  0.5007769  0.53146116  0.5396802
  # sigmaEpsSq       1.36035637 0.08822314  1.1960919  1.2498013  1.47417042  1.5404201
  # 
  # $Mdm
  # Est          SD       Q2.5        Q10          Q90      Q97.5
  # (Int)           -1.52384488 0.089736498 -1.6990790 -1.6383775 -1.409769043 -1.3500377
  # urb              0.87864330 0.047206337  0.7869515  0.8184040  0.939141815  0.9716784
  # access          -0.09049749 0.077034204 -0.2383977 -0.1897603  0.008674919  0.0591374
  # elev             0.05105080 0.253198063 -0.4409549 -0.2751215  0.374459752  0.5434192
  # distRiversLakes  0.03451236 0.426822880 -0.7988400 -0.5076707  0.581463763  0.8744702
  # pop              0.52947219 0.051857242  0.4283591  0.4629830  0.596604883  0.6295748
  # sigmaSq          0.39957936 0.037013459  0.3328260  0.3533029  0.448288244  0.4768059
  # phi              0.51884753 0.009880103  0.4996800  0.5061710  0.531656748  0.5381797
  # sigmaEpsSq       1.85097677 0.092683319  1.6757986  1.7347862  1.970900723  2.0369814
  # 
  # $MDM
  #                         Est          SD       Q2.5         Q10        Q90      Q97.5
  # (Int)           -2.29139454 0.100884525 -2.4888623 -2.41969663 -2.1618440 -2.0951868
  # urb              0.78463597 0.045132599  0.6983759  0.72695004  0.8426546  0.8729954
  # access           0.07691087 0.103974703 -0.1251478 -0.05516393  0.2106283  0.2765126
  # elev            -0.02118382 0.315581474 -0.6415125 -0.42754692  0.3835422  0.5902366
  # distRiversLakes  0.04003146 0.480813411 -0.8960995 -0.56906068  0.6573468  0.9755319
  # pop              1.26075480 0.062704648  1.1392759  1.17982723  1.3408474  1.3817880
  # sigmaSq          0.37362761 0.036815456  0.3071483  0.32805320  0.4224336  0.4500976
  # phi              0.82253770 0.005424617  0.8117663  0.81553814  0.8294229  0.8329414
  # sigmaEpsSq       1.41296094 0.065694832  1.2892284  1.32994917  1.4973411  1.5449716
  
  # new results: ----
}

# function for getting combined MICS and DHS direct estimates at the admin1 level
# clustDatDHS: a subset of edVal
# clustDatMICS: a subset of edMICSval
getCombinedDirectEsts = function(clustDatDHS=edVal, clustDatMICS=edMICSval, 
                                 divideWeight=TRUE, signifs=c(.5, .8, .9, .95), 
                                 leftOutFolds=6:10) {
  
  clustDatDHS = clustDatDHS[clustDatDHS$fold %in% leftOutFolds,]
  clustDatMICS = clustDatMICS[clustDatMICS$fold %in% leftOutFolds,]
  
  # convert DHS data to the correct format
  # clustDatDHS$area = NULL
  # clustDatDHS$subarea = NULL only remove this because it could be confusing
  # names(clustDatDHS)[grepl("Stratum", names(clustDatDHS))] = "area"
  names(clustDatDHS)[grepl("clusterID", names(clustDatDHS))] = "clustID"
  
  estsDHS = getDirectEsts(clustDatDHS, divideWeight=divideWeight, signifs=signifs)
  
  names(clustDatMICS)[grepl("ys", names(clustDatMICS))] = "y"
  names(clustDatMICS)[grepl("ns", names(clustDatMICS))] = "n"
  names(clustDatMICS)[grepl("Area", names(clustDatMICS))] = "area"
  
  estsMICS = getDirectEsts(clustDatMICS, divideWeight=divideWeight, 
                           signifs=signifs, customStratVarName="Stratum")
  
  # take precision weighted average of admin1 estimates
  precsDHS = 1/estsDHS$logit.var
  precsMICS = 1/estsMICS$logit.var
  wsDHS = precsDHS/(precsDHS + precsMICS)
  wsMICS = 1 - wsDHS
  ests = estsDHS
  ests$logit.est = estsDHS$logit.est * wsDHS + estsMICS$logit.est * wsMICS
  ests$est = logitNormMeanSimple(cbind(ests$logit.est, sqrt(ests$logit.var)))
  ests$logit.var = estsDHS$logit.var * wsDHS^2 + estsMICS$logit.var * wsMICS^2
  ests$var = logitNormVarSimple(cbind(ests$logit.est, sqrt(ests$logit.var)))
  
  # calculate CI quantiles of the new estimator. Convergence is only true if both converged
  lowerQuants = (1-signifs)/2
  upperQuants = 1 - (1-signifs)/2
  ests[,6:(5+length(signifs))] = sweep(t(sapply(sqrt(ests$logit.var), qnorm, mean=0, p=lowerQuants)), 1, ests$logit.est, "+")
  ests[,(6+length(signifs)):(5+2*length(signifs))] = sweep(t(sapply(sqrt(ests$logit.var), qnorm, mean=0, p=upperQuants)), 1, ests$logit.est, "+")
  ests$converge = apply(cbind(estsDHS$converge, estsMICS$converge), 1, max)
  
  save(ests, estsDHS, estsMICS, file="savedOutput/validation/directEsts.RData")
  
  list(ests=ests, estsDHS=estsDHS, estsMICS=estsMICS)
}