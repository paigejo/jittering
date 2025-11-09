# run motivating simulation study

runSimStudy1I = function(i, mod=c("M_M", "M_DM", "M_D", "Md"), regenData=FALSE, signif=c(.8, .95)) {
  mod = match.arg(mod)
  
  # load data
  out = load("savedOutput/simStudy1/simPopsSurveys.RData")
  # subareaPops, areaPops, stratumPops, surveysDHS, surveysMICS
  
  thisDHS = surveysDHS[[i]]
  thisMICS = surveysMICS[[i]]
  
  thisSubareaPop = subareaPops[,i]
  thisStratumPop = stratumPops[,i]
  thisAreaPop = areaPops[,i]
  
  # check if we've computed the DHS integration points for this run
  # dhsIntFile = paste0("savedOutput/simStudy1/intPtsDHS_sim", i, ".RData")
  dhsIntFile = paste0("savedOutput/simStudy1/intPtsDHS_simStudy1_", i, ".RData")
  if(file.exists(dhsIntFile)) {
    out = load(dhsIntFile)
  } else {
    KDHSurb = 11 # 3 rings of 5 each
    JInnerUrban = 3
    KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
    JInnerRural = 3
    JOuterRural = 1
    intPtsDHS = makeAllIntegrationPointsDHS(cbind(thisDHS$east, thisDHS$north), thisDHS$urban, 
                                            areaNames=thisDHS$subarea, popPrior=TRUE, 
                                            numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                            JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                            JOuterRural=JOuterRural, adminMap=adm2Full, saveOutput=FALSE)
    
    save(intPtsDHS, file=dhsIntFile)
  }
  
  scoreFile = paste0("savedOutput/simStudy1/scores", mod, "_SepRepar_", i, ".RData")
  if(regenData || !file.exists(scoreFile)) {
    # fit TMB model
    if(mod == "M_M") {
      out = fitMM(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, repar=TRUE)
    } else if(mod == "M_DM") {
      out = fitMDM(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, repar=TRUE)
    } else if(mod == "M_D") {
      out = fitMD(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, repar=TRUE)
    } else if(mod == "Md") {
      out = fitMd(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, repar=TRUE)
    }
    
    # generate prediction grid
    system.time(gridPreds <- predGrid(out$TMBsd, popMat=popMatNGAThresh, nsim=5000, admLevel="stratMICS", 
                                      quantiles=c(0.025, 0.1, 0.9, 0.975), sep=TRUE))[3]
    # 228.571  seconds
    
    # get aggregated predictions and parameter summary table
    stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
    admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
    admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
    
    parTab = summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, 
                            gridPreds=gridPreds)
    
    # save results
    # save(stratPreds, admin1Preds, admin2Preds, parTab, parMat, file=predFile)
    # save(stratPreds, admin1Preds, parTab, parMat, file=predFile)
    
    # calculate scores for predictions
    scoresAdm2 = getScores(thisSubareaPop, estMat=admin2Preds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
    scoresAdm1 = getScores(thisAreaPop, estMat=admin1Preds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
    scoresStratum = getScores(thisStratumPop[match(row.names(stratPreds$aggregationResults$p), names(thisStratumPop))], estMat=stratPreds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
    
    # calculate scores for parameters
    parMat = rbind(gridPreds$alphaDraws, 
                   gridPreds$betaDraws, 
                   gridPreds$sigmaSqDraws, 
                   gridPreds$sigmaEpsSqDraws)
    
    truePar = c(-1.25, 1, 0, 0, 0, 0.5, 0.5, 1.5)
    scoresPar = getScores(truePar, estMat=parMat, significance=signif, doFuzzyReject=FALSE, getAverage=FALSE)
    row.names(scoresPar) = c("(Int)", "Urban", "Healthcare inaccessibility", 
                             "Elevation", "Dist. rivers & lakes", "Population", 
                             "sigmaSq", "sigmaEpsSq")
    
    save(scoresAdm2, scoresAdm1, scoresStratum, parTab, scoresPar, file=scoreFile)
  } else {
    invisible(NULL)
  }
}

runSimStudy1IPar = function(i, mod=c("M_M", "M_DM", "M_D", "Md"), regenData=FALSE) {
  mod = match.arg(mod)
  
  tmplogfile <- paste0("savedOutput/simStudy1/simRun_", mod, "_", i, "_tmp.txt")
  sink(tmplogfile)
  cat("Started simRun job i =", i, ":\n")
  sink()
  
  tryCatch(runSimStudy1I(i, mod=mod, regenData=regenData), 
           error = function(e) {
             logfile <- paste0("savedOutput/simStudy1/simRun_", mod, "_", i, "_err.txt")
             sink(logfile)
             cat("Error at i =", i, ":\n")
             cat(paste("Call stack:\n", paste(deparse(sys.calls()), collapse = "\n")), "\n")
             cat(conditionMessage(e), "\n")
             sink()
             stop("error")
           })
  
  # remove the tmp log file to indicate the job is complete
  system(paste0("rm ", tmplogfile))
  
}

runSimStudy1All = function(doPar=FALSE, nCores=16, regenData=FALSE, mod=c("M_M", "M_DM", "M_D", "Md")) {
  mod = match.arg(mod)
  
  is = 1:100
  
  if(doPar) {
    # start parallel cluster
    cl = makeCluster(nCores)
    clusterEvalQ(cl, source("code/setup.R"))
    
    # generate well data in parallel
    tmp = parLapply(cl, is, runSimStudy1IPar, mod=mod, regenData=regenData)
    
    # remember to stop the cluster
    stopCluster(cl)
  } else {
    for(i in is) {
      runSimStudy1I(i, mod=mod, regenData=regenData)
    }
  }
  
}


# aggregate/summarize scores and parameter estimates
simStudyScores = function(mod=c("M_M", "M_DM", "M_D", "Md")) {
  require(abind)
  mod = match.arg(mod)
  
  allScores2 = c() # admin2
  allScores1 = c() # admin1
  allScoresS = c() # MICS stratum
  allScoresP = c() # parameters
  # allParEsts = c()
  meanParStats = c()
  totalFiles = 0
  for(i in 1:100) {
    # load scores
    scoreFile = paste0("savedOutput/simStudy1/scores", mod, "_SepRepar_", i, ".RData")
    
    if(file.exists(scoreFile)) {
      totalFiles = totalFiles + 1
      out = load(scoreFile)
      # scoresAdm2, scoresAdm1, scoresStratum, parTab
      
      # concatenate
      allScores2 = rbind(allScores2, scoresAdm2)
      allScores1 = rbind(allScores1, scoresAdm1)
      allScoresS = rbind(allScoresS, scoresStratum)
      allScoresP = abind(allScoresP, scoresPar, along=3)
      
      # allParEsts = cbind(allParEsts, parTab[,1])
      
      # calculate running average of parameter summary statistics
      if(i == 1) {
        meanParStats = parTab
      } else {
        meanParStats = ((totalFiles-1) * meanParStats + parTab) * (1/totalFiles)
      }
    }
  }
  
  # average scores
  meanScoresS = colMeans(allScoresS)
  meanScores1 = colMeans(allScores1)
  meanScores2 = colMeans(allScores2)
  meanScoresP = apply(allScoresP, c(1,2), mean)
  
  # calculate score
  
  # make tables
  savedColI = c(1, 4, 5, 9, 13, 17, 18)
  savedColI = c(1, 5, 9, 18)
  savedColI = c(1:5, 7, 9, 11)
  digits = c(rep(3, 6), 2, 3)
  digitsPar = rep(2, 9)
  
  # population prediction scores
  scoreTab = rbind(allScoresS, 
                   allScores1, 
                   allScores2)
  scoreTab = cbind("Areal level"=rep(c("MICS stratum", "Admin1", "Admin2"), each=totalFiles), scoreTab)
  
  allMeanScores = rbind(meanScoresS, 
                        meanScores1, 
                        meanScores2)
  row.names(allMeanScores) = c("MICS Stratum scores", "Admin1 scores", "Admin2 scores")
  
  print(xtable(roundCols(allMeanScores[,savedColI], digits=digits), digits=c(0, digits)))
  # print(xtable(roundCols(t(as.matrix(meanScores1[savedColI])), digits=digits), digits=c(0, digits)))
  # print(xtable(roundCols(t(as.matrix(meanScores2[savedColI])), digits=digits), digits=c(0, digits)))
  # print(xtable(roundCols(t(as.matrix(meanScoresS[savedColI])), digits=digits), digits=c(0, digits)))
  
  # parameter scores
  # print(xtable(roundCols(meanParStats, digits=digits), digits=c(0, digits)))
  
  print(xtable(roundCols(meanScoresP[,savedColI], digits=digits), digits=c(0, digits)))
  
  # make plots?
  browser()
  truePar = c(-1.25, 1, 0, 0, 0, 0.5, 0.5, NA, 1.5)
  
  # M_M (full)
  # Browse[1]> allMeanScores
  #                           Bias          Var         MSE       RMSE       CRPS
  # MICS Stratum scores 0.01873885 0.0070810369 0.007468301 0.08476590 0.05080257
  # Admin1 scores       0.01866146 0.0009599265 0.001341765 0.03635591 0.02097000
  # Admin2 scores       0.01713520 0.0160948908 0.016425490 0.12673605 0.08723693
  #                     IntervalScore80 IntervalScore95 Coverage80 Coverage95
  # MICS Stratum scores       0.4065014       1.0572593  0.4541463  0.6073171
  # Admin1 scores             0.1337058       0.1863101  0.6710811  0.8600000
  # Admin2 scores             0.7488141       2.2177870  0.2406727  0.3593920
  #                        Width80   Width95
  # MICS Stratum scores 0.07818489 0.1193567
  # Admin1 scores       0.07626589 0.1164021
  # Admin2 scores       0.07945773 0.1212731
  # 
  # M_DM (partial)
  # Browse[1]> allMeanScores
  #                           Bias          Var          MSE       RMSE       CRPS
  # MICS Stratum scores 0.01835489 0.0063093591 0.0066695277 0.08026978 0.04836411
  # Admin1 scores       0.01833716 0.0005664941 0.0009257887 0.03018491 0.01772875
  # Admin2 scores       0.01699506 0.0155902382 0.0159050976 0.12472942 0.08841043
  #                     IntervalScore80 IntervalScore95 Coverage80 Coverage95
  # MICS Stratum scores       0.3999818       1.1159749  0.4039634  0.5556402
  # Admin1 scores             0.1144049       0.1637545  0.6283784  0.8201014
  # Admin2 scores             0.7848603       2.4929046  0.1915023  0.2880417
  #                        Width80    Width95
  # MICS Stratum scores 0.06149922 0.09385563
  # Admin1 scores       0.05904221 0.09012877
  # Admin2 scores       0.06180419 0.09435781
  
}

# check correlation between mean population density and education rate mean for 
# sim study 1
testPopDensCor = function() {
  
  # get pop density means per MICS stratum, urban and rural
  # out = load("savedOutput/global/popMatNGAedThresh.RData")
  out = load("savedOutput/global/intPtsMICS_100.RData")
  intPtsMICS = straightenMICS(intPtsMICS)
  XRur = intPtsMICS$XRur[,names(intPtsMICS$XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  XUrb = intPtsMICS$XUrb[,names(intPtsMICS$XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  wUrban = intPtsMICS$wUrban
  wRural = intPtsMICS$wRural
  
  getStratMean = function(i, j=7, urb=TRUE) {
    if(urb) {
      X = XUrb
      w = wUrban
    } else {
      X = XRur
      w = wRural
    }
    is = rep(F, 41)
    is[i] = TRUE
    is = rep(is, 100)
    sum(X[is,j] * w[i,])
  }
  getStratArea = function(i, urb=TRUE) {
    5 * sum((popMatNGAedThresh$stratumMICS == intPtsMICS$strataMICS[i]) & (popMatNGAedThresh$urban == urb))
  }
  getStratUrbWeight = function(i, urb=TRUE) {
    thisPops = popMatNGAedThresh$pop[(popMatNGAedThresh$stratumMICS == intPtsMICS$strataMICS[i])]
    thisUrbs = popMatNGAedThresh$urban[(popMatNGAedThresh$stratumMICS == intPtsMICS$strataMICS[i])]
    sum(thisPops[thisUrbs])/sum(thisPops)
  }
  normPopsUrb = sapply(1:41, getStratMean, urb=TRUE)
  normPopsRur = sapply(1:41, getStratMean, urb=FALSE)
  
  truePopUrb = poppStratMICSThresh$popUrb
  truePopRur = poppStratMICSThresh$popRur
  trueAreaUrb = sapply(1:41, getStratArea, urb=TRUE)
  trueAreaRur = sapply(1:41, getStratArea, urb=FALSE)
  
  if(FALSE) {
    cbind(normPopsUrb, truePopUrb/trueAreaUrb)
    plot(normPopsUrb, log(truePopUrb/trueAreaUrb))
    plot(normPopsRur, log(truePopRur/trueAreaRur))
  }
  
  out = load("savedOutput/simStudy1/simPopsSurveys.RData")
  aggWeightUrb = sapply(1:41, getStratUrbWeight)
  aggWeightRur = 1-aggWeightUrb
  
  cors = apply(stratumPops, 2, function(x) {
    cor(x, aggWeightUrb * normPopsUrb + aggWeightRur * normPopsRur)
    })
  mean(cors)
  hist(cors)
  abline(v=mean(cors), col="blue")
  
  plot(rep(aggWeightUrb * normPopsUrb + aggWeightRur * normPopsRur, 100), c(stratumPops))
}







