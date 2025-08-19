# run motivating simulation study

runSimStudy1I = function(i, mod=c("M_M", "M_DM"), regenData=FALSE, signif=c(.8, .95)) {
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
  dhsIntFile = paste0("savedOutput/simStudy1/intPtsDHS_sim", i, ".RData")
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
    scoresStratum = getScores(thisStratumPop, estMat=stratPreds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
    
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

runSimStudy1IPar = function(i, mod=c("M_M", "M_DM"), regenData=FALSE) {
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

runSimStudy1All = function(doPar=FALSE, nCores=16, regenData=FALSE, mod=c("M_M", "M_DM")) {
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
simStudyScores = function(mod=c("M_M", "M_DM")) {
  require(abind)
  mod = match.arg(mod)
  
  allScores2 = c() # admin2
  allScores1 = c() # admin1
  allScoresS = c() # MICS stratum
  allScoresP = c() # parameters
  # allParEsts = c()
  meanParStats = c()
  for(i in 1:100) {
    # load scores
    out = load(paste0("savedOutput/simStudy1/scores", mod, "_SepRepar_", i, ".RData"))
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
      meanParStats = ((i-1)/100) * meanParStats + (.01) * parTab
    }
  }
  
  # average scores
  meanScoresS = colMeans(allScoresS)
  meanScores1 = colMeans(allScores1)
  meanScores2 = colMeans(allScores2)
  meanScoresP = apply(allScoresP, c(1,2), mean)
  
  # calculate score
  
  # make tables
  browser()
  savedColI = c(1, 4, 5, 9, 13, 17, 18)
  savedColI = c(1, 5, 9, 18)
  savedColI = c(1:5, 7, 9, 11)
  digits = c(rep(3, 6), 2, 3)
  digitsPar = rep(2, 9)
  
  # population prediction scores
  scoreTab = rbind(allScoresS, 
                   allScores1, 
                   allScores2)
  scoreTab = cbind("Areal level"=rep(c("MICS stratum", "Admin1", "Admin2"), each=100), scoreTab)
  
  allMeanScores = rbind(meanScoresS, 
                        meanScores1, 
                        meanScores2)
  row.names(allMeanScores) = c("MICS Stratum scores", "Admin1 scores", "Admin2 scores")
  
  print(xtable(roundCols(allMeanScores[,savedColI], digits=digits), digits=c(0, digits)))
  print(xtable(roundCols(meanScores1[,savedColI], digits=digits), digits=c(0, digits)))
  print(xtable(roundCols(meanScores2[,savedColI], digits=digits), digits=c(0, digits)))
  print(xtable(roundCols(meanScoresS[,savedColI], digits=digits), digits=c(0, digits)))
  
  # parameter scores
  print(xtable(roundCols(meanParStats, digits=digits), digits=c(0, digits)))
  
  print(xtable(roundCols(allScoresP[,savedColI], digits=digits), digits=c(0, digits)))
  
  # make plots?
  browser()
  truePar = c(-1.25, 1, 0, 0, 0, 0.5, 0.5, NA, 1.5)
  
  
  
}







