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
  
  predFile = paste0("savedOutput/simStudy1/allPreds", mod, "_SepRepar_", i, ".RData")
  if(regenData || !file.exists(predFile)) {
    # fit TMB model
    if(mod == "M_M") {
      out = fitMM(datDHS=thisDHS, datMICS=thisMICS, repar=TRUE)
    } else if(mod == "M_DM") {
      out = fitMDM(datDHS=thisDHS, datMICS=thisMICS, repar=TRUE)
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
    save(stratPreds, admin1Preds, admin2Preds, parTab, file=predFile)
  } else {
    load(predFile)
  }
  
  # calculate scores
  scoresAdm2 = getScores(thisSubareaPop, estMat=admin2Preds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
  scoresAdm1 = getScores(thisAreaPop, estMat=admin1Preds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
  scoresStratum = getScores(thisStratumPop, estMat=stratPreds$aggregationResults$p, significance=signif, doFuzzyReject=FALSE)
  
  scoreFile = paste0("savedOutput/simStudy1/scores", mod, "_SepRepar_", i, ".RData")
  save(scoresAdm2, scoresAdm1, scoresStratum, parTab)
}


# aggregate/summarize scores and parameter estimates








