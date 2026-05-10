source("code/setup.R")
# DHS-only FE+nugget GH model: fit, predictions and plots

# load datasets ----
load("savedOutput/global/ed.RData")
load("savedOutput/global/edMICS.RData")
edMICS = sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

# fit model ----
Qgh   <- 10
KDHSu <- 16
KDHSr <- 21
fitCacheFile = "savedOutput/global/fitFED_final.RData"
if(file.exists(fitCacheFile)) {
    load(fitCacheFile)
} else {
    res = fitFED(datDHS=ed, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
    if(!dir.exists("savedOutput/global")) dir.create("savedOutput/global", recursive=TRUE)
    save(res, file=fitCacheFile)
}
SD0 = res$TMBsd
obj  = res$TMBobj

# grid predictions ----
# FE+nugget model: no spatial component (epsilon=0).
# predGrid auto-detects noSpatial and produces alpha + X*beta predictions only.
gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=5000, obj=obj,
                     admLevel="stratMICS",
                     quantiles=c(0.025, 0.1, 0.9, 0.975))

# areal predictions ----
stratPreds  = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area",        orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea",     orderedAreas=adm2@data$NAME_2)

# parameter summary table ----
summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, gridPreds=gridPreds)

# plots ----
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=stratPreds,
          plotNameRoot="edFE_D_final", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=admin1Preds,
          plotNameRoot="edFE_D_final", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=admin2Preds,
          plotNameRoot="edFE_D_final", plotNameRootAreal="Admin2")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=NULL,
          plotNameRoot="edFE_D_final")
