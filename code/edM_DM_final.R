source("code/setup.R")
# DHS+MICS fusion BYM2 GH model: fit, predictions and plots

# load datasets ----
out = load("savedOutput/global/ed.RData")
out = load("savedOutput/global/edMICS.RData")
edMICS = sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

# fit model ----
Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
res = fitMDM(datDHS=ed, datMICS=edMICS, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
SD0 = res$TMBsd
obj  = res$TMBobj

# grid predictions ----
# New GH BYM2 uses w_bym2Free (explicit alpha, no random effects).
# predGrid auto-detects this and treats it as constr2n2 parameterization.
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
          plotNameRoot="edM_DM_final", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=admin1Preds,
          plotNameRoot="edM_DM_final", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=admin2Preds,
          plotNameRoot="edM_DM_final", plotNameRootAreal="Admin2")
plotPreds(SD0, obj, popMat=popMatNGAThresh,
          gridPreds=gridPreds, arealPreds=NULL,
          plotNameRoot="edM_DM_final")
