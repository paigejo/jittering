# test script for modM_DSep.R, ... (TMB models)

# test modM_DSepRepar case:


modI = 1
mods = c("Md", "M_D", "M_M", "M_DM")
modName = mods[modI]
if(modName == "Md") {
  out = fitMd()
} else if(modName == "M_D") {
  out = fitMD()
} else if(modName == "M_M") {
  out = fitMM()
} else if(modName == "M_DM") {
  out = fitMDM()
}

system.time(gridPreds <- predGrid(out$TMBsd, popMat=popMatNGAThresh, nsim=5000, admLevel="stratMICS", 
                     quantiles=c(0.025, 0.1, 0.9, 0.975), sep=TRUE))[3]
# 228.571  seconds

# save(gridPreds, file=paste0("savedOutput/ed/gridPreds", modName, "SepRepar.RData"))
# out = load(paste0("savedOutput/ed/gridPreds", modName, "SepRepar.RData"))

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
# save(stratPreds, file=paste0("savedOutput/ed/stratPreds", modName, "SepRepar.RData"))
# save(admin1Preds, file=paste0("savedOutput/ed/admin1Preds", modName, "SepRepar.RData"))
# save(admin2Preds, file=paste0("savedOutput/ed/admin2Preds", modName, "SepRepar.RData"))
# out = load(paste0("savedOutput/ed/stratPreds", modName, "SepRepar.RData"))
# out = load(paste0("savedOutput/ed/admin1Preds", modName, "SepRepar.RData"))
# out = load(paste0("savedOutput/ed/admin2Preds", modName, "SepRepar.RData"))

summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, 
               gridPreds=gridPreds)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrr}
# \hline
# & Est & Q0.025 & Q0.975 \\ 
# \hline
# X.Int. & -1.68 & -1.93 & -1.35 \\ 
# beta & 0.47 & 0.24 & 0.70 \\ 
# beta.1 & -0.14 & -0.23 & -0.05 \\ 
# beta.2 & 0.14 & 0.02 & 0.27 \\ 
# beta.3 & 0.09 & -0.03 & 0.21 \\ 
# beta.4 & 0.82 & 0.63 & 1.02 \\ 
# sigmaSq & 0.60 & 0.35 & 0.95 \\ 
# phi & 0.84 & 0.46 & 0.98 \\ 
# sigmaEpsSq & 1.33 & 1.15 & 1.53 \\ 
# \hline
# \end{tabular}
# \end{table}
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=NULL, 
          plotNameRoot=paste0("TMBmodTest", modName, "SepRepar"))
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot=paste0("TMBmodTest", modName, "SepRepar"), plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot=paste0("TMBmodTest", modName, "SepRepar"), plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot=paste0("TMBmodTest", modName, "SepRepar"), plotNameRootAreal="Admin2")