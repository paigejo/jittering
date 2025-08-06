# test script for modM_DSep.R, ... (TMB models)

# test modM_DSepRepar case:

out = fitMD()
system.time(gridPreds <- predGrid(out$TMBsd, popMat=popMatNGAThresh, nsim=5000, admLevel="stratMICS", 
                     quantiles=c(0.025, 0.1, 0.9, 0.975), sep=TRUE))[3]
# 228.571  seconds

save(gridPreds, file="savedOutput/ed/gridPredsM_DSepRepar.RData")
out = load("savedOutput/ed/gridPredsM_DSepRepar.RData")

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPreds, file="savedOutput/ed/stratPredsM_DSepRepar.RData")
save(admin1Preds, file="savedOutput/ed/admin1PredsM_DSepRepar.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsM_DSepRepar.RData")
out = load("savedOutput/ed/stratPredsM_DSepRepar.RData")
out = load("savedOutput/ed/admin1PredsM_DSepRepar.RData")
out = load("savedOutput/ed/admin2PredsM_DSepRepar.RData")

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
          plotNameRoot="edFusionM_DSepRepar")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edFusionM_DSepRepar", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edFusionM_DSepRepar", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edFusionM_DSepRepar", plotNameRootAreal="Admin2")