#!/usr/bin/env Rscript
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
cat("Loading saved DHS intPts and building TMB inputs...\n")

load("savedOutput/global/intPtsDHS_16_21.RData") # loads intPtsDHS

inputsMDM <- makeInputsMDM(datDHS=ed, datMICS=edMICS, intPtsDHS=intPtsDHS, KMICS=100,
                           KDHSurb=16, JInnerUrban=4, KDHSrur=21, JInnerRural=4, saveNewIntPts=FALSE)

cat("Saving TMB inputs to savedOutput/global/edInputs_16_21.RData...\n")
list2env(inputsMDM, envir=.GlobalEnv)
save(AUrbMICS, ARurMICS, AUrbDHS, ARurDHS, intPtsDHS, intPtsMICS, 
     areaidxlocUrbanMICS, areaidxlocRuralMICS, areaidxlocUrbanDHS, areaidxlocRuralDHS,
     ysUrbMICS, nsUrbMICS, ysRurMICS, nsRurMICS,
     ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS,
     file="savedOutput/global/edInputs_16_21.RData")

cat("Saved edInputs_16_21.RData\n")
