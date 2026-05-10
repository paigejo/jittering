#!/usr/bin/env Rscript
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
cat("Generating DHS integration points (KDHSurb=16, KDHSrur=21)...\n")

intPtsDHS <- makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban,
                                         areaNames=ed$subarea, popPrior=TRUE,
                                         numPointsUrban=16, numPointsRural=21,
                                         JInnerUrban=4, JInnerRural=4, JOuterRural=1,
                                         adminMap=adm2Full,
                                         outFile="savedOutput/global/intPtsDHS_16_21.RData",
                                         saveOutput=TRUE)

cat("Saved intPtsDHS to savedOutput/global/intPtsDHS_16_21.RData\n")
